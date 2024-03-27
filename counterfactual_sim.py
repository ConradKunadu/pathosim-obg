'''
Specify counterfactual simulations
'''

import numpy as np
import sciris as sc
import datetime as dt
from . import sim as inf
from . import misc as cvm
from . import utils as cvu
from . import base as cvb
from . import version as cvv
from . import defaults as cvd
from . import parameters as cvpar
from collections import defaultdict
import tempfile
import pandas as pd

__all__ = ['CounterfactualSim', 'CounterfactualMultiSim']

class CounterfactualSim(cvb.ParsObj):

    def __init__(self, pars=None, datafile=None, label=None, simfile=None,
                popfile=None, people=None, version=None, intervention_packages=None, **kwargs):
       # Set attributes
        self.label         = label    # The label/name of the simulation
        self.created       = None     # The datetime the sim was created
        self.simfile       = simfile  # The filename of the sim
        self.datafile      = datafile # The name of the data file
        self.popfile       = popfile  # The population file
        self.data          = None     # The actual data
        self.popdict       = people   # The population dictionary
        self._default_ver  = version  # Default version of parameters used
        self.initialized   = False    # Whether or not initialization is complete
        self.people_file_temp = None  # Temporary file for people

        # ensure that intervention packages are a list of lists
        if intervention_packages is not None:
            if not isinstance(intervention_packages, dict):
                raise ValueError(f'intervention_packages must be a dict, not {type(intervention_packages)}.')
            for pkg_key, pkg in intervention_packages.items():
                if not isinstance(pkg, list):
                    raise ValueError(f'Error for package {pkg_key}: an intervention package must be a list of interventions.')
        self.intervention_packages = intervention_packages # The intervention packages

        # sim objects
        self.sim_baseline = None
        self.pars_baseline = None
        self.sims_counterfactual = None
        
        # summary objects
        self.summary_baseline = None
        self.summaries_counterfactual = None

        default_pars = cvpar.make_pars(version=version) # Start with default pars
        super().__init__(default_pars) # Initialize and set the parameters as attributes

        self.set_metadata(simfile) # Set the metadata
        self.update_sim_pars(pars, **kwargs) # Update the parameters, if provided

    def update_sim_pars(self, pars=None, create=False, **kwargs):
        pars = sc.mergedicts(pars, kwargs)
        super().update_pars(pars=pars, create=create)

    def set_metadata(self, simfile, **kwargs):
        ''' Set the metadata for the simulation -- creation time and filename '''
        self.created = kwargs.get('created', sc.now())
        self.version = kwargs.get('version', cvv.__version__)
        self.git_info = kwargs.get('git_info', cvm.git_info())
        if simfile is None:
            self.simfile = 'covasim.sim'
        return
    
    def initialize(self):
        # initialize baseline sim
        self.sim_baseline = inf.Sim(self.pars)
        self.sim_baseline.initialize()
        self.pars_baseline = self.sim_baseline._orig_pars
        
        # store pop file
        # create temporary file using tempfile
        with tempfile.NamedTemporaryFile(delete=False) as fp:
            self.sim_baseline.people.save(fp.name)
            self.people_file_temp = fp.name

        # initialize object for storing counterfactual sims
        if self.intervention_packages is not None:
            self.sims_counterfactual = sc.odict()
            self.summaries_counterfactual = sc.odict()
            for k in self.intervention_packages.keys():
                self.sims_counterfactual[k] = dict()
                self.summaries_counterfactual[k] = dict()

        self.initialized = True
        return

    def run_baseline(self, verbose = False):
        if not self.initialized:
            self.initialize()
        
        if verbose:
            print(f'Running baseline (seed {self.pars_baseline["rand_seed"]}) simulation.')
        self.sim_baseline.check_already_run()
        
        # run the baseline simulation
        self.sim_baseline.run()

        # store the baseline summary
        summary_df = pd.DataFrame((self.sim_baseline.summary[p] for p in range(self.sim_baseline.n_pathogens)))
        summary_df.insert(0, "pathogen", summary_df.index.values)
        self.summary_baseline = summary_df
        return
    
    def run_counterfactual(self, intervention_package_key, detection_times, store_sims = True, verbose = False):
        for detection_time in detection_times:
            if verbose:
                print(f'Running counterfactual (seed {self.pars_baseline["rand_seed"]}) for intervention package "{intervention_package_key}" with detection time {detection_time}.')
            
            # create counterfactual simulation safely (using deep copy where necessary)
            sim_cf = inf.Sim(
                sc.dcp(self.pars_baseline), # deepcopy is necessary because Sim.run() modifies the parameters
                popfile=self.people_file_temp,
                interventions = sc.dcp(self.intervention_packages[intervention_package_key])
                )
            
            # run the counterfactual simulation
            sim_cf.run()

            # store results
            summary_df = pd.DataFrame((sim_cf.summary[p] for p in range(sim_cf.n_pathogens)))
            summary_df.insert(0, "pathogen", summary_df.index.values)
            self.summaries_counterfactual[intervention_package_key][detection_time] = summary_df
            
            if store_sims:
                self.sims_counterfactual[intervention_package_key][detection_time] = sim_cf
        return
    
    def scan_detection_range(self, pathogen_index = 0, intervention_package_keys = None, n_steps = None, store_sims = True, verbose = False):
        if not self.initialized:
            self.initialize()

        if not self.sim_baseline.results_ready:
            self.run_baseline()

        d_range = self.sim_baseline.get_detection_ranges()[pathogen_index]
        if n_steps is None:
            n_steps = d_range["upper"] - d_range["lower"] + 1
        d_times = np.round(np.linspace(d_range['lower'], d_range['upper'], n_steps)).astype(int)

        if intervention_package_keys is None:
            intervention_package_keys = list(self.intervention_packages.keys())
        if not isinstance(intervention_package_keys, list):
            intervention_package_keys = [intervention_package_keys]

        for package_key in intervention_package_keys:
            self.run_counterfactual(package_key, detection_times = d_times, store_sims = store_sims, verbose = verbose)

        return
    
    def get_summaries_df(self):
        def label_df(df, seed, intervention_package, delay):
            df = df.copy()
            df.insert(0, "seed", seed)
            df.insert(1, "intervention_package", intervention_package)
            df.insert(2, "delay", delay)
            return df

        summaries_df = pd.concat(
            [
                label_df(self.summary_baseline, seed = self.pars_baseline["rand_seed"], intervention_package = "baseline", delay = 0),
                *[label_df(df, seed = self.pars_baseline["rand_seed"], intervention_package = interv, delay = d) for interv, dictres in self.summaries_counterfactual.items() for d, df in dictres.items()]
            ]
        )
        return summaries_df
    
class CounterfactualMultiSim(cvb.ParsObj):

    def __init__(self, pars=None, datafile=None, label=None, simfile=None,
                popfile=None, people=None, version=None, intervention_packages=None, n_sims=1, maxcpu = 0.9, maxmem = 0.9, interval = 0.5, **kwargs):
         # Set attributes
        self.pars          = pars
        self.datafile      = datafile # The name of the data file
        self.label         = label    # The label/name of the simulation
        self.simfile       = simfile  # The filename of the sim
        self.popfile       = popfile  # The population file
        self.people        = people  # The population dictionary
        self.version       = version  # version
        self.intervention_packages = intervention_packages
        self.n_sims        = n_sims # number of simulations
        self.sims          = list()
        self.kwargs        = kwargs

        self.set_opts_parallelization(maxcpu, maxmem, interval)

        self.initialized   = False    # Whether or not initialization is complete
        return
    
    def initialize(self, n_sims = None):
        if not self.initialized:
            if n_sims is not None:
                self.n_sims = n_sims
            for i in range(self.n_sims):
                simpars = sc.dcp(self.pars)
                simpars['rand_seed'] = simpars['rand_seed'] + i
                sim = CounterfactualSim(pars = simpars, datafile=self.datafile, label=self.label, simfile=self.simfile, popfile=self.popfile, people=sc.dcp(self.people), version=self.version, intervention_packages=sc.dcp(self.intervention_packages), **self.kwargs)
                sim.initialize()
                self.sims.append(sim)
            self.initialized = True
        return
    
    def set_opts_parallelization(self, maxcpu = 0.9, maxmem = 0.9, interval = 0.5):
        self.maxcpu = maxcpu
        self.maxmem = maxmem
        self.interval = interval
        return

    def run_baseline(self, maxcpu = None, maxmem = None, verbose = False, interval = None):
        if not self.initialized:
            self.initialize()
        
        def run_indi_sim(indi_sim, verbose):
            indi_sim.run_baseline(verbose = verbose)
            return indi_sim
        
        if maxcpu is None:
            maxcpu = self.maxcpu
        if maxmem is None:
            maxmem = self.maxmem
        if interval is None:
            interval = self.interval

        self.sims = sc.parallelize(run_indi_sim, self.sims, verbose = verbose, maxcpu = maxcpu, maxmem = maxmem, interval = interval)
        return
    
    def run_counterfactual(self, intervention_package_key, detection_times, store_sims = True, verbose = False, maxcpu = None, maxmem = None, interval = None):
        if not self.initialized:
            self.initialize()

        def run_indi_sim(indi_sim, intervention_package_key, detection_times, store_sims, verbose):
            indi_sim.run_counterfactual(intervention_package_key, detection_times, store_sims = store_sims, verbose = verbose)
            return indi_sim
        
        if maxcpu is None:
            maxcpu = self.maxcpu
        if maxmem is None:
            maxmem = self.maxmem
        if interval is None:
            interval = self.interval
        
        self.sims = sc.parallelize(run_indi_sim, self.sims, intervention_package_key = intervention_package_key, detection_times = detection_times, store_sims = store_sims, verbose = verbose, maxcpu = maxcpu, maxmem = maxmem, interval = interval)
        return
    
    def scan_detection_range(self, pathogen_index = 0, intervention_package_keys = None, n_steps = None, store_sims = True, verbose = False, maxcpu = None, maxmem = None, interval = None):
        if not self.initialized:
            self.initialize()

        def run_indi_sim(indi_sim, pathogen_index, intervention_package_keys, n_steps, store_sims, verbose):
            indi_sim.scan_detection_range(pathogen_index = pathogen_index, intervention_package_keys = intervention_package_keys, n_steps = n_steps, store_sims = store_sims, verbose = verbose)
            return indi_sim
        
        if maxcpu is None:
            maxcpu = self.maxcpu
        if maxmem is None:
            maxmem = self.maxmem
        if interval is None:
            interval = self.interval
        
        self.sims = sc.parallelize(run_indi_sim, self.sims, pathogen_index = pathogen_index, intervention_package_keys = intervention_package_keys, n_steps = n_steps, store_sims = store_sims, verbose = verbose, maxcpu = maxcpu, maxmem = maxmem, interval = interval)
        return
    
    def get_summaries_df(self):
        summaries_df = pd.concat([sim.get_summaries_df() for sim in self.sims])
        return summaries_df