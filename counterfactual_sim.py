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
            for k in self.intervention_packages.keys():
                self.sims_counterfactual[k] = dict()

        self.initialized = True
        return

    def run_baseline(self):
        if not self.initialized:
            self.initialize()
        
        self.sim_baseline.check_already_run()
        
        # run the baseline simulation
        self.sim_baseline.run()
        return
    
    def run_counterfactual(self, intervention_package_key, detection_times, verbose = False):
        for detection_time in detection_times:
            if verbose:
                print(f'Running counterfactual (seed {self.pars_baseline["rand_seed"]}) for intervention package "{intervention_package_key}" with detection time {detection_time}.')
            sim_cf = inf.Sim(
                sc.dcp(self.pars_baseline), # deepcopy is necessary because Sim.run() modifies the parameters
                popfile=self.people_file_temp,
                interventions = sc.dcp(self.intervention_packages[intervention_package_key])
                )
            sim_cf.run()
            self.sims_counterfactual[intervention_package_key][detection_time] = sim_cf
        return
    
    def scan_detection_range(self, pathogen_index = 0, intervention_package_keys = None, n_steps = None, verbose = False):
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
            self.run_counterfactual(package_key, detection_times = d_times, verbose = verbose)

        return
    
class CounterfactualMultiSim(cvb.ParsObj):

    def __init__(self, pars=None, datafile=None, label=None, simfile=None,
                popfile=None, people=None, version=None, intervention_packages=None, n_sims=1, **kwargs):
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

    def run_baseline(self):
        if not self.initialized:
            self.initialize()
        for sim in self.sims:
            sim.run_baseline()
        return
    
    def run_counterfactual(self, intervention_package_key, detection_times, verbose = False):
        if not self.initialized:
            self.initialize()
        for sim in self.sims:
            sim.run_counterfactual(intervention_package_key, detection_times, verbose = verbose)
        return
    
    def scan_detection_range(self, pathogen_index = 0, intervention_package_keys = None, n_steps = None, verbose = False):
        if not self.initialized:
            self.initialize()
        for sim in self.sims:
            sim.scan_detection_range(pathogen_index = pathogen_index, intervention_package_keys = intervention_package_keys, n_steps = n_steps, verbose = verbose)
        return