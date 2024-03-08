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

__all__ = ['CounterfactualSim']

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
        self.intervention_packages = intervention_packages # The intervention packages
        self.initialized   = False    # Whether or not initialization is complete
        self.people_file_temp = None  # Temporary file for people

        # sim objects
        self.sim_baseline = None

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
        self.sim_baseline.init_people()
        
        # store pop file
        # create temporary file using tempfile
        with tempfile.NamedTemporaryFile(delete=False) as fp:
            self.sim_baseline.people.save(fp.name)
            self.people_file_temp = fp.name
        return

    def run_baseline(self):
        # run the baseline simulation
        self.sim_baseline.run()
        return
