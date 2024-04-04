'''Measure the burden of disease'''

import numpy as np

__all__ = ['Burden', 'LifeExpectancy']

class Burden:

    def __init__(self, life_expectancy):
        self.results = None 
        self.result_keys = ["yll"]
        self.life_expectancy = life_expectancy
        self.initialized = False
        return
    
    def initialize(self, sim):
        if self.initialized:
            raise ValueError("Burden object already initialized.")
        self.n_days = sim['n_days']

        # burden result metrics
        self.results = {
            "new_yll": np.zeros(self.n_days+1), # years of life lost
            "proj_yll": np.zeros(self.n_days+1) # projected years of life lost
        }
        
        self.initialized = True

    def update(self, sim):
        if not self.initialized:
            raise ValueError("Burden object must be initialized before updating.")
        
        self.compute_results(sim)
        self.compute_projections(sim)

        # compute cumulative results
        for key in self.result_keys:
            self.results["cum_" + key] = np.cumsum(self.results["new_" + key]) + self.results["proj_" + key]

        return
    
    def compute_results(self, sim):
        """
        compute burden for today
        """
        died_today = sim.people.date_dead == sim.t
        self.results["new_yll"][sim.t] = self.compute_yll(sim, died_today)

        return
    
    def compute_projections(self, sim):
        """
        compute projected burden after today
        """

        # determine latest date of projections
        if np.all(np.isnan(sim.people.date_dead)):
            last_t = sim["n_days"]
        else:    
            last_t = int(np.nanmax(sim.people.date_dead))

        # initialize projections
        yll_dead = 0

        # compute projections
        for t in range(sim.t+1, last_t+1):
            died_today = sim.people.date_dead == t
            yll_dead += self.compute_yll(sim, died_today)

        self.results["proj_yll"][sim.t] = yll_dead
        return

    def compute_yll(self, sim, inds):
        if np.all(inds == False):
            return 0
        ages = sim.people.age[inds]
        sexes = sim.people.sex[inds]
        yll_dead = [self.life_expectancy.get_remaining_years(age, sex) for age, sex in zip(ages, sexes)]
        return sum(yll_dead)
    
class LifeExpectancy:

    def __init__(self, life_remaining_males, life_remaining_females):
        """
        Initialize the life expectancy object

        Parameters
        ----------
        life_remaining_males: numpy array with expected remaining male life years for each age, starting from age 0
        life_remaining_females: numpy array with expected remaining male life years for each age, starting from age 0
        """
        self.life_remaining_males = life_remaining_males
        self.life_remaining_females = life_remaining_females
        return
    
    def get_life_expectancy(self, age, sex):
        """
        Get the life expectancy for a given age and sex
        """
        
        if sex == "male":
            if age < len(self.life_remaining_males):
                return age + self.life_remaining_males[int(age)]
            else:
                return age + self.life_remaining_males[-1]
        elif sex == "female":
            if age < len(self.life_remaining_females):
                return age + self.life_remaining_females[int(age)]
            else:
                return age + self.life_remaining_females[-1]
        else:
            raise ValueError("Sex is binary here, must be either 'male' or 'female'.")
        
    def get_remaining_years(self, age, sex):
        """
        Get the expected number of remaining life years for a given age and sex
        """
        
        if sex == 0: # male
            if age < len(self.life_remaining_males):
                return self.life_remaining_males[int(age)]
            else:
                return self.life_remaining_males[-1]
        elif sex == 1: # female
            if age < len(self.life_remaining_females):
                return self.life_remaining_females[int(age)]
            else:
                return self.life_remaining_females[-1]
        else:
            raise ValueError("Sex is binary here, must be either 0 (for male) or 1 (for female).")