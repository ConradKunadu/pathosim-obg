'''Measure the burden of disease'''

import numpy as np

__all__ = ['Burden', 'LifeExpectancy']

class Burden:

    def __init__(self, life_expectancy):
        self.yll = None # years of life lost
        self.life_expectancy = life_expectancy
        self.initialized = False
        return
    
    def initialize(self, sim):
        self.n_days = sim['n_days']
        self.yll = np.zeros(self.n_days+1)
        self.initialized = True

    def update(self, sim):
        if not self.initialized:
            raise ValueError("Burden object must be initialized before updating.")
        ages = sim.people.age[sim.people.dead]
        sexes = sim.people.sex[sim.people.dead]
        yll_dead = [self.life_expectancy.get_remaining_years(age, sex) for age, sex in zip(ages, sexes)]
        self.yll[sim.t] = sum(yll_dead)
        return
    
    def get_results(self):
        if not self.initialized:
            raise ValueError("Burden object must be initialized before getting results.")
        results = {
            "new_yll": self.yll,
            "cum_yll": np.cumsum(self.yll)
        }
        return results
    
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