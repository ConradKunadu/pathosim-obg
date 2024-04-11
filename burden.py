'''Measure the burden of disease'''

import numpy as np

__all__ = ['Burden', 'LifeExpectancy']

class Burden:

    def __init__(self, life_expectancy = None):
        self.results = dict()
        self.result_keys = ["yll", "hrql_loss"]
        self.life_expectancy = life_expectancy
        self.initialized = False
        return

    def initialize(self, sim):
        if self.initialized:
            raise ValueError("Burden object already initialized.")
        self.n_days = sim['n_days']

        # burden result metrics
        self.results = {p: {
            "new_yll": np.full(self.n_days+1, np.nan), # years of life lost
            "new_hrql_loss": np.full(self.n_days+1, np.nan) # health-related quality of life loss
        } for p in range(len(sim.pathogens))}

        self.projections = {p: {
            "proj_yll": np.nan, # projected years of life lost
            "proj_hrql_loss": np.nan # projected health-related quality of life loss
        } for p in range(len(sim.pathogens))}

        self.initialized = True

    def update(self, sim):
        if not self.initialized:
            raise ValueError("Burden object must be initialized before updating.")
        
        self.compute_results(sim)

        # compute cumulative results
        for p in range(len(sim.pathogens)):
            for key in self.result_keys:
                self.results[p]["cum_" + key] = np.cumsum(self.results[p]["new_" + key])

        return
    
    def compute_results(self, sim):
        """
        compute burden for today
        """
        for p in range(len(sim.pathogens)):
            self.results[p]["new_yll"][sim.t] = self.compute_yll(sim, sim.t, p)
            self.results[p]["new_hrql_loss"][sim.t] = self.compute_hrql_loss(sim, sim.t, p)

        return
    
    def compute_projections(self, sim):
        """
        compute projected burden after today
        """
        for p in range(len(sim.pathogens)):
            # determine latest date of projections
            if np.all(np.isnan(sim.people.date_p_recovered[p])):
                last_t_recovered = sim["n_days"]
            else:    
                last_t_recovered = np.nanmax(sim.people.date_p_recovered[p])

            if np.all(np.isnan(sim.people.date_p_dead[p])):
                last_t_dead = sim["n_days"]
            else:
                last_t_dead = np.nanmax(sim.people.date_p_dead[p])

            last_t = int(np.nanmax([last_t_recovered, last_t_dead]))

            # initialize projections
            yll_dead = 0
            hrql_loss = 0

            # compute projections
            for t in range(sim.t+1, last_t+1):
                yll_dead += self.compute_yll(sim, t, p)
                hrql_loss += self.compute_hrql_loss(sim, t, p)

            self.projections[p]["proj_yll"] = yll_dead
            self.projections[p]["proj_hrql_loss"] = hrql_loss
        return

    def compute_yll(self, sim, t, p = 0):

        if self.life_expectancy is None:
            return np.nan

        inds = sim.people.date_p_dead[p] == t
        if np.all(inds == False):
            return 0
        ages = sim.people.age[inds]
        sexes = sim.people.sex[inds]
        yll_dead = [self.life_expectancy.get_remaining_years(age, sex) for age, sex in zip(ages, sexes)]
        return sum(yll_dead)
    
    def compute_hrql_loss(self, sim, t, p):

        if sim.pathogens[p].hrql_loss is None or any([v is None for v in sim.pathogens[p].hrql_loss.values()]):
            return np.nan

        before_symptomatic = (t < sim.people.date_p_symptomatic[p])
        before_severe = (t < sim.people.date_p_severe[p])
        before_critical = (t < sim.people.date_p_critical[p])
        before_recovered = (t < sim.people.date_p_recovered[p])
        before_dead = (t < sim.people.date_p_dead[p])
        not_symptomatic = np.isnan(sim.people.date_p_symptomatic[p])
        not_severe = np.isnan(sim.people.date_p_severe[p])
        not_critical = np.isnan(sim.people.date_p_critical[p])
        not_recovered = np.isnan(sim.people.date_p_recovered[p])
        not_dead = np.isnan(sim.people.date_p_dead[p])

        symptomatic_now = ~(before_symptomatic | not_symptomatic) & (before_severe | not_severe) & (before_critical | not_critical) & (before_dead | not_dead) & (before_recovered | not_recovered)
        hrql_loss_symptomatic = np.sum(symptomatic_now) * sim.pathogens[p].hrql_loss["symptomatic"]
        
        severe_now = ~(before_severe | not_severe) & (before_critical | not_critical) & (before_dead | not_dead) & (before_recovered | not_recovered)
        hrql_loss_severe = np.sum(severe_now) * sim.pathogens[p].hrql_loss["severe"]
        
        critical_now = ~(before_critical | not_critical) & (before_dead | not_dead) & (before_recovered | not_recovered)
        hrql_loss_critical = np.sum(critical_now) * sim.pathogens[p].hrql_loss["critical"]
    
        return(hrql_loss_symptomatic + hrql_loss_severe + hrql_loss_critical)

class LifeExpectancy:

    def __init__(self, years_remaining_males, years_remaining_females):
        """
        Initialize the life expectancy object

        Parameters
        ----------
        years_remaining_males: numpy array with expected remaining male life years for each age, starting from age 0
        years_remaining_females: numpy array with expected remaining male life years for each age, starting from age 0
        """
        self.years_remaining_males = years_remaining_males
        self.years_remaining_females = years_remaining_females
        return
    
    def get_life_expectancy(self, age, sex):
        """
        Get the life expectancy for a given age and sex
        """
        
        if sex == "male":
            if age < len(self.years_remaining_males):
                return age + self.years_remaining_males[int(age)]
            else:
                return age + self.years_remaining_males[-1]
        elif sex == "female":
            if age < len(self.years_remaining_females):
                return age + self.years_remaining_females[int(age)]
            else:
                return age + self.years_remaining_females[-1]
        else:
            raise ValueError("Sex is binary here, must be either 'male' or 'female'.")
        
    def get_remaining_years(self, age, sex):
        """
        Get the expected number of remaining life years for a given age and sex
        """
        
        if sex == 0: # male
            if age < len(self.years_remaining_males):
                return self.years_remaining_males[int(age)]
            else:
                return self.years_remaining_males[-1]
        elif sex == 1: # female
            if age < len(self.years_remaining_females):
                return self.years_remaining_females[int(age)]
            else:
                return self.years_remaining_females[-1]
        else:
            raise ValueError("Sex is binary here, must be either 0 (for male) or 1 (for female).")