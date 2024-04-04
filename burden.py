'''Measure the burden of disease'''

import numpy as np
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