import behaviour as bh
import pathosim as inf 
import numpy as np
 
# ---- DEMO SHOWCASING HOW TO CONFIGURE A GENERIC PATHOGEN'S IMMUNITY (USING THE GENERALIZED IMMUNITY SYSTEM) ---- #

pop_size = 20000 
  
#Create the circulating pathogen
MyNewPathogen = inf.SARS_COV_2(50)   

#Configure the immunity of that pathogen. 4 parameters have to be passed in this function. See Pathosim Documentation.
MyNewPathogen.configure_generalized_immunity(0.35, 0.95, 180, 14)  #Here, we pass parameters such that the spread of MyNewPathogen ressembles COVID-19's spread.

my_intervention = inf.intervention_bucket()
detection_event = inf.FixedEvent('detection', 30)
delayed_event1 = inf.DelayedTrigger('mask mandate', start_delay=0, stop_delay=50, trigger_event='detection')
delayed_event2 = inf.DelayedTrigger('quarantine', start_delay=10, stop_delay=20, trigger_event='detection')
ratio_event = inf.RatioEvent('social distancing', start_delay=1, stop_delay=7, start_threshold=0.02, stop_threshold=0.01)

all_events = [delayed_event1, delayed_event2, ratio_event]

#Initialize and run the simulation
sim = inf.Sim(pop_size = pop_size, pathogens = [MyNewPathogen], interventions = [my_intervention], events = all_events, detection_events = [detection_event], n_days = 365)

#sim.events = {40:['detection']}

sim.run() 
sim.plot() 

print(sim.events)
