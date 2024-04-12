'''
Lorem ipsum
'''

"""
import numpy as np
import pandas as pd
import scipy as sp
import pylab as pl
import sciris as sc
import inspect
import json
import datetime as dt
from . import misc as cvm
from . import utils as cvu
from . import base as cvb
from . import defaults as cvd
from . import parameters as cvpar
from . import immunity as cvi
from . import pathogens as pat 
from collections import defaultdict
from . import stratify as strat
"""
from abc import ABC, abstractmethod


__all__ = ['FixedEvent', 'DelayedTrigger', 'RatioTrigger', 'RateTrigger']


# TODO consider making factories for events


class Event(ABC):
    """
    Lorem ipsum
    """

    @abstractmethod
    def __init__(self, name):
        self.event_name = name

    def get_event_name(self):
        return self.event_name

    def add_event(self, event_list, t):
        if t in range(len(event_list)):
            event_list[t].append(self.get_event_name())

    @abstractmethod
    def update(self):
        pass



class StartStopEvent(Event, ABC):
    """
    Lorem ipsum

    TODO: multiple pathogens?
    """

    @abstractmethod
    def __init__(self, name, start_delay=0, stop_delay=0):
        super().__init__(name)
        self.start_token = 'START ' # TODO can these maybe be made static instead of having to be instantiated for each object?
        self.stop_token = 'STOP '
        self.start_delay = start_delay
        self.stop_delay = stop_delay

    def get_start_event_name(self):
        return self.start_token + self.event_name

    def get_stop_event_name(self):
        return self.stop_token + self.event_name

    def add_start_event(self, event_list, t):
        if t in range(len(event_list)):
            event_list[t].append(self.get_start_event_name())

    def add_stop_event(self, event_list, t):
        if t in range(len(event_list)):
            event_list[t].append(self.get_stop_event_name())



class FixedEvent(Event):
    """
    Event that is set for a fixed, pre-defined time
    """

    def __init__(self, name, t_event):
        super().__init__(name)
        self.t_event = t_event
        self.initialized = False

    def update(self, sim):
        if not self.initialized:
            self.add_event(sim.events, self.t_event)
            self.initialized = True


class DelayedTrigger(StartStopEvent):
    """
    Triggers a start (and stop) event at a certain time after some other event
    """

    def __init__(self, name, start_delay=0, stop_delay=None, trigger_event=None):
        super().__init__(name, start_delay, stop_delay)
        self.trigger_event = trigger_event

    def update(self, sim):
        if self.trigger_event in sim.events[sim.t]:
            self.add_start_event(sim.events, sim.t+self.start_delay)
            if self.stop_delay is not None:
                self.add_stop_event(sim.events, sim.t+self.stop_delay)



class RatioTrigger(StartStopEvent):
    """
    Triggers a start (and stop) event at a certain ratio of infections of the alive population

    TODO: ratio by region
    TODO: allow to adjust if counting e.g. infectious, exposed or symptomatic
    TODO: allow for multiple pathogens
    """

    def __init__(self, name, start_delay=0, stop_delay=0, start_threshold=1, stop_threshold=None, active=False):
        super().__init__(name, start_delay, stop_delay)
        self.start_threshold = start_threshold
        self.stop_threshold = stop_threshold
        self.active = active
        if stop_threshold > start_threshold:
            raise Exception("Stop threshold has been defined above start threshold")

    # TODO one should probably make sure that a stop event is not placed before a start event by 
    # having conditions on start_delay and stop_delay
    def update(self, sim):
        n_susc = sim.results[0]['n_susceptible'][sim.t-1]
        if not self.active:
            if n_susc == 0:
                ratio = 0
            else:
                ratio = sim.results[0]['n_infectious'][sim.t-1] / n_susc
            if ratio > self.start_threshold:
                self.add_stop_event(sim.events, sim.t+self.stop_delay)
                self.active = True
        else:
            if n_susc == 0:
                ratio = 1
            else:
                ratio = sim.results[0]['n_infectious'][sim.t-1] / n_susc
            if ratio < self.stop_threshold:
                self.add_stop_event(sim.events, sim.t+self.stop_delay)
                self.active = False
        
        

class RateTrigger(StartStopEvent):
    """
    Triggers a start (and stop) event at a certain rate of increase in cases

    TODO: rate by region
    TODO: allow for multiple pathogens
    """

    def __init__(self, name, start_delay=0, stop_delay=0, start_threshold=1, stop_threshold=None, active=False):
        super().__init__(name, start_delay, stop_delay)
        self.start_threshold = start_threshold
        self.stop_threshold = stop_threshold
        self.active = active
        if stop_threshold > start_threshold:
            raise Exception("Stop threshold has been defined above start threshold")

    # TODO one should probably make sure that a stop event is not placed before a start event by 
    # having conditions on start_delay and stop_delay
    def update(self, sim):
        n_susc = sim.results[0]['n_susceptible'][sim.t-1]
        if not self.active:
            if n_susc == 0:
                ratio = 0
            else:
                ratio = sim.results[0]['new_infections'][sim.t - 1] / n_susc
            if ratio > self.start_threshold:
                self.add_stop_event(sim.events, sim.t+self.stop_delay)
                self.active = True
        else:
            if n_susc == 0:
                ratio = 1
            else:
                ratio = sim.results[0]['new_infections'][sim.t - 1] / n_susc
            if ratio < self.stop_threshold:
                self.add_stop_event(sim.events, sim.t+self.stop_delay)
                self.active = False