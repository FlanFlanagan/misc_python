import numpy as np
import matplotlib.pyplot as plt
import rx_data

from scipy import linalg
from ReactorSlice import ReactorSlice

class Reactor:
    def __init__(self, slices, reactor_size, enrichment, nu, rxdata, simulation_time, dt, power):
        self.reactor_size = reactor_size
        self.dx = float(reactor_size/(slices-1))
        self.enrich = enrichment
        self.slices = slices
        self.time = 0
        self.rxdata = rxdata
        self.microxs = rxdata[0]
        self.f_yields = rxdata[1]
        self.nucs = list(self.microxs.keys())
        self.nu = nu
        self.rxslices = []
        self.flux = []
        self.simulation_time = simulation_time
        self.dt = dt
        self.power = power
        self.build_reactor_slices()
        self.scaling_constant = 0.
        self.k = []

    def build_reactor_slices(self):
        #initializes and stores the reactors slices

    def build_constant_array(self, time):
        #builds the arrays used for the finite differencing

    def build_diffusion_matrix(self, d, sigma_a):
        #builds the diff_matrix

    def build_fission_matrix(self, sigma_f):
        #builds the fission matrix

    def calculate_scaled_flux(self, flux):
        #scales the flux to match power demand

    def calculate_flux(self, time):
        #performs the required operations to solve the diffusion matrix and scale the flux

    def calculate_criticality(self):
        #calculates K

    def burn_reactor(self):
        #pushes the reactor forward in time
