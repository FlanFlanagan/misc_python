import numpy as np

class ReactorSlice:
    def __init__(self, dx, compositions, rxdata, nuclides):
        self.evolution_composition = {}
        self.evolution_composition[0] = compositions
        self.dx = dx
        self.microxs = rxdata[0]
        self.yields = rxdata[1]
        self.rrm = self.compute_reaction_rate_matrix
        self.nuclides = nuclides
        self.d = None
        self.sigma_a = None
        self.sigma_f = None
        self.flux = 0.
        
    def calculate_d(self, time):


    def calculate_sigma_a(self, time):


    def calculate_sigma_f(self, time):


    def calculate_constants(self, time):


    def compute_reaction_rate_matrix(self, xs, fy):


    def evolution_solver(self, dt, time):

    

