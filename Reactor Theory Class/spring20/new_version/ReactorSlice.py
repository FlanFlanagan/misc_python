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
        E_S = 0.
        for i in range(len(self.evolution_composition[time])):
            E_S += self.evolution_composition[time][i]*self.microxs[self.nuclides[i]]['xs_scattering']*1.E-24
        d = 1/(3*E_S)
        self.d = d

    def calculate_sigma_a(self, time):
        sigma_a = 0.
        for i in range(len(self.evolution_composition[time])):
            sigma_a += self.evolution_composition[time][i]*self.microxs[self.nuclides[i]]['xs_absorption']*1.E-24
        self.sigma_a = sigma_a

    def calculate_sigma_f(self, time):
        sigma_f = 0.
        for i in range(len(self.evolution_composition[time])):
            sigma_f += self.evolution_composition[time][i]*self.microxs[self.nuclides[i]]['xs_fission']*1.E-24
        self.sigma_f = sigma_f

    def calculate_constants(self, time):
        self.calculate_d(time)
        self.calculate_sigma_a(time)
        self.calculate_sigma_f(time)
        return [self.d, self.sigma_a, self.sigma_f]

    def compute_reaction_rate_matrix(self, xs, fy):
        rrm = np.zeros((6, 6))
        # U235
        rrm[0, 0] = -self.flux * (xs['U235']['xs_absorption'] + xs['U235']['xs_fission']) * 1.E-24
        # U238
        rrm[1, 1] = -self.flux * xs['U238']['xs_absorption'] * 1.E-24
        # Pu239
        rrm[2, 2] = -self.flux * (xs['Pu239']['xs_absorption'] + xs['Pu239']['xs_fission']) * 1.E-24
        rrm[2, 1] = self.flux * xs['U238']['xs_absorption'] * 1.E-24
        # Xe135
        rrm[3, 3] = -self.flux * xs['Xe135']['xs_absorption'] * 1.E-24 - xs['Xe135']['lambda']
        rrm[3, 0] = self.flux * fy['Xe135'] * xs['U235']['xs_fission'] * 1.E-24
        rrm[3, 2] = self.flux * fy['Xe135'] * xs['Pu239']['xs_fission'] * 1.E-24
        # Lumped fission products
        rrm[4, 0] = self.flux * (2 - fy['Xe135']) * xs['U235']['xs_fission'] * 1.E-24
        rrm[4, 2] = self.flux * (2 - fy['Xe135']) * xs['Pu239']['xs_fission'] * 1.E-24
        rrm[4, 3] = self.flux * xs['Xe135']['xs_absorption'] * 1.E-24 + xs['Xe135']['lambda']
        self.rrm = rrm

    def evolution_solver(self, dt, time):
        rrm_id = np.identity(len(self.evolution_composition[time])) - (dt * self.rrm)
        # Ax = B --> A : rrm_id, B : previous step composition
        compositions = np.linalg.solve(rrm_id, self.evolution_composition[time])
        self.evolution_composition[time+dt] = compositions

    

