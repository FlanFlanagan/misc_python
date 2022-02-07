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
        self.initial_comp = self.initialize_EN()
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
        for _ in range(self.slices):
            self.rxslices.append(ReactorSlice(self.dx, self.initial_comp, self.rxdata, self.nucs))

    def initialize_EN(self):
        N_U = 18.65/ ((1-self.enrich)*238.0 + (self.enrich)*235.0) * 6.022E23
        N_U235 = self.enrich * N_U * 0.5
        N_U238 = (1.- self.enrich) * N_U * 0.5
        N_H20 = 0.5/18.0 * 6.022E23
        return [N_U235, N_U238, 0, 0, 0, N_H20]

    def build_constant_array(self, time):
        d = [0] * self.slices
        sigma_a = [0] * self.slices
        sigma_f = [0] * self.slices
        i = 0
        for rx in self.rxslices:
            constants = rx.calculate_constants(time)
            d[i] = constants[0]
            sigma_a[i] = constants[1]
            sigma_f[i] = constants[2]
            i+=1
        d.insert(0, d[0])
        d.append(d[-1])
        return d, sigma_a, sigma_f

    def build_diffusion_matrix(self, d, sigma_a):
        matrix = []
        for i in range(1, self.slices+1):
            row = [0] * self.slices
            if i == 1:
                row[i-1] = - (1/(2*self.dx**2)) * (d[i-1] + 2 * d[i] + d[i+1]) - sigma_a[i-1]
                row[i] = (1/(2* self.dx**2)) * (d[i] + d[i+1])
            elif i == self.slices:
                row[i-2] = (1/(2*self.dx**2)) * (d[i-1] + d[i])
                row[i-1] = -(1/(2*self.dx**2)) * (d[i-1] + 2 * d[i] + d[i+1]) - sigma_a[i-1]
            else:
                row[i - 2] = (1 / (2 * self.dx ** 2)) * (d[i-1] + d[i])
                row[i - 1] = -(1 / (2 * self.dx ** 2)) * (d[i-1] + 2 * d[i] + d[i + 1]) - sigma_a[i-1]
                row[i] = (1 / (2 * self.dx ** 2)) * (d[i] + d[i + 1])
            matrix.append(row)
        matrix = np.array(matrix)
        return matrix

    def build_fission_matrix(self, sigma_f):
        sigma_f_m = np.identity(self.slices) * sigma_f * self.nu
        return sigma_f_m

    def calculate_scaled_flux(self, flux):
        if self.scaling_constant == 0: 
            denom = 0.
            for s in self.rxslices:
                denom += s.flux * s.sigma_f
            self.scaling_constant = self.power/(3.204E-17*self.dx*denom)
        else:
            for s in self.rxslices:
                s.flux = s.flux*self.scaling_constant
        return flux*self.scaling_constant

    def calculate_flux(self, time):
        d, sigma_a, sigma_f = self.build_constant_array(time)
        diff_matrix = self.build_diffusion_matrix(d, sigma_a)
        fission_matrix = self.build_fission_matrix(sigma_f)
        w, v = linalg.eig(diff_matrix, b=fission_matrix)
        flux = np.abs(v[:,np.argmax(w)])
        for i in range(self.slices):
            self.rxslices[i].flux = flux[i]
        flux = self.calculate_scaled_flux(flux)
        self.flux.append(flux)

    def calculate_criticality(self):
        neutrons_produced = 0.0
        neutrons_absorbed = 0.0
        for s_slice in self.rxslices:
            neutrons_produced += s_slice.flux * (s_slice.sigma_f * self.nu)
            neutrons_absorbed += s_slice.flux * (s_slice.sigma_f + s_slice.sigma_a)
        criticality = neutrons_produced/neutrons_absorbed
        return criticality

    def burn_reactor(self):
        for t in np.arange(0, self.simulation_time, self.dt):
            self.calculate_flux(t)
            for s in self.rxslices:
                s.compute_reaction_rate_matrix(self.microxs, self.f_yields)
                s.evolution_solver(self.dt, t)  
            self.k.append(self.calculate_criticality())
            if self.calculate_criticality() <= 1.0:
                return
    
    def plot_flux(self):
        i = 0
        for _ in self.rxslices:
            flux = []
            for t in range(len(self.flux)):
                flux.append(self.flux[t][i])
            label = 'flux of slice '+str(i)
            plt.plot(flux, label=label)
            i += 1
        plt.title('Flux by Slice')
        plt.legend()
        plt.show()

    def plot_iso(self, isotope):
        i = 0
        for s_slice in self.rxslices:
            iso_array = []
            for t in range(len(s_slice.evolution_composition)):
                iso_array.append(s_slice.evolution_composition[t*self.dt][isotope])
            label = str(self.nucs[isotope]) + ' number density for slice ' + str(i)
            plt.plot(iso_array, label=label)
            i+=1
        plt.title(self.nucs[isotope]+' number density')
        plt.legend()
        plt.show()

    def plot_flux_by_time(self):
        l = len(self.flux)
        plt.plot(self.flux[0], label='Start')
        plt.plot(self.flux[int(l/5)], label='1L')
        plt.plot(self.flux[int(2*l/5)], label='2L')
        plt.plot(self.flux[int(3*l/5)], label='3L')
        plt.plot(self.flux[int(4*l/5)], label='4L')
        plt.plot(self.flux[-1], label='End')
        #plt.yscale('log')
        plt.title('Flux evolution')
        plt.legend()
        plt.show()


def reactor():
    rxdata = [rx_data.microscopic_cross_sections, rx_data.fission_yield]
    x = Reactor(20, 1, 0.10, 2.4, rxdata, 1000*86400, 86400, 0.001)
    x.burn_reactor()
    return x

#import cProfile, pstats
#profiler = cProfile.Profile()
#profiler.enable()
x = reactor()
#profiler.disable()
#stats = pstats.Stats(profiler).sort_stats('cumtime')
#stats.print_stats()



plt.plot(x.k)
plt.title('Criticality')
plt.show()
x.plot_flux()
x.plot_iso(0)
x.plot_iso(1)
x.plot_iso(2)
x.plot_iso(3)
x.plot_iso(4)

#print(x.initial_comp)
x.plot_flux_by_time()
