import numpy as np
import matplotlib.pyplot as plt
import rx_data
import math

from scipy import linalg
from ReactorSlice import ReactorSlice

class Reactor:
    def __init__(self, slices, reactor_size, enrichment, nu, rxdata, simulation_time, dt, power, shield=False):
        self.reactor_size = reactor_size
        self.dx = float(reactor_size/(slices-1))
        self.enrich = enrichment
        self.slices = slices
        self.time = 0
        self.rxdata = rxdata
        self.microxs = rxdata[0]
        self.f_yields = rxdata[1]
        self.shielding = shield
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
        self.d = []
        self.siga = []
        self.sigf = []

    def build_reactor_slices(self):
        if self.shielding is False:
            for _ in range(self.slices):
                self.rxslices.append(ReactorSlice(self.dx, self.initial_comp, self.rxdata, self.nucs))
        if self.shielding is True:
            l = math.floor(self.slices/2)
            for _ in range(l):
                self.rxslices.append(ReactorSlice(self.dx, [0, 0, 0, 0, 0, 1/18.0 * 6.022E23], self.rxdata, self.nucs))     
            for _ in range(self.slices):
                self.rxslices.append(ReactorSlice(self.dx, self.initial_comp, self.rxdata, self.nucs))
            for _ in range(l):
                self.rxslices.append(ReactorSlice(self.dx, [0, 0, 0, 0, 0, 1/18.0 * 6.022E23], self.rxdata, self.nucs))                 

    def initialize_EN(self):
        N_U = 18.65/ ((1-self.enrich)*238.0 + (self.enrich)*235.0) * 6.022E23
        N_U235 = self.enrich * N_U * 0.5
        N_U238 = (1.- self.enrich) * N_U * 0.5
        N_H20 = 0.5/18.0 * 6.022E23
        return [N_U235, N_U238, 0, 0, 0, N_H20]

    def build_constant_array(self, time):
        d = [0] * len(self.rxslices)
        sigma_a = [0] * len(self.rxslices)
        sigma_f = [0] * len(self.rxslices)
        i = 0
        for rx in self.rxslices:
            constants = rx.calculate_constants(time)
            d[i] = constants[0]
            sigma_a[i] = constants[1]
            sigma_f[i] = constants[2]
            i+=1
        d.insert(0,0)
        d.append(0)
        self.d = np.array(d)
        self.siga = np.array(sigma_a)
        self.sigf = np.array(sigma_f)
        return d, sigma_a, sigma_f

    def build_diffusion_matrix(self, d, sigma_a):
        matrix = []
        for i in range(1, len(self.rxslices)+1):
            row = [0] * len(self.rxslices)
            if i == 1:
                row[i-1] = - (1/(2*self.dx**2)) * (d[i-1] + 2 * d[i] + d[i+1]) - sigma_a[i-1]
                row[i] = (1/(2* self.dx**2)) * (d[i] + d[i+1])
            elif i == len(self.rxslices):
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
        sigma_f_m = np.identity(len(self.rxslices)) * sigma_f * self.nu
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
        w = w.real
        w[w > 0 ] = -1e299
        flux = np.abs(v[:,np.argmax(w)])
        for i in range(len(self.rxslices)):
            self.rxslices[i].flux = flux[i]
        flux = self.calculate_scaled_flux(flux)
        self.flux.append(flux)

    def calculate_criticality(self, time):
        if self.shielding is True:
            l = math.floor(self.slices/2)
            neutrons_produced = sum(self.flux[time][l:l+self.slices] * (self.sigf[l:l+self.slices] * self.nu))*self.dx**2
            neutrons_absorbed = sum(self.flux[time][l:l+self.slices] * (self.siga[l:l+self.slices]))*self.dx**2
            n_p = (self.flux[time][l] + self.flux[time][l+1])*2*self.dx*self.d[l+1]
        else:
            neutrons_produced = self.flux * (self.sigf * self.nu)
            neutrons_absorbed = self.flux * (self.siga)
            n_p = (self.flux[time][-1] + 0)*2*self.dx*self.d[-1]
        print(neutrons_absorbed)
        print(neutrons_produced)
        print(n_p)
        criticality = neutrons_produced - neutrons_absorbed-n_p
        self.k.append(criticality)
        return criticality

    def burn_reactor(self):
        i = 0
        for t in np.arange(0, self.simulation_time, self.dt):
            self.calculate_flux(t)
            for s in self.rxslices:
                s.compute_reaction_rate_matrix(self.microxs, self.f_yields)
                s.evolution_solver(self.dt, t)  
            self.k.append(self.calculate_criticality(i))
            i+=1
    
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
    x = Reactor(200, 1, 0.0373023571410201, 2.4, rxdata, 1*86400, 86400, 0.001, shield=True)
    #x = Reactor(10, 1, 0.007, 2.4, rxdata, 500*86400, 86400, 0.001, shield=True)
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
#x.plot_flux()
#x.plot_iso(0)
#x.plot_iso(1)
#x.plot_iso(2)
#x.plot_iso(3)
#x.plot_iso(4)

#print(x.initial_comp)
#x.plot_flux_by_time()
