from scipy import linalg
from matplotlib import pyplot as plt

import eig_stub as es
import rx_data as rx

''' Define the geometry. Set up the width and stepsize of your geometry.'''
system_size = None
deltax = None
geo = [size, deltax]

def read_xs():
    microxs = rx.microscopic_cross_sections
    fiss_yield = rx.fission_yield
    return microxs, fiss_yield

def main():
    #read in information
    microxs, fiss_yield = read_xs()
    #calculate the constants
    consts = es.calculate_constants()
    #calculate the flux (the eigen vector from the correct eigen value)
    eigen_values, eigen_vectors = es.solve_flux(geo, microxs, fiss_yield, consts)
    #plot the correct eigen vector using plt.plot()
    

    '''
    This time loop will be used later on to facilitate the burnup section of the code. The homework that will
    put this together will have more instruction for doing this. 
    delta_t = None
    simulation_time = None
    time_steps = [0]
    evolution_flux = []
    evolution_comp = []
    for t in time_steps:
        # initial values
        ### solve_flux() this is the eigen value code
        ### burn_up() this is your burnup code (when it's done)
        ### update number density
    '''

if __name__ == "__main__":
    main()