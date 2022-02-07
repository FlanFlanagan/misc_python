import numpy as np
from scipy import linalg
from matplotlib import pyplot as plt

def calculate_constants():
    '''You will need to update this section of code. To calculate this you will need to create homogeneous
    values for each of the values based on the reactor composition'''
    # Define constants assuming the reactor is homoegeneous these are of the following units
    # d = (1/cm), siga (1/cm)
    con_dict = {}
    num_den = density calculator()
    num_den['U235'] = 
    # Solve for D
    con_dict['d'] = None #replace the none with your calculation

    # Solve for sigma_a 
    con_dict['siga'] = None #replace the none with your calculation

    # assume a given nu of 2.4
    con_dict['nu'] = 2.4

    # solve for sigma_f
    con_dict['sigf'] = None #replace the none with your calculation
    return con_dict

def solve_flux(geo, microxs, fiss_yield, constants):
    '''Also update this section. You will need to build the matrix.'''
    # Build the matrix. You will need to build each row of the matrix and then add 
    # them to this matrix array. Then convert it to a np.array. 
    matrix = []
        #build matrix here
    matrix = np.array(matrix)

    # Construct the identiy matrix for fission. Recall this is our source of neutrons for each 'slice' of the reactor.  
    
    idsigf = None #replace the none with your calculation

    #use the scipi linear algebra package to solve for the eigenvalues and eigenvectors. 
    w, v = linalg.eig(matrix, b=idsigf)
    return w, v