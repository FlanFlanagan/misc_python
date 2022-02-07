# Define the geometry --> This comes from the eigenvalue solver
system_size = None
deltax = None

# Define the initial material --> This comes from the burnup solver 
N = None  # number densities  --> N = [n1, n2, n3, n4, n5]
xs = None  # microscopic cross section
lambd = None  # lambda --> library

# Define time steps and period --> This comes from the burnup solver
delta_t = None
simulation_time = None
time_steps = []

# Now, we have our initial parameters --> We need to start the time loop
evolution_flux = [] # stores the evolution of the flux over time. This should be an array of arrays. 
evolution_comp = [] # stores the evolution of each isotope over time, this should be an array of arrays (or dictionaries)
for t in time_steps:
    # Calculate the flux:
    # Get material properties (diffusion coefficient, macroscopic cross sections)
    diffusion_coefficient = None
    macroscopic_cross_sections = None
    # Using what you've done in the eigenvalue solver: Build the matrix
    matrix_diffusion_equation = None
    right_hand_side = None
    # Solve the eigenvalue problem to get the flux
    w, v = np.linalg.eig(matrix_diffusion_equation, b=right_hand_side)
    # Scale the flux relative to a given power, 1 MW/m3
    # Your flux is not homogeneous, it varies depending on which node you are on, where in your system you're looking
    # Calculate the burnup:
    # Build the reaction rate matrix --> This comes from the burnup solver (compute_reaction_rate_matrix function)
    # Save the results for each time step in a list

#plot flux, and U235 composition.


