microscopic_cross_sections = {
    'U235': {'xs_fission': 585.1,
             'xs_absorption': 673.9,
             'xs_scattering': 15.12},
    'U238': {'xs_fission': 0.,
             'xs_absorption': 2.6, 
             'xs_scattering': 9.3},
    'Pu239': {'xs_fission': 747.4,
              'xs_absorption': 1020.0,
              'xs_scattering': 8.813},
    'Xe135': {'xs_fission': 0.,
              'xs_absorption': 24,
              'xs_scattering': 4.3,
              'lambda': 3e-5},
    'Lumped': {'xs_fission': 0.,
              'xs_absorption': 0.,
              'xs_scattering': 0.,
              'lambda': (1-3e-5)},
    'H2O': {'xs_fission': 0.,
            'xs_absorption': 0.664,
            'xs_scattering': 103.0}
}

fission_yield = {'Xe135': 0.06}