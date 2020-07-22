# Copyright (c) 2009-2019 The Regents of the University of Michigan
# This file is part of the HOOMD-blue project, released under the BSD 3-Clause License.

# this simple python interface just activates the c++ ExampleUpdater from cppmodule
# Check out any of the python code in lib/hoomd-python-module/hoomd_script for moreexamples

# First, we need to import the C++ module. It has the same name as this module (example_plugin) but with an underscore
# in front


# Next, since we are extending an updater, we need to bring in the base class updater and some other parts from
# hoomd_script
import hoomd
from hoomd.md.pair import pair
from hoomd.quincke import _quincke

## Zeroes all particle velocities
#
# Every \a period time steps, particle velocities are modified so that they are all zero
#
class quincke_force(pair):
    ## Implement quincke forces
    #
    #
    #
    # \b Examples:
    # \code
    # import hoomd, hoomd.quincke
    # qactive = hoomd.quincke.compute.quincke_force(
    #    group = hoomd.group.all(), nlist = hoomd.md.nlist.cell(), 
    #    Dpassive = 0.01, 
    #    Dactive = 1,
    #    Ee = 1,
    #    Ecut = 1,
    #    rcut = 5,
    #    sigma21 = -0.5,
    #    H = H,
    #    epsilon = 1,
    #    seed = 2
    #    )
    # \endcode
    #
    # \a period can be a function: see \ref variable_period_docs for details
    def __init__(self, nlist, group, Dpassive, Dactive, Ee, Ecut, rcut, sigma21, H, epsilon, seed, name=None):
        hoomd.util.print_status_line();

        # initialize base class
        # hoomd.md.force._force.__init__(self, name);
        pair.__init__(self, rcut, nlist, name);

        # # setup the neighbor list
        # self.nlist = nlist
        # self.nlist.subscribe(lambda:rcut)
        # self.nlist.update_rcut()

        # initialize the reflected c++ class
        if not hoomd.context.exec_conf.isCUDAEnabled():
            self.cpp_force = _quincke.QuinckeForceCompute(hoomd.context.current.system_definition, nlist.cpp_nlist, group.cpp_group, Dpassive, Dactive, Ee, Ecut, rcut, sigma21, H, epsilon, seed);
        else:
            self.nlist.cpp_nlist.setStorageMode(_md.NeighborList.storageMode.full);
            self.cpp_force = _quincke.QuinckeForceComputeGPU(hoomd.context.current.system_definition, nlist.cpp_nlist, group.cpp_group, Dpassive, Dactive, Ee, Ecut, rcut, sigma21, H, epsilon, seed);

        # store metadata
        self.metadata_fields = ['group', 'Dpassive', 'Dactive', 'Ee', 'Ecut', 'rcut', 'sigma21', 'H', 'epsilon', 'seed']
        
        self.group = group
        self.Dpassive = Dpassive
        self.Dactive = Dactive
        self.Ee = Ee
        self.Ecut = Ecut
        self.rcut = rcut
        self.sigma21 = sigma21
        self.H = H
        self.epsilon = epsilon
        self.seed = seed

        hoomd.context.current.system.addCompute(self.cpp_force, self.force_name);



    # there are no coeffs to update for now
    def update_coeffs(self):
        pass



