# 3D slab box for Quincky particles

import hoomd, hoomd.md, hoomd.quincke
import numpy as np


hoomd.context.initialize('--mode=cpu --notice-level=0')


# Setting up simulation box dimensions
Lx = 12
Ly = 12
Lz = 12


# Number of particles
N = 300

# Intializing random positions
position = np.random.random((N,3))*[Lx,Ly,0.1] - [Lx/2,Ly/2,0.1/2]


orientation     = [[1,0,0,0]]*N
types           = ['A']
typeid          = np.ones(N)*0
diameter        = np.ones(N)
mass            = np.ones(N)
body            = np.ones(N)*-1
moment_inertia  = [[1,1,1]]*N

s = hoomd.data.make_snapshot(
    N = N, 
    box = hoomd.data.boxdim(Lx = Lx, Ly = Ly, Lz = Lz, dimensions=3))

if hoomd.comm.get_rank() == 0:
    s.particles.types = types
    s.particles.typeid[:] = typeid
    s.particles.position[:] = position
    s.particles.orientation[:] = orientation
    s.particles.diameter[:] = diameter
    s.particles.mass[:] = mass
    s.particles.body[:] = body


# system = hoomd.init.read_snapshot(s)
system = hoomd.init.read_gsd('init.gsd', frame=-1)

# must define moment of inertia
for i in range(N):
    system.particles[i].moment_inertia = moment_inertia[i]



# create cell list
nl = hoomd.md.nlist.cell()
nl2 = hoomd.md.nlist.cell()

# Use dpd conservative force to resolve overlap of randomly placed particles
# generate restart file using this
# dpdc = hoomd.md.pair.dpd_conservative(r_cut=1.0, nlist=nl)
# dpdc.pair_coeff.set(types, types, A=10.0, r_cut = 1.5)


# set lj force
lj = hoomd.md.pair.lj(r_cut=1.0, nlist=nl)
lj.pair_coeff.set(types, types, 
    epsilon = 0.1, 
    sigma   = 1/2.0**(1.0/6.0), 
    alpha   = 1.0, 
    r_cut   = 1)
lj.set_params(mode = "shift")
# lj.disable()



# Quincke force
# Activates when local field is beyond Ecut
H=2.5
qactive = hoomd.quincke.compute.quincke_force(
    group = hoomd.group.all(),
    nlist = nl2,
    Dpassive = 0.0,
    Dactive = 0.1,
    Ee = 0.5,
    Ecut = 1,
    rcut = 5,
    sigma21 = -0.5,
    H = H,
    epsilon = 1,
    seed = 2
    )


# place a wall if confinement in z direction is required
walls=hoomd.md.wall.group(
    [hoomd.md.wall.plane(
    origin=[0,0,H/2], 
    normal=(0.0, 0.0, -1.0), inside=True),
    hoomd.md.wall.plane(
    origin=[0,0,-H/2], 
    normal=(0.0, 0.0, 1.0), inside=True)])
wall_force_lj=hoomd.md.wall.lj(walls, r_cut=3.0)
wall_force_lj.force_coeff.set(
    types, 
    epsilon=0.1, 
    alpha=1,
    sigma=0.5/2.0**(1.0/6.0),
    r_cut=0.5)



integrator_mode = hoomd.md.integrate.mode_standard(dt=1e-4, aniso=False)
all = hoomd.group.all()

seed = int(np.random.randint(1, 100 + 1))
kT = 0.0 # doesn't matter with quincke
bd = hoomd.md.integrate.brownian(
    group=hoomd.group.all(), 
    kT=kT, seed=seed,
    noiseless_t=False,
    noiseless_r=True)
bd.set_gamma(types, gamma = 1.0)


hoomd.dump.gsd(filename="dump.gsd", overwrite=True, period=1e3, group=hoomd.group.all(), phase=0)

hoomd.run(1e5)




