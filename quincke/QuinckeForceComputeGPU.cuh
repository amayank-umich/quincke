// Copyright (c) 2009-2019 The Regents of the University of Michigan
// This file is part of the HOOMD-blue project, released under the BSD 3-Clause License.

#include "hoomd/HOOMDMath.h"
#include "hoomd/ParticleData.cuh"

#ifndef _QUINCKEFORCECOMPUTE_CUH_
#define _QUINCKEFORCECOMPUTE_CUH_


/*! \file ExampleUpdater.cuh
    \brief Declaration of CUDA kernels for ExampleUpdater
*/

cudaError_t gpu_compute_quincke_force_set_forces(const unsigned int group_size,
                                                unsigned int *d_rtag,
                                                unsigned int *d_groupTags,
                                                Scalar4 *d_force,
                                                Scalar4 *d_torque,
                                                Scalar4 *d_orientation,
                                                Scalar Dpassive,
                                                Scalar Dactive,
                                                Scalar Ee,
                                                Scalar Ecut,
                                                Scalar rcut,
                                                Scalar sigma21,
                                                Scalar H,
                                                Scalar epsilon,
                                                const unsigned int N,
                                                unsigned int block_size);

#endif // _QUINCKEFORCECOMPUTE_CUH_
