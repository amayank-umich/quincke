// Copyright (c) 2009-2019 The Regents of the University of Michigan
// This file is part of the HOOMD-blue project, released under the BSD 3-Clause License.

#include "QuinckeForceComputeGPU.cuh"
#include "hoomd/RandomNumbers.h"
#include "hoomd/RNGIdentifiers.h"
using namespace hoomd;

#include <assert.h>

/*! \file QuinckeForceCompute.cu
    \brief CUDA kernels for QuinckeForceCompute
*/

// First, the kernel code for zeroing the velocities on the GPU
//! Kernel that zeroes velocities on the GPU
/*! \param d_vel Velocity-mass array from the ParticleData
    \param N Number of particles

    This kernel executes one thread per particle and zeros the velocity of each. It can be run with any 1D block size
    as long as block_size * num_blocks is >= the number of particles.
*/

__global__ void gpu_compute_active_force_set_forces_kernel(const unsigned int group_size,
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
                                                    const unsigned int N)
    {
    unsigned int group_idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (group_idx >= group_size)
        return;

    unsigned int tag = d_groupTags[group_idx];
    unsigned int idx = d_rtag[tag];


    // for now do nothing
    //quat<Scalar> quati(h_orientation.data[idx]);
    //Scalar3 f = make_scalar3(0,0,m_params);
    //vec3<Scalar> fi = rotate(quati, vec3<Scalar>(f));
    //d_force[idx].x = fi.x;
    //d_force[idx].y = fi.y;
    //d_force[idx].z = fi.z;  

    
    }



/*! \param d_vel Velocity-mass array from the ParticleData
    \param N Number of particles
    This is just a driver for gpu_zero_velocities_kernel(), see it for the details
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
                                                    unsigned int block_size)
    {
    // setup the grid to run the kernel
    dim3 grid( group_size / block_size + 1, 1, 1);
    dim3 threads(block_size, 1, 1);

    // run the kernel
    cudaMemset(d_force, 0, sizeof(Scalar4)*N);
    gpu_compute_quincke_force_set_forces_kernel<<< grid, threads>>>( group_size,
                                                                    d_rtag,
                                                                    d_groupTags,
                                                                    d_force,
                                                                    d_torque,
                                                                    d_orientation,
                                                                    Dpassive,
                                                                    Dactive,
                                                                    Ee,
                                                                    Ecut,
                                                                    rcut,
                                                                    sigma21,
                                                                    H,
                                                                    epsilon,
                                                                    N);
    return cudaSuccess;
    }