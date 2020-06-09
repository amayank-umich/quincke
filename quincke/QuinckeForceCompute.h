// Copyright (c) 2009-2019 The Regents of the University of Michigan
// This file is part of the HOOMD-blue project, released under the BSD 3-Clause License.

// **********************
// This is a simple example code written for no function purpose other then to demonstrate the steps needed to write a
// c++ source code plugin for HOOMD-Blue. This example includes an example Updater class, but it can just as easily be
// replaced with a ForceCompute, Integrator, or any other C++ code at all.

// inclusion guard
#ifndef _QUINCKEFORCECOMPUTE_H_
#define _QUINCKEFORCECOMPUTE_H_

/*! \file QuinckeForceCompute.h
    \brief Declaration of class for computing Quincke Force
*/

#include "hoomd/ForceCompute.h"
#include "hoomd/ParticleGroup.h"
#include <memory>
#include <hoomd/extern/pybind/include/pybind11/pybind11.h>
#include "hoomd/extern/pybind/include/pybind11/numpy.h"
#include "hoomd/HOOMDMath.h"
#include "hoomd/VectorMath.h"
#include "hoomd/GlobalArray.h"
#include "hoomd/ForceCompute.h"
#include "hoomd/md/NeighborList.h"
#include "hoomd/RandomNumbers.h"
#include "hoomd/RNGIdentifiers.h"

// pybind11 is used to create the python bindings to the C++ object,
// but not if we are compiling GPU kernels
#ifdef NVCC
#error This header cannot be compiled by nvcc
#endif

// (if you really don't want to include the whole hoomd.h, you can include individual files IF AND ONLY IF
// hoomd_config.h is included first)
// For example:
// #include <hoomd/Updater.h>

// Second, we need to declare the class. One could just as easily use any class in HOOMD as a template here, there are
// no restrictions on what a template can do

//! A nonsense particle updater written to demonstrate how to write a plugin
/*! This updater simply sets all of the particle's velocities to 0 when update() is called.
*/
class PYBIND11_EXPORT QuinckeForceCompute : public ForceCompute
    {
    public:
        //! Constructor
        QuinckeForceCompute(std::shared_ptr<SystemDefinition> sysdef,
                            std::shared_ptr<NeighborList> nlist,
                            std::shared_ptr<ParticleGroup> group,
                            Scalar Dpassive,
                            Scalar Dactive,
                            Scalar Ee,
                            Scalar Ecut,
                            Scalar rcut,
                            Scalar sigma21,
                            Scalar H,
                            Scalar epsilon,
                            int seed);

        //! Destructor
        ~QuinckeForceCompute();

        protected:
        //! Actually compute the forces
        virtual void computeForces(unsigned int timestep);

        //! Set forces for particles
        virtual void setForces(unsigned int timestep);

        std::shared_ptr<NeighborList> m_nlist;
        std::shared_ptr<ParticleGroup> m_group;   //!< Group of particles on which this force is applied
        Scalar m_Dpassive;
        Scalar m_Dactive;
        Scalar m_Ee;
        Scalar m_Ecut;
        Scalar m_rcut;
        Scalar m_sigma21;
        Scalar m_H;
        Scalar m_epsilon;
        int m_seed;

        unsigned int last_computed;
    };

//! Export the QuinckeForceCompute class to python
void export_QuinckeForceCompute(pybind11::module& m);

// Third, this class offers a GPU accelerated method in order to demonstrate how to include CUDA code in pluins
// we need to declare a separate class for that (but only if ENABLE_CUDA is set)

// #ifdef ENABLE_CUDA

// //! A GPU accelerated nonsense particle updater written to demonstrate how to write a plugin w/ CUDA code
// /*! This updater simply sets all of the particle's velocities to 0 (on the GPU) when update() is called.
// */
// class ExampleUpdaterGPU : public ExampleUpdater
//     {
//     public:
//         //! Constructor
//         ExampleUpdaterGPU(std::shared_ptr<SystemDefinition> sysdef);

//         //! Take one timestep forward
//         virtual void update(unsigned int timestep);
//     };

// //! Export the ExampleUpdaterGPU class to python
// void export_ExampleUpdaterGPU(pybind11::module& m);

// #endif // ENABLE_CUDA

#endif // _QUINCKEFORCECOMPUTE_H_
