// Copyright (c) 2009-2019 The Regents of the University of Michigan
// This file is part of the HOOMD-blue project, released under the BSD 3-Clause License.

// Include the defined classes that are to be exported to python
#include "QuinckeForceCompute.h"


#include <hoomd/extern/pybind/include/pybind11/pybind11.h>

// specify the python module. Note that the name must explicitly match the PROJECT() name provided in CMakeLists
// (with an underscore in front)
PYBIND11_MODULE(_quincke, m)
    {
    export_QuinckeForceCompute(m);

    #ifdef ENABLE_CUDA
    export_QuinckeForceComputeGPU(m);
    #endif
    }