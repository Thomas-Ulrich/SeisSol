#! /usr/bin/python
##
# @file
# This file is part of SeisSol.
#
# @author Sebastian Rettenberger (sebastian.rettenberger @ tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
#
# @section LICENSE
# Copyright (c) 2012, SeisSol Group
# All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
# 
# 1. Redistributions of source code must retain the above copyright notice,
#    this list of conditions and the following disclaimer.
# 
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
# 
# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from this
#    software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
#
# @section DESCRIPTION
# Build parameters for SuperMUC Phase 2 and Hornet @ HLRS
#

# build options
compileMode                 = 'release'
#compileMode                 = 'relWithDebInfo'
#compileMode                 = 'debug'
parallelization             = 'hybrid'
#parallelization             = 'mpi'
generatedKernels            = 'yes'
measureNodeLevelPerformance = 'none'
useExecutionEnvironment     = 'yes'
order = 3
equations='elastic'
# machine dependent options
cppCompiler          = 'mpiCC'
fortranCompiler      = 'mpif90'

netcdf='yes'
hdf5='yes'

netcdfDir='/lrz/sys/libraries/netcdf/4.3.3/intel/ibmmpi_poe1.4_1505'
hdf5Dir='/lrz/sys/libraries/hdf5/1.8.14/ibmmpi_poe1.4_15.0.5'


phase=1
if phase==1:
   arch = 'dsnb'
else:
   arch = 'dhsw'
   #commThread ='yes'

plasticity='yes'

#logLevel                    = 'warning'
logLevel                    = 'warning'
logLevel0                   = 'info'
