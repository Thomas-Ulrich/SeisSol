#! /usr/bin/python
##
# @file
# This file is part of SeisSol.
#
# @author Alex Breuer (breuer AT mytum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
# @author Sebastian Rettenberger (sebastian.rettenberger AT tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
# @author Alexander Heinecke (heinecke AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Alexander_Heinecke,_M.Sc.,_M.Sc._with_honors)
#
# @section LICENSE
# Copyright (c) 2012-2015, SeisSol Group
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
# Builds the SeisSol code with several options.
#

# operation system (required for exectuion environment)
import os
import sys

# print the welcome message
print '********************************************'
print '** Welcome to the build script of SeisSol **'
print '********************************************'
print 'Copyright (c) 2012-2015, SeisSol Group'

# Check if we the user wants to show help only
if '-h' in sys.argv or '--help' in sys.argv:
  helpMode = True
else:
  helpMode = False
  
def ConfigurationError(msg):
    """Print the error message and exit. Continue only
    if the user wants to show the help message"""
    
    if not helpMode:
        print msg
        Exit(1) 

#
# set possible variables
#
vars = Variables()

# read parameters from a file if given
vars.AddVariables(
  PathVariable( 'buildVariablesFile', 'location of the python file, which contains the build variables', None, PathVariable.PathIsFile )
)
env = Environment(variables=vars)
if 'buildVariablesFile' in env:
  vars = Variables(env['buildVariablesFile'])

# SeisSol specific variables
vars.AddVariables(
  EnumVariable( 'equations',
                'system of PDEs that will be solved',
                'elastic',
                allowed_values=('elastic', 'viscoelastic')
              ),
  

  EnumVariable( 'order',
                'convergence order of the ADER-DG method',
                'none',
                allowed_values=('none', '2', '3', '4', '5', '6', '7', '8')
              ),

  ( 'numberOfMechanisms', 'Number of anelastic mechanisms (needs to be set if equations=viscoelastic).', '0' ),

  ( 'libxsmmGenerator', 'Path to code generator from libxsmm (needs to be set if equations=viscoelastic).' ),

  ( 'programName', 'name of the executable', 'none' ),

  PathVariable( 'buildDir', 'where to build the code', 'build', PathVariable.PathIsDirCreate ),
  
  EnumVariable( 'compileMode', 'mode of the compilation', 'release',
                allowed_values=('debug', 'release', 'relWithDebInfo')
              ),

  EnumVariable( 'parallelization', 'level of parallelization', 'none',
                allowed_values=('none', 'omp', 'mpi', 'hybrid')
              ),

  BoolVariable( 'generatedKernels', 'use generated kernels', False ),
  
  BoolVariable( 'vecReport', 'print verbose vectorization report when using Intel Compiler suite', False ),
  
  BoolVariable( 'hdf5', 'use hdf5 library for data output', False ),
  
  BoolVariable( 'netcdf', 'use netcdf library for mesh input', False ),
  
  BoolVariable( 'sionlib', 'use sion library for checkpointing', False ),
  
  BoolVariable( 'memkind', 'use memkind library for hbw memory support', False ),
  
  EnumVariable( 'unitTests', 'builds additional unit tests',
                'none',
                allowed_values=('none', 'fast', 'all') ),

  EnumVariable( 'logLevel',
                'logging level. \'debug\' runs assertations and prints all information available, \'info\' prints information at runtime (time step, plot number), \'warning\' prints warnings during runtime, \'error\' is most basic and prints errors only',
                'info',
                allowed_values=('debug', 'info', 'warning', 'error')
              ),
                
  EnumVariable( 'logLevel0',
                'logging level for rank 0. Default is same as logLevel',
                'none',
                allowed_values=('none', 'debug', 'info', 'warning', 'error')
              ),

  EnumVariable( 'numberOfTemporalIntegrationPoints',
                'number of temporal integration points for the dynamic rupture boundary integration.; \'auto\' uses the number of temporal integration points required to reach formal convergence order.',
                'auto',
                allowed_values=('auto', '1', '2', '3', '4', '5', '6')
              ),

  BoolVariable( 'commThread', 'use communication thread for MPI progression (option has no effect when not compiling hybrid target)', False ),

  BoolVariable( 'plasticity', 'enable plasticity (generated kernels only)', False )
)

# external variables
vars.AddVariables(
  PathVariable( 'memkindDir',
                'memkind installation directory',
                None,
                PathVariable.PathAccept ),
                  
  PathVariable( 'netcdfDir',
                'NetCDF installation directory',
                None,
                PathVariable.PathAccept ),
                  
  PathVariable( 'hdf5Dir',
                'HDF5 installation directory',
                None,
                PathVariable.PathAccept ),
                  
  PathVariable( 'zlibDir',
                'zlib installation directory',
                None,
                PathVariable.PathAccept ),
                  
  EnumVariable( 'compiler',
                'Select the compiler (default: intel)',
                'intel',
                allowed_values=('intel', 'gcc')),
                
  BoolVariable( 'useExecutionEnvironment',
                'set variables set in the execution environment',
                True ),

  EnumVariable( 'arch',
                'precision -- s for single- and d for double precision -- and architecture used. Warning: \'noarch\' calls the fall-back code and is outperformed by architecture-specific optimizations (if available) greatly.',
                'dnoarch',
                allowed_values=( 'snoarch', 'dnoarch', 'swsm', 'dwsm', 'ssnb', 'dsnb', 'sknc', 'dknc', 'shsw', 'dhsw', 'sskx', 'dskx', 'sknl', 'dknl' )
              ),

  EnumVariable( 'scalasca', 'instruments code with scalasca. \n \'default\': instruments only outer loops. \n'+\
                                                              ' \'kernels\': additionally instruments inner kernels.\n'+\
                                                              ' \'default_2.x\': outer loops with Scalasca version 2.x\n'+\
                                                              ' \'kernels_2.x\': loops and kernels with Scalasca version 2.x\n',
                'none',
                allowed_values=('none', 'default', 'kernels', 'default_2.x', 'kernels_2.x')
              ),
)

env.Tool('MPITool', vars=vars, toolpath=['build/scons/Tools'])

# set environment
env = Environment(variables=vars)

if env['useExecutionEnvironment']:
    env['ENV'] = os.environ

# generate help text
Help(vars.GenerateHelpText(env))

# handle unknown, maybe misspelled variables
unknownVariables = vars.UnknownVariables()

# remove the buildVariablesFile from the list of unknown variables (used before)
if 'buildVariablesFile' in unknownVariables:
  unknownVariables.pop('buildVariablesFile')

# exit in the case of unknown variables
if unknownVariables:
  ConfigurationError("*** The following build variables are unknown: " + str(unknownVariables.keys()))
  
if env['order'] == 'none':
  ConfigurationError("*** Convergence order not set.")
  
if env['equations'] == 'viscoelastic':
  if env['numberOfMechanisms'] == '0':
    ConfigurationError("*** Number of mechanisms not set.")
  if not os.path.exists(os.path.expanduser(env['libxsmmGenerator'])):
    ConfigurationError("*** Path to libxsmm code generator not set.")
  
# check for architecture
if env['arch'] == 'snoarch' or env['arch'] == 'dnoarch':
  print "*** Warning: Using fallback code for unknown architecture. Performance will suffer greatly if used by mistake and an architecture-specific implementation is available."
  
if not env['generatedKernels'] and ( env['parallelization'] == 'omp' or env['parallelization'] == 'hybrid' ):
  ConfigurationError("*** Classic version does not support hybrid parallelization")

#
# preprocessor, compiler and linker
#

# Basic compiler setting
if env['compiler'] == 'intel':
    env['CC'] = 'icc'
    env['CXX'] = 'icpc'
    env['F90'] = 'ifort'
elif env['compiler'] == 'gcc':
    env['CC'] = 'gcc'
    env['CXX'] = 'g++'
    env['F90'] = 'gfortran'
else:
    assert(false)
    
# Parallel compiler required?
if env['parallelization'] in ['mpi', 'hybrid']:
    env.Tool('MPITool', toolpath=['build/scons/Tools'])
    
    # Do not include C++ MPI Bindings
    env.Append(CPPDEFINES=['OMPI_SKIP_MPICXX'])

# Include preprocessor in all Fortran builds
env['F90COM'] = env['F90PPCOM']

# Use Fortran for linking
env['LINK'] = env['F90']

#
# Scalasca
#

# instrument code with scalasca
if env['scalasca'] in ['default', 'kernels']:
  for mode in ['CC', 'CXX', 'F90', 'LINK']:
    l_scalascaPrelude = 'scalasca -instrument -comp=none -user '

    if env['parallelization'] in ['mpi', 'hybrid']:
      l_scalascaPrelude = l_scalascaPrelude + '-mode=MPI '

    env[mode] = l_scalascaPrelude + env[mode]

if env['scalasca'] in ['default_2.x', 'kernels_2.x']:
  l_scorepArguments = " --noonline-access --nocompiler --user "
  if env['parallelization'] == 'none':
    l_scorepArguments = l_scorepArguments + ' --mpp=none '
  if env['parallelization'] in ['mpi', 'hybrid']:
    l_scorepArguments = l_scorepArguments + ' --mpp=mpi '

  if env['parallelization'] in ['mpi', 'none']:
    l_scorepCxxArguments = l_scorepArguments + ' --thread=none '
  else:
    if env['commThread']:
      # Seems to work with "RF_output_on = 0"
      l_scorepCxxArguments = l_scorepArguments + ' --thread=pthread '
    else:
      l_scorepCxxArguments = l_scorepArguments + ' --thread=omp '

  for mode in ['F90']:
    env[mode] = 'scorep' + l_scorepArguments + ' --thread=none ' + env[mode]
  for mode in ['CC', 'CXX', 'LINK']:
    env[mode] = 'scorep' + l_scorepCxxArguments + env[mode]
    
# kernel instrumentation with scalasca
if env['scalasca'] == 'kernels_2.x':
  env.Append(CPPDEFINES=['INSTRUMENT_KERNELS'])
    
#
# Common settings
#

# enforce restrictive C/C++-Code
env.Append(CFLAGS   = ['-Wall', '-Werror', '-ansi'],
           CXXFLAGS = ['-Wall', '-Werror', '-ansi'])
if env['compiler'] == 'intel':
    env.Append(CXXFLGAS = ['-wd13379'])
elif env['compiler'] == 'gcc':
    # TODO Fix kernel generation
    env.Append(CXXFLAGS = ['-Wno-error=unknown-pragmas'])

# generate vector report (only if requested)
if env['vecReport']:
  env.Append(CXXFLAGS = ['-vec-report3'])
  env.Append(CFLAGS   = ['-vec-report3'])

# run preprocessor before compiling
env.Append(F90FLAGS=['-cpp'])

# Use complete line
if env['compiler'] == 'gcc':
    env.Append(F90FLAGS=['-ffree-line-length-none'])
  
# enforce 8 byte precision for reals (required in SeisSol) and compile time boundary check
if env['compiler'] == 'intel':
    env.Append(F90FLAGS=['-r8', '-WB'])
elif env['compiler'] == 'gcc':
    env.Append(F90FLAGS=['-fdefault-real-8'])
    
# Align structs and arrays
if env['compiler'] == 'intel':
    # TODO Check if Fortran alignment is still necessary in the latest version 
    env.Append(F90LFAGS=['-align', '-align', 'array64byte'])

# Add  Linker-flags  for cross-compiling
if env['compiler'] == 'intel':
    env.Append(LINKFLAGS=['-nofor-main', '-cxxlib']) #Add -ldmalloc for ddt
elif env['compiler'] == 'gcc':
    env.Append(LIBS=['stdc++'])
    
#
# Architecture dependent settings
#

# set vector instruction set
if env['arch'] in ['snoarch', 'dnoarch']:
  env['alignment'] = 16
elif env['arch'] in ['ssnb', 'dsnb']:
  env['alignment'] = 32
  env.Append( CFLAGS    = ['-mavx'],
              CXXFLAGS  = ['-mavx'],
              F90FLAGS  = ['-mavx'],
              LINKFLAGS = ['-mavx'] )
elif env['arch'] in ['swsm', 'dwsm']:
  env['alignment'] = 16
  env.Append( CFLAGS    = ['-msse3'],
              CXXFLAGS  = ['-msse3'],
              F90FLAGS  = ['-msse3'],
              LINKFLAGS = ['-msse3']  )
elif env['arch'] in ['shsw', 'dhsw']:
  env['alignment'] = 32
  if env['compiler'] == 'intel':
    env.Append( CFLAGS    = ['-xCORE-AVX2', '-fma'],
                CXXFLAGS  = ['-xCORE-AVX2', '-fma'],
                F90FLAGS  = ['-xCORE-AVX2', '-fma'],
                LINKFLAGS = ['-xCORE-AVX2', '-fma']  )
  else:
    env.Append( CFLAGS    = ['-mavx2', '-mfma'],
                CXXFLAGS  = ['-mavx2', '-mfma'],
                F90FLAGS  = ['-mavx2', '-mfma'],
                LINKFLAGS = ['-mavx2', '-mfma']  )
elif env['arch'] in ['sskx', 'dskx']:
  env['alignment'] = 64
  if env['compiler'] == 'intel':
    env.Append( CFLAGS    = ['-xCORE-AVX512', '-fma'],
                CXXFLAGS  = ['-xCORE-AVX512', '-fma'],
                F90FLAGS  = ['-xCORE-AVX512', '-fma'],
                LINKFLAGS = ['-xCORE-AVX512', '-fma'] ) 
  else:
    env.Append( CFLAGS    = ['-mavx512f', '-mavx512cd', '-mfma'],
                CXXFLAGS  = ['-mavx512f', '-mavx512cd', '-mfma'],
                F90FLAGS  = ['-mavx512f', '-mavx512cd', '-mfma'],
                LINKFLAGS = ['-mavx512f', '-mavx512cd', '-mfma']  )
elif env['arch'] in ['sknc', 'dknc']:
  env['alignment'] = 64
  env.Append( CFLAGS    = ['-mmic', '-fma'],
              CXXFLAGS  = ['-mmic', '-fma'],
              F90FLAGS  = ['-mmic', '-fma'],
              LINKFLAGS = ['-mmic', '-fma'] )
elif env['arch'] in ['sknl', 'dknl']:
  env['alignment'] = 64
  if env['compiler'] == 'intel':
    env.Append( CFLAGS    = ['-xMIC-AVX512', '-fma', '-DENABLE_MATRIX_PREFETCH', '-DENABLE_STREAM_MATRIX_PREFETCH'],
                CXXFLAGS  = ['-xMIC-AVX512', '-fma', '-DENABLE_MATRIX_PREFETCH', '-DENABLE_STREAM_MATRIX_PREFETCH'],
                F90FLAGS  = ['-xMIC-AVX512', '-fma', '-DENABLE_MATRIX_PREFETCH', '-DENABLE_STREAM_MATRIX_PREFETCH'],
                LINKFLAGS = ['-xMIC-AVX512', '-fma', '-DENABLE_MATRIX_PREFETCH', '-DENABLE_STREAM_MATRIX_PREFETCH'] ) 
  else:
    env.Append( CFLAGS    = ['-mavx512f', '-mavx512cd', '-mavx512pf', '-mavx512er', '-mfma', '-DENABLE_MATRIX_PREFETCH', '-DENABLE_STREAM_MATRIX_PREFETCH'],
                CXXFLAGS  = ['-mavx512f', '-mavx512cd', '-mavx512pf', '-mavx512er', '-mfma', '-DENABLE_MATRIX_PREFETCH', '-DENABLE_STREAM_MATRIX_PREFETCH'],
                F90FLAGS  = ['-mavx512f', '-mavx512cd', '-mavx512pf', '-mavx512er', '-mfma', '-DENABLE_MATRIX_PREFETCH', '-DENABLE_STREAM_MATRIX_PREFETCH'],
                LINKFLAGS = ['-mavx512f', '-mavx512cd', '-mavx512pf', '-mavx512er', '-mfma', '-DENABLE_MATRIX_PREFETCH', '-DENABLE_STREAM_MATRIX_PREFETCH']  )
else:
  #assert(env['compileMode'] == 'debug')
  pass

env.Append(CPPDEFINES=['ALIGNMENT='+str(env['alignment']),
                       str(env['arch']).upper()])

# enable interproc. opts for small cores
if env['arch'] in ['sknc', 'dknc', 'sknl', 'dknl']:
  env.Append( F90FLAGS  = ['-ip', '-ipo'],
              CFLAGS    = ['-ip', '-ipo'],
              CXXFLAGS  = ['-ip', '-ipo'],
              LINKFLAGS = ['-ip', '-ipo'] )
  
#
# Compile mode settings
#

# set (pre-)compiler flags for the compile modes
if env['compileMode'] == 'debug':
  env.Append(F90FLAGS = ['-O0'],
             CLFGAGS  = ['-O0'],
             CXXFLAGS = ['-O0'])
  if env['compiler'] == 'intel':
      env.Append(F90FLAGS = ['-shared-intel', '-check'],
                 CLFGAGS  = ['-shared-intel'],
                 CXXFLAGS = ['-shared-intel'])
  else:
      env.Append(F90FLAGS = ['-fcheck=all'])
if env['compileMode'] in ['debug', 'relWithDebInfo']:
  env.Append(F90FLAGS  = ['-g'],
             CFLAGS    = ['-g'],
             CXXFLAGS  = ['-g'],
             LINKFLAGS = ['-g', '-rdynamic'])
  if env['compiler'] == 'intel':
      env.Append(F90FLAGS  = ['-traceback'],
                 CFLAGS    = ['-traceback'],
                 CXXFLAGS  = ['-traceback'])
  else:
      env.Append(F90FLAGS  = ['-fbacktrace'])

if env['compileMode'] in ['relWithDebInfo', 'release']:
    env.Append(CPPDEFINES = ['NDEBUG'])
    env.Append(F90FLAGS = ['-O2'],
               CFLAGS   = ['-O2'], 
               CXXFLAGS = ['-O2'])
    if env['compiler'] == 'intel':
        env.Append(F90FLAGS = ['-fno-alias'])
    
#
# Basic preprocessor defines
#

# set precompiler mode for the number of quantities and basis functions
env.Append(CPPDEFINES=['CONVERGENCE_ORDER='+env['order']])

numberOfQuantities = { 'elastic' : 9, 'viscoelastic' : 9 + int(env['numberOfMechanisms']) * 6 }
env.Append(CPPDEFINES=['NUMBER_OF_QUANTITIES=' + str(numberOfQuantities[ env['equations'] ]), 'NUMBER_OF_RELAXATION_MECHANISMS=' + str(env['numberOfMechanisms'])])

# set number of temporal integration points for dynamic ruputure boundary conditions
if( env['numberOfTemporalIntegrationPoints'] != 'auto' ):
  env.Append(CPPDEFINES=['NUMBER_OF_TEMPORAL_INTEGRATION_POINTS='+env['numberOfTemporalIntegrationPoints']])

# add parallel flag for mpi
if env['parallelization'] in ['mpi', 'hybrid']:
    # TODO rename PARALLEL to USE_MPI in the code
    env.Append(CPPDEFINES=['PARALLEL', 'USE_MPI'])

# add OpenMP flags
if env['parallelization'] in ['omp', 'hybrid']:
    env.Append(CPPDEFINES=['OMP'])
    env.Append(CFLAGS    = ['-fopenmp'],
               CXXFLAGS  = ['-fopenmp'],
               F90FLAGS  = ['-fopenmp'],
               LINKFLAGS = ['-fopenmp'])

if( env['plasticity'] ):
  env.Append(CPPDEFINES=['USE_PLASTICITY'])

# set pre compiler flags for matrix optimizations
if env['generatedKernels']:
  env.Append(CPPDEFINES=['GENERATEDKERNELS', 'CLUSTERED_LTS'])

# set pre compiler flags and link flags for commuincation thread
if env['commThread']:
  env.Append(CPPDEFINES=['USE_COMM_THREAD'])
  env.Append(LINKFLAGS=['-lpthread'] )
 
# Default log level for rank 0 is same as logLevel
if env['logLevel0'] == 'none':
  env['logLevel0'] = env['logLevel']

# set level of logger for fortran
if env['logLevel'] == 'debug':
  env.Append(CPPDEFINES=['LOGLEVEL=3'])
elif env['logLevel'] == 'info':
  env.Append(CPPDEFINES=['LOGLEVEL=2'])
elif env['logLevel'] == 'warning':
  env.Append(CPPDEFINES=['LOGLEVEL=1'])
elif env['logLevel'] == 'error':
  env.Append(CPPDEFINES=['LOGLEVEL=0'])
else:
  assert(false)
  
# set level of logger for rank 0 and C++
if env['logLevel0'] == 'debug':
  env.Append(CPPDEFINES=['LOGLEVEL0=3', 'LOG_LEVEL=3'])
elif env['logLevel0'] == 'info':
  env.Append(CPPDEFINES=['LOGLEVEL0=2', 'LOG_LEVEL=2'])
elif env['logLevel0'] == 'warning':
  env.Append(CPPDEFINES=['LOGLEVEL0=1', 'LOG_LEVEL=1'])
elif env['logLevel0'] == 'error':
  env.Append(CPPDEFINES=['LOGLEVEL0=0', 'LOG_LEVEL=0'])
else:
  assert(false)

# add include path for submodules
env.Append( CPPPATH=['#/submodules'] )

#
# add libraries
#

# Library pathes
env.Tool('DirTool', fortran=True, toolpath=['build/scons/Tools'])

# HDF5
if env['hdf5']:
    env.Tool('Hdf5Tool', required=(not helpMode), toolpath=['build/scons/Tools'])
    env.Append(CPPDEFINES=['USE_HDF'])

# memkind
if env['memkind']:
  env.Tool('MemkindTool', toolpath=['build/scons/Tools'])
  env.Append(CPPDEFINES=['USE_MEMKIND'])

# netCDF
if env['netcdf']:
    env.Tool('NetcdfTool', parallel=(env['parallelization'] in ['hybrid', 'mpi']), required=(not helpMode), toolpath=['build/scons/Tools'])
    env.Append(CPPDEFINES=['USE_NETCDF'])
    
# sionlib still need to create a Tool for autoconfiguration
if env['sionlib'] and not ( env['parallelization'] in ['mpi', 'hybrid']):
  env['sionlib'] = False
  print "conflict: non mpi or hybrid build but sionlib requires mpi:"
  print "          ... deactivating sionlib: env['sionlib'] =",env['sionlib']
  print "          during runtime: sionlib will use posix as fallback solution"

if env['sionlib']:
  env.Append(F90FLAGS=['-DUSE_SIONLIB'])
  env.Append(CPPDEFINES=['USE_SIONLIB'])
  env.Append(CXXFLAGS=['-I/lrz/sys/libraries/sionlib/1.5p1/intel/ibmmpi/include'])
  env.Append(LIBPATH=['/lrz/sys/libraries/sionlib/1.5p1/intel/ibmmpi/lib'])
  env.Append(LINKFLAGS=['-lsionmpi_cxx_64'])
    #    env.Append(LINKFLAGS=['-lsionser_cxx_64'])
  env.Append(LIBS=['sionmpi_64'])
  env.Append(LIBS=['sionser_64'])
  env.Append(LIBS=['sioncom_64'])
  env.Append(LIBS=['sioncom_64_lock_none'])
  env.Append(LIBS=['rt'])

# add pathname to the list of directories wich are search for include
env.Append(F90FLAGS=['-Isrc'])
env.Append(CPPPATH=['#/src', '#/src/Equations/' + env['equations'], '#/src/Equations/' + env['equations'] + '/generated_code'])

#
# setup the program name and the build directory
#
if env['programName'] == 'none':
  program_suffix = '%s_%s_%s_%s_%s_%s_%s' %(
    env['compileMode'],
    'generatedKernels' if env['generatedKernels'] else 'classic',
    env['arch'],
    env['parallelization'],
    'scalasca' if env['scalasca'] != 'none' else 'none',
    numberOfQuantities[ env['equations'] ],
    env['order']
  )
  env['programFile'] = '%s/SeisSol_%s' %(
    env['buildDir'],
    program_suffix
  )
else:
  program_suffix = env['programName']
  env['programFile'] = env['buildDir']+'/'+env['programName']

# build directory

env['buildDir'] = '%s/build_%s' %(env['buildDir'], program_suffix)

# set sub directories (important for scons tree)
buildDirectories = ['Checkpoint', 'Monitoring', 'Reader', 'Physics', 'Geometry', 'Numerical_aux', 'Initializer', 'Solver', 'ResultWriter']

for buildDir in range(len(buildDirectories)):
  buildDirectories[buildDir] = '#/'+env['buildDir'] + '/' + buildDirectories[buildDir]
env.AppendUnique(F90PATH=buildDirectories)

# set module path
if env['compiler'] == 'intel':
    env.Append(F90FLAGS='-module ${TARGET.dir}')
elif env['compiler'] == 'gcc':
    env.Append(F90FLAGS='-J ${TARGET.dir}')

# get the source files
env.sourceFiles = []
env.testSourceFiles = []

Export('env')
SConscript('generated_code/SConscript', variant_dir='#/'+env['buildDir'], src_dir='#/', duplicate=0)
SConscript('src/SConscript', variant_dir='#/'+env['buildDir'], src_dir='#/', duplicate=0)
SConscript('submodules/SConscript', variant_dir='#/'+env['buildDir']+'/submodules', duplicate=0)
Import('env')

# remove .mod entries for the linker
sourceFiles = []
for sourceFile in env.sourceFiles:
  sourceFiles.append(sourceFile[0])

#print env.Dump()

# build standard version
env.Program('#/'+env['programFile'], sourceFiles)

# build unit tests
if env['unitTests'] != 'none':
  # Anything done here should only affect tests
  env = env.Clone()
    
  # define location of cxxtest
  env['CXXTEST'] = 'submodules/cxxtest'
  
  # Continue testing if tests fail
  env['CXXTEST_SKIP_ERRORS'] = True
  
  # Parallel tests?
  if env['parallelization'] in ['mpi', 'hybrid']:
      env['CXXTEST_COMMAND'] = 'mpiexec -np 3 %t'
      
  # Fail on error (as we can't see OK messages in the output)
  env.Append(CPPDEFINES=['CXXTEST_HAVE_EH', 'CXXTEST_ABORT_TEST_ON_FAIL'])
  
  # add cxxtest-tool
  env.Tool('cxxtest', toolpath=[env['CXXTEST']+'/build_tools/SCons'])
  
  # Get test source files
  env.sourceFiles = []
  
  Export('env')
  SConscript('src/tests/SConscript', variant_dir='#/'+env['buildDir']+'/tests', src_dir='#/')
  Import('env')

  # Remove main() to avoid double definition
  sourceFiles = filter(lambda sf: os.path.basename(str(sf)) != 'main.o', sourceFiles)
  # Remove .mod files from additional Fortran files
  for sourceFile in env.sourceFiles:
    sourceFiles.append(sourceFile[0])

  if env.testSourceFiles:
    # build unit tests
    env.CxxTest(target='#/'+env['buildDir']+'/tests/cxxtest_runner',
                source=sourceFiles+env.testSourceFiles)