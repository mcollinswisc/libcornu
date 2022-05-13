from __future__ import print_function
import os, os.path
import multiprocessing

# Set up environment
env = Environment(ENV = os.environ)

env.Tool('warnings')

env.AppendUnique(LIBS = ['m'])

#env.AppendUnique(CFLAGS = ['-g'], MEXFLAGS = ['-g'], CXXFLAGS = ['-g'], LINKFLAGS = ['-g'])
env.AppendUnique(CFLAGS = ['-O3'], LINKFLAGS = ['-O3'], CXXFLAGS = ['-O3'])
env.AppendUnique(CFLAGS = ['-flto'], LINKFLAGS = ['-flto'], CXXFLAGS = ['-flto'])
#env.Append(CFLAGS = ['-g', '-pg'], CXXFLAGS = ['-g', '-pg'], LINKFLAGS = ['-g', '-pg'])

# Assemble core sources
core_src = ['cornu.c', 'integration.c', 'newton.c',
            'fresnel.c', 'fresnel_exact.c', 'fleckner.c']
core_src = [os.path.join('src', f) for f in core_src]
core_obj = [env.Object(f) for f in core_src]

util_obj = env.Object('src/timing.c') + env.Object('src/progress.c')

# Generate the precomputed Fresnel integral table
table_env = env.Clone()

table_env.Tool('gmp')
table_env.Tool('mpfr')
table_env.AppendUnique(LIBS = ['pthread'])

table_size = Value(1024)
nthreads = Value(multiprocessing.cpu_count())

table_prog = table_env.Program('bin/fresnel-table', ['src/fresnel_table.c'] + core_obj[:-1] + util_obj)

env.Command('src/fresnel_table.h', table_prog + [table_size, nthreads],
            '$SOURCES $TARGET')

# Compile library
static_lib = env.StaticLibrary('lib/cornu', core_src)
shared_lib = env.SharedLibrary('lib/cornu', core_src)

# Add include files
env.Install('include', 'src/cornu.h')

# Compile tests
env.Tool('cairo')
env.Tool('vxl')

common_test_code = util_obj
common_test_code += static_lib

if env['WITH_CAIRO']:
    env.Program('bin/test-cornu', ['src/test_cornu.c'] + common_test_code)

if env['WITH_CAIRO'] and env['WITH_VXL']:
    kimia = SConscript('external/kimia/SConscript', exports = 'env')
    env.Program('bin/cornu-speed', ['src/cornu_speed.cc'] + kimia + common_test_code)
