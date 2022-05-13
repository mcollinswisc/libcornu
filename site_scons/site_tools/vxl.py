from __future__ import print_function
from SCons.Script import *
import os, os.path

def generate(env):
    if not os.path.isdir('/usr/include/vxl/vcl/'):
        env['WITH_VXL'] = False
        return
    env.AppendUnique(CPPPATH = '/usr/include/vxl/vcl/')

    conf = Configure(env)

    env['WITH_VXL'] = True

    needed_headers = ['vcl_iostream.h', 'vcl_fstream.h', 'vcl_complex.h', 'vcl_cmath.h', 'vcl_algorithm.h', 'vcl_vector.h']
    for hdr in needed_headers:
        if not conf.CheckCXXHeader(hdr):
            env['WITH_VXL'] = False
            break

    conf.Finish()


def exists(env):
    return True
