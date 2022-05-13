from __future__ import print_function
from SCons.Script import *

def generate(env):
    conf = Configure(env)

    if not conf.CheckLibWithHeader('gmp', 'gmp.h', 'c'):
        print('GNU Multiple Precision Arithmetic Library (GMP) required')
        Exit(1)

    conf.Finish()

def exists(env):
    return True
