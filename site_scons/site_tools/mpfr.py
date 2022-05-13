from __future__ import print_function
from SCons.Script import *

def generate(env):
    conf = Configure(env)

    if not conf.CheckLibWithHeader('mpfr', 'mpfr.h', 'c'):
        print('GNU MPFR required')
        Exit(1)

    conf.Finish()

def exists(env):
    return True
