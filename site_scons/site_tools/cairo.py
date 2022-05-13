from __future__ import print_function
from SCons.Script import *

def generate(env):
    try:
        env.ParseConfig('pkg-config --cflags --libs cairo')
    except:
        env['WITH_CAIRO'] = False
        return

    env['WITH_CAIRO'] = True

def exists(env):
    return True
