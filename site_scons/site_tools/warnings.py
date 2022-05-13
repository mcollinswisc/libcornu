
gcc_warnings = ['implicit',
                'format',
                'missing-braces',
                'missing-include-dirs',
                'sequence-point',
                'switch',
                'parentheses',
                #'redundant-decls',
                'bad-function-cast',
                'uninitialized',
                'return-type',
                'sign-compare',
                'ignored-qualifiers',
                'error']

def generate(env):
    for warning in gcc_warnings:
        flag = '-W' + warning
        env.AppendUnique(CFLAGS = [flag])

def exists(env):
    # TODO: check that this is gcc > 4.0
    return True
