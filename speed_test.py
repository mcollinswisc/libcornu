from __future__ import print_function
import re
import subprocess
import numpy

re_perspiral = '([0-9]*(\.?[0-9]*))([num]?s) per spiral'
re_perspiral = re.compile(re_perspiral)

units_map = {'ns': 1e-9, 'us': 1e-6, 'ms': 1e-3, 's': 1}

ntrials = 100
kimia = []
walton = []

print('Running %d trials of Cornu spiral speed test' % ntrials)

for i in range(ntrials):
    proc = subprocess.Popen('bin/cornu-speed', stdout = subprocess.PIPE)

    nextline = 'walton'

    for line in proc.stdout:
        print(line[:-1])
        m = re_perspiral.search(line)
        if m == None:
            continue

        time = float(m.groups()[0])
        units = m.groups()[2]

        time *= units_map[units]

        if nextline == 'walton':
            walton.append(time)
            nextline = 'kimia'
        elif nextline == 'kimia':
            kimia.append(time)
            nextline = None
        else:
            print('Huh?')
            exit(1)

walton_mu = numpy.mean(walton)
walton_sigma = numpy.std(walton)

kimia_mu = numpy.mean(kimia)
kimia_sigma = numpy.std(kimia)

print(u'Walton: %e \xB1 %e' % (walton_mu, walton_sigma))
print(u'Kimia: %e \xB1 %e' % (kimia_mu, kimia_sigma))



