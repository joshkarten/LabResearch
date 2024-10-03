import numpy as np
from Commands import dec2int
#from line_profiler import LineProfiler
import os
import random
#os.environ["LINE_PROFILER"] = "1"
# int(binary, base=2) for bin to int


def boltzmann_probability(initial, final, beta):
    return np.exp(-beta * (final - initial))/2

def glauber_energy(initial, final, beta):
    return (np.exp(-beta*((final-initial)))/(1+np.exp(-beta*((final-initial)))))

def energy(position_x, position_y, array, length_x, length_y, k, two_chains):
    if two_chains: # needs to be fixed for later. bad imp
        return ((-1) **((array[position_x]&(1<<position_y))>>(position_y)) * k *((3-2*( ((array[position_x]&(1 << (position_y-1)%length_y))>>(position_y-1)%length_y) + ((array[position_x]&(1<<(position_y+1)%length_y))>>(position_y+1)%length_y) + 
                                                                       ((array[1-position_x]&(1<<position_y))>>(position_y))))))


    else:
        return ((-1) **(array[position_x]&(1<<position_y))) * k * (4-2*(array[position_x]&(1 << (position_y-1)%length_y) + array[position_x]&(1<<(position_y+1)%length_y) # perhaps change to (0b111^num).count('1')
                                                    + array[1-position_x]&(1<<position_y)+array[1+position_x]&(1<<position_y)))

def OneDEnergy(position_y, chain,length_y, k):
    return ((-1) **((chain&(1<<position_y))>>(position_y)) * k *
            ((2-2*( ((chain&(1 << (position_y-1)%length_y))>>(position_y-1)%length_y) 
                   + ((chain&(1<<(position_y+1)%length_y))>>(position_y+1)%length_y) ))))
# chain&(1<<position_y) removes less significant bits, >>position_y removes the more significant bits, in effect this selects a bit at a position from an integer representation
# the %length_y is done to account for position_y=0
# multiplying by k accounts for ferromagnetism or antiferromagnetism being favorable (k = -1 or 1)
'''
length = 10
K = -1
number1 = random.random()
number2 = random.random()
#number3 = random.random()
#number = tests[k]
rep1 = dec2int(number1, length)
rep2 = dec2int(number2, length)
lattice = [format(rep1, '0'+str(length)+'b'), format(rep2, '0'+str(length)+'b')]
lp = LineProfiler()
lp_wrapper = lp(energy)
lp_wrapper(1, 0, lattice, 3, length, K, True)
lp.print_stats()
'''
