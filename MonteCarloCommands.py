import numpy as np
from Commands import dec2int
from line_profiler import LineProfiler
import os
import random

#os.environ["LINE_PROFILER"] = "1"
# int(binary, base=2) for bin to int


def boltzmann_probability(initial, final, beta):
    return np.exp(-beta * (final - initial))/2



def energy(position_x, position_y, array, length_x, length_y, k, two_chains):
    if two_chains: # needs to be fixed for later. bad imp
        return ((-1) **((array[position_x]&(1<<position_y))>>(position_y)) * k *((3-2*( ((array[position_x]&(1 << (position_y-1)%length_y))>>(position_y-1)%length_y) + ((array[position_x]&(1<<(position_y+1)%length_y))>>(position_y+1)%length_y) + 
                                                                       ((array[1-position_x]&(1<<position_y))>>(position_y))))))


    else:
        return ((-1) **(array[position_x]&(1<<position_y))) * k * (4-2*(array[position_x]&(1 << (position_y-1)%length_y) + array[position_x]&(1<<(position_y+1)%length_y) # perhaps change to (0b111^num).count('1')
                                                    + array[1-position_x]&(1<<position_y)+array[1+position_x]&(1<<position_y)))

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
