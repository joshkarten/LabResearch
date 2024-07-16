import numpy as np
from Commands import dec2int
from line_profiler import LineProfiler
import os
import random

#os.environ["LINE_PROFILER"] = "1"
# int(binary, base=2) for bin to int


def boltzmann_probability(initial, final, beta):
    return np.exp(-beta * (final - initial))



def energy(position_x, position_y, array, length_x, length_y, k, two_chains):
    num = array[position_x]
    if two_chains: # needs to be fixed for later. bad imp, 2nd if is redundant, all are missing y edge conditions
        if position_y==length_y-1:
            return ((-1) **(~int(num[position_y], base=2))) * k * (3 - 2*(num[slice(0, position_y, position_y-1)].count('0') # perhaps change to (0b111^num).count('1')
                                                    + array[1-position_x][position_y].count('0')))
        if position_y==0:
            return ((-1 )** (~int(num[position_y], base=2))) * k * (3 - 2*(num[slice(1, length_y, length_y - 2)].count('0')
                                                        + array[1 - position_x][position_y].count('0')))
        else:
            return ((-1) ** (~int(num[position_y], base=2))) * k*(3 - 2*(num[slice(position_y-1, position_y+2, 2)].count('0')
                                                    + array[1-position_x][position_y].count('0')))

    elif position_x == 0:
        if position_y == length_y-1:
            return ((-1) ** ~int(num[position_y], base=2)) *k* (4 - 2*(num[slice(0, position_y, length_y - 2)].count('0')
                                                         + array[position_x - 1][position_y].count('0')
                                                         + array[length_x-1][position_y].count('0')))
        if position_y==0:
            return ((-1) ** ~int(num[position_y], base=2)) * k * (4 - 2 * (num[slice(1, length_y, length_y - 2)].count('0')
                                                                 + array[position_x + 1][position_y].count('0')
                                                                 + array[length_x - 1][position_y].count('0')))
        else:
            return ((-1) ** ~int(num[position_y], base=2)) *k* (4 - 2*(num[slice(position_y - 1, position_y + 2, 2)].count('0')
                                                         + array[position_x + 1][position_y].count('0')
                                                         + array[length_x-1][position_y].count('0')))
    else:
        if position_y == length_y-1:
            return ((-1) ** ~int(num[position_y], base=2)) *k* (4 - 2* (num[slice(0, length_y, length_y-2)].count('0')
                                                         + array[position_x + 1][position_y].count('0')
                                                         + array[position_x - 1][position_y].count('0')))
        if position_y == 0:
            return ((-1) ** ~int(num[position_y], base=2)) * k * (4 - 2 * (num[slice(1, length_y, length_y - 2)].count('0')
                                                                 + array[position_x + 1][position_y].count('0')
                                                                 + array[position_x - 1][position_y].count('0')))
        else:
            return ((-1) ** ~int(num[position_y], base=2)) *k* (4 - 2*(num[slice(position_y - 1, position_y + 2, 2)].count('0')
                                                         + array[position_x + 1][position_y].count('0')
                                                         + array[position_x - 1][position_y].count('0')))

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
