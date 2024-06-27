import numpy as np
from Commands import dec2int
# int(binary, base=2) for bin to int


def boltzmann_probability(initial, final, beta):
    return np.exp(-beta * (final - initial))


def energy(position_x, position_y, array, length_x, length_y, k, two_chains):
    num = array[position_x]
    # likely can save a lot of time here, but time saves are likely correlated to the data stored in lattice
    if two_chains:
        if position_y == length_y - 1:
            return (-1 ** (~int(num[position_y]))) * k * (3 - (num[slice(0, position_y, position_y - 1)].count('0')
                                                               + array[1 - position_x][position_y].count('0')))
        if position_y == 0:
            return (-1 ** (~int(num[position_y]))) * k * (3 - (num[slice(1, length_y, length_y - 2)].count('0')
                                                               + array[1 - position_x][position_y].count('0')))
        else:
            return (-1 ** (~int(num[position_y]))) * k * (3 - (num[slice(position_y - 1, position_y + 2, 2)].count('0')
                                                               + array[1 - position_x][position_y].count('0')))

    elif position_x == 0:
        if position_y == length_y - 1:
            return (-1 ** ~int(num[position_y])) * k * (4 - 2 * (num[slice(0, position_y, length_y - 2)].count('0')
                                                                 + array[position_x - 1][position_y].count('0')
                                                                 + array[length_x - 1][position_y].count('0')))
        if position_y == 0:
            return (-1 ** ~int(num[position_y])) * k * (4 - 2 * (num[slice(1, length_y, length_y - 2)].count('0')
                                                                 + array[position_x + 1][position_y].count('0')
                                                                 + array[length_x - 1][position_y].count('0')))
        else:
            return (-1 ** ~int(num[position_y])) * k * (
                        4 - 2 * (num[slice(position_y - 1, position_y + 2, 2)].count('0')
                                 + array[position_x + 1][position_y].count('0')
                                 + array[length_x - 1][position_y].count('0')))
    else:
        if position_y == length_y - 1:
            return (-1 ** ~int(num[position_y])) * k * (4 - 2 * (num[slice(0, length_y, length_y - 2)].count('0')
                                                                 + array[position_x + 1][position_y].count('0')
                                                                 + array[position_x - 1][position_y].count('0')))
        if position_y == 0:
            return (-1 ** ~int(num[position_y])) * k * (4 - 2 * (num[slice(1, length_y, length_y - 2)].count('0')
                                                                 + array[position_x + 1][position_y].count('0')
                                                                 + array[position_x - 1][position_y].count('0')))
        else:
            return (-1 ** ~int(num[position_y])) * k * (
                        4 - 2 * (num[slice(position_y - 1, position_y + 2, 2)].count('0')
                                 + array[position_x + 1][position_y].count('0')
                                 + array[position_x - 1][position_y].count('0')))
