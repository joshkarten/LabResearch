import numpy as np
import matplotlib as plt
from Commands import dec2int
from Mon
from fractions import Fraction
import random


def boltzmann_probability(initial, final, beta):
    return np.exp(-beta*(final-initial))


length = 40
square = False
chains = 1
# initialize ising lattice
if square:
    lattice = np.zeros(length, Fraction)
    for i in range(length):
        lattice[i] = format(random.getrandbits(length), '0'+str(length)+'b')
        # initializes array with binary representations

else:
    # ladders
    lattice = np.zeros(chains, Fraction)
    for i in range(chains):
        lattice[i] = format(random.getrandbits(length), '0'+str(length)+'b')
# choosing position of flip
if square:
    x_pos = random.randrange(0, length)
    y_pos = random.randrange(0, length)
    E_i =
    lattice[x_pos] = lattice[x_pos]

else:
    x_pos = random.randrange(0, chains)
    y_pos = random.randrange(0, length)
    lattice[x_pos] = lattice[x_pos]

