import numpy as np
import matplotlib as plt
from Commands import dec2int
from MonteCarloCommands import boltzmann_probability, energy
from fractions import Fraction
import random





length = 4
square = True
chains = 1
# initialize ising lattice
if square:
    lattice = []
    for i in range(length):
        lattice.append(format(dec2int(Fraction(1,3), length), '0'+str(length)+'b'))
        print(lattice[i], format(dec2int(Fraction(1,3), length), '0'+str(length)+'b'))
        # initializes array with binary representations

else:
    # ladders
    lattice = np.zeros(chains)
    for i in range(chains):
        lattice[i] = format(random.getrandbits(length), '0'+str(length)+'b')
# choosing position of flip
if square:
    print(lattice)
    x_pos = random.randrange(0, length)
    y_pos = random.randrange(0, length)
    E_i = energy(x_pos,y_pos,lattice, length, length, 1, False)
    old = lattice[x_pos]
    # noinspection PyTypeChecker
    lattice[x_pos]= format(int(lattice[x_pos], base=2) ^ (1<<(length-y_pos-1)),'0'+str(length)+'b')
    E_f = energy(x_pos,y_pos,lattice,length,length,1,False)
    if random.random() <= boltzmann_probability(E_i,E_f,0):
        pass
    else:
        lattice[x_pos] = old
    print(lattice, E_i, E_f)
else:
    x_pos = random.randrange(0, chains)
    y_pos = random.randrange(0, length)
    E_i = energy(x_pos,y_pos,lattice, 2, length, 1, True)
    lattice[x_pos][y_pos] = lattice[x_pos]

