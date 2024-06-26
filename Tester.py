import numpy as np
import random
from Commands import dec2int
from fractions import Fraction
from MonteCarloCommands import boltzmann_probability, energy

def add(X):
    Temp = X
    for k in range(len(X[0,:])):
        for l in range(len(X[:,0])):
            if l == len(X[0,:])-1:
                Temp[l, k] = X[l, k] + X[0, k]
            else:
                Temp[l, k] = X[l, k] + X[l + 1, k]
    return Temp
random.seed(10)
#X = np.zeros((4, 1))
'''for i in range(X.shape[0]):
    for j in range(X.shape[1]):
        X[i, j] = random.random()'''
a = 10
n = 2
lattice = []
length = 50
number1 = random.random()
number2 = random.random()
number3 = random.random()
'''
for i in range(4):
    lattice.append(format(dec2int(Fraction(1,3),length), '0'+str(length)+'b'))
num = lattice[0]
position_y = 2
print(num)
print(num[slice(1, length, length-2)].count('0'))
'''
rep1 = dec2int(number1, length)
rep2 = dec2int(number2, length)
rep3 = dec2int(number3, length)
lattice = [format(rep1, '0'+str(length)+'b'), format(rep2, '0'+str(length)+'b'),
                       format(rep3, '0'+str(length)+'b')]
for i in range(10000):
    for x_pos in range(0, 3):
        for y_pos in range(0, length):
            E_i = energy(x_pos, y_pos, lattice, 3, length, 1, False)
            old = lattice[x_pos]
            # noinspection PyTypeChecker
            lattice[x_pos] = format(int(lattice[x_pos], base=2) ^ (1 << (length - y_pos - 1)), '0' + str(length) + 'b')
            E_f = energy(x_pos, y_pos, lattice, 3, length, 1, False)
            if random.random() <= boltzmann_probability(E_i, E_f, 0):
                pass
            else:
                lattice[x_pos] = old
print(lattice)
'''
print(bin(0))
print( a^(0<<n), bin(a^(0<<n)))
print(f'{a:b}', f'{a:b}'[1])
c = format(a, '0'+str(a)+'b')
print(c, 0&(~0b111), 'c')
print(int(c, base=2))
print(c[slice(2,9,6)], 'test')
print((-1)**~int('0',base=2), 'an')
print(len(c))
n=3
print((a>>n) % 2)
print(np.full((2,2), random.random()))
print(1/(1-np.tanh(.5))*(.5-np.tanh(.5)))
rep1 = format(a,'b')
print(format(a,'b'))
temp1 = str(rep1[0:2]) + str(rep1[4-1])
print(temp1)
print(range)
print(10^(1<<(3-2)), bin(10), bin(10^(1<<(3-2))))
print([format(dec2int(Fraction(1,3), 10), '0'+str(10)+'b')])
'''
