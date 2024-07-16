import numpy as np
import random
from Commands import dec2int, order_parameter, bernoulli, control
from fractions import Fraction
from MonteCarloCommands import boltzmann_probability, energy
import time
import timeit

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
'''
length = 200
prob = 50
K = -1
Beta = 1
number1 = random.random()
number2 = random.random()
#number3 = random.random()
#number = tests[k]
rep1 = dec2int(number1, length)
rep2 = dec2int(number2, length)
#rep3 = dec2int(number3, length)
#print(k, rep, dec2int(Fraction(2, 3), length), dec2int(Fraction(1, 3), length))
for i in range((length**2)//2):
    if random.random() > (float(prob) / 100):
        rep1 = bernoulli(rep1, length)
    else:
        rep1 = control(rep1, length)
    if random.random() > (float(prob)/100):
        rep2 = bernoulli(rep2, length)
    else:
        rep2 = control(rep2, length)
#if random.random() > (float(prob)/100):
#    rep3 = bernoulli(rep3, length)
#else:
#    rep3 = control(rep3, length)

lattice = [format(rep1, '0'+str(length)+'b'), format(rep2, '0'+str(length)+'b')] # each instance of format takes about .3E-6 s to run
print(timeit.timeit(stmt="format(rep1, '0'+str(length)+'b')", setup =
import random
from Commands import dec2int
number1 = random.random()
length = 100
rep1 = dec2int(number1, length)
                    , number=1))
            #format(rep3, '0'+str(length)+'b')]
begin = time.time_ns()
for x_pos in range(0, 2):
    for y_pos in range(0, length):
        E_i = energy(x_pos, y_pos, lattice, 2, length, K, True)
        old = lattice[x_pos]
        # noinspection PyTypeChecker
        lattice[x_pos] = format(int(lattice[x_pos], base=2) ^ (1 << (length - y_pos - 1)), '0' + str(length) + 'b')
        E_f = energy(x_pos, y_pos, lattice, 2, length, K, True)
        if random.random() <= boltzmann_probability(E_i, E_f, Beta):
            pass
        else:
            lattice[x_pos] = old
end = time.time_ns()
print(end-begin)
num= [rep1]
print(timeit.timeit(stmt="format(rep1, '0'+str(length)+'b')", setup =
import random
from Commands import dec2int
number1 = random.random()
length = 100
rep1 = dec2int(number1, length)
                    ))
#X = np.zeros((4, 1))
for i in range(X.shape[0]):
    for j in range(X.shape[1]):
        X[i, j] = random.random()
a = 10
n = 2
lattice = []
length = 50
testar = np.zeros((4,4))
for i in range(4):
    for j in range(4):
        testar[i,j] = i*j
print(testar[:,0])
#print(bin(10)[10])
test = ['0110','0100','0101']
print(test[1][1]+test[0][1], int(test[1][1]+test[0][1],base=2), order_parameter(int(test[1][1]+test[0][1],base=2),2)
      , order_parameter(int('000',base=2),3))
print(2^4)

for i in range(4):
    lattice.append(format(dec2int(Fraction(1,3),length), '0'+str(length)+'b'))
num = lattice[0]
position_y = 2
print(num)
print(num[slice(1, length, length-2)].count('0'))
'''

K = -1
length = 100
rep1 = dec2int(random.random(), length)
rep2 = dec2int(random.random(), length)
Beta=100
lattice = [format(rep1, '0'+str(length)+'b'), format(rep2, '0'+str(length)+'b')]
print(lattice)
for i in range(1):
    for x_pos in range(0, 2):
        for y_pos in range(0, length):
            E_i = energy(x_pos, y_pos, lattice, 2, length, K, True)
            old = lattice[x_pos]
            # noinspection PyTypeChecker
            lattice[x_pos] = format(int(lattice[x_pos], base=2) ^ (1 << (length - y_pos - 1)), '0' + str(length) + 'b')
            E_f = energy(x_pos, y_pos, lattice, 2, length, K, True)
            if random.random() <= boltzmann_probability(E_i, E_f, Beta):
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
