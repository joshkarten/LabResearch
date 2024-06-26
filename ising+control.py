import numpy as np
import random
from fractions import Fraction
from Commands import dec2int, bernoulli, control, order_parameter
from MonteCarloCommands import energy, boltzmann_probability
import matplotlib.pyplot as plt

length = 50
random.seed(10)
times = 100
low_prob = 00
high_prob = 100 + 1
record = np.zeros([high_prob-low_prob, 1], dtype=float)
for k in range(0, times):
    j = 0
    print(k)
    for prob in range(low_prob, high_prob, 1):
        #print(prob)
        number1 = random.random()
        number2 = random.random()
        number3 = random.random()
        #number = tests[k]
        rep1 = dec2int(number1, length)
        rep2 = dec2int(number2, length)
        rep3 = dec2int(number3, length)
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
            if random.random() > (float(prob)/100):
                rep3 = bernoulli(rep3, length)
            else:
                rep3 = control(rep3, length)
            lattice = [format(rep1, '0'+str(length)+'b'), format(rep2, '0'+str(length)+'b'),
                       format(rep3, '0'+str(length)+'b')]
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
        record[j] += (order_parameter(rep1, length)+order_parameter(rep2,length)+order_parameter(rep3, length))/(3*times)
        j += 1
        #print(j, number, rep, record[j])
plt.rcParams.update({
    "text.usetex": False,
})
print(record[99])
fig, ax = plt.subplots()
ax.plot( np.linspace(0,1,high_prob-low_prob), record, marker='o')
ax.set_xlabel(r'Probability of Selecting the Control Map')
ax.set_ylabel(r'$ <\hat{O}>$')
ax.set_title(r'Results of a Simulation with Bitstring length '+str(length))
#fig.save("200 Length BitString, Full Probability Spectrum")
plt.show()
fig, ax = plt.subplots()
ax.plot( np.linspace(.3,.6,30), record[30:60], marker='o')
ax.set_xlabel(r'Probability of Selecting the Control Map')
ax.set_ylabel(r'$ <\hat{O}>$')
ax.set_title(r'Results of a Simulation with Bitstring length '+str(length))
#fig.save("200 Length BitString, Full Probability Spectrum")
plt.show()
