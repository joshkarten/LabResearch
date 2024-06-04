import numpy as np
import random
import matplotlib.pyplot as plt
from fractions import Fraction
"""
num << 1 for leftshift 
num >> 1 for rightshift 
x ^ y for XOR 
Changing the control map to work for a coupled bernoulli map 
Initial step is to reproduce the ising transition in the second part of 
https://academic.oup.com/ptp/article/80/1/7/1873647?login=false

"""

def _init_( num, Length):
    num = num
    Length = Length


def to_bin_array( num, Length):
    array = np.empty(Length, bool)
    for i in range(Length, -1, 1):
        if num % 2 == 1:
            array[i] = 1
            num = num // 2
        else:
            array[i] = 0
            num = num / 2
    return array


def dec2int( num, length):
    # Converts a decimal to an equivalent integer representation to work with binary operations
    assert 0 <= num < 1, length > 0
    return int(num * (2 ** length))


def reset(bin_num):
    """
    :type bin_num: int
    odd numbers end with b_L =1
    num - 1 for set b_L = 0 if odd
    """
    if bin_num % 2 == 1:
        bin_num = bin_num - 1
    return bin_num


def left(num, length):
    """

    :param num: int
    :param length: int
    :return: int
    """
    if num >= 2 << (length - 2):
        num = (num << 1) ^ 1
    else:
        num = num << 1
    return num % (2 ** (length))


def right(num):
    """

    :param num: int
    :return: int
    """
    num = reset(num)
    # if num % 2 ==1 and (num//2) % 2 == 0 and (num//4)%2 == 1:
    #    num = num + 2
    num = num // 2
    return num


def order_parameter(num, length):
    param = float(0)
    binrep ="{0:b}".format(left(num, length)^num)
    count = 2*binrep.count('1')/length-1
    #for i in range(length):
    #    if i == length-1:
    #        param += -((-1) ** (((num >> i) % 2) + 1)) * (-1) ** ((num % 2) + 1)
        #if i == 0:
        #    param += ((-1) ** ((num % 2) ^ 1)) * (-1) ** (((num >> 1) % 2) ^ 1)
    #    else:
    #        param += -((-1) ** (((num >> i) % 2) + 1)) * (-1) ** (((num >> (i + 1)) % 2) + 1)
    return count


def bernoulli(num, length):
    num = left(num, length)
    num = num & (~0b111)
    return num + random.randint(0, 8)


def control(num, length):

    if num < 2 ** (length-1):
        xj = dec2int(Fraction(1, 3), length)
    else:
        xj = dec2int(Fraction(2, 3), length)
    if num == xj:
        return num
    else:
        temp = right(num) + right(xj)
    if temp == num:
        return num + 1
    return temp


def S(chain_prev, chain_current):
    if chain_current < chain_prev:
        return -1
    else:
        return 1


def Delta(chain_prev, chain_current, J):
    Del = np.zeros(chain_prev.shape)
    imax = chain_prev.shape[0]-1
    jmax = chain_prev.shape[1]-1
    for i in range(imax+1):
        if i == 0:
            for j in range(jmax+1):
                if j == 0:
                    Del[i, j] = np.tanh(J / 4 * (S(chain_prev[imax, j], chain_current[imax, j])
                                                 + S(chain_prev[i + 1, j], chain_current[i + 1, j])
                                                 + S(chain_prev[i, jmax], chain_current[i, jmax])
                                                 + S(chain_prev[i, j + 1], chain_current[i, j + 1])))
                elif j == chain_prev.shape[1] - 1:
                    Del[i, j] = np.tanh(J / 4 * (S(chain_prev[imax, j], chain_current[imax, j])
                                                 + S(chain_prev[i + 1, j], chain_current[i + 1, j])
                                                 + S(chain_prev[i, j - 1], chain_current[i, j - 1])
                                                 + S(chain_prev[i, 0], chain_current[i, 0])))
                else:
                    Del[i, j] = np.tanh(J / 4 * (S(chain_prev[imax, j], chain_current[imax, j])
                                                 + S(chain_prev[i + 1, j], chain_current[i + 1, j])
                                                 + S(chain_prev[i, j - 1], chain_current[i, j - 1])
                                                 + S(chain_prev[i, j + 1], chain_current[i, j + 1])))
        elif i == chain_prev.shape[0]-1:
            for j in range(jmax +1):
                if j == 0:
                    Del[i, j] = np.tanh(J / 4 * (S(chain_prev[i - 1, j], chain_current[i - 1, j])
                                                 + S(chain_prev[0, j], chain_current[0, j])
                                                 + S(chain_prev[i, jmax], chain_current[i, jmax])
                                                 + S(chain_prev[i, j + 1], chain_current[i, j + 1])))
                elif j == chain_prev.shape[1]-1:
                    Del[i, j] = np.tanh(J / 4 * (S(chain_prev[i - 1, j], chain_current[i - 1, j])
                                                 + S(chain_prev[0, j], chain_current[0, j])
                                                 + S(chain_prev[i, j - 1], chain_current[i, j - 1])
                                                 + S(chain_prev[i, 0], chain_current[i, 0])))
                else:
                    Del[i, j] = np.tanh(J / 4 * (S(chain_prev[i - 1, j], chain_current[i - 1, j])
                                                 + S(chain_prev[0, j], chain_current[0, j])
                                                 + S(chain_prev[i, j - 1], chain_current[i, j - 1])
                                                 + S(chain_prev[i, j + 1], chain_current[i, j + 1])))
        else:
            for j in range(jmax +1):
                if j == 0:
                    Del[i, j] = np.tanh(J / 4 * (S(chain_prev[i - 1, j], chain_current[i - 1, j])
                                                 + S(chain_prev[i + 1, j], chain_current[i + 1, j])
                                                 + S(chain_prev[i, jmax], chain_current[i, jmax])
                                                 + S(chain_prev[i, j + 1], chain_current[i, j + 1])))
                elif j == chain_prev.shape[1]-1:
                    Del[i, j] = np.tanh(J / 4 * (S(chain_prev[i - 1, j], chain_current[i - 1, j])
                                                 + S(chain_prev[i + 1, j], chain_current[i + 1, j])
                                                 + S(chain_prev[i, j - 1], chain_current[i, j - 1])
                                                 + S(chain_prev[i, 0], chain_current[i, 0])))
                else:
                    Del[i, j] = np.tanh(J / 4 * (S(chain_prev[i - 1, j], chain_current[i - 1, j])
                                                 + S(chain_prev[i + 1, j], chain_current[i + 1, j])
                                                 + S(chain_prev[i, j - 1], chain_current[i, j - 1])
                                                 + S(chain_prev[i, j + 1], chain_current[i, j + 1])))
    return Del


def D_Bernoulli(chain_prev, chain_current, J):
    # takes in an initial array and returns a modified array by applying the bernoulli rule
    chain_next = np.zeros(chain_prev.shape)
    delta = Delta(chain_prev, chain_current, J)
    for i in range(chain_next.shape[0]):
        for j in range(chain_next.shape[1]):
            if chain_prev[i, j] < delta[i, j]:
                chain_next[i, j] = (1/delta[i, j])*chain_prev[i, j]
            else:
                chain_next[i, j] = (1/(1-delta[i, j]))*(chain_prev[i, j]- delta[i, j])
    return chain_next


def D_Order(chain_prev, chain_current):
    param = 0
    size = np.size(chain_prev)
    for i in range(np.shape(chain_prev)[0]):
        for j in range(np.shape(chain_prev)[1]):
            param += (1/size) * S(chain_prev[i, j], chain_current[i, j])
    return param

def to_dec( integer, length):
    ans = 0
    i = length
    while integer != 0:
        if integer % 2 == 1:
            ans += 2**(-i)
        i-=1
        integer = integer//2
    return ans

np.random.seed(10)
length = 10
# Testing the sign function
# print(S(1,1), S(1,0), S(0,0), S(0,1))
# Testing the Delta function
'''
D = Delta(np.array([[1,1],[1,1]]),np.array([[1,1],[1,1]]),1)
print(np.tanh(1),"\n", D)
'''
# Delta function should work, needs more extensive testing

#chain0 = np.random.rand(length)
'''
chain0 = np.array([[.5,.5],[.5,.5]])
chain1 = chain0
chain2 = D_Bernoulli(chain0, chain1, .5)
print(Delta(chain0, chain1, .5))
print(chain0, "\n", chain1, "\n","\n", chain2)
'''
param = np.zeros(201)
n = 0
for Jtemp in range(100, 301):
    J = Jtemp / 100
    print(J)
    for loops in range(100):
        #print(loops)
        chain0 = np.random.rand(length,length)
        #chain0 = np.array([chain0, chain0])
        chain1 = chain0
        for steps in range((length**2)//2):
            chain2 = D_Bernoulli(chain0, chain1, J)
            chain0 = chain1
            chain1=chain2
        param[n] += D_Order(chain0, chain1)/100
    n += 1
fig, ax = plt.subplots()
ax.plot( np.linspace(1.00,3.01, 201), param, marker='o')
ax.set_xlabel(r'Coupling')
ax.set_ylabel(r'$ <\hat{O}>$')
ax.set_title(r'Results of a Simulation with Chain length '+str(length))
#fig.save("200 Length BitString, Full Probability Spectrum")
plt.show()


"""
# Tests the conversion and control
test = dec2int(0.5, 4)
print(test)
print(right(test))
print(right(control(test, 4)))
"""
tests = [0, Fraction(3,4), Fraction(2, 3),  Fraction(1, 3),Fraction(1, 2), Fraction(1, 5)]
'''
# Tests iterations of the bernoulli and control map
for i in range(6):
    num = dec2int(tests[i], 4)
    print("new", num)
    for loops in range(10):
        if True:
            num = bernoulli(num, 4)
            print(num)
        else:
            num = control(num, 4)
    print(num)
'''
'''
# Tests the order parameter
for i in range(6):
    num = dec2int(tests[i],  200)
    print(tests[i], order_parameter(num, 200))
'''
# Tests the full thing with random ints over all probabilities
'''
length = 0
random.seed(10)
times = 100
low_prob = 100
high_prob = 100 + 1
record = np.zeros([high_prob-low_prob, 1], dtype=float)
for k in range(0, times):
    j = 0
    print(k)
    for prob in range(low_prob, high_prob, 1):
        #print(prob)
        number = random.random()
        #number = tests[k]
        length = 200
        rep = dec2int(number, length)
        #print(k, rep, dec2int(Fraction(2, 3), length), dec2int(Fraction(1, 3), length))
        for i in range((length**2)//2):
            if random.random() > (float(prob)/100):
                rep = bernoulli(rep, length)
            else:
                rep = control(rep, length)
            #rep = control(rep, length)
            #print(rep)
        binary = to_bin_array(rep, length)
        record[j] += order_parameter(rep, length)/times
        j += 1
        #print(j, number, rep, record[j])
plt.rcParams.update({
    "text.usetex": False,
})
print(record[0])
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
'''
