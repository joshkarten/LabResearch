import numpy as np
import random
import matplotlib.pyplot as plt
from fractions import Fraction
"""
num << 1 for leftshift
num >> 1 for rightshift
x ^ y for XOR
application of the control map fully works along with bernoulli map
have tested with [0, 1/4, 1/3, 1/2, 2/3, 0.9]
and length 4, 5, 20
Issue: Order Parameter fails for large bin lengths
Solution: Change to fraction representations
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
    """

    :type num: Fraction
    :type length: int
    """
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
    for i in range(length):
        if i == length-1:
            param += -((-1) ** (((num >> i) % 2) + 1)) * (-1) ** ((num % 2) + 1)
        #if i == 0:
        #    param += ((-1) ** ((num % 2) ^ 1)) * (-1) ** (((num >> 1) % 2) ^ 1)
        else:
            param += -((-1) ** (((num >> i) % 2) + 1)) * (-1) ** (((num >> (i + 1)) % 2) + 1)
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


def to_dec( integer, length):
    ans = 0
    i = length
    while integer != 0:
        if integer % 2 == 1:
            ans += 2**(-i)
        i-=1
        integer = integer//2
    return ans


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

