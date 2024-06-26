import numpy as np
import random
import matplotlib.pyplot as plt
from fractions import Fraction
from Commands import dec2int, order_parameter, left, right, reset, control, bernoulli

K = 1
length = 100
random.seed(10)
times = 100
low_prob = 0
high_prob = 100 + 1
record1 = np.zeros([high_prob-low_prob, 1], dtype=float)
record2 = record1
for k in range(0, times):
    j = 0
    print(k)
    for prob in range(low_prob, high_prob, 1):
        #print(prob)
        number1 = random.random()
        number2 = random.random()
        #number = tests[k]
        length = 10
        rep1 = dec2int(number1, length)
        rep2 = dec2int(number2, length)
        #print(k, rep, dec2int(Fraction(2, 3), length), dec2int(Fraction(1, 3), length))
        for i in range((length**2)//2):
            temp1 = 0
            temp2 = 0
            if random.random() > (float(prob)/100):
                rep1 = bernoulli(rep1, length)
            else:
                rep1 = control(rep1, length)
            if random.random() > (float(prob)/100):
                rep2 = bernoulli(rep2, length)
            else:
                rep2 = control(rep2, length)
            # apply transfer matrix
            for j in range(length):
                if j == 0:
                    a1 = int(format(rep1, '0'+str(length)+'b')[length-1])
                    a2 =int(format(rep1, '0'+str(length)+'b')[0])
                    a3 =int(format(rep1, '0'+str(length)+'b')[1])
                    b1 = int(format(rep2, '0'+str(length)+'b')[length-1])
                    b2 = int(format(rep2, '0'+str(length)+'b')[0])
                    b3 =int(format(rep1, '0'+str(length)+'b')[1])
                if j == length-1:
                    a1 = int(format(rep1, '0' + str(length) + 'b')[length - 2])
                    a2 = int(format(rep1, '0' + str(length) + 'b')[length -1])
                    a3 = int(format(rep1, '0' + str(length) + 'b')[0])
                    b1 = int(format(rep2, '0' + str(length) + 'b')[length - 2])
                    b2 = int(format(rep2, '0' + str(length) + 'b')[length -1])
                    b3 = int(format(rep1, '0' + str(length) + 'b')[0])
                else:
                    a1 = int(format(rep1, '0' + str(length) + 'b')[j-1])
                    a2 = int(format(rep1, '0' + str(length) + 'b')[j])
                    a3 = int(format(rep1, '0' + str(length) + 'b')[j+1])
                    b1 = int(format(rep2, '0' + str(length) + 'b')[j- 1])
                    b2 = int(format(rep2, '0' + str(length) + 'b')[j])
                    b3 = int(format(rep1, '0' + str(length) + 'b')[j+1])
                if a1 == 0:
                    a1 = -1
                if a2 == 0:
                    a2 = -1
                if a3 == 0:
                    a3 = -1
                if b1 == 0:
                    b1 = -1
                if b2 == 0:
                    b2 = -1
                if b3 == 0:
                    b3 = -1

                transfer1 = np.tanh(K*(a1*a2 + a2*a3 +a2*b2))
                transfer2 = np.tanh(K*(b1*b2 + b2 * b3 + a2 * b2))
                b2 = transfer2 * b2
                a2 = transfer1 * a2
                if b2 < 0:
                    b2 = 0
                else:
                    b2 = 1
                if a2<0:
                    a2 = 0
                else:
                    a2 = 1
                temp1 = temp1 + (a2<<j)
                temp2 = temp2 + (b2<<j)
            rep1 = temp1
            rep2 = temp2
            #rep = control(rep, length)
            #print(rep)
        record1[j] += order_parameter(rep1, length) / times
        record2[j] += order_parameter(rep2, length)/times
        j += 1
    print(bin(rep1), bin(rep2))

        #print(j, number, rep, record[j])
plt.rcParams.update({
    "text.usetex": False,
})
print(record1[-1])
fig, ax = plt.subplots()
ax.plot( np.linspace(0,1,high_prob-low_prob), record1-record2, marker='o')
ax.set_xlabel(r'Probability of Selecting the Control Map')
ax.set_ylabel(r'$ <\hat{O}>$')
ax.set_title(r'Results of two chains subtracted  '+str(length))
#fig.save("200 Length BitString, Full Probability Spectrum")
plt.show()
fig, ax = plt.subplots()
ax.plot( np.linspace(0,1,high_prob-low_prob), record1, marker='o')
ax.set_xlabel(r'Probability of Selecting the Control Map')
ax.set_ylabel(r'$ <\hat{O}>$')
ax.set_title(r'Results of a Simulation with Bitstring length '+str(length))
#fig.save("200 Length BitString, Full Probability Spectrum")
plt.show()
fig, ax = plt.subplots()
ax.plot( np.linspace(0,1,high_prob-low_prob), record2, marker='o')
ax.set_xlabel(r'Probability of Selecting the Control Map')
ax.set_ylabel(r'$ <\hat{O}>$')
ax.set_title(r'Results of a Simulation with Bitstring length '+str(length))
#fig.save("200 Length BitString, Full Probability Spectrum")
plt.show()
