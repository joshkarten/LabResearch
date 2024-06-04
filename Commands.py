import random

import numpy as np
import scipy


def to_bin_array(num, Length):
    array = np.empty(Length, bool)
    assert 0 <= num < 1
    for i in range(Length):
        if num * 2 >= 1:
            array[i] = 1
            num = num * 2 - 1
        else:
            array[i] = 0
            num = num * 2
    return array


def left_shift(binArray):
    store = binArray[0]
    for i in range(len(binArray) - 1):
        binArray[i] = binArray[i + 1]
    binArray[len(binArray) - 1] = store
    return binArray


def right_shift(binArray):
    for i in range(len(binArray) - 1, 0, -1):
        binArray[i] = binArray[i - 1]
    binArray[0] = 0
    return binArray


def apply_control(binArray):
    if binArray[0] == 1:
        return to_bin_array(2 / 3, len(binArray))
    else:
        return to_bin_array(1 / 3, len(binArray))


def adder(bin1, bin2):
    assert len(bin1) == len(bin2)
    print(bin1)
    store = False
    bin3=np.empty(len(bin1),bool)
    for i in range(len(bin1) - 1, -1, -1):
        #print("bin3", bin3, "\nstore",store, "\n", (bin2[i] and store) or (store and bin1[i]) or (bin1[i] and bin2[i]))
        bin3[i] = bin1[i] ^ bin2[i] ^ store
        #print("bin3", bin3, "\nbin1", bin1, "\nbin2", bin2)
        if (bin2[i] and store) or (store and bin1[i]) or (bin1[i] and bin2[i]): store = True
        else: store = False


    return bin3


def bernoulli(bin):
    bin = left_shift(bin)
    # Then set last 3 digits of bin to 0
    # After that, generate a random number between 0 and 8

    return bin


def control(bin):
    if bin[0]:
        return adder(right_shift(bin), right_shift(to_bin_array(2 / 3)))
    else:
        return adder(right_shift(bin), right_shift(to_bin_array(1 / 3)))


binArray = to_bin_array(2 / 3, 5)
print("binarray", binArray)
left = left_shift(binArray)
print("leftshift", left)
right = right_shift(binArray)
print("rightshift", right)
add = adder(binArray, left)
print(add)


'''def dec2bin(xj, L, base=2):
    return int(xj*(1<<L))


def _initialize_binary(L, xj):
        binary_xj = {xj / 2: dec2bin(xj / 2, L) for xj in xj}
        return binary_xj


x = np.binary_repr(10, 5)
print(x)
vec=dec2bin(5,5)
vv = _initialize_binary(5,5)
print(vec)'''
