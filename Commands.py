import numpy as np
import random
import matplotlib.pyplot as plt
from fractions import Fraction


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
    '''for i in range(length):
        if i == length-1:
            param += -((-1) ** (((num >> i) % 2) + 1)) * (-1) ** ((num % 2) + 1)
        #if i == 0:
        #    param += ((-1) ** ((num % 2) ^ 1)) * (-1) ** (((num >> 1) % 2) ^ 1)
        else:
            param += -((-1) ** (((num >> i) % 2) + 1)) * (-1) ** (((num >> (i + 1)) % 2) + 1)
            '''
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