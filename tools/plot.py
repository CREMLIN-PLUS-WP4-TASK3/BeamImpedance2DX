#!/usr/bin/env python3

from sys import exit
from enum import Enum, auto
import pylab as plt
import pandas as pd
import numpy as np
from sys import argv


if len(argv) != 2:
    exit("Usage:\t./plot data_file.csv")


class DataSign(Enum):
    POSITIVE = auto()
    NEGATIVE = auto()
    ZERO = auto()
    MIXED = auto()


def get_sign(array):
    if np.all(array > 0):
        return DataSign.POSITIVE
    elif np.all(array < 0):
        return DataSign.NEGATIVE
    elif np.all(array == 0.0):
        return DataSign.ZERO
    else:
        return DataSign.MIXED


data = pd.read_csv(argv[1], delimiter=" ", names=["f", "re", "im"])


for i, col, name in [(1, data.re, r"$\Re(z)$"), (2, data.im, r"$\Im(z)$")]:
    sign = get_sign(col)
    plt.subplot(2, 1, i)
    if sign == DataSign.NEGATIVE:
        col *= -1
    plt.plot(data.f, col)
    if sign != DataSign.MIXED and sign != DataSign.ZERO:
        plt.yscale("log")
    plt.xscale("log")
    plt.xlabel("f")
    plt.ylabel(f"-{name}" if sign == DataSign.NEGATIVE else name)
    plt.grid()

plt.tight_layout()
plt.show()
