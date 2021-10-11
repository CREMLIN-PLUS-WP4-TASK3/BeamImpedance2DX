#!/usr/bin/env python3

import pylab as plt
import pandas as pd
import numpy as np
from sys import argv

data = pd.read_csv(argv[1], delimiter=" ", names=["f", "re", "im"])

plt.subplot(2, 1, 1)
plt.plot(data.f, data.re)
plt.xscale("log")
plt.yscale("log")
plt.grid()
plt.subplot(2, 1, 2)
plt.plot(data.f, -data.im)
plt.xscale("log")
plt.yscale("log")
plt.grid()
plt.show()
