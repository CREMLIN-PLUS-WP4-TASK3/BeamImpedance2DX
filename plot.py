#!/usr/bin/env python3

import pylab as plt
import pandas as pd

data = pd.read_csv("out.csv", delimiter=" ", names=["f", "re", "im"])

plt.subplot(2, 1, 1)
plt.plot(data.f, data.re)
plt.xscale("log")
# plt.yscale("log")
plt.grid()
plt.subplot(2, 1, 2)
plt.plot(data.f, -data.im)
plt.xscale("log")
plt.yscale("log")
plt.grid()
plt.show()
