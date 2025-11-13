#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

# Simple test plot
x = np.linspace(-5, 5, 500)
y = np.sin(x)

plt.plot(x, y)
plt.title("Test sine wave")
plt.savefig("test_plot.png", dpi=150)
print("Done. Plot saved as test_plot.png")
