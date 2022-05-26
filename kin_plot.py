import numpy as np
import matplotlib.pyplot as plt

x = np.arange(0, 10, 1)
y = np.zeros(10)

for i in range(10):
    if i < 2:
        y[i] = 1.0/3.0
    if i >= 2 and i < 4:
        y[i] = (np.sin(90*(i-3)/180*np.pi) + 2)/3.0
    if i >= 4:
        y[i] = 1.0

plt.plot(x, y)
plt.show()
