#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt

mu, sigma = 0, 1 # mean and standard deviation
f = np.random.normal(mu, sigma, 1000) # generate feature-vector with normal distribution

# plot the histogram - check the distribution
count, bins, ignored = plt.hist(f, 30, normed=True)

plt.plot(bins, 1/(sigma * np.sqrt(2 * np.pi)) *
                np.exp( - (bins - mu)**2 / (2 * sigma**2) ),
        linewidth=2, color='r')
plt.xlabel('Values')
plt.ylabel('Probability')
plt.title('Histogram')
plt.text(60, .025, r'$\mu=0,\ \sigma=1$')
plt.axis([-0.4, 0.3, 0, 5])
plt.grid(True)
plt.show()
