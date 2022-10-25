import numpy as np
import matplotlib.pyplot as plt

arr = np.random.random((100, 100))

# plot CDF of arr
count, bins_count = np.histogram(arr.flatten(), bins=100)

# finding the PDF of the histogram using count values
pdf = count / sum(count)

# using numpy np.cumsum to calculate the CDF
# We can also find using the PDF values by looping and adding
cdf = np.cumsum(pdf)

# plotting PDF and CDF
plt.plot(bins_count[1:], pdf, label="PDF")
plt.plot(bins_count[1:], cdf, label="CDF")
plt.legend()
plt.show()