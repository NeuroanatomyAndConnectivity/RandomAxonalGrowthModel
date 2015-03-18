import numpy
import matplotlib.pyplot as plt

data = numpy.load("chemical_gradient_field.npy")

sample = data[100, :, :, 3]
fig, ax = plt.subplots()
im = plt.pcolor(sample)
fig.colorbar(im, ax=ax)
plt.show()
