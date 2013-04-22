import matplotlib.pyplot as plt
import numpy as np

fig = plt.figure(figsize=(8.5,11))

x = np.arange(100)
fig.add_subplot(311)
plt.plot(x,np.random.random(100))
plt.legend('Random Values')
fig.add_subplot(312)
plt.plot(x,np.random.random(100))
plt.legend('Random Values')
fig.add_subplot(313)
plt.plot(x,np.random.random(100))
plt.legend('Random Values')

plt.savefig('test_figsize.pdf')
