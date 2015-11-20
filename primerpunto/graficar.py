import numpy as np
import matplotlib.pyplot as plt
import matplotlib

archivo = np.loadtxt("data.dat")

a = archivo[:,0] 
d = archivo[:,1]
b = archivo[:,2]
c = archivo[:,3]
l = archivo[:,4]

fig = plt.figure(figsize=(15,15))

ax1 = plt.subplot(5,5,1)
plt.title("a vs d")
ax1.scatter(a,d)

ax1 = plt.subplot(5,5,6)
plt.title("a vs d")
plt.scatter(a,b)

ax2 = plt.subplot(5,5,7)
plt.title("a vs c")
plt.scatter(a,c)

ax3 = plt.subplot(5,5,11)
plt.title("a vs l")
plt.scatter(a,l)

ax4 = plt.subplot(5,5,12)
plt.title("d vs b")
plt.scatter(d,b)

ax5 = plt.subplot(5,5,13)
plt.title("d vs c")
plt.scatter(d,c)

ax6 = plt.subplot(5,5,16)
plt.title("d vs l")
plt.scatter(d,l)

ax7 = plt.subplot(5,5,17)
plt.title("b vs c")
plt.scatter(b,c)

ax8 = plt.subplot(5,5,18)
plt.title("b vs l")
plt.scatter(b,l)

plt.show()
