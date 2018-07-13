import matplotlib.pyplot as plt
import numpy as np
import os

main_dir = "/Users/kaimibk/Documents/ASTR_301/A1/"
SED_dir = main_dir+"Galaxy SEDs and Filter curves/GalaxySEDs/"
Filter_dir = main_dir+"Galaxy SEDs and Filter curves/FilterCurves/"
out_dir = main_dir+"output/"

plt.clf()

#'''
f, ax = plt.subplots()
ax.set_xscale("log", nonposx='clip')
ax.set_yscale("log", nonposy='clip')

for filename in os.listdir(SED_dir):
    if filename.endswith(".sed"):
        label = filename[0:-5]
        X,Y = np.genfromtxt(SED_dir+filename, unpack=True)
        ax.plot(X, Y, label=label)

ax.set_xlabel("Wavelength [$\AA$]")
ax.set_ylabel("Relative Flux Density")
ax.legend(bbox_to_anchor=(1.01, 0.7), loc=2, borderaxespad=0.)
#plt.show()
plt.savefig(out_dir+"SEDs.png", bbox_inches="tight")

#'''

'''
plt.clf()
## too lazy to auto do the colors,
colors = ("blue","orange", "g", "violet", "red")

i = 0
for filename in os.listdir(Filter_dir):
    if filename.endswith(".txt"):
        label = filename[0:-4]
        X,Y = np.genfromtxt(Filter_dir+filename, unpack=True)
        plt.plot(X, Y, label=label, color=colors[i])
        plt.fill_between(X, np.zeros(len(X)), Y, alpha=0.4, color=colors[i])

        i += 1

#plt.xlabel("log$_{10}$($\ lambda$)")
#plt.ylabel("log$_{10}$($f$)")
plt.legend(bbox_to_anchor=(1.01, 0.7), loc=2, borderaxespad=0.)


plt.xlabel("$\lambda$ [nm]")

plt.ylim(0, 100)
plt.ylabel("Stuff [units]")

plt.grid()
plt.savefig(out_dir+"filters.png", bbox_inches="tight")
'''
