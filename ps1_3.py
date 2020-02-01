#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np




class Map:
    def __init__(self, map_type):
        self.description = map_type

    def transition(self, x, r):
        pass

    def transition_range(self, x, r, n):
        xn = x
        for i in range(n):
            xn = self.transition(x, r)
        return xn


class Logistic(Map):
    def __init__(self):
        super().__init__("Logistic Map")

    def transition(self, x, r):
        return r * x * (1 - x)


class Henon(Map):
    def __init__(self, y0, b):
        super().__init__("Henon Map")
        self.y = y0
        self.b = b

    def transition(self, x, a):
        xn = self.y + 1 - a * x**2
        self.y = self.b * xn
        return xn









##
# @brief Displays a bifurcation plot of m iterates
# within the 
#
# @param rmin Minimum r value
# @param rmax Maximum r value
# @param iterates Number of iterations per R value
# @param cut_off Starting iteration for x, algorithm starts at x0, plots cut_off number of
# @param points, then plots iterates number of points (so a total of cut_off + iterates distinct iterations are
# carried out)
# @param x0 Initial x value
# @param dr R increment
# @param map_func The mapping function should take in one number jand spit out another number
def biffurcation(rmin, rmax, iterates, x0, mapping: Map, cut_off=0, dr=0.01):
    
    r_vals = np.arange(rmin, rmax, dr, dtype=np.float64)

    # Array to hold outputs (only sized iterates)
    xarr = np.zeros(iterates, dtype=np.float64)
    
    

    for i in range(r_vals.shape[0]):

        x = x0
        for j in range(cut_off):
            x = mapping.transition(x, r_vals[i])


        for j in range(iterates):
            x = mapping.transition(x, r_vals[i])
            xarr[j] = x

        plt.scatter([r_vals[i] for r in range(iterates)], xarr, c="black", s=0.0003)

    plt.xlabel('R Value')
    plt.ylabel(f'x_n (iterations (m), cut_off (l) = {iterates, cut_off}')
    plt.title(f'Bifuracation Diagram for the {mapping.description} x0 = {x0}')

def feigenbaum(rmin, rmax, iterates, x0, mapping: Map, cut_off=0, dr=0.01, niterates = 1):

    # Array to hold outputs (only sized iterates)
    xarr = np.zeros(iterates, dtype=np.float64)
    r_vals = np.arange(rmin, rmax, dr, dtype=np.float64)

    std_prev = 0
    stds = np.zeros(r_vals.shape[0])


    for i in range(r_vals.shape[0]):

        x = x0
        for j in range(cut_off):
            x = mapping.transition(x, r_vals[i])


        for j in range(iterates):
            x = mapping.transition(x, r_vals[i])
            xarr[j] = x

        std_temp = np.std(xarr[0::niterates])
        print(std_temp)

        stds[i] = std_temp - std_prev

        std_prev = std_temp

    plt.plot(r_vals, stds, label=f'Cycle = {niterates}')







            









if __name__ == '__main__':
    
    log = Logistic()
    henon = Henon(0.3, 0.3)

    #  biffurcation(rmin=2.8, rmax=4.0, x0=0.5, cut_off=500, iterates=600, dr=0.0005, mapping = log)
    #  biffurcation(rmin=0, rmax=2, x0=0.5, cut_off=500, iterates=600, dr=0.0005, mapping = henon)
    #
    #
    #
    max_r = 2
    min_r = 0
    x0 = 0.5
    cut_off = 500
    iterates = 600
    dr = 0.000001
    mapping = log

    #  max_r = 3.57
    #  min_r = 2.8
    #  x0 = 0.5
    #  cut_off = 500
    #  iterates = 600
    #  dr = 0.0001
    #  mapping = log

    feigenbaum(rmin=2.9994, rmax=3.0006, x0=x0, cut_off=cut_off, iterates=iterates, dr=dr, mapping = mapping, niterates=1)
    #  feigenbaum(rmin=min_r, rmax=max_r, x0=x0, cut_off=cut_off, iterates=iterates, dr=dr, mapping = mapping, niterates=2)
    #  feigenbaum(rmin=min_r, rmax=max_r, x0=x0, cut_off=cut_off, iterates=iterates, dr=dr, mapping = mapping, niterates=4)
    #  feigenbaum(rmin=min_r, rmax=max_r, x0=x0, cut_off=cut_off, iterates=iterates, dr=dr, mapping = mapping, niterates=8)

    print("here")
    #  plt.xticks( np.arange(min_r, max_r, dr) )


    plt.xlabel('R Value')
    plt.ylabel('Standard Deviation of N cycles')
    plt.title('Analysis of n cycles')
    plt.legend()
    plt.grid()
    plt.show()










