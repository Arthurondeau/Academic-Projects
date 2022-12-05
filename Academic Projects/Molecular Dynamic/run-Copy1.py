from simul import Simul
from animatesimul import AnimateSimul
import numpy as np
import matplotlib.pyplot as plt

def main():
    sample_time=0.5
    np.random.seed(10)  # set random numbers to be always the same
    simulation = Simul(sample_time, sigma=0.2 , L=4,Np=15) #  sigma particle radius # L box size
    print(simulation.__doc__)  # print the documentation from the class

    animate = AnimateSimul(simulation)
    animate.go(nframes=100)  # number of animation steps
    print(simulation)  #  print last configuration to screen

    
    
if __name__ == '__main__':
    main()
