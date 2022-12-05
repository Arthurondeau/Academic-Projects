from simul_osmos import Simul
from animatesimul import AnimateSimul
import numpy as np
import matplotlib.pyplot as plt

def main():
    sample_time=0.5 #Time of each iteration
    np.random.seed(10)  # set random numbers to be always the same
    simulation = Simul(sample_time,sigma_m = 0.1, sigma_p=0.3 ,omega = 0.1, L=4,Nb=10,No=2) 
    #  sigma_p particle radius # sigma_m particle membrane radius # L box size # omega distance between each membrane particle # Nb number of blue particle # No number of orange particle
    print(simulation.__doc__)  # print the documentation from the class

    animate = AnimateSimul(simulation)
    animate.go(nframes=100)  # number of animation steps
    print(simulation)  #  print last configuration to screen

    
    
if __name__ == '__main__':
    main()
