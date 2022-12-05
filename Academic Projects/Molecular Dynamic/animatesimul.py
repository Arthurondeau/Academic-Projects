import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.collections import EllipseCollection
from matplotlib.patches import Rectangle
from simul_osmos import Simul
import time


class AnimateSimul:
    def __init__(self, simulation):
        print('AnimateSimul')
        self.simulation = simulation
        self.fig, self.ax = plt.subplots(figsize=(5, 5))  # initialise  graphics
        self.circles = EllipseCollection(widths=2*simulation.sigma_all, heights=2*simulation.sigma_all,
        angles=0,units='x',offsets=simulation.position_all, transOffset=self.ax.transData)  # circles at pos
        #self.circles_membrane = EllipseCollection(widths=2*simulation.sigma_m, heights=2*simulation.sigma_m, angles=0, units='x',         offsets=simulation.position_m, transOffset=self.ax.transData)  # circles at pos
        self.ax.add_collection(self.circles)
        #self.ax.add_collection(self.circles_membrane) # problème ici ! 
       
        rect = Rectangle((0, 0), simulation.L[0], simulation.L[1], ec='black', facecolor='none')   #   enclosing box
        self.ax.set_xlim(left=-.5, right=simulation.L[0]+.5)  # viewport
        self.ax.set_ylim(bottom=-.5, top=simulation.L[1]+.5)
        self.ax.add_patch(rect)  # draw the box
        self.ax.set_aspect(1)
        Nm = simulation.Nm
        No = simulation.No
        Nb = simulation.Nb
        color=[]
        #listcolor = ['blue','orange','black'] # liste avec deux types de particules
        listcolor = ['blue','blue','black']
        for i,N in enumerate([Nb,No,Nm]) :
            for k in range(N) :
                color.append(listcolor[i])
        self.circles.set_color(color)
        

    def _anim_step(self, m):  # m is the number of calls that have occurred to this function
        print('anim_step m = ', m)
        if m == 0 :
            time.sleep(1)

        self.simulation.md_step()  # perform simulation step
        self.circles.set_offsets(self.simulation.position_all)

    def go(self, nframes):
        print('go')
        self._ani = animation.FuncAnimation(self.fig, func=self._anim_step, frames=nframes,
                                      repeat=False, interval=20)  # run animation
        plt.show()
