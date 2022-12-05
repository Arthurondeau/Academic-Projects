import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.collections import EllipseCollection
from matplotlib.patches import Rectangle
from simul import Simul
import time


class AnimateSimul:
    def __init__(self, simulation):
        print('AnimateSimul')
        self.simulation = simulation
        self.fig, self.ax = plt.subplots(figsize=(5, 5))  # initialise  graphics
        self.circles = EllipseCollection(widths=2*simulation.sigma, heights=2*simulation.sigma, angles=0, units='x',
                                         offsets=simulation.position, transOffset=self.ax.transData)  # circles at pos
        self.ax.add_collection(self.circles)
        rect = Rectangle((0, 0), simulation.L, simulation.L, ec='black', facecolor='none')   #   enclosing box
        self.ax.set_xlim(left=-.5, right=simulation.L+.5)  # viewport
        self.ax.set_ylim(bottom=-.5, top=simulation.L+.5)
        self.ax.add_patch(rect)  # draw the box
        self.ax.set_aspect(1)

    def _anim_step(self, m):  # m is the number of calls that have occurred to this function
        print('anim_step m = ', m)
        if m == 0 :
            time.sleep(1)

        self.simulation.md_step()  # perform simulation step
        self.circles.set_offsets(self.simulation.position) 

    def go(self, nframes):
        print('go')
        self._ani = animation.FuncAnimation(self.fig, func=self._anim_step, frames=nframes,
                                      repeat=False, interval=20)  # run animation
        plt.show()
