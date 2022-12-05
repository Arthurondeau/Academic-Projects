import numpy as np


class Simul:
    """ 
    This is the prototype of the simulation code
    It moves the particles with at _velocity, using a vector notation: numpy should be used.
    """
    def __init__(self, sample_time, sigma_m,sigma_p,omega,L,Nb,No):
        np.seterr(all='ignore')  # remove errors in where statements
        print("Simul init")
        #self.position = np.array([[0.25, 0.25], [0.75, 0.25], [0.25, 0.75], [0.75, 0.75]])  # starting positions
        #self.position = np.array([[1., 1.], [3., 1.], [1., 3.], [3., 3.]])  # starting positions
        self.No = No
        self.Nb = Nb
        self.Np = Nb+No
        self.sigma = sigma_p # blue particle radius
        self.sigma_m = sigma_m
        #sigma_o = sigma_p*f # orange particle radius, must be lower than the distance betwee membrane particule and where is the ratio of radii particles and f² is the ratio of particle masses
        
        sigma_o = sigma_p
        f = sigma_o/sigma_p
        self._sample_time = sample_time
        self.L = np.array([3*L,L]) # boite rectangulaire
        self.Lgrid = [self.L for i in range(self.Np)]
        
        self.omega = omega # distance between membrane particle
        self.Nm = int(L//(self.omega+3*sigma_m))+1 # Number of membrane particule
        self.Lgridall = [self.L for i in range(self.Np+self.Nm)]
        positionY_m = np.linspace(sigma_m,self.L[1]-sigma_m,self.Nm) # on positionne les particules à une distance minimale égale à leur rayon et à une distance maximale égale à longueur Ly moins leur rayon.
        self.position_m = np.zeros((self.Nm,2))
        self.position_m[:,0] = 3*L/2
        self.position_m[:,1] = positionY_m
        self.position = np.random.random_sample((Nb,2))*(self.L-2*sigma_p) + sigma_p # position of blue particle
        self.position_o = np.random.random_sample((No,2))*(self.L-2*sigma_o) + sigma_o # position of orange particle
        self._velocity = np.random.normal(size=self.position.shape)  # random velocities of blue particle
        self._velocity_o = np.random.normal(size=self.position_o.shape)  # random velocities of blue particle
        self._i, self._j = np.triu_indices(self.position.shape[0], k=1)  # all pairs of indices between particles
        self.mass_b = [1 for i in range(self.Nb)]
        self.mass_o = [f**2 for i in range(self.No)]
        self.mass_m = [10**(9) for i in range(self.Nm)]
        
        #Creation of lists for the function wall_time,pair_time,md_step
        #The most convenient way to implement two types of particules is to create an array by concatenating the two position array of each type of particle (same for their radius). Thanks to that definition we can keep the same formalism than before.
        self.sigma_p = np.zeros(Nb+No)
        self.sigma_p[:Nb] = sigma_p
        self.sigma_p[Nb::No] = sigma_o
        self.sigma_pgrid = np.zeros((self.Np,2))
        self.sigma_pgrid[:,0] = self.sigma_p
        self.sigma_pgrid[:,1] = self.sigma_p
        
        self.position_p = np.concatenate((self.position,self.position_o),axis=0)
        self._velocity_p = np.random.normal(size=self.position_p.shape)
        
        #creation of the lists for the animation
        self.position_all = np.concatenate((self.position_p,self.position_m), axis=0)
        self._velocity_m = np.zeros((self.Nm,2))
        self._velocity_all = np.concatenate((self._velocity_p,self._velocity_m), axis=0)
        
        self.sigma_all = np.zeros(Nb+No+self.Nm)
        self.sigma_all[:Nb] = sigma_p
        self.sigma_all[Nb:No+Nb] = sigma_o
        self.sigma_all[No+Nb:] = sigma_m
        self.mass_all = np.concatenate((self.mass_b,self.mass_o,self.mass_m), axis=0)
        self.sigma_pgridall = np.zeros((self.Np+self.Nm,2))
        self.sigma_pgridall[:,0] = self.sigma_all
        self.sigma_pgridall[:,1] = self.sigma_all
        #self.pressure = []

    def _wall_time(self):
	#Index_positive = np.where(velocity >0)
	#Index_negative = np.where(velocity < 0)
	#time_part = L - position[Index_positive]/velocity[Index_positive] + position[negative]/velocity[Index_negative]
        #J'ai trouvé le problème ! Dans la fonction md_step on met a jour les coordonnées _all mais à aucun moment les coordonées 
        # _p sont elles mises à jour ! 
        velocity = self._velocity_all[:self.Np]
        position = self.position_all[:self.Np]
        time_part = np.where(velocity > 0, (self.Lgrid - self.sigma_pgrid- position) / velocity, (self.sigma_pgrid-position / velocity))
        #find the disk and direction corresponding to smallest v0
        coord = np.unravel_index(time_part.argmin(), time_part.shape)
        first_collision_time = np.min(time_part)
        particle = coord[0]
        direction = coord[1]
        
        return first_collision_time, particle, direction
	# calculate time of first collision, particle involved and direction
        

    def _pair_time(self):
        i, j = np.triu_indices(self.position_all.shape[0], k=1)
        rij = self.position_all[i]-self.position_all[j]  # set of all 6 separation vectors
        vij = self._velocity_all[i]-self._velocity_all[j]
        vij_2 = (vij**2).sum(1)
        rij_2 = (rij**2).sum(1) 
        b = (rij*vij).sum(1)
        c = (rij_2-(self.sigma_all[i]+self.sigma_all[j])**2)
        Delta = 4*((rij*vij).sum(1))**2 - 4*vij_2*c
        
        # Calcul du discriminant pour le calcul des racines (temps de collision)
        time_collision = np.where((Delta>0) & (b<0)& (c>0), (-2*b-np.sqrt(Delta))/(2*vij_2),np.infty)
        #print('time collision matrix',time_collision)
        #print('time_collision',time_collision)
        # index pour les couples ij correspondant aux racines réelles et la condition sur b correspond aux racines positives
        #Il faut maintenant prendre le minimum des temps positifs
        
        couple_collision_coord = np.unravel_index(time_collision.argmin(), time_collision.shape)
        time_collision = np.min(time_collision)
        coord_particule = [i[couple_collision_coord],j[couple_collision_coord]]
        return time_collision,coord_particule,c 
    
                                  
    def md_step(self):
        print('Simul::md_step')
        
        pressure = 0 
        current_time = 0 
        w_time, particle, direction = self._wall_time()
        choc_time,particule, c = self._pair_time()
        bool_c = np.where(c>0,True,False) # liste testant si des couples ij de particules se superposent ( vrai si pas de superposition faux sinon) 
        bool_c = np.all(bool_c) # renvoie vrai si toutes les particules sont non superposées
        #print("arguments pair time", choc_time,particule)

        while (w_time < self._sample_time -current_time) | (choc_time < self._sample_time -current_time) :   # think about this
            #Update position but with two differents cases
            if(w_time<choc_time) :
                
                self.position_all += w_time* self._velocity_all
                #Update the pressure
                pressure += 2*abs(self._velocity_all[particle][direction])/(self._sample_time*self.L[direction])
                self._velocity_all[particle][direction] = - self._velocity_all[particle][direction]

                #Update current time 
                current_time += w_time
                w_time, particle, direction = self._wall_time()  # update collisions times
                choc_time,particule, c = self._pair_time()
            else :
                self.position_all += choc_time* self._velocity_all
                r = self.position_all[particule[0]] - self.position_all[particule[1]]
                r = r/np.linalg.norm(r)
                m1,m2 = self.mass_all[particule[0]], self.mass_all[particule[1]]
                
                delta_v = self._velocity_all[particule[0]]-self._velocity_all[particule[1]] # vecteur des différences des vitesses entre les particules
                self._velocity_all[particule[0]] += - (2*m2/(m1+m2))*r*((r*delta_v).sum())
                self._velocity_all[particule[1]] += + (2*m1/(m1+m2))*r*((r*delta_v).sum())
                
                #Update current time 
                current_time += choc_time
                choc_time,particule, c  = self._pair_time()  # update collisions times
                w_time, particle, direction = self._wall_time()  # update collisions times
                
            
            
            
        #Update position after the last collision
        self.position_all += (self._sample_time-current_time) * self._velocity_all
        
        #self.pressure.append(pressure)
        return pressure, bool_c

    def __str__(self):   # this is used to print the position and velocity of the particles

        p = np.array2string(self.position_all)
        v = np.array2string(self._velocity_all)
        return 'pos= '+p+'\n'+'vel= '+v+'\n'
