#!/usr/bin/python
## Objected oriented Ising model for Chem 160 ##
import random  #Need for randum number generation
import math    #Need for exponental

# Define spin class
class spin:
    def __init__(self):
        self.spin=1
    def randomize(self):
        self.spin=1
        if random.random()<0.5:
            self.spin=-1
    def myspin(self):
        return(self.spin)
    def flip(self):
        self.spin=-self.spin
    def char(self):
        if self.spin==1:
            return "U"
        else:
            return "D"

# Define cell class
class cell:
    def __init__(self, i, j, n):
        self.spin=spin()
        self.left=((i-1+n)%n,j)
        self.right=((i+1)%n,j)
        self.up=(i,(j+1)%n)
        self.down=(i,(j-1+n)%n)
        self.propcount=0
    def cellspin(self):
        return self.spin.myspin()
    def energy(self,cells):
        energy=0
        energy-=cells[self.left[0]][self.left[1]].cellspin()
        energy-=cells[self.right[0]][self.right[1]].cellspin()
        energy-=cells[self.up[0]][self.up[1]].cellspin()
        energy-=cells[self.down[0]][self.down[1]].cellspin()
        energy*=self.cellspin()
        return energy

# Define ising class
class ising:
    def __init__(self, Temp, n):
        self.Temp=Temp
        self.n=n
        self.n2=n*n
        #Create lattice of cells
        self.cells=[]
        for i in range(n):
            self.cells.append([])
            for j in range(n):
                acell=cell(i,j,n)
                self.cells[i].append(acell)
        self.energy=0.0
        self.magnetism=0.0
    def randomize(self):
        for i in range(self.n):
            for j in range(self.n):
                self.cells[i][j].spin.randomize()
    def changeTemp(self,newTemp):
        self.Temp=newTemp
    def printsys(self):
        for i in range(self.n):
            for j in range(self.n):
                print("%2s"%(self.cells[i][j].spin.char()),end=" ")
            print("")
    def cellenergy(self,i,j):
        return self.cells[i][j].energy(self.cells)
    def ave_energy(self):
        energy=0
        for i in range(self.n):
            for j in range(self.n):
                energy+=self.cellenergy(i,j)
        return energy/(self.n*self.n)
    def ave_magnetization(self):
        magnetization=0
        for i in range(self.n):
            for j in range(self.n):
                magnetization+=self.cells[i][j].cellspin()
        return abs(magnetization)/(self.n*self.n)
    def fliptry(self):
        i=int(self.n*random.random())
        j=int(self.n*random.random())
        deltaE=-2.*self.cellenergy(i,j)
        deltaM=-2.*self.cells[i][j].cellspin()
        boltzmann=-deltaE/self.Temp
        #print(i, j)
        #print(self.cellenergy(i,j), boltzmann, math.exp(boltzmann))
        if math.exp(boltzmann)>random.random():
            self.cells[i][j].spin.flip()
            return deltaE,deltaM
        return 0.,0.
    def trial(self):
        deltaE,deltaM=self.fliptry()
        self.energy+=2*deltaE/self.n2
        self.magnetism+=deltaM/self.n2
    def trials(self,ntrials):
        for i in range(ntrials):
            self.trial()
    def addprops(self):
        if self.energy_flag:
            self.energy=self.ave_energy()
            self.magnetism=self.ave_magnetization()
            self.energy_flag=False
        energy=self.energy
        magnetism=self.magnetism
        self.sum_energy+=energy
        self.sum_energy2+=energy**2
        self.sum_magnetism+=magnetism
        self.sum_magnetism2+=magnetism**2
        self.propcount+=1
    def resetprops(self):
        self.energy_flag=True
        self.sum_energy=0
        self.sum_energy2=0
        self.sum_magnetism=0
        self.sum_magnetism2=0
        self.propcount=0
    def calcprops(self):
        energy_ave=self.sum_energy/self.propcount
        magnetism_ave=self.sum_magnetism/self.propcount
        energy2_ave=self.sum_energy2/self.propcount
        magnetism2_ave=self.sum_magnetism2/self.propcount
        Cv=1./self.Temp**2*(energy2_ave-energy_ave**2)
        magnetic_susceptibility=1./self.Temp*(magnetism2_ave-magnetism_ave**2)
        print("%8.4f %10.6f %10.6f %10.6f %10.6f"%\
            (self.Temp, energy_ave, Cv, magnetism_ave, magnetic_susceptibility))
