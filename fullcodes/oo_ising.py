from ising_class import *
Temp=2.3
n=20
ntrials=500000
nequil=100000

ising1=ising(Temp,n)
ising1.randomize()
ising1.trials(nequil)
ising1.resetprops()
for i in range(ntrials):
    ising1.trial()
    ising1.addprops()
ising1.calcprops()
ising1.printsys()

