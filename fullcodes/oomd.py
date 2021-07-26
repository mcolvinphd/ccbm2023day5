#!/usr/bin/python
### Soft Sphere Molecular Dynamics Program ###
########### Prof. Michael Colvin #############
######## School of Natural Sciences ##########
################# UC Merced ##################

#Computational core and some parameter names based closely on mdss.f program in
#Molecular Dynamics Simulation by J.M. Haile, John Wiley & Sons, 1992

###### Run Simulations #######
from oomd_class import *
import time

#Set system state--used to create the system object
NAtom=108            #Number of atoms
# Set density and temperatures (gm/cm^3 and K)
Dens=.1
Temp=150.

#Set number of equilibration and simulation time steps
MaxEqb=100
MaxStep=1000

# Convert to dimensionless units
Density=Dens/1.68247         #Density  (1.0=1.68247 gm/cm^3)
Temperature=Temp/119.4      #Temperature (1.0=119.4K )


###### Start timing ######
start_time=time.process_time()

###### Set random seed ######
random.seed(1)

###### Setup Simulations #######       
asystem=system(NAtom, Density, Temperature, MaxEqb, MaxStep) #Create system object
asystem.setfcc()          #Fill with NAtom atoms on FCC lattice
asystem.setvel()          #Set atom velocities
asystem.printparam()      #Print simulation parameters
asystem.evalforce()       #Initialize force calculation
asystem.scaleforces()     #Scale forces for initial step
asystem.printproplabel()  #Print header for property output

###### Run Equilibration #######       
for step in range(1,asystem.MaxEqb+1):
    asystem.setstep(step)
    asystem.update()      #Update positions and velocites
    asystem.properties()  #Calculate the properties every step during equil
    asystem.printproperties(1) #Print out properties
    asystem.escale()           #Rescale energy

###### Run Simulation #######       
print('Starting Production run')
equil_time=time.process_time()
asystem.resetprop()
asystem.printproplabel()

for step in range(1,asystem.MaxStep+1):
    asystem.setstep(step)
    asystem.update()
    asystem.propupdate()

#End of loop over production timesteps
print('Average temperature=%lf (%lf K)'%(asystem.TempAve,asystem.TempAve*119.4))
print('Average pressure=%lf'%(asystem.PressAve))

if asystem.Kdiffuse:
    asystem.diffuseave()
if asystem.Kvelauto:
    asystem.velautocorr()
    asystem.printvelautocorr()
    print('Diffusion constant from vel autocorr=%lf'%(asystem.autodiffconst))
    diffusion_cgs=asystem.autodiffconst*0.1344*asystem.Step*100000.
    print('Diffusion constant from vel autocorr=%lf (x10-5 cm^2/sec)'%(diffusion_cgs))
                                                    
#Print timing info
stop_time=time.process_time()
print("Equilibration time=%.3f minutes"%((equil_time-start_time)/60.))
print("   Production time=%.3f minutes"%((stop_time-equil_time)/60.))
print("    Total run time=%.3f minutes"%((stop_time-start_time)/60.))
