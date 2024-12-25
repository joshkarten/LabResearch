from ClassesSim import MCCTClassical
import numpy as np
import os
import line_profiler
import sys
os.environ["LINE_PROFILE"]='1'
betarange = [0]#np.array([0,0.2,0.3,0.4,0.44,(2/np.log(1+np.sqrt(2)))**(-1),0.45,0.5,0.6,0.7,1])
probrange= [0]#np.array([0,10,30,45,50,55,70,90,100])#np.array(sys.argv[4],dtype = int)
times = 2
chains=1
length=100
cf=1
mf=1
timeSim=MCCTClassical(length,chains,10,length*20,betarange,probrange,2,-1,cf,mf)
timeSim.Simulation()
print('hi')
    

