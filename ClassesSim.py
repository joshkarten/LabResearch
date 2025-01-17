from fractions import Fraction
from itertools import repeat
import random
import numpy as np
import os
from line_profiler import profile

class MCCTClassical:
    def __init__(self,Length, numChains,sampleFreq,totSteps, betarange, probrange, iterations, coupling,ControlFreq=1, MonteCarloFreq=1,boltzmannMod=1/2, perLatSweep=True  ):
        self.Length = Length # how long is each bit string
        self.sampleFreq = sampleFreq #after how many update steps are the parameters calculated
        self.numChains = numChains # how many bit strings are there
        self.ControlFreq=ControlFreq # how many times per step is the stochasitic control performed
        self.iterations = iterations # how many times is the simulation repeated for each set a variables
        self.K=coupling # how strong the nearest neighbor spin interactions are
        self.MonteCarloFreq=MonteCarloFreq # how many times per step is the monte carlo step performed
        self.betarange=betarange # over what thermodynamic betas are we sampling
        self.probrange=probrange # over what probability of selecting the control map (inclusive) are we sampling
        self.allones = (1<<(self.Length))
        self.totSteps=totSteps # total number of update steps
        self.xj = np.array([Fraction(1/3),Fraction(2/3)])# fixed points for the control map
        self.xjint =np.array([self.dec2int(self.xj[0]),self.dec2int(self.xj[1])])
        self.twoChains=True if numChains==2 else False
        self.latticeSize=Length*numChains
        self.boltzmannMod=boltzmannMod
        if perLatSweep:
            self.mcrep = self.latticeSize
        else:
            self.mcrep=1
        # number of data points along the time axis
        # we add 50 additional points to cluster our samples close to the start
        if not (self.totSteps%self.sampleFreq):
            self.time_size = (self.totSteps)//self.sampleFreq+1
            
        else: 
            self.time_size = (self.totSteps)//self.sampleFreq+2
        self.setDimensionality() # chooses functions to used based on number of chains
        self.createSamplingArrays() # these arrays store the calculated parameters
        if self.MonteCarloFreq==0:
            betarange=[0]
        if self.ControlFreq==0:
            probrange=[0]    
        
    def dec2int(self, num):
        """
        Conversion of a decimal between 0 and 1 into an integer between 0 and 2**Length-1
        :type num: Fraction
        :type length: int
        """
        assert 0 <= num < 1, self.Length > 0
        return int(num * (2 ** self.Length))


    def reset(self,bin_num):
        """
        :type bin_num: int
        odd numbers end with b_L =1
        num - 1 for set b_L = 0 if odd
        """
        #if bin_num % 2 == 1:
        #    bin_num = bin_num - 1
        return bin_num&(~1)


    def left(self,num):
        """

        :param num: int
        :param length: int
        :return: int
        """
        #if num >= 2 << (length - 2):
        #    num = (num << 1) ^ 1
        #else:
        #    num = num << 1
        return (num>>self.Length-1)^(num<<1)&((1<<self.Length)-1)#num % (2 ** (length))


    def right(self,num):
        """

        :param num: int
        :return: int
        """
        num = self.reset(num)
        # if num % 2 ==1 and (num//2) % 2 == 0 and (num//4)%2 == 1:
        #    num = num + 2
        num = num // 2
        return num

    @profile
    def order_parameter(self,num):
        '''
        Calculates sigma_i sigma_j in an integer represenation of a bitstring.'''
        binrep =bin(self.left(num)^num) # \sigma_{i+1}*\sigma_{i}
        count = 2*binrep.count('1')/self.Length-1

        return count

    @profile
    def bernoulli(self,num):
        num = self.left(num)
        num = num & (~0b111)
        return num + random.getrandbits(3)

    @profile
    def control(self,num):

        if num < 2 ** (self.Length-1):
            xj = self.dec2int(Fraction(1, 3))
        else:
            xj = self.dec2int(Fraction(2, 3))
        if num == xj:
            return num
        else:
            temp = self.right(num) + self.right(xj)
        if temp == num:
            return num + 1
        return temp
    @profile
    def boltzmann_probability(self,initial, final, beta):
        '''
        Calculates a value from a modified Boltzmann distribution. The choice of 1/2 (which modifies the distribution) is to create bit scrambling when Beta=0 '''
        return np.exp(-beta * (final - initial))*self.boltzmannMod
        #boltzmann weight. we take 1/2 such that the state is continuously randomized at infinite
    @profile
    def BitConfiguration1d(self,position_x,position_y):
    #problem: individual calculations are incredibly inefficient
    #solution: Make lookup tables for everything
    #so, from the bits, we can just determine the energy without additional calculation.   
        '''        
        Returns the targeted bit and the nearest neighbor bits (up,down, left) in an integer.'''
        mainBit =(self.lattice[position_x]&(1<<position_y))>>(position_y)
        ubit = (self.lattice[position_x]&(1 << (position_y-1)%self.Length))>>(position_y-1)%self.Length
        dbit= ((self.lattice[position_x]&(1<<(position_y+1)%self.Length))>>(position_y+1)%self.Length)
        return mainBit,ubit,dbit

    def BitConfigurationLad(self,position_x,position_y):
    #problem: individual calculations are incredibly inefficient
    #solution: Make lookup tables for everything
    #so, from the bits, we can just determine the energy without additional calculation.   
        '''        
        Returns the targeted bit and the nearest neighbor bits (up,down, left) in an integer.'''
        mainBit =(self.lattice[position_x]&(1<<position_y))>>(position_y)
        ubit = (self.lattice[position_x]&(1 << (position_y-1)%self.Length))>>(position_y-1)%self.Length
        dbit= ((self.lattice[position_x]&(1<<(position_y+1)%self.Length))>>(position_y+1)%self.Length)
        lbit =((self.lattice[position_x-1]&(1<<position_y))>>(position_y))
        return mainBit,ubit,dbit, lbit

    def BitConfiguration2d(self,position_x,position_y):
    #problem: individual calculations are incredibly inefficient
    #solution: Make lookup tables for everything
    #so, from the bits, we can just determine the energy without additional calculation.   
        '''        
        Returns the targeted bit and the nearest neighbor bits (up,down, left, right) in an integer.'''
        mainBit =(self.lattice[position_x]&(1<<position_y))>>(position_y)
        ubit = (self.lattice[position_x]&(1 << (position_y-1)%self.Length))>>(position_y-1)%self.Length
        dbit= ((self.lattice[position_x]&(1<<(position_y+1)%self.Length))>>(position_y+1)%self.Length)
        lbit =((self.lattice[position_x-1]&(1<<position_y))>>(position_y))
        rbit=(self.lattice[(1+position_x)%self.numChains]&(1<<position_y)>>position_y)
        return mainBit,ubit,dbit,lbit,rbit
   
    # The lookup dictionary is designed around taking in as many bits as there are nearest neighbors.
    #  Thus this imp has a bunch of different tables for each dimensionality
    @profile
    def LookupEnergy1d(self,K,Beta):
        self.energyDict={}
        self.energyDict[0,0,0] =min(self.boltzmannMod,self.boltzmann_probability(2*K,-2*K,Beta))
        self.energyDict[1,1,1]=self.energyDict[0,0,0]

        self.energyDict[0,1,0]=min(self.boltzmannMod,self.boltzmann_probability(0,0,Beta))
        self.energyDict[0,0,1]=self.energyDict[0,1,0]
        self.energyDict[1,0,1]=self.energyDict[0,1,0]
        self.energyDict[1,1,0]=self.energyDict[0,1,0]

        self.energyDict[0,1,1]=min(self.boltzmannMod,self.boltzmann_probability(-2*K,2*K,Beta))
        self.energyDict[1,0,0]=self.energyDict[0,1,1]
        return self.energyDict
    
    def LookupEnergyLad2d(self,K,Beta):
        self.energyDict = {}
        
        self.energyDict[0,0,0,1] =  min(self.boltzmannMod,self.boltzmann_probability(K,-K,Beta))
        self.energyDict[0,0,1,0]=self.energyDict[0,0,0,1]
        self.energyDict[0,1,0,0]=self.energyDict[0,0,0,1]

        self.energyDict[1,0,0,1]=min(self.boltzmannMod, self.boltzmann_probability(-K,K,Beta))
        self.energyDict[1,0,1,0]=self.energyDict[1,0,0,1]
        self.energyDict[1,1,0,0]=self.energyDict[1,0,0,1]

        self.energyDict[0,0,0,0]=min(self.boltzmannMod,self.boltzmann_probability(K*3,K*-3,Beta))
        self.energyDict[1,0,0,0]=min(self.boltzmannMod,self.boltzmann_probability(K*-3,K*3,Beta))
        self.energyDict[0,1,1,1]=self.energyDict[1,0,0,0]#8
        self.energyDict[1,1,1,1]=self.energyDict[0,0,0,0] 

        self.energyDict[1,0,1,1]=self.energyDict[0,0,0,1] 
        self.energyDict[1,1,0,1]=self.energyDict[0,1,0,0] 
        self.energyDict[1,1,1,0]=self.energyDict[0,1,0,0] 

        self.energyDict[0,1,1,0]=self.energyDict[1,0,0,1]
        self.energyDict[0,0,1,1]=self.energyDict[1,0,0,1]
        self.energyDict[0,1,0,1]=self.energyDict[1,0,0,1]
        return self.energyDict
    
    def LookupEnergy2d(self,K,Beta):
        self.energyDict = {}
        self.energyDict[0,0,0,1,1] = min(self.boltzmannMod,self.boltzmann_probability(0,0,Beta))
        self.energyDict[0,0,1,0,1]=self.energyDict[0,0,0,1,1]
        self.energyDict[0,1,0,0,1]=self.energyDict[0,0,0,1,1]
        self.energyDict[0,1,1,0,0]=self.energyDict[0,0,0,1,1] #4
        self.energyDict[0,1,0,1,0]=self.energyDict[0,0,0,1,1] #5
        self.energyDict[0,0,1,1,0]=self.energyDict[0,0,0,1,1] #6

        self.energyDict[1,0,0,1,1]=self.energyDict[0,0,0,1,1]
        self.energyDict[1,0,1,0,1]=self.energyDict[0,0,0,1,1]
        self.energyDict[1,1,0,0,1]=self.energyDict[0,0,0,1,1]
        self.energyDict[1,1,1,0,0]=self.energyDict[0,0,0,1,1]
        self.energyDict[1,1,0,1,0]=self.energyDict[0,0,0,1,1] #5
        self.energyDict[1,0,1,1,0]=self.energyDict[0,0,0,1,1] #6

        self.energyDict[0,0,0,0,0]=min(self.boltzmannMod,self.boltzmann_probability(K*4,K*-4,Beta))
        self.energyDict[1,0,0,0,0]=min(self.boltzmannMod,self.boltzmann_probability(K*-4,K*4,Beta))
        self.energyDict[0,1,1,1,1]=self.energyDict[1,0,0,0,0]#8
        self.energyDict[1,1,1,1,1]=self.energyDict[0,0,0,0,0] 

        self.energyDict[0,1,0,0,0]=min(self.boltzmannMod,self.boltzmann_probability(K*2,K*-2,Beta))
        self.energyDict[0,0,1,0,0]=self.energyDict[0,1,0,0,0] 
        self.energyDict[0,0,0,1,0]=self.energyDict[0,1,0,0,0] 
        self.energyDict[0,0,0,0,1]=self.energyDict[0,1,0,0,0]  #12

        self.energyDict[1,0,1,1,1]=self.energyDict[0,1,0,0,0] 
        self.energyDict[1,1,0,1,1]=self.energyDict[0,1,0,0,0] 
        self.energyDict[1,1,1,0,1]=self.energyDict[0,1,0,0,0] 
        self.energyDict[1,1,1,1,0]=self.energyDict[0,1,0,0,0]  #12

        self.energyDict[1,1,0,0,0]=min(self.boltzmannMod,self.boltzmann_probability(K*-2,K*2,Beta))
        self.energyDict[1,0,1,0,0]=self.energyDict[1,1,0,0,0] 
        self.energyDict[1,0,0,1,0]=self.energyDict[1,1,0,0,0] 
        self.energyDict[1,0,0,0,1]=self.energyDict[1,1,0,0,0]  #16
        self.energyDict[0,0,1,1,1]=min(self.boltzmannMod,self.boltzmann_probability(K*-2,K*2,Beta))
        self.energyDict[0,1,0,1,1]=self.energyDict[1,1,0,0,0] 
        self.energyDict[0,1,1,0,1]=self.energyDict[1,1,0,0,0] 
        self.energyDict[0,1,1,1,0]=self.energyDict[1,1,0,0,0]  #16
        return self.energyDict
    # check this for proper iplementation
    # verified for lattices [0,0],[15,15],[0,15],[15,0],[5,10],[0,10],[14,15]
    def LatticeOrderParameterLad(self):
        param = 0
        for i in self.lattice:
            param += (bin(self.left(i)^i).count('1')) # \sigma_(i,j)*\sigma_(i,j+1)
        param= (2*param/self.Length-2 +2*bin(self.lattice[0]^self.lattice[1]).count('1')/self.Length-1 )/3# \sigma_(i,j)*\sigma(i+j,j)
        return param
    
    def LatticeOrderParameter2d(self):
        param,param2 = 0,0
        for i in self.lattice:
            param += (bin(self.left(i,)^i).count('1')) # \sigma_(i,j)*\sigma_(i,j+1)
        param = (2*param/(self.latticeSize)-1)/2
        for i in range(self.numChains):
            param2+= (bin(self.lattice[i]^self.lattice[i-1]).count('1')) # \sigma_(i,j)*\sigma(i-1,j)
        param2 = (2*param2/(self.latticeSize)-1)/2
        return param+param2
#works as intended
    @profile
    def Magnetization(self):
        param =0
        for i in self.lattice:
            param+= bin(i).count('1')
        param = 1- 2*param/self.latticeSize
        return param
#works as intended
    @profile
    def StaggeredMagnetization(self):
        param =0
        for i in range(0,self.numChains,2):
            param+= bin(self.lattice[i]^self.xjint[1]).count('1')
        for i in range(1, self.numChains,2):
            param += bin(self.lattice[i]^self.xjint[0]).count('1')
        param = 1-2*param/self.latticeSize
        return param
    
    # These are done just for some small optimizations and convenience
    def setDimensionality(self):
        if self.numChains==1:
            self.LookupEnergy = self.LookupEnergy1d
            self.BitConfiguration = self.BitConfiguration1d
            self.stochasticControl=self.stochasticControl1d
            self.monteCarlo=self.monteCarlo1d
            self.LatticeOrderParameter = self.order_parameter
            self.Step=self.Step1D
        elif self.twoChains:
            self.LookupEnergy = self.LookupEnergyLad2d
            self.BitConfiguration = self.BitConfigurationLad
            self.stochasticControl=self.stochasticControl2d
            self.monteCarlo=self.monteCarloLadder
            self.LatticeOrderParameter=self.LatticeOrderParameterLad
            self.Step = self.StepLad
        else:
            self.LookupEnergy=self.LookupEnergy2d
            self.BitConfiguration = self.BitConfiguration2d
            self.stochasticControl=self.stochasticControl2d
            self.monteCarlo=self.monteCarlo2d
            self.LatticeOrderParameter=self.LatticeOrderParameter2d
            self.Step = self.Step2D
        return
    @profile
    def createLattice(self):
        lat = []
        for i in repeat(None,self.numChains):
            lat.append(random.getrandbits(self.Length))
        self.lattice=np.array(lat)
        return
    
    def createSamplingArrays(self):
        self.acceptance = np.zeros((self.iterations,self.time_size,len(self.betarange),len(self.probrange))) # keeps track of monte carlo acceptance per step
        self.time_size=2*self.time_size #take two samples each time 

        if self.numChains==1:
            self.record1= np.zeros((self.iterations,self.time_size,len(self.betarange),len(self.probrange))) #\sigma_i\sigma_j in a chain, effectively our energy
            self.recordMag =  np.zeros((self.iterations,self.time_size,len(self.betarange),len(self.probrange)))#Magnetization
            self.recordMagS =  np.zeros((self.iterations,self.time_size,len(self.betarange),len(self.probrange))) #staggered magnetization
            self.time_size=self.time_size//2

            return
        if self.twoChains:
            self.record1= np.zeros((self.iterations,self.time_size,len(self.betarange),len(self.probrange)))#1st chan
            self.record2= np.zeros((self.iterations,self.time_size,len(self.betarange),len(self.probrange)))# 2nd chain
            self.recordlong= np.zeros((self.iterations,self.time_size,len(self.betarange),len(self.probrange)))  # between chains
            
            self.recordMag =  np.zeros((self.iterations,self.time_size,len(self.betarange),len(self.probrange)))
            self.recordMagS =  np.zeros((self.iterations,self.time_size,len(self.betarange),len(self.probrange)))
            self.time_size=self.time_size//2

            return
        self.recordtot= np.zeros((self.iterations,self.time_size,len(self.betarange),len(self.probrange))) # over entire lattice
        self.recordMag =  np.zeros((self.iterations,self.time_size,len(self.betarange),len(self.probrange)))
        self.recordMagS =  np.zeros((self.iterations,self.time_size,len(self.betarange),len(self.probrange)))
        self.time_size=self.time_size//2
        return
    
    def stochasticControl2d(self,prob):
        for i in range(self.numChains):
            if random.random()>(prob):
                self.lattice[i]=self.bernoulli(self.lattice[i])
            else:
                self.lattice[i]=self.control(self.lattice[i])
        return
    # the correct phase transition is maintained in 1d
    @profile
    def stochasticControl1d(self,prob):
        if random.random()>(prob):
            self.lattice[0]=self.bernoulli(self.lattice[0])
        else:
            self.lattice[0]=self.control(self.lattice[0])
        return
    @profile
    def monteCarlo1d(self,time,betaNum,probNum,itt):
        if not (time%self.sampleFreq):
            for nr in repeat(None,self.mcrep):
                y_pos = random.randrange(self.Length)
                mbit,ubit,dbit = self.BitConfiguration(0,y_pos)
                if self.energyDict[mbit,ubit,dbit]>random.random(): 
                    self.acceptance[itt,time//self.sampleFreq,betaNum,probNum] += 1
                    self.lattice[0] = self.lattice[0]^(0b1<<y_pos)

        else:
            for nr in repeat(None,self.mcrep):
                y_pos = random.randrange(self.Length)
                mbit,ubit,dbit, = self.BitConfiguration(0, y_pos)
                if self.energyDict[mbit,ubit,dbit]>random.random(): 
                    self.lattice[0] = self.lattice[0]^(0b1<<y_pos)

        return
    
    def monteCarloLadder(self,time,betaNum,probNum,itt):
        if not (time%self.sampleFreq):
            for nr in repeat(None,self.mcrep):
                x_pos =  random.randrange(self.numChains)
                y_pos =  random.randrange(self.Length)

                mbit,ubit,dbit,lbit = self.BitConfiguration(x_pos,y_pos)

                if self.energyDict[mbit,ubit,dbit,lbit]>random.random(): 
                    self.acceptance[itt,time//self.sampleFreq,betaNum,probNum] += 1
                    self.lattice[x_pos] = self.lattice[x_pos]^(0b1<<y_pos)

        else:
            for nr in repeat(None,self.mcrep):
                x_pos =  random.randrange(self.numChains)
                y_pos = random.randrange(self.Length)
 
                mbit,ubit,dbit,lbit = self.BitConfiguration(x_pos,y_pos)

                if self.energyDict[mbit,ubit,dbit,lbit]>random.random(): 
                    self.lattice[x_pos] = self.lattice[x_pos]^(0b1<<y_pos)

        return   
    def monteCarloLadder2(self,time,betaNum,probNum1,probNum2,itt):
        if not (time%self.sampleFreq):
            for nr in repeat(None,self.mcrep):
                x_pos =  random.randrange(self.numChains)
                y_pos =  random.randrange(self.Length)

                mbit,ubit,dbit,lbit = self.BitConfiguration(x_pos,y_pos)

                if self.energyDict[mbit,ubit,dbit,lbit]>random.random(): 
                    self.acceptance[itt,time//self.sampleFreq,betaNum,probNum1,probNum2] += 1
                    self.lattice[x_pos] = self.lattice[x_pos]^(0b1<<y_pos)

        else:
            for nr in repeat(None,self.mcrep):
                x_pos =  random.randrange(self.numChains)
                y_pos = random.randrange(self.Length)
 
                mbit,ubit,dbit,lbit = self.BitConfiguration(x_pos,y_pos)

                if self.energyDict[mbit,ubit,dbit,lbit]>random.random(): 
                    self.lattice[x_pos] = self.lattice[x_pos]^(0b1<<y_pos)

        return             
    
    def monteCarlo2d(self,time,betaNum,probNum,itt):
        if not (time%self.mcrep):
            for nr in repeat(None,self.latticeSize):
                x_pos =  random.randrange(self.numChains)
                y_pos =  random.randrange(self.Length)

                mbit,ubit,dbit,lbit,rbit = self.BitConfiguration(x_pos,y_pos)

                if self.energyDict[mbit,ubit,dbit,lbit,rbit]>random.random(): 
                    self.acceptance[itt,time//self.sampleFreq,betaNum,probNum] += 1
                    self.lattice[x_pos] = self.lattice[x_pos]^(0b1<<y_pos)

        else:
            for nr in repeat(None,self.mcrep):
                x_pos =  random.randrange(self.numChains)
                y_pos =  random.randrange(self.Length)

                mbit,ubit,dbit,lbit,rbit = self.BitConfiguration(x_pos,y_pos)

                if self.energyDict[mbit,ubit,dbit,lbit,rbit]>random.random(): 
                    self.lattice[x_pos] = self.lattice[x_pos]^(0b1<<y_pos)
        return
    @profile
    def Step1D(self,time,prob,b,p,itt):
        if not (time%self.sampleFreq):
            self.record1[itt,2*time//self.sampleFreq,b,p]=self.order_parameter(self.lattice[0])
            self.recordMag[itt,2*time//self.sampleFreq,b,p]=self.Magnetization()
            self.recordMagS[itt,2*time//self.sampleFreq,b,p]=self.StaggeredMagnetization() 
            for controlfreq in repeat(None,self.ControlFreq):
                self.stochasticControl(prob)
            self.record1[itt,2*time//self.sampleFreq+1,b,p]=self.order_parameter(self.lattice[0])
            self.recordMag[itt,2*time//self.sampleFreq+1,b,p]=self.Magnetization()
            self.recordMagS[itt,2*time//self.sampleFreq+1,b,p]=self.StaggeredMagnetization()    
            for mcfreq in repeat(None,self.MonteCarloFreq):
               self.monteCarlo(time,b,p,itt)
        else:
            for controlfreq in repeat(None,self.ControlFreq):
                self.stochasticControl(prob)
            for mcfreq in repeat(None,self.MonteCarloFreq):
               self.monteCarlo(time,b,p,itt)
        return
    
    def StepLad(self,time,prob,b,p,itt):
        if not (time%self.sampleFreq):
            self.record1[itt,2*time//self.sampleFreq,b,p]=self.order_parameter(self.lattice[0])
            self.record2[itt,2*time//self.sampleFreq,b,p]=self.order_parameter(self.lattice[1])
            self.recordlong[itt,2*time//self.sampleFreq,b,p]=2*(bin(self.lattice[0]^self.lattice[1]).count('1'))/self.Length-1
            self.recordMag[itt,2*time//self.sampleFreq,b,p]=self.Magnetization()
            self.recordMagS[itt,2*time//self.sampleFreq,b,p]=self.StaggeredMagnetization() 
            for controlfreq in repeat(None,self.ControlFreq):
                self.stochasticControl(prob)
            self.record1[itt,2*time//self.sampleFreq+1,b,p]=self.order_parameter(self.lattice[0])
            self.record2[itt,2*time//self.sampleFreq+1,b,p]=self.order_parameter(self.lattice[1])
            self.recordlong[itt,2*time//self.sampleFreq+1,b,p]=2*(bin(self.lattice[0]^self.lattice[1]).count('1'))/self.Length-1
            self.recordMag[itt,2*time//self.sampleFreq+1,b,p]=self.Magnetization()
            self.recordMagS[itt,2*time//self.sampleFreq+1,b,p]=self.StaggeredMagnetization() 
            for mcfreq in repeat(None,self.MonteCarloFreq):
               self.monteCarlo(time,b,p,itt)
        else:
            for controlfreq in repeat(None,self.ControlFreq):
                self.stochasticControl(prob)
            for mcfreq in repeat(None,self.MonteCarloFreq):
                self.monteCarlo(time,b,p,itt)
        return
    
    def Step2D(self,time,prob,b,p,itt):
        if not (time%self.sampleFreq):
            self.recordtot[itt,2*time//self.sampleFreq,b,p] = self.LatticeOrderParameter2d()
            self.recordMag[itt,2*time//self.sampleFreq,b,p]=self.Magnetization()
            self.recordMagS[itt,2*time//self.sampleFreq,b,p]=self.StaggeredMagnetization() 
            for controlfreq in repeat(None,self.ControlFreq):
                self.stochasticControl(prob)
            self.recordtot[itt,2*time//self.sampleFreq+1,b,p] = self.LatticeOrderParameter2d()
            self.recordMag[itt,2*time//self.sampleFreq+1,b,p]=self.Magnetization()
            self.recordMagS[itt,2*time//self.sampleFreq+1,b,p]=self.StaggeredMagnetization() 
            for mcfreq in repeat(None,self.MonteCarloFreq):
               self.monteCarlo(time,b,p,itt)
        else:
            for controlfreq in repeat(None,self.ControlFreq):
                self.stochasticControl(prob)
            for mcfreq in repeat(None,self.MonteCarloFreq):
                self.monteCarlo(time,b,p,itt)
        return
    
    @profile
    def Simulation(self):
        b=0
        for beta in self.betarange:
            p=0
            self.LookupEnergy(self.K,beta)
            for prob in self.probrange:
                
                for itt in range(self.iterations):
                    self.createLattice()
                    for time in range(self.totSteps):
                        self.Step(time,prob,b,p,itt)
                p+=1
            b+=1
        if self.twoChains:
            self.recordtot = (self.record1+self.record2+self.recordlong)/3
        return
        
    def FullyStochasticSimulation(self,pstoch,pctrl):
        '''
        This is for a simulation where you replace the individual steps of the simulation with probabilities to select a certain step. 
        To make this easier to compute, this choice is divided into two steps: MC vs Stochastic -> (if Stochastic) Control vs Bernoulli.
        
        Currently only works for a two-chain system'''
        self.pstoch = pstoch
        self.pctrl = pctrl
        self.acceptance = np.zeros((self.iterations,self.time_size,len(self.betarange),len(self.pstoch),len(self.pctrl))) # keeps track of monte carlo acceptance per step

        self.record1= np.zeros((self.iterations,self.time_size,len(self.betarange),len(self.pstoch),len(pctrl)))#1st chan
        self.record2= np.zeros((self.iterations,self.time_size,len(self.betarange),len(self.pstoch),len(pctrl)))# 2nd chain
        self.recordlong= np.zeros((self.iterations,self.time_size,len(self.betarange),len(self.pstoch),len(pctrl)))  # between chains
        
        self.recordMag =  np.zeros((self.iterations,self.time_size,len(self.betarange),len(self.pstoch),len(pctrl)))
        self.recordMagS =  np.zeros((self.iterations,self.time_size,len(self.betarange),len(self.pstoch),len(pctrl)))
        b=0
        for beta in self.betarange:
            p=0
            self.LookupEnergy(self.K,beta)
            for probS in self.pstoch:
                p2=0
                for probC in self.pctrl:
                    for itt in range(self.iterations):
                        self.createLattice()
                        for time in range(self.totSteps):
                            self.StepStochLadder(time,probS,probC,b,p,p2,itt)
                        self.record1[itt,time//self.sampleFreq,b,p,p2]=self.order_parameter(self.lattice[0])
                        self.record2[itt,time//self.sampleFreq,b,p,p2]=self.order_parameter(self.lattice[1])
                        self.recordlong[itt,time//self.sampleFreq,b,p,p2]=2*(bin(self.lattice[0]^self.lattice[1]).count('1'))/self.Length-1
                        self.recordMag[itt,time//self.sampleFreq,b,p,p2]=self.Magnetization()
                        self.recordMagS[itt,time//self.sampleFreq,b,p,p2]=self.StaggeredMagnetization() 
                    p2+=1
                p+=1
            b+=1
        if self.twoChains:
            self.recordtot = (self.record1+self.record2+self.recordlong)/3
        return
    def StepStochLadder(self,time,probS,probC,b,p1,p2,itt):
        if not (time%self.sampleFreq):
            self.record1[itt,time//self.sampleFreq,b,p1,p2]=self.order_parameter(self.lattice[0])
            self.record2[itt,time//self.sampleFreq,b,p1,p2]=self.order_parameter(self.lattice[1])
            self.recordlong[itt,time//self.sampleFreq,b,p1,p2]=2*(bin(self.lattice[0]^self.lattice[1]).count('1'))/self.Length-1
            self.recordMag[itt,time//self.sampleFreq,b,p1,p2]=self.Magnetization()
            self.recordMagS[itt,time//self.sampleFreq,b,p1,p2]=self.StaggeredMagnetization() 
        if random.random()< probS:
            self.stochasticControl(probC)
        else:
            self.monteCarloLadder2(time,b,p1,p2,itt)
        return
