from fractions import Fraction
from itertools import repeat
import random
import numpy as np
import os
from line_profiler import profile

class MCCTClassical:
    def __init__(self,Length, numChains,sampleFreq,totSteps, betarange, probrange, iterations, coupling,ControlFreq=1, MonteCarloFreq=1,boltzmannMod=1/2, perLatSweep=False, fullStoch=True, setXJ=False, xj= np.array([Fraction(1,3),Fraction(2,3)]), randomControl=True, magField=0):
        self.Length = Length # how long is each bit string
        self.sampleFreq = sampleFreq #after how many update steps are the parameters calculated
        self.numChains = numChains # how many bit strings are there
        self.ControlFreq=ControlFreq # how many times per step is the stochasitic control performed
        self.iterations = iterations # how many times is the simulation repeated for each set a variables
        self.K=coupling # how strong the nearest neighbor spin interactions are
        self.MonteCarloFreq=MonteCarloFreq # how many times per step is the monte carlo step performed
        self.betarange=betarange # over what thermodynamic betas are we sampling
        self.probrange=probrange # over what probability of selecting the control map (inclusive) are we sampling
        self.allones = (1<<(self.Length))-1
        self.totSteps=totSteps # total number of update steps
        self.neel=np.array((self.dec2int(Fraction(1,3)),self.dec2int(Fraction(2,3))))
        self.xj = xj# fixed points for the control map
        self.xjint =np.array([self.dec2int(i) for i in xj],dtype=object)
        self.twoChains=True if numChains==2 else False
        self.latticeSize=Length*numChains
        self.boltzmannMod=boltzmannMod
        self.fullStoch=fullStoch
        self.Lm=self.Length-1
        self.Lmm = self.Length-2
        self.setXJ=setXJ
        self.randomControl=randomControl
        self.magField=magField
        if perLatSweep:
            self.mcrep = self.latticeSize
        else:
            self.mcrep=1
        # number of data points along the time axis
        if not (self.totSteps%self.sampleFreq):
            self.time_size = (self.totSteps)//self.sampleFreq+1
            
        else: 
            self.time_size = (self.totSteps)//self.sampleFreq+2
        self.setDimensionality() # chooses functions to used based on number of chains, gets around python's lack of overloading
        self.createSamplingArrays() # these arrays store the calculated parameters
        if self.MonteCarloFreq==0:
            betarange=[0]
        if self.ControlFreq==0:
            probrange=[0]    
        
    def dec2int(self, num):
        """
        Conversion of a decimal between 0 and 1 into an integer between 0 and 2**Length-1
        :type num: Fraction
        """
        assert 0 <= num < 1, self.Length > 0
        return int(num * (2 ** self.Length))


    def reset(self,bin_num):
        """
        odd numbers end with b_L =1
        num - 1 for set b_L = 0 if odd
        """
        #if bin_num % 2 == 1:
        #    bin_num = bin_num - 1
        return bin_num&(~1)


    def left(self,num):
        """

  
        """
        #if num >= 2 << (length - 2):
        #    num = (num << 1) ^ 1
        #else:
        #    num = num << 1
        return (num>>self.Length-1)^(num<<1)&((1<<self.Length)-1)#num % (2 ** (length))


    def right(self,num):
        """

        """
        num = self.reset(num)
        # if num % 2 ==1 and (num//2) % 2 == 0 and (num//4)%2 == 1:
        #    num = num + 2
        num = num // 2
        return num

    
    def order_parameter(self,num):
        '''
        Calculates sigma_i sigma_j in an integer represenation of a bitstring.'''
        binrep =bin(self.left(num)^num) # \sigma_{i+1}*\sigma_{i}
        count = 2*binrep.count('1')/self.Length-1

        return count

    
    def bernoulli(self,num):
        num = self.left(num)
        num = num & (~0b111)
        return num + random.getrandbits(3)

    def controlSet(self,num,xj):
        if num == xj:
            return num
        else:
            temp = self.right(num) + self.right(xj)
        if temp == num:
            return num + 1
        return temp
    def controlNorm(self,num,xj):

        if num < 2 ** (self.Length-1):
            xj = self.xjint[0]
        else:
            xj = self.xjint[1]
        if num == xj:
            return num
        else:
            temp = self.right(num) + self.right(xj)
        if temp == num:
            return num + 1
        return temp
    
    def boltzmann_probability(self,initial, final, beta):
        '''
        Calculates a value from a modified Boltzmann distribution. The choice of 1/2 (which modifies the distribution) is to create bit scrambling when Beta=0 '''
        return np.exp(-beta * (final - initial))*self.boltzmannMod
        #boltzmann weight. we take 1/2 such that the state is continuously randomized at infinite
    
    def BitConfiguration1d(self,position_x,position_y):
    #problem: individual calculations are incredibly inefficient
    #solution: Make lookup tables for everything
    #so, from the bits, we can just determine the energy without additional calculation.   
        '''        
        Returns the targeted bit and the nearest neighbor bits (up,down, left) in an integer.'''
        if position_y==self.Lm: #select bits centered on position L-1 on the bitstring x
           chainbits = ((self.lattice[position_x]&1)<<2)+(((self.lattice[position_x]&((3<<(self.Lmm))))>>(self.Lmm)))
        elif position_y: 
           chainbits = (self.lattice[position_x]&(7<<(position_y-1)))>>(position_y-1)
        else: #select bits centered on position 0 on the bitstring x
            #dbit = (self.lattice[position_x]&(1<<self.Length))>>self.Length
            #twobits= (self.lattice[position_x]&3)>>1
            chainbits = ((self.lattice[position_x]&3)<<1)+(((self.lattice[position_x])&(1<<self.Lm))>>(self.Lm))
        #mainBit =(self.lattice[position_x]&(1<<position_y))>>(position_y)
        #ubit = (self.lattice[position_x]&(1 << (position_y-1)%self.Length))>>(position_y-1)%self.Length
        #dbit= ((self.lattice[position_x]&(1<<(position_y+1)%self.Length))>>(position_y+1)%self.Length)
        return chainbits

    def BitConfigurationLad(self,position_x,position_y):
    #problem: individual calculations are incredibly inefficient
    #solution: Make lookup tables for everything
    #so, from the bits, we can just determine the energy without additional calculation.   
        '''        
        Returns the targeted bit and the nearest neighbor bits (up,down, left) in an integer.'''
        if position_y==self.Lm: #select bits centered on position L-1 on the bitstring x
           chainbits = ((self.lattice[position_x]&1)<<2)+(((self.lattice[position_x]&((3<<(self.Lmm))))>>(self.Lmm)))
        elif position_y: 
           chainbits = (self.lattice[position_x]&(7<<(position_y-1)))>>(position_y-1)
        else: #select bits centered on position 0 on the bitstring x
            #dbit = (self.lattice[position_x]&(1<<self.Length))>>self.Length
            #twobits= (self.lattice[position_x]&3)>>1
            chainbits = ((self.lattice[position_x]&3)<<1)+(((self.lattice[position_x])&(1<<self.Lm))>>(self.Lm))
        #mainBit =(self.lattice[position_x]&(1<<position_y))>>(position_y)
        #ubit = (self.lattice[position_x]&(1 << (position_y-1)%self.Length))>>(position_y-1)%self.Length
        #dbit= ((self.lattice[position_x]&(1<<(position_y+1)%self.Length))>>(position_y+1)%self.Length)
        lbit =((self.lattice[position_x-1]&(1<<position_y))>>(position_y))
        #return mainBit,ubit,dbit, lbit
        return chainbits,lbit
    def BitConfiguration2d(self,position_x,position_y):
    #problem: individual calculations are incredibly inefficient
    #solution: Make lookup tables for everything
    #so, from the bits, we can just determine the energy without additional calculation.   
        '''        
        Returns the targeted bit and the nearest neighbor bits (up,down, left, right) in an integer.'''
       
        if position_y==self.Lm: #select bits centered on position L-1 on the bitstring x
           chainbits = ((self.lattice[position_x]&1)<<2)+(((self.lattice[position_x]&((3<<(self.Lmm))))>>(self.Lmm)))
        elif position_y: 
           chainbits = (self.lattice[position_x]&(7<<(position_y-1)))>>(position_y-1)
        else: #select bits centered on position 0 on the bitstring x
            #dbit = (self.lattice[position_x]&(1<<self.Length))>>self.Length
            #twobits= (self.lattice[position_x]&3)>>1
            chainbits = ((self.lattice[position_x]&3)<<1)+(((self.lattice[position_x])&(1<<self.Lm))>>(self.Lm))


        #mainBit =(self.lattice[position_x]&(1<<position_y))>>(position_y)
        #ubit = (self.lattice[position_x]&(1 << (position_y-1)%self.Length))>>(position_y-1)%self.Length
        #dbit= ((self.lattice[position_x]&(1<<(position_y+1)%self.Length))>>(position_y+1)%self.Length)
        lbit =((self.lattice[position_x-1]&(1<<position_y))>>(position_y))
        rbit=(((self.lattice[(1+position_x)%self.numChains])&(1<<position_y))>>position_y)
        #return mainBit,ubit,dbit,lbit,rbit
        return chainbits,lbit,rbit
    # The lookup dictionary is designed around taking in as many bits as there are nearest neighbors.
    #  Thus this imp has a bunch of different tables for each dimensionality
    
    def LookupEnergy1d(self,K,Beta):
        self.energyDict={} #more efficient imp has not been verified explicitly
        self.energyDict[0b000] =min(self.boltzmannMod,self.boltzmann_probability(2*K,-2*K,Beta))
        self.energyDict[0b111]=self.energyDict[0b000]

        self.energyDict[0b011]=min(self.boltzmannMod,self.boltzmann_probability(0,0,Beta))
        self.energyDict[0b001]=self.energyDict[0b011]
        self.energyDict[0b100]=self.energyDict[0b011]
        self.energyDict[0b110]=self.energyDict[0b011]

        self.energyDict[0b010]=min(self.boltzmannMod,self.boltzmann_probability(-2*K,2*K,Beta))
        self.energyDict[0b101]=self.energyDict[0b010]
        #
        #self.energyDict[0,0,0] =min(self.boltzmannMod,self.boltzmann_probability(2*K,-2*K,Beta))
        #self.energyDict[1,1,1]=self.energyDict[0,0,0]

        #self.energyDict[0,1,0]=min(self.boltzmannMod,self.boltzmann_probability(0,0,Beta))
        #self.energyDict[0,0,1]=self.energyDict[0,1,0]
        #self.energyDict[1,0,1]=self.energyDict[0,1,0]
        #self.energyDict[1,1,0]=self.energyDict[0,1,0]

        #self.energyDict[0,1,1]=min(self.boltzmannMod,self.boltzmann_probability(-2*K,2*K,Beta))
        #self.energyDict[1,0,0]=self.energyDict[0,1,1]
        return self.energyDict
    
    def LookupEnergyLad2d(self,K,Beta):
        self.energyDict = {}
        # One diff, up
        self.energyDict[0b000,1] =  min(self.boltzmannMod,self.boltzmann_probability(K-self.magField,-K+self.magField,Beta))
        self.energyDict[0b001,0]=self.energyDict[0b000,1]
        self.energyDict[0b100,0]=self.energyDict[0b000,1]
        # One same, up
        self.energyDict[0b010,1]=min(self.boltzmannMod, self.boltzmann_probability(-K+self.magField,K-self.magField,Beta))
        self.energyDict[0b011,0]=self.energyDict[0b010,1]
        self.energyDict[0b110,0]=self.energyDict[0b010,1]
        # All/None Same
        self.energyDict[0b000,0]=min(self.boltzmannMod,self.boltzmann_probability(K*3-self.magField,K*-3+self.magField,Beta))
        self.energyDict[0b010,0]=min(self.boltzmannMod,self.boltzmann_probability(K*-3+self.magField,K*3-self.magField,Beta))
        self.energyDict[0b101,1]=min(self.boltzmannMod, self.boltzmann_probability(-3*K-self.magField,K*3+self.magField,Beta))#8
        self.energyDict[0b111,1]= min(self.boltzmannMod,self.boltzmann_probability(K*3+self.magField,-K*3-self.magField,Beta)) 
        # One Diff, down
        self.energyDict[0b011,1]=min(self.boltzmannMod,self.boltzmann_probability(K+self.magField,-K-self.magField,Beta))
        self.energyDict[0b110,1]=self.energyDict[0b011,1] 
        self.energyDict[0b111,0]=self.energyDict[0b011,1] 
        # One Same, down
        self.energyDict[0b101,0]=min(self.boltzmannMod,self.boltzmann_probability(K*-3-self.magField,K*3+self.magField,Beta))
        self.energyDict[0b001,1]=self.energyDict[0b101,0]
        self.energyDict[0b100,1]=self.energyDict[0b101,0]
        #
        #self.energyDict[0,0,0,1] =  min(self.boltzmannMod,self.boltzmann_probability(K,-K,Beta))
        #self.energyDict[0,0,1,0]=self.energyDict[0,0,0,1]
        #self.energyDict[0,1,0,0]=self.energyDict[0,0,0,1]

        #self.energyDict[1,0,0,1]=min(self.boltzmannMod, self.boltzmann_probability(-K,K,Beta))
        #self.energyDict[1,0,1,0]=self.energyDict[1,0,0,1]
        #self.energyDict[1,1,0,0]=self.energyDict[1,0,0,1]

        #self.energyDict[0,0,0,0]=min(self.boltzmannMod,self.boltzmann_probability(K*3,K*-3,Beta))
        #self.energyDict[1,0,0,0]=min(self.boltzmannMod,self.boltzmann_probability(K*-3,K*3,Beta))
        #self.energyDict[0,1,1,1]=self.energyDict[1,0,0,0]#8
        #self.energyDict[1,1,1,1]=self.energyDict[0,0,0,0] 

        #self.energyDict[1,0,1,1]=self.energyDict[0,0,0,1] 
        #self.energyDict[1,1,0,1]=self.energyDict[0,1,0,0] 
        #self.energyDict[1,1,1,0]=self.energyDict[0,1,0,0] 

        #self.energyDict[0,1,1,0]=self.energyDict[1,0,0,1]
        #self.energyDict[0,0,1,1]=self.energyDict[1,0,0,1]
        #self.energyDict[0,1,0,1]=self.energyDict[1,0,0,1]
        return self.energyDict
    
    def LookupEnergy2d(self,K,Beta):
        self.energyDict = {}
        # This can also be a loop without significant performance decrease. Just make sure that, when sampling over parameters that Beta changes, and thus 
        # updates to this energy dictionary, are in the outer most loop
        #2**5=32
        self.energyDict[0b111,1,1] = min(self.boltzmannMod,self.boltzmann_probability(K*4+self.magField,K*-4-self.magField,Beta))
        self.energyDict[0b000,0,0] = min(self.boltzmannMod,self.boltzmann_probability(K*4-self.magField,K*-4+self.magField,Beta))
        #2 All same
        self.energyDict[0b101,1,1] = min(self.boltzmannMod,self.boltzmann_probability(K*(-4)-self.magField,K*4+self.magField,Beta))
        self.energyDict[0b010,0,0] =  min(self.boltzmannMod,self.boltzmann_probability(K*(-4)+self.magField,K*4-self.magField,Beta))
        #2 All diff
        self.energyDict[0b001,1,1] = min(self.boltzmannMod,self.boltzmann_probability(K*-2-self.magField,K*2+self.magField,Beta))
        self.energyDict[0b100,1,1] = self.energyDict[0b001,1,1]
        self.energyDict[0b101,0,1] = self.energyDict[0b001,1,1]
        self.energyDict[0b101,1,0] = self.energyDict[0b001,1,1]
        #4 One same
        self.energyDict[0b110,0,0] = min(self.boltzmannMod,self.boltzmann_probability(K*-2+self.magField,K*2-self.magField,Beta))
        self.energyDict[0b011,0,0] = self.energyDict[0b110,0,0]
        self.energyDict[0b010,1,0] = self.energyDict[0b110,0,0]
        self.energyDict[0b010,0,1] = self.energyDict[0b110,0,0]
        #4 One Same
        self.energyDict[0b011,1,1] = min(self.boltzmannMod,self.boltzmann_probability(K*2+self.magField,K*-2-self.magField,Beta))
        self.energyDict[0b110,1,1] = self.energyDict[0b011,1,1]
        self.energyDict[0b111,0,1] = self.energyDict[0b011,1,1]
        self.energyDict[0b111,1,0] = self.energyDict[0b011,1,1]
        #4 One diff
        self.energyDict[0b100,0,0] = min(self.boltzmannMod,self.boltzmann_probability(K*2-self.magField,K*-2+self.magField,Beta))
        self.energyDict[0b001,0,0] =  self.energyDict[0b100,0,0]
        self.energyDict[0b000,1,0] =  self.energyDict[0b100,0,0]
        self.energyDict[0b000,0,1] =  self.energyDict[0b100,0,0]
        #4 One diff
        self.energyDict[0b010,1,1] = min(self.boltzmannMod,self.boltzmann_probability(self.magField,-self.magField,Beta))
        self.energyDict[0b101,0,0] = self.energyDict[0b010,1,1]
        self.energyDict[0b100,1,0] =  min(self.boltzmannMod,self.boltzmann_probability(-self.magField,self.magField,Beta))
        self.energyDict[0b100,0,1] = self.energyDict[0b100,1,0]
        self.energyDict[0b001,1,0] = self.energyDict[0b100,1,0]
        self.energyDict[0b001,0,1] = self.energyDict[0b100,1,0]
        #6 Eq
        self.energyDict[0b111,0,0] = self.energyDict[0b010,1,1]
        self.energyDict[0b110,1,0] = self.energyDict[0b010,1,1]
        self.energyDict[0b110,0,1] = self.energyDict[0b010,1,1]
        self.energyDict[0b011,1,0] = self.energyDict[0b010,1,1]
        self.energyDict[0b011,0,1] = self.energyDict[0b010,1,1]
        self.energyDict[0b000,1,1] = self.energyDict[0b100,1,0]
       
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
    def Magnetization(self):
        param =0 
        for i in self.lattice:
            param+= bin(i).count('1')
        param = 1- 2*param/self.latticeSize #0 is positive, 1 is negative
        return param
#works as intended
    
    def StaggeredMagnetization(self):
        param =0
        for i in range(0,self.numChains,2):
            param+= bin(self.lattice[i]^self.neel[1]).count('1')
        for i in range(1, self.numChains,2):
            param += bin(self.lattice[i]^self.neel[0]).count('1')
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
            self.Measure=self.Measure1d
            if self.fullStoch:
                self.monteCarlo=self.monteCarlo1ds
                self.Simulation=self.FullyStochasticSimulation
                if self.randomControl:
                    self.stochasticControl=self.stochasticControlRand 
        elif self.twoChains:
            self.LookupEnergy = self.LookupEnergyLad2d
            self.BitConfiguration = self.BitConfigurationLad
            self.stochasticControl=self.stochasticControl2d
            self.monteCarlo=self.monteCarloLadder
            self.LatticeOrderParameter=self.LatticeOrderParameterLad
            self.Step = self.StepLad
            self.Measure=self.MeasureLadder

            if self.fullStoch:
                self.monteCarlo=self.monteCarloLadder2
                self.Simulation=self.FullyStochasticSimulation
                if self.randomControl:
                    self.stochasticControl=self.stochasticControlRand 
        else:
            self.LookupEnergy=self.LookupEnergy2d
            self.BitConfiguration = self.BitConfiguration2d
            self.stochasticControl=self.stochasticControl2d
            self.monteCarlo=self.monteCarlo2d
            self.LatticeOrderParameter=self.LatticeOrderParameter2d
            self.Step = self.Step2D
            self.Measure=self.Measure2d

            if self.fullStoch:
                self.monteCarlo=self.monteCarlo2ds
                self.Simulation=self.FullyStochasticSimulation
                if self.randomControl:
                    self.stochasticControl=self.stochasticControlRand 
        if self.setXJ and self.numChains!=1:
            self.control=self.controlSet
        else:
            self.control=self.controlNorm
        return
    
    def createLattice(self):
        lat = []
        for i in repeat(None,self.numChains):
            lat.append(random.getrandbits(self.Length))
        self.lattice=np.array(lat,dtype=object)
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
        self.recordChain= np.zeros((self.iterations,self.time_size,len(self.betarange),len(self.probrange))) # over all chains
        self.recordMag =  np.zeros((self.iterations,self.time_size,len(self.betarange),len(self.probrange)))
        self.recordMagS =  np.zeros((self.iterations,self.time_size,len(self.betarange),len(self.probrange)))
        self.time_size=self.time_size//2
        return
    
    def stochasticControl2d(self,prob):
        for i in range(self.numChains):
            if random.random()>(prob):
                self.lattice[i]=self.bernoulli(self.lattice[i])
            else:
                self.lattice[i]=self.control(self.lattice[i],self.xjint[i%2])
        return
    def stochasticControlRand(self,prob):
        i = random.randrange(0,self.numChains)
        if random.random()>prob:   
            self.lattice[i]=self.bernoulli(self.lattice[i])
        else:
            self.lattice[i]=self.control(self.lattice[i],self.xjint[i%2])
    # the correct phase transition is maintained in 1d
    
    def stochasticControl1d(self,prob):
        if random.random()>(prob):
            self.lattice[0]=self.bernoulli(self.lattice[0])
        else:
            self.lattice[0]=self.control(self.lattice[0],None) #none is just a placeholder for convenience
        return
    
    def monteCarlo1d(self,time,betaNum,probNum,itt):
        if not (time%self.sampleFreq):
            for nr in repeat(None,self.mcrep):
                y_pos = random.randrange(self.Length)
                chainbits = self.BitConfiguration(0,y_pos)
                if self.energyDict[chainbits]>random.random(): 
                    self.acceptance[itt,time//self.sampleFreq,betaNum,probNum] += 1
                    self.lattice[0] = self.lattice[0]^(0b1<<y_pos)

        else:
            for nr in repeat(None,self.mcrep):
                y_pos = random.randrange(self.Length)
                chainbits= self.BitConfiguration(0, y_pos)
                if self.energyDict[chainbits]>random.random(): 
                    self.lattice[0] = self.lattice[0]^(0b1<<y_pos)

        return
    def monteCarlo1ds(self,time,betaNum,probNum,probNum2,itt):
        #if not (time%self.sampleFreq):
        #    for nr in repeat(None,self.mcrep):
        #        y_pos = random.randrange(self.Length)
        #        chainbits = self.BitConfiguration(0,y_pos)
        #        if self.energyDict[chainbits]>random.random(): 
        #            self.acceptance[itt,time//self.sampleFreq,betaNum,probNum,probNum2] += 1
        #            self.lattice[0] = self.lattice[0]^(0b1<<y_pos)

        #else:
        for nr in repeat(None,self.mcrep):
            y_pos = random.randrange(self.Length)
            chainbits= self.BitConfiguration(0, y_pos)
            if self.energyDict[chainbits]>random.random(): 
                self.lattice[0] = self.lattice[0]^(0b1<<y_pos)

        return
    def monteCarloLadder(self,time,betaNum,probNum,itt):
        if not (time%self.sampleFreq):
            for nr in repeat(None,self.mcrep):
                x_pos =  random.randrange(self.numChains)
                y_pos =  random.randrange(self.Length)

                chainbits,lbit = self.BitConfiguration(x_pos,y_pos)
                #cbits,lbit = self.BitConfiguration(x_pos,y_pos)

                if self.energyDict[chainbits,lbit]>random.random(): 
                    self.acceptance[itt,time//self.sampleFreq,betaNum,probNum] += 1
                    self.lattice[x_pos] = self.lattice[x_pos]^(0b1<<y_pos)

        else:
            for nr in repeat(None,self.mcrep):
                x_pos =  random.randrange(self.numChains)
                y_pos = random.randrange(self.Length)
 
                chainbits,lbit = self.BitConfiguration(x_pos,y_pos)
                #cbits,lbit = self.BitConfiguration(x_pos,y_pos)

                if self.energyDict[chainbits,lbit]>random.random(): 
                    self.lattice[x_pos] = self.lattice[x_pos]^(0b1<<y_pos)

        return   
    def monteCarloLadder2(self,time,betaNum,probNum1,probNum2,itt):
        #if not (time%self.sampleFreq): # can do massive speedup here by not calculating acceptance ratio
        #    for nr in repeat(None,self.mcrep):
        #        x_pos =  random.randrange(self.numChains)
        #        y_pos =  random.randrange(self.Length)

         #       chainbits,lbit = self.BitConfiguration(x_pos,y_pos)
                #cbits,lbit = self.BitConfiguration(x_pos,y_pos)

        #        if self.energyDict[chainbits,lbit]>random.random(): 
        #            self.acceptance[itt,time//self.sampleFreq,betaNum,probNum1,probNum2] += 1
        #            self.lattice[x_pos] = self.lattice[x_pos]^(0b1<<y_pos)

        #else:
        for nr in repeat(None,self.mcrep):
            x_pos =  random.randrange(self.numChains)
            y_pos = random.randrange(self.Length)

            chainbits,lbit = self.BitConfiguration(x_pos,y_pos)
            #cbits,lbit = self.BitConfiguration(x_pos,y_pos)

            if self.energyDict[chainbits,lbit]>random.random(): 
                self.lattice[x_pos] = self.lattice[x_pos]^(0b1<<y_pos)

        return             
    
    def monteCarlo2d(self,time,betaNum,probNum,itt):
        #if not (time%self.mcrep):
        #    for nr in repeat(None,self.mcrep):
        #        x_pos =  random.randrange(self.numChains)
        #        y_pos =  random.randrange(self.Length)

        #        chainbits,lbit,rbit = self.BitConfiguration(x_pos,y_pos)
                #cbits,lbit,rbit = self.BitConfiguration(x_pos,y_pos)
        #        if self.energyDict[chainbits,lbit,rbit]>random.random(): 
        #            self.acceptance[itt,time//self.sampleFreq,betaNum,probNum] += 1
        #            self.lattice[x_pos] = self.lattice[x_pos]^(0b1<<y_pos)

        #else:
        for nr in repeat(None,self.mcrep):
            x_pos =  random.randrange(self.numChains)
            y_pos =  random.randrange(self.Length)

            chainbits,lbit,rbit = self.BitConfiguration(x_pos,y_pos)
            #cbits,lbit,rbit = self.BitConfiguration(x_pos,y_pos)

            if self.energyDict[chainbits,lbit,rbit]>random.random(): 
                self.lattice[x_pos] = self.lattice[x_pos]^(0b1<<y_pos)
        return
    def monteCarlo2ds(self,time,betaNum,probNum,probNum2,itt):
        #if not (time%self.mcrep):
        #    for nr in repeat(None,self.mcrep):
        #        x_pos =  random.randrange(self.numChains)
        #        y_pos =  random.randrange(self.Length)

        #        chainbits,lbit,rbit = self.BitConfiguration(x_pos,y_pos)
                #cbits,lbit,rbit = self.BitConfiguration(x_pos,y_pos)
        #        if self.energyDict[chainbits,lbit,rbit]>random.random(): 
        #            self.acceptance[itt,time//self.sampleFreq,betaNum,probNum,probNum2] += 1
        #            self.lattice[x_pos] = self.lattice[x_pos]^(0b1<<y_pos)

        #else:
        
        for nr in repeat(None,self.mcrep):
            x_pos =  random.randrange(self.numChains)
            y_pos =  random.randrange(self.Length)

            chainbits,lbit,rbit = self.BitConfiguration(x_pos,y_pos)

            if self.energyDict[chainbits,lbit,rbit]>random.random(): 
                self.lattice[x_pos] = self.lattice[x_pos]^(0b1<<y_pos)
        return
    def monteCarlo2dns(self,time,betaNum,probNum,probNum2,itt):
     
        
        for nr in repeat(None,self.mcrep):
            x_pos =  random.randrange(self.numChains)
            y_pos =  random.randrange(self.Length)

            chainbits,lbit,rbit = self.BitConfiguration(x_pos,y_pos)

            if self.energyDict[chainbits,lbit,rbit]>random.random(): 
                self.lattice[x_pos] = self.lattice[x_pos]^(0b1<<y_pos)
                self.acceptance[itt,time,betaNum,probNum,probNum2] += 1
            else:
                self.stochasticControl(self.pctrl[probNum])

        return
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
        '''
        self.pstoch = pstoch
        self.pctrl = pctrl
        self.acceptance = np.zeros((self.iterations,self.time_size,len(self.betarange),len(self.pstoch),len(self.pctrl))) # keeps track of monte carlo acceptance per step

        self.record1= np.zeros((self.iterations,self.time_size,len(self.betarange),len(self.pstoch),len(pctrl)))#1st chan
        self.record2= np.zeros((self.iterations,self.time_size,len(self.betarange),len(self.pstoch),len(pctrl)))# 2nd chain
        self.recordlong= np.zeros((self.iterations,self.time_size,len(self.betarange),len(self.pstoch),len(pctrl)))  # between chains
        self.recordtot=np.zeros((self.iterations,self.time_size,len(self.betarange),len(self.pstoch),len(pctrl))) 
        self.recordChain=np.zeros((self.iterations,self.time_size,len(self.betarange),len(self.pstoch),len(pctrl))) 
        self.recordMag =  np.zeros((self.iterations,self.time_size,len(self.betarange),len(self.pstoch),len(pctrl)))
        self.recordMagS =  np.zeros((self.iterations,self.time_size,len(self.betarange),len(self.pstoch),len(pctrl)))
        self.finalState=np.zeros((self.iterations,len(self.betarange),len(self.pstoch),len(pctrl)),dtype=object)

        b=0
        for beta in self.betarange:
            p=0
            self.LookupEnergy(self.K,beta)
            for probS in self.pstoch:
                p2=0
                for probC in self.pctrl:
                    for itt in range(self.iterations):
                        self.createLattice()
                       
                        for time in range(self.time_size):
                            if probS==1:
                                self.Measure(time,b,p,p2,itt)
                                for i in repeat(None,self.sampleFreq):
                                    self.stochasticControl(probC)
                            elif probS==0:
                                self.Measure(time,b,p,p2,itt)
                                for i in repeat(None,self.sampleFreq):
                                    self.monteCarlo(time,b,p,p2,itt)
                            else:
                                self.StepStoch(time,probS,probC,b,p,p2,itt)
                        self.finalState[itt,b,p,p2]=self.lattice

                        self.record1[itt,-1,b,p,p2]=self.order_parameter(self.lattice[0])
                        self.record2[itt,-1,b,p,p2]=self.order_parameter(self.lattice[1])
                        self.recordlong[itt,-1,b,p,p2]=2*(bin(self.lattice[0]^self.lattice[1]).count('1'))/self.Length-1
                        self.recordMag[itt,-1,b,p,p2]=self.Magnetization()
                        self.recordMagS[itt,-1,b,p,p2]=self.StaggeredMagnetization() 
                    p2+=1
                p+=1
            b+=1
        if self.twoChains:
            self.recordtot = (self.record1+self.record2+self.recordlong)/3
        return
    def StepStoch(self,time,probS,probC,b,p1,p2,itt):
        #if not (time%self.sampleFreq):
        self.Measure(time,b,p1,p2,itt)
        for i in repeat(None,self.sampleFreq):
            if random.random()< probS:
                self.stochasticControl(probC)
            else:
                self.monteCarlo(time,b,p1,p2,itt)
        return
    
    def SimBetaLink(self,pctrl):
        self.pstoch=[0]
        self.pctrl = pctrl
        self.acceptance = np.zeros((self.iterations,self.time_size,len(self.betarange),len(self.pstoch),len(self.pctrl))) # keeps track of monte carlo acceptance per step

        self.record1= np.zeros((self.iterations,self.time_size,len(self.betarange),len(self.pstoch),len(pctrl)))#1st chan
        self.record2= np.zeros((self.iterations,self.time_size,len(self.betarange),len(self.pstoch),len(pctrl)))# 2nd chain
        self.recordlong= np.zeros((self.iterations,self.time_size,len(self.betarange),len(self.pstoch),len(pctrl)))  # between chains
        self.recordtot=np.zeros((self.iterations,self.time_size,len(self.betarange),len(self.pstoch),len(pctrl))) 
        self.recordChain=np.zeros((self.iterations,self.time_size,len(self.betarange),len(self.pstoch),len(pctrl))) 
        self.recordMag =  np.zeros((self.iterations,self.time_size,len(self.betarange),len(self.pstoch),len(pctrl)))
        self.recordMagS =  np.zeros((self.iterations,self.time_size,len(self.betarange),len(self.pstoch),len(pctrl)))
        self.finalState=np.zeros((self.iterations,len(self.betarange),len(self.pstoch),len(pctrl)),dtype=object)
        p=0
        b=0
        for beta in self.betarange:
            p2=0
            self.LookupEnergy(self.K,beta)
            for probC in self.pctrl:
                for itt in range(self.iterations):
                    self.createLattice()
                    
                    for time in range(self.time_size):
                    
                        self.Measure(time,b,p,p2,itt)
                        for i in repeat(None,self.sampleFreq):
                            self.monteCarlo2dns(time,b,p,p2,itt)
                        self.acceptance[itt,time,b,p,p2] =self.acceptance[itt,time,b,p,p2]/(self.mcrep*self.sampleFreq)
                    self.finalState[itt,b,p,p2]=self.lattice
                    self.Measure(-1,b,p,p2,itt)
                p2+=1
            
            b+=1

    def Measure1d(self,time,b,p1,p2,itt):
        self.record1[itt,time,b,p1,p2]=self.order_parameter(self.lattice[0])
        self.recordMag[itt,time,b,p1,p2]=self.Magnetization()
        self.recordMagS[itt,time,b,p1,p2]=self.StaggeredMagnetization() 
        return
    def MeasureLadder(self,time,b,p1,p2,itt):
        self.record1[itt,time,b,p1,p2]=self.order_parameter(self.lattice[0])
        self.record2[itt,time,b,p1,p2]=self.order_parameter(self.lattice[1])
        self.recordlong[itt,time,b,p1,p2]=2*(bin(self.lattice[0]^self.lattice[1]).count('1'))/self.Length-1

        self.recordMag[itt,time,b,p1,p2]=self.Magnetization()
        self.recordMagS[itt,time,b,p1,p2]=self.StaggeredMagnetization() 
        return
    def Measure2d(self,time,b,p1,p2,itt):
        self.record1[itt,time,b,p1,p2]=self.order_parameter(self.lattice[0]) # just to see if this equilibriates

        self.recordtot[itt,time,b,p1,p2]=self.LatticeOrderParameter()
        self.recordMag[itt,time,b,p1,p2]=self.Magnetization()
        self.recordMagS[itt,time,b,p1,p2]=self.StaggeredMagnetization() 
        param = 0
        for i in self.lattice:
            param += (bin(self.left(i,)^i).count('1')) # \sigma_(i,j)*\sigma_(i,j+1)
        self.recordChain[itt,time,b,p1,p2] = (2*param/(self.latticeSize)-1)
        return
