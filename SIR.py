import numpy as np #Numpy for everything
from scipy.special import stdtrit #Gives critical values in t-disturbution
import matplotlib.pyplot as plt #For plotting

'''Tools for solving the taks in problem 1'''
 
class SIR:
    '''Using class mostly for fun here'''
    def __init__(self,alpha = 0.005,beta = 0.01,gamma = 0.10, T = 300, vax = 0):
        '''Class for simulation the SIR model
            input:
            alpha : prob of R -> S
            beta: prob of S -> I
            gamma: prob of I ->R
            T: endtime of simulations
            vax: number of people vacinated
        '''
        self.alpha = alpha #Intitialize object
        self.beta = beta
        self.gamma = gamma
        self.T = T
        self.vax = vax

    def simInd(self, X0):
        '''simulateSIR(alpha, beta , gamma, T)
        function which simulates a single individual over time T

        input:
        SIR-object
        T: End time in days
        X0: Starting state
        
        output:
        update the objects X value
        X : array of states which X_i is in to time i
        '''
        P = np.array([[1 - self.beta, self.beta, 0], [0 , 1 - self.gamma, self.gamma], [self.alpha, 0, 1- self.alpha]]) # Transition probabillity matrix
        X = np.zeros(self.T) #Array of state to time i
        X[0] = X0

        for i in range(1,X.size): #Iterate through time
            r = np.random.choice(3, 1, p = P[int(X[i-1])]) #Pick out random state from disturbution based on P
            X[i] = r

        self.X = X #Assign the individual states to the objects state array

    def calculate95CI(self,other = [False],N = 1000):
        '''calculate95CI(X,N)
        Calculate a 95% CI of the mean of a array x by using student t disturbution.

        input:
        N: Number of simulations/tests default value of 1000
        Other: If we want to use something else than self.X

        output:
        CI: 95% confidence interval of the mean of X 
        
        '''
        if other[0] != False:
            meanX = np.mean(other) #Find the mean
            varX = np.var(other) #Find the estimated variance
            t_alphaN = stdtrit(N - 1, 0.025) #Critcal value of student t-disturbution t_alphaN is such that P(T>t_alphaN) = alpha
            CI = np.array([meanX + np.sqrt(varX/N)*t_alphaN, meanX - np.sqrt(varX/N)*t_alphaN])
        
        else:
            meanX = np.mean(self.X) #Find the mean
            varX = np.var(self.X) #Find the estimated variance
            t_alphaN = stdtrit(N - 1, 0.025) #Critcal value of student t-disturbution t_alphaN is such that P(T>t_alphaN) = alpha
            CI = np.array([meanX + np.sqrt(varX/N)*t_alphaN, meanX - np.sqrt(varX/N)*t_alphaN])


        return CI #Return the CIs

    def simMult(self,Y0):
        '''simMult(self,Y0)
        Function for simulation of multiple individuals over time T, where the probabillity of becoming infected
        is dependent of the number of infected people

        input:
        Y0: Starting state on form (S_0, I_0, R_0) where S_0 + I_0 + R_0 = population size

        output:
        updates the objects state array X value to Y
        Y: array with values of (S_n, I_n, R_n) at each time step.
        '''
        N = np.sum(Y0) #Number of individuals
        Y = np.zeros((self.T, 3)) #array of (S_n, I_n, R_n) at each time step given at form Y[time, state]
        Y[0] = Y0
        Y[0,2] += self.vax #Adding to vaccinated to "recoverd" since they have immunity
        Y[0,0] -= self.vax #Subtract from suspictible

        for i in range(1,self.T):
            self.beta = 0.5 * Y[i-1,1]/N #S -> I
            P = np.array([self.alpha, self.beta , self.gamma]) #Probabillity of changing state, i.e "sucsess"
            #Using the binomial disturbution to find how many that changes state
            #Vacinated people cannot lose immunity therefore they do not take part in the chain and must be removed from the sample
            r = np.array([np.random.binomial(Y[i-1,2] - self.vax,P[0]),np.random.binomial(Y[i-1,0] ,P[1]),np.random.binomial(Y[i-1,1] ,P[2])])#New disturbution of people in each state
            Y[i] = np.array([Y[i-1,0] + r[0] - r[1], Y[i-1,1] + r[1] - r[2], Y[i-1,2] + r[2] - r[0]])

        self.X = Y #Update state-array
    
def excpetedValueCI(SIRobj,N,y0,cases = np.array([0])):
    '''excpetedValueCI(SIRobj,N,y0,cases)
    function for calculating 95%CI for expected value of maximum number of infected and time of maximum number of infected

    input:
    SIRobj: an SIR object
    N: Number of simulations
    y0: Starting state
    cases: array of with number of vacinated people, default value of 1dimensial array with only one 0 element in it

    output:
    CIEImax: 95%-CIs of maximum number of infected for each case
    CIETmax: 95%-CIs of time of maxiumum number of infected.
    
    '''
    EImax = np.zeros((N,cases.size)) #Array for collecting of maximum number of infected for each case and for each simulation
    ETmax = np.zeros((N,cases.size)) #Array for collecting time of maxiumum number of infected for each case and for each simulation.
    
    CIEImax = np.zeros((cases.size,2)) #Array for collecting CIs of max number of infected for each case
    CIETmax = np.zeros((cases.size,2)) #Array for collecting CIs of time of max number of infected for each case
    
    for j in range(cases.size): #Iterate through each case
        for i in range(N): #Run N simulations
            SIRobj.vax = cases[j]
            SIRobj.simMult(y0) #Use simulation function
            
            EImax[i,j] = np.max(SIRobj.X[:,1]) #Find maxiumum
            ETmax[i,j] = np.argmax(SIRobj.X[:,1]) #Find time of maxiumum
        
        CIEImax[j] = SIRobj.calculate95CI(other  = EImax[:,j],N = N) #Use CI function to calculate CIs for each case
        CIETmax[j] = SIRobj.calculate95CI(other  = ETmax[:,j],N = N)

    return CIEImax, CIETmax #Return CIs
    

        