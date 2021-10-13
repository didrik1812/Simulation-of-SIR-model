import SIR #Out Simulation tools
import numpy as np #Numpy for everything
from scipy.special import stdtrit #Gives critical values in studen t disturbution
import matplotlib.pyplot as plt #For plotting

'''Solve the different tasks in problem 1'''

def task1c(taskFile):
    '''Solve task 1c , writes data to file taskFile'''

    singleSIR = SIR.SIR(T = 7300) #Create SIR object for simulation of single individuals
    N = 30 #Number of individuals
    pi = np.zeros((N, 3)) #Array for collecting approximate pi-values
    T2 = 7300//2 #Half the time 
    
    for i in range(N): #Iterate through all individuals
        singleSIR.simInd(0) #simulate individual starting state in S
        pi[i] = np.array([np.count_nonzero(singleSIR.X[T2:] == 0),np.count_nonzero(singleSIR.X[T2:] == 1),np.count_nonzero(singleSIR.X[T2:] == 2)])/(T2) #Find pi-value for this individual

    CIpi = np.array([singleSIR.calculate95CI(other = pi[:,0],N = N), singleSIR.calculate95CI(other  = pi[:,1],N = N), singleSIR.calculate95CI(other  = pi[:,2],N = N)]) #Calculate confidence interval

    taskFile.write('Confidence intervals for limiting disturbution:\n') #Wrtie data to file
    np.savetxt(taskFile,CIpi,fmt='%1.3f')

def task1e():
    '''Solve task 1e'''

    y0 = np.array([950,50,0]) #Intial states
    multSIR = SIR.SIR() #create SIR object for simulation
    multSIR.simMult(y0) #Simulate over many individuals
    t = np.arange(0,multSIR.T,1) #time array for plotting
    
    '''Realization of S_n, I_n, R_n'''
    plt.figure() #Plotting the devolopment of S_n, I_n, and R_n over time
    plt.step(t,multSIR.X[:,0], label = 'S')
    plt.step(t,multSIR.X[:,1], label = 'I')
    plt.step(t,multSIR.X[:,2], label = 'R')
    plt.ylim((0,1000)) #Make plot look nice
    plt.xlim((0,multSIR.T))
    plt.xlabel('Days')
    plt.ylabel('Number of people')
    plt.legend()
    plt.title('Visualization of SIR-model in population of size N = 1000 over 300 days')
    plt.savefig('problem1djpeg.jpeg') #Save plot
    plt.savefig('problem1dpdf.pdf')

def task1fg(taskFile):
    '''Solve task 1f and 1g and write data to the file taskFile'''
    y0 = np.array([950,50,0]) #Initial state
    multSIR = SIR.SIR()#create SIR object for simulation
    t = np.arange(0,multSIR.T,1) #time array for plotting
    vaxCase = np.array([0,100,600,800]) #Different vacination scenarios
    Ycases = np.array([np.zeros((multSIR.T,3)),np.zeros((multSIR.T,3)),np.zeros((multSIR.T,3)),np.zeros((multSIR.T,3))]) #multdimensional array for collecting the outcome of the different cases 

    i = 0 #Index
    plt.figure() #Plot visualizations of the different cases
    for case in vaxCase:
        multSIR.vax = case #set number of vaccinated to the case we are in
        multSIR.simMult(y0) #simulate over many individuals
        Ycases[i] = multSIR.X #Simulation of SIR with case number i

        plt.step(t,Ycases[i,:,1], label = f'{case} vacinated') #plotting visualization of case number i
        i+=1 #Next case

    plt.legend() #Make the plots look nice
    plt.ylabel('Number of Infected individuals')
    plt.xlabel('Time in days')
    plt.ylim((0))
    plt.xlim((0,300)) 
    plt.title('Visualization of number of infected people in population of N = 1000')
    plt.savefig('problem1gjpeg.jpeg') #Save plots
    plt.savefig('problem1gpdf.pdf')

    '''Calculate confidence intervals for expected maximum and time of expected maximum'''

    N = 1000 #number of simulations
    CIEImax,CIETmax = SIR.excpetedValueCI(multSIR,N,y0,cases = vaxCase) #Get expected values

    taskFile.write('\n\nRUNNING SIR SIMULATION N = 1000 TIMES WITH RESPECTIVELY (0,100,600,800) PEOPLE VACINATED') #Write to file
    taskFile.write('\n\nConfidence intervals for maximum number of infected: \n')
    np.savetxt(taskFile,CIEImax, fmt='%1.3f')

    taskFile.write('\n\nConfidence intervals for time of maximum number of infected: \n')
    np.savetxt(taskFile,CIETmax, fmt='%1.3f')


