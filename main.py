import SIR #Our simulation tools
import numpy as np #numpy for everything
from scipy.special import stdtrit #gives critical values in student t disturbution
import matplotlib.pyplot as plt #For plotting
from tasks import * #functions which solves the individual tasks

def main():
    '''Solving problem 1 by simulation of the SIR-model'''

    print('Starting simulation\n')
    taskFile = open('problem1.txt','w') #file which date is to be saved in
    taskFile.write("DATA FROM PROJECT 1\n\nSIMULATION OF SINGLE INDIVIDUAL N = 30 TIMES\n")
    '''Single individual'''
    task1c(taskFile)

    '''Introducing multiple people and dependence of number of infected'''
    task1e()

    '''Introducing vacination'''
    task1fg(taskFile)

    taskFile.close()
    print('Simulation complete')
    return 0
main()