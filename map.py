from __future__ import division
import ConfigParser
import numpy as np
from scipy.optimize import minimize
import learndnakinetics
import myenums
import sys
import os
sys.path.insert(0,os.path.realpath('../reactions'))
import parent

def map ( ):
	""" This function runs the MAP approach with  the Nelder-Mead optimization technique!"""
	learndnakinetics.METHOD = myenums.MethodName.MAPNAME.value # letting the software know to use the MAP approach
	learndnakinetics.set_configuration() #setting some configuations
	
	#Initializing the parameters for Nelder-Mead
	if parent.rate_method == myenums.ModelName.ARRHENIUSMODELNAME.value :
		theta = [  13.0580, 3,  13.0580, 3,   13.0580, 3,  13.0580 , 3,   13.0580, 3,  13.0580,  3,   13.0580 , 3,    0.0402 ]
	elif parent.rate_method == myenums.ModelName.METROPOLISMODELNAME.value:
		theta = [8.2 *  (10 **6), 3.3  * (10**5) ]
	else:
		raise ValueError('Error: Please specify rate_method to be Arrhenius or Metropolis!')
	theta = theta + [1]
	thetaOptmized = minimize(learndnakinetics.objective_function, theta, method='nelder-mead',options={ 'disp': True})
	print(thetaOptmized.x)

if __name__ == "__main__":
	map()
