import numpy as np
import numpy as np
import emcee
import learndnakinetics
import cPickle as pickle
from copy import deepcopy
import sys
import parent
import learndnakinetics
import myenums

def mcmc_ensemble():
	"""This function runs the MCMC ensemble method with the emcee software!"""
	
	# use 100 walkers, where each walkers takes 50 * 15 = 750 steps. Results are saved every 50 (N_STEPS) steps
	FINAL_STEP =15
	N_WALKERS =100
	N_STEPS =50
	N_INITIAL_WALKERS = 100 * N_WALKERS
	# An initial parameter set might lead to lnprobability of -np.inf ( equivalently probability 0 or error npinf) and the walker might not be able to update to a higher probability during the next iterations.
	# Therefore we initialize N_INITAL_WALKERS walkers, but only  use N_WALKERS walkers which do not have an lnprobability of -np.inf
	
	
	learndnakinetics.METHOD = myenums.MethodName.MCMCNAME.value  # letting the software know to use the MCMC ensemble approach
	learndnakinetics.set_configuration() # setting some configurations
	
	#Initializing the parameters for each walker"""
	if parent.rate_method ==myenums.ModelName.METROPOLISMODELNAME.value :
		# The Metropolis model has 2 parameters """
		n_dim = 2
		max_range_initial = 10 ** 8
		p1 = [np.random.uniform(0, max_range_initial, n_dim) for i in xrange(N_INITIAL_WALKERS)]
	elif parent.rate_method == myenums.ModelName.ARRHENIUSMODELNAME.value :
		#The Arrhenius model has 15 parameters"""
		n_dim = 15
		max_range_initial = 10
		p1= np.zeros((N_INITIAL_WALKERS , n_dim))
		p1A =  [np.random.uniform(0, 9.2103, 14 ) for i in xrange(N_INITIAL_WALKERS)]
		p1E =  [np.random.uniform(0, 6 , 14 ) for i in xrange(N_INITIAL_WALKERS)]
		p1alpha =  [np.random.uniform(0, 10 , 1 ) for i in xrange(N_INITIAL_WALKERS)]
		for i in range( 15):
			for j in range (N_INITIAL_WALKERS) :
				if i == 14:
					p1 [j , i] = p1alpha [j][ 0]
				elif i % 2 ==  0 :
					p1 [ j, i] = p1A [ j ][ i]
				elif i %2 == 1:
					p1 [ j , i] = p1E [j ][  i ]
	else :
		raise ValueError('Error: specify rate_method to be Arrhenius or Metropolis!')
	n_dim +=1
	p3 = [np.random.uniform( 0 , 1, 1) for i in xrange(N_INITIAL_WALKERS)] # Sample the noise parameter  sigma
	p0 = [0 for i in xrange(N_INITIAL_WALKERS)]
	for i in range(N_INITIAL_WALKERS) :
		p0[i] = list( p1[i] )+  list (p3[i])
	pos = []
	
	
	# An initial parameter set might lead to lnprobability of -np.inf ( equivalently probability 0 or error npinf) and the walker might not be able to update to a higher probability during the next iterations.
	# Therefore we initialize N_INITAL_WALKERS walkers, but only  use N_WALKERS walkers which do not have an lnprobability of -np.inf
	for i in range (N_INITIAL_WALKERS) :
		lnprob  =learndnakinetics.objective_function (p0[i])
		if  lnprob != -np.inf:
			pos.append(p0[i])
			print  " lnprob added" +str (lnprob ) + "\n"
		else :
			print "lnprob not added " + str(lnprob) + "\n"
		if len (pos) == N_WALKERS :
			break
	if len(pos)  < N_WALKERS:
		# This means it  was not able to find N_WALKERS which did not have lnprobability -np.inf. Therefore have to increase N_INITIAL_WALKERS
		raise ValueError('Increase N_INITIAL_WALKERS')
	
	
	
	for i in range (0, FINAL_STEP) :
		#Each walker will take FINAL_STEP * N_STEPS iterations. The results are saved every N_STEPS steps """
		sampler = emcee.EnsembleSampler(N_WALKERS, n_dim, learndnakinetics.objective_function)
		pos, prob, state = sampler.run_mcmc(pos, N_STEPS)
		f_out = open(learndnakinetics.parameterFolder + '/' +learndnakinetics.parameterFolder  + 'burn' +str(i)+'.pkl', 'wb')
		pickle.dump(n_dim, f_out)
		pickle.dump(N_WALKERS, f_out)
		pickle.dump(N_STEPS * i, f_out)
		pickle.dump(N_STEPS, f_out)
		pickle.dump(p0, f_out)
		pickle.dump(sampler, f_out)
		pickle.dump(prob,  f_out)
		pickle.dump(pos, f_out)
		pickle.dump(state, f_out)
		f_out.close()
		if i <  FINAL_STEP -1 :
			sampler.reset()
	print("Mean acceptance fraction:", np.mean(sampler.acceptance_fraction))
	print("Autocorrelation time:", sampler.get_autocorr_time())


if __name__ == "__main__":
	mcmc_ensemble()
