from __future__ import division
import parent
from parent import *

class HairpinComplex(ParentComplex):
	"""Hairpin opening or closing reaction"""
	def __init__(self,  myPickles, theta,beta, gamma , betaprime,  zip, hairpin, m , L ,  T, concentration, sodium, magnesium,  dataset_name, docID ):
		ParentComplex.__init__(self , myPickles, theta ,  T, concentration, sodium, magnesium,  dataset_name, docID )
		self.beta = beta
		self.gamma = gamma
		self.betaprime= betaprime
		self.hairpin = hairpin
		self.m = m
		self.L = L
		self.zip = zip
	

	def dot_paren(self, state):
		"""Returns the structure of the complex in dot-paren notation"""
		m, L = self.m, self.L
		i, j = state
		dotPar = '.' * i + '('* (j-i) + '.' * (m - j ) + '.' *  (L - 2*m)  + '.' *( m - j ) + ')'* (j-i) + '.' * i
		return dotPar
	

	def dot_paren_modify(self, state):
		"""Insert 	`*' at the start and end of the dot-paren notation, before and after all `+' signs, and also before and after every space"""
		m, L = self.m, self.L
		i, j = state
		dotPar = '*'+ '.' * i + '('* (j-i) + '.' * (m - j ) + '.' *  (L - 2*m)  + '.' *( m - j ) + ')'* (j-i) + '.' * i   + '*'
		return dotPar
	
	def sequence(self, state):
		"""Returns the sequence of the complex as NUPACK expects. The
		   first line is the number of independent strands, and the last
		   line determines how many times each strand appears."""
		return ('1\n' + self.hairpin.sequence  +  '\n1\n')
	
	def num_complex(self) :
		"""Counts the number of complexes in each state """
		self.n_complex = dict()
		for state in self.statespace:
			self.n_complex[state]  = 1
	  
	def calculate_energy(self):
		ParentComplex.calculate_energy( self, AUXILIARY_NUPACK_FILE+"/hairpin")
	 
	def allowed_state(self, state):
		"""Check that a state is allowed."""
		i, j = state
		return 0 <= i <= j <= self.m
	
	def possible_states(self, state ):
		"""Returns the neighbors of state"""
		i, j = state
		if (i == j):
			states = [(n, n + 1) for n in range(0, self.m)]
		else:
			states = [(i - 1, j),
					  (i + 1, j),
					  (i, j - 1),
					  (i, j + 1)]
		removed = False
		removeList = []
		for s in states :
			if s[0] == s[1] and 0  < s[0]  <=   self.m  :
				removeList.append((s[0],s[1]))
				removed= True
		for s in removeList:
			states.remove(s )
		if removed == True :
			states.append((0,0))
		return filter(self.allowed_state, states)
		
	def initial_final_state_config(self ):
  		"""sets the initial and final state for hairpin closing (zip == True) and hairpin opening (zip == False ) """
		if self.zip == True :
			initialStateConfig = (0, 0 )
			finalStateConfig  = (0, self.m)
		if self.zip == False:
			initialStateConfig = (0, self.m )
			finalStateConfig  = (0, 0)
		return [initialStateConfig, finalStateConfig]
	
def main( myPickles,real_Rate, theta, beta, gamma, zip, T, concentration, sodium, magnesium  ,  dataset_name, docID ):
	
	betaprime = ''.join(list(reversed(beta))).translate(
		TRANSLATION_TABLE)
	hairpinString = beta+gamma+betaprime

	hairpin = MyStrand(hairpinString)
	hairpin_complex = HairpinComplex( myPickles, theta, beta, gamma , betaprime, zip, hairpin, len(beta) , len(hairpin) ,   T, concentration, sodium, magnesium ,  dataset_name, docID )
	bimolTransition  = False
	return hairpin_complex.find_answers(concentration, real_Rate, bimolTransition  )
   
	
if __name__ == "__main__":
	main()
