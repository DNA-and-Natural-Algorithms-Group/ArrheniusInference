from __future__ import division
import parent
from parent import *

class HelixComplex(ParentComplex):
	"""Helix association or disocciation reaction"""
	def __init__(self , myPickles, dangleleft, dangleright, theta, zip, strand1, strand2,T, concentration, sodium, magnesium,  dataset_name, docID, name ):
		ParentComplex.__init__(self,myPickles,  theta ,  T, concentration, sodium, magnesium, dataset_name, docID )
		self.name = name
		assert strand1.complement == strand2
		self.strand1 = strand1
		self.strand2 = strand2
		self.L = len(strand1)
		assert self.L == len(strand2)
		self.zip = zip
		self.dangleleft = dangleleft
		self.dangleright = dangleright
	
		
	def dot_paren(self, state):
		"""Returns the structure of the complex in dot-paren notation"""
		L = self.L
		i, j = state
		strand2 = '.' * len(self.dangleleft)+'.' * (L - j) + '(' * (j - i) + '.' * i+ '.' * len(self.dangleright)
		strand1 =    '.' * i       + ')' * (j - i) + '.' * (L-j)
		if i == j:
			return strand1
		else:
			return strand2 + '+' + strand1
	
	def dot_paren_modify(self, state):
		"""Insert 	`*' at the start and end of the dot-paren notation, before and after all `+' signs, and also before and after every space"""
		L = self.L
		i, j = state
		strand2 = '.' * len(self.dangleleft)+'.' * (L - j) + '(' * (j - i) + '.' * i+ '.' * len(self.dangleright)
		strand1 ='.' * i       + ')' * (j - i) + '.' * (L-j)
		return '*' + strand2 + '*' +'+' +'*' + strand1  +'*'

	def sequence(self, state):
		"""Returns the sequence of the complex as NUPACK expects. The
		   first line is the number of independent strands, and the last
		   line determines how many times each strand appears."""
		i, j = state
		if i == j:
			return ('1\n' +   self.strand1.sequence +  '\n1\n')
		else:
			return ('2\n' + self.dangleleft +self.strand2.sequence +self.dangleright+ '\n' +   self.strand1.sequence  + '\n1 2 \n')
	
	def num_complex(self) :
		"""counts the number of complexes in each state """
		self.n_complex = dict()
		for state in self.statespace:
			i, j = state
			self.n_complex[state]  =  2 if i == j else 1
	
	def calculate_energy(self):
		ParentComplex.calculate_energy( self, AUXILIARY_NUPACK_FILE + "/helix"+str(self.name))
	
	def allowed_state(self, state):
		"""Check that a state is allowed."""
		i, j= state
		return 0 <= i <= j <=  self.L
	
	def possible_states(self, state):
		"""Returns the neighbors of state"""
		i, j= state
		if (i == j):
			states = [(n, n + 1) for n in range(0, self.L)]
		else:
			states = [(i - 1, j),
					  (i + 1, j),
					  (i, j - 1),
					  (i, j + 1)]
		removed = False
		removeList = []
		for s in states :
			if s[0] == s[1] and 0 < s[0] <= self.L  :
				removeList.append((s[0],s[1]))
				removed= True
		for s in removeList:
			states.remove(s )
		if removed == True :
			states.append((0,0))
		return filter(self.allowed_state, states)

	def initial_final_state_config(self ):
		"""sets the initial and final state for helix association (zip == True) and helix dissociation (zip == False ) """
		if self.zip == True :
			initialStateConfig = (0, 0 )
			finalStateConfig  = (0, self.L)
		if self.zip == False:
			initialStateConfig = (0, self.L )
			finalStateConfig  = (0, 0)
		return [initialStateConfig, finalStateConfig]
  
def main(myPickles,  real_rate, theta, beta, zip, T, concentration, sodium, magnesium,  dangle, dataset_name, docID ,name  ):
	if name ==myenums.DatasetName.REYNALDODISSOCIATE.value :
		dangleleft = "GAA"
		dangleright =dangle [len(dangleleft) + len(beta):   ]
	else :
		dangleleft = ""
		dangleright=""
	strand = MyStrand(beta)
	strandComplement = strand.complement
	helix_complex = HelixComplex( myPickles, dangleleft, dangleright, theta, zip, strand, strandComplement, T, concentration, sodium, magnesium , dataset_name, docID, name)
	return helix_complex.find_answers(concentration,real_rate, zip  )
 
if __name__ == "__main__":
	main()