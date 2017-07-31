from __future__ import division
import sys
import os
sys.path.insert( 0 , os.path.realpath('../parent'))
import parent
from parent import *

class   ThreeWayBranchComplex(ParentComplex)  :
	"""Three-way strand displacement reaction.
	For the Machinek paper considers mismatches"""
	def __init__( self,myPickles ,beta, gamma,toehold_length, hasincumbentDangle, incumbentDangle,  dangleleft, dangleright, theta, substrate, attacker, incumbent,T, concentration, sodium, magnesium , dataset_name, docID , name  , substrateDangle, mismatchPosition):
		ParentComplex.__init__(self, myPickles, theta ,  T, concentration, sodium, magnesium, dataset_name, docID )
		self.substrateDangle = substrateDangle
		self.mismatchPosition = mismatchPosition
		self.beta = beta
		self.gamma = gamma
		self.hasincumbentDangle = hasincumbentDangle
		self.incumbentDangle = incumbentDangle
		self.name = name
		self.dangleleft = dangleleft
		self.dangleright = dangleright
		self.toehold_length = toehold_length
		self.substrate = substrate
		self.attacker = attacker
		self.incumbent = incumbent
		self.L = len(substrate)
		assert self.L == len(attacker)
		self.n = len(incumbent)
		self.m = self.L - self.n
		
	def dot_paren( self, state):
		"""Returns the structure of the complex in dot-paren notation"""
		# Note that if i == j, then the attacker strand is
		# not present. If k == l, the incumbent is not present.
		m, L = self.m, self.L
		i, j, k, l = state
		
		attacker = '.' * (L - j) + '(' * (j - i) + '.' * i
		if self.name ==myenums.DatasetName.ZHANG.value or self.name ==myenums.DatasetName.MACHINEK.value:
			substrate =  ('.' * i       + ')' * (j - i) + '.' * (k - j) +
					  '(' * (l - k) + '.' * (L - l))
		if self.name ==myenums.DatasetName.REYNALDOSEQUENTIAL.value:
			substrate =  (  '.'*len(self.dangleleft )+ '.' * i       + ')' * (j - i) + '.' * (k - j) +
					  '(' * (l - k) + '.' * (L - l)  + '.'*len(self.dangleright))
		
		if self.hasincumbentDangle == True  :
			incumbent = '.'*  len(self.incumbentDangle) +  '.' * (L - l) + ')' * (l - k) + '.' * (k - m)
		else :
			incumbent = '.' * (L - l) + ')' * (l - k) + '.' * (k - m)
		if self.name == myenums.DatasetName.MACHINEK.value :
			if self.mismatchPosition != -1:
				attacker  =  attacker [: len (substrate ) - self.mismatchPosition   -1 ] + '.' + attacker[len (substrate ) - self.mismatchPosition    : ]
				
				if   substrate[self.mismatchPosition    ] ==')' :

					substrate = substrate[:self.mismatchPosition] + '.' + substrate[ self.mismatchPosition  + 1   : ]
			substrate =  len (self.substrateDangle ) * '.' + substrate
		if i == j:
			return substrate + '+' + incumbent
		elif k == l:
			return attacker + '+' + substrate
		else:
			return attacker + '+' + substrate + '+' + incumbent
	
	def dot_paren_modify(self, state):
		"""Insert 	`*' at the start and end of the dot-paren notation, before and after all `+' signs, and also before and after every space"""
		m, L = self.m, self.L
		i, j, k, l = state
		attacker = '.' * (L - j) + '(' * (j - i) + '.' * i
		if self.name == myenums.DatasetName.ZHANG.value or self.name == myenums.DatasetName.MACHINEK.value:
			substrate = ('.' * i + ')' * (j - i) + '.' * (k - j) +
						 '(' * (l - k) + '.' * (L - l))
		if self.name == myenums.DatasetName.REYNALDOSEQUENTIAL.value:
			substrate = ('.' * len(self.dangleleft) + '.' * i + ')' * (j - i) + '.' * (k - j) +
						 '(' * (l - k) + '.' * (L - l) + '.' * len(self.dangleright))
		if self.hasincumbentDangle == True:
			incumbent = '.' * len(self.incumbentDangle) + '.' * (L - l) + ')' * (l - k) + '.' * (k - m)
		else:
			incumbent = '.' * (L - l) + ')' * (l - k) + '.' * (k - m)
		if self.name == myenums.DatasetName.MACHINEK.value:
			if self.mismatchPosition != -1:
				attacker = attacker[: len(substrate) - self.mismatchPosition - 1] + '.' + attacker[len(
					substrate) - self.mismatchPosition:]
				if substrate[self.mismatchPosition] == ')':
					substrate = substrate[:self.mismatchPosition] + '.' + substrate[self.mismatchPosition + 1:]
			substrate = len(self.substrateDangle) * '.' + substrate
		return '*' + attacker + '*' + '+' + '*' + substrate + '*' + '+' + '*' + incumbent + '*'
	
	def sequence(self,  state):
		"""Returns the sequence of the complex as NUPACK expects. The
		   first line is the number of independent strands, and the last
		   line determines how many times each strand appears."""
		# Note that if i == j, then the attacker strand is not present.
		# If k == l, the incumbent is not present.
		i, j, k, l = state
		realincumbent= self.incumbent.sequence
		if self.hasincumbentDangle == True :
			realincumbent = self.incumbentDangle + realincumbent
		realsubstrate = self.substrate.sequence
		if self.name == myenums.DatasetName.MACHINEK.value :
			realsubstrate =self.substrateDangle +  realsubstrate

		if self.name == myenums.DatasetName.REYNALDOSEQUENTIAL.value :
			realsubstrate= self.dangleleft+realsubstrate+ self.dangleright
		if i == j:
			return ('2\n' + realsubstrate + '\n' +  realincumbent + '\n1 2\n')
		elif k == l:
			return ('2\n' + self.attacker.sequence + '\n' +
					realsubstrate + '\n1 2\n')
		else:
			return ('3\n' + self.attacker.sequence + '\n' +
					realsubstrate + '\n' +
					realincumbent+ '\n1 2 3\n')

	def num_complex(self) :
		"""counts the number of complexes in each state """
	
		self.n_complex = dict()
		
		for state in self.statespace:
			i, j, k, l = state
			self.n_complex[state]= 2 if i == j or k == l else 1
			

	def calculate_energy(self):
		ParentComplex.calculate_energy( self, AUXILIARY_NUPACK_FILE+"/three_waystranddisplacement" + str(self.name))
	def allowed_state(self, state):
		"""Checks that a state is allowed."""
		i, j, k, l = state
		allow =  0 <= i <= j <= k <= l <= self.L and k >= self.m
		#Further prune the statespace to make computations tractable
		if self.name != myenums.DatasetName.MACHINEK.value :
			if  ( self.toehold_length >0 )  and ( k > self.m  ) and ( i != 0 or j < k- 1) :
				allow =False
			if ( (i == j and k == l)  or abs ( j - k  ) > self.toehold_length  + 2):
				allow = False
		elif self.name ==myenums.DatasetName.MACHINEK.value:
			if  ( self.toehold_length >0 )  and ( k > self.m  ) and ( i != 0 or j < k-5) :
				allow =False
			if ( (i == j and k == l)  or abs ( j - k  ) > self.toehold_length  + 4):
				allow = False
		return allow

	def possible_states( self, state):
		"""Returns the neighbors of state"""

		i, j, k, l = state
		if (i == j):
			states = [ ]
			if self.name ==myenums.DatasetName.REYNALDOSEQUENTIAL.value :
				states += [(i-1, j, k, l) ,  (i+1 , j,  k , l) , (i, j-1 ,  k , l), (i, j+1 ,  k , l), (i, j,  k- 1 , l) ,(i, j,  k+ 1 , l), (i, j,  k , l-1) , (i, j,  k, l +1)   ]
			if self.name != myenums.DatasetName.MACHINEK.value  :
				states += [(n, n + 1, k, l) for n in range(0, k)]
			 
			if self.name ==myenums.DatasetName.MACHINEK.value:
				states += [(n, n + 1, k, l) for n in range(0, k) if n !=  self.mismatchPosition  ]
		else:
			if self.name != myenums.DatasetName.MACHINEK.value :
				states = [(i - 1, j, k, l),
						  (i + 1, j, k, l),
						  (i, j - 1, k, l),
						  (i, j + 1, k, l),
						  (i, j, k - 1, l),
						  (i, j, k + 1, l),
						  (i, j, k, l - 1),
						  (i, j, k, l + 1)]
			if self.name == myenums.DatasetName.MACHINEK.value :
				
				states =  [ (i, j, k - 1, l),
						  (i, j, k + 1, l),
						  (i, j, k, l - 1),
						  (i, j, k, l + 1)]
				if self.mismatchPosition != -1:
					if i -1!=self.mismatchPosition    :
						states+= [(i - 1, j, k, l)]
					else :
						states+= [(i - 2, j, k, l)]
					if i + 1 != self.mismatchPosition     :
						states+=  [ (i  +  1, j, k, l) ]
					else :
						states += [(i  +  2, j, k, l) ]
						
					if j -1 -1 !=self.mismatchPosition    :
						states+= [(i, j -1 , k, l)]
					else :
						states+= [(i , j -2 , k, l)]
					if j+  1  - 1 != self.mismatchPosition     :
						states+=  [ (i , j + 1 , k, l) ]
					else :
						states += [(i , j +2   , k, l) ]
				else  :
					states += [(i - 1, j, k, l),  (i + 1, j, k, l),         (i, j - 1, k, l),            (i, j + 1, k, l) ]
			  
	
		removed = False
		removeList = []
		for s in states :
			if s[2] == s[3] and  0 <= s[2] < self.L :
				removeList.append((s[0],s[1], s[2] , s[3] ))
				removed= True
		for s in removeList:
			states.remove(s )
		if removed == True :
			states.append((s[0], s[1], self.L , self.L  ))
		
		
		removed = False
		removeList = []
		for s in states :
			if s[0] == s[1] and  0 < s[2] <= self.L :
				removeList.append((s[0],s[1], s[2] , s[3] ))
				removed= True
		for s in removeList:
			states.remove(s )
		if removed == True :
			states.append(( 0 , 0 , s[2],s[3]  ))

		return filter(self.allowed_state, states)

	def initial_final_state_config(self ):
		"""sets the initial and final state for three-way strand displacement """
		initialStateConfig = (0, 0 , self.m, self.L )
		finalStateConfig = (0, self.L , self.L , self.L )
		return [initialStateConfig, finalStateConfig]
	
def main(myPickles, hasincumbentDangle,  incumbentDangle,real_rate, toehold_length,beta,gamma,  theta, T, concentration, sodium, magnesium   , dangle,  dataset_name, docID, name , incumbentMachinek, targetMachinek, invaderMachinek , mismatchMachinek ):
	
	beta= beta.encode()
	gamma = gamma.encode()
	if name ==myenums.DatasetName.MACHINEK.value :
		invaderMachinek = invaderMachinek.encode()
		incumbentMachinek = incumbentMachinek.encode()
		targetMachinek = str(targetMachinek.encode() )
		aa =  len (targetMachinek) - len( invaderMachinek)

		attacker = MyStrand(invaderMachinek)
		targetMachinekReal = targetMachinek [aa : ]
		substrateDangle = targetMachinek[0: aa]

	 
		substrate = MyStrand(targetMachinekReal)
		incumbent = MyStrand( incumbentMachinek)
		if mismatchMachinek != -1 :
		 
			mismatchPosition = len (attacker) -  (  toehold_length + mismatchMachinek )
			mismatchPosition = len (attacker ) - mismatchPosition - 1
		else :
			mismatchPosition =-1
	  
		dangleleft =""
		dangleright=""
	if name == myenums.DatasetName.ZHANG.value :
		attacker = MyStrand(beta + gamma[:toehold_length])

		substrate = attacker.complement
		incumbent = MyStrand(attacker.sequence[:-toehold_length])
		dangleleft =""
		dangleright=""
		substrateDangle =""
		mismatchPosition =""
	if name == myenums.DatasetName.REYNALDOSEQUENTIAL.value :
		attacker = MyStrand(beta)
		substrate = attacker.complement
		incumbent = MyStrand(beta)
		dangleleft = "GAA"
		dangleright = dangle[len(dangleleft) + len(beta):   ]
 
		substrateDangle =""
		mismatchPosition =""
	threewaybranch_complex = ThreeWayBranchComplex(myPickles,  beta, gamma, toehold_length ,  hasincumbentDangle, incumbentDangle, dangleleft, dangleright,  theta,substrate,
		attacker, incumbent, T, concentration, sodium, magnesium , dataset_name, docID  , name , substrateDangle, mismatchPosition)
	bimolTransition  = True
 
	return threewaybranch_complex.find_answers(concentration, real_rate, True )
   
  
   
if __name__ == "__main__":
	main()
