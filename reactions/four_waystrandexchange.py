from __future__ import division
import parent
from parent import *


class  FourWayBranchComplex(ParentComplex):
	"""Four-way strand exchange reaction"""
	def __init__( self, myPickles, ms, ns ,   theta,complex1, complex2, reporter1, reporter2, lenX, lenm, lenn  ,T, concentration, sodium, magnesium, dataset_name, docID):
		ParentComplex.__init__(self, myPickles, theta ,  T, concentration, sodium, magnesium, dataset_name, docID )
		self.complex1  = complex1
		self.complex2 = complex2
		self.reporter1   =  reporter1
		self.reporter2 = reporter2
		self.T = T
		self.concentration = concentration
		self.sodium = sodium
		self.magnesium = magnesium
		self.X  = lenX
		self.m = lenm
		self.n =  lenn
		self.ms = ms
		self.ns = ns

	def dot_paren( self, state):
		"""Returns the structure of the complex in dot-paren notation"""
		i, j, ip, jp , k , l= state
		reporter1= '(' * k                      + '.' * ( i - k  )     +'('*( ip - i )  + '.'*( self.X + self.ms - ip )
		complex1 = '.' * (self.X +  self.m - ip ) + ')' * (ip- i )  + '.' * (   i - l   )  + '(' * l
		complex2 = ')' *  l                     + '.' *  (  j -l  )    + '(' *(  jp - j  )    + '.' * (self.X+ self.n - jp )
		reporter2= '.' * (self.X+ self.ns - jp ) + ')' * ( jp  - j)  + '.' *  (  j -k )     + ')' * k
		if ( i == ip and j == jp ):
			return [complex1  + '+' + complex2 , reporter1 + '+'  +  reporter2 ]
		elif k == 0  and l ==0 :
			return [reporter1   + '+' + complex1 , complex2   + '+'+ reporter2 ]

		else:
			return [reporter1 +'+'+  complex1+'+' + complex2+'+' + reporter2, -1  ]

	def dot_paren_modify( self,state):
		"""Insert 	`*' at the start and end of the dot-paren notation, before and after all `+' signs, and also before and after every space"""
		i, j, ip, jp , k , l= state
		reporter1= '(' * k                      + '.' * ( i - k  )     +'('*( ip - i )  + '.'*( self.X + self.ms - ip )
		complex1 = '.' * (self.X +  self.m - ip ) + ')' * (ip- i )  + '.' * (   i - l   )  + '(' * l
		complex2 = ')' *  l                     + '.' *  (  j -l  )    + '(' *(  jp - j  )    + '.' * (self.X+ self.n - jp )
		reporter2= '.' * (self.X+ self.ns - jp ) + ')' * ( jp  - j)  + '.' *  (  j -k )     + ')' * k
		return  '*' + reporter1 + '*' + '+'+   '*' + complex1+ '*' + '+' + '*' +  complex2+ '*' + '+' + '*' +  reporter2 +  '*'

	def sequence( self, state):
		"""Returns the sequence of the complex as NUPACK expects. The
		   first line is the number of independent strands, and the last
		   line determines how many times each strand appears."""
		i, j, ip, jp,  k , l = state
		if ( i == ip   and j == jp  ):
			#return [complex1  + '+' + complex2 , reporter1 + '+'  +  reporter2 ]
			return [ ('2\n' + self.complex1.sequence + '\n' +  self.complex2.sequence+ '\n1 2\n'),
			('2\n' + self.reporter1.sequence + '\n' +  self.reporter2.sequence + '\n1 2\n') ]
		elif k == 0  and l == 0 :
			return [ ('2\n' + self.reporter1.sequence + '\n' +  self.complex1.sequence+ '\n1 2\n'),
			('2\n' + self.complex2.sequence + '\n' +  self.reporter2.sequence + '\n1 2\n') ]
		else:
			return [   ('4\n' + self.reporter1.sequence + '\n' +
					self.complex1.sequence + '\n' +
					self.complex2.sequence + '\n' +
					self.reporter2.sequence + '\n1 2 3 4\n'), -1 ]


	def num_complex(self) :
		"""counts the number of complexes in each state """
		self.n_complex = dict()
		for state in self.statespace:
			dot_parenstate1 = self.dot_paren(state)
			self.n_complex[state]  = 2 if  dot_parenstate1[1]  != -1  else 1

	def calculate_energy(self):
		try:
			fs =open(self.dataset_name+"/"+myenums.Permanent_Folder.ENERGY.value+ "/"  + myenums.Permanent_Folder.ENERGY.value+ str(self.docID)   , "rb" )
			self.energies  = pickle.load(fs )
			fs.close()
		except:
			self.energies = dict()
			print myenums.Permanent_Folder.ENERGY.value +" creating RESUBLE files in the simplifiedstatespace folder!.   Please wait untill you see a message showing the computed lnprobability and iteration time! These messages will not appear again in the next iteration!"
			shell = Popen("bash", stdin=PIPE, stdout=PIPE)
			for state in self.statespace:

				sequence1, sequence2   = self.sequence(state)
				dot_paren1 , dot_paren2 = self.dot_paren(state)
				file_contents1 = sequence1 + dot_paren1
				filename1 = AUXILIARY_NUPACK_FILE+ "/fourway"+ str(self.docID) +"a"
				energy1  =self.nupack_energy( shell, filename1, file_contents1 )
				energy2  = 0
				if dot_paren2 != -1 :
					filename2 = AUXILIARY_NUPACK_FILE + "/four_waystrandexhcange"  + str(self.docID) +"a"
					file_contents2 = sequence2 + dot_paren2
					energy2 =self.nupack_energy(shell, filename2, file_contents2  )
				energy = energy1 + energy2
				"""NUPACK uses mole fraction units when computing free energy. We correct to use molar concentration units:  DeltaG_molar = DeltaG - (N-1)RT*log(molarity(H2O)) where N is the number of interacting complexes"""
				complex_count = self.n_complex[state]
				energy -= ((complex_count - 1) * R * (self.T + 273.15) *
						   np.log(MOLARITY_OF_WATER))
				self.energies[state] =  energy
			shell.stdout.close()
			shell.stdin.close()
			pickle.dump(        self.energies ,  open(self.dataset_name+"/"+myenums.Permanent_Folder.ENERGY.value+ "/"  + myenums.Permanent_Folder.ENERGY.value+str(self.docID) , "wb") )

	def allowed_state(self , state):
		"""Checks that a state is allowed."""
		i, j, ip, jp ,  k , l = state
		allow =  (l <= i  and l <= j and k <= i and k <= j and  0 <=k <= self.X and 0 <= l <= self.X  and  0 <= i  <= ip and ip <= self.X+ self.m  and  0 <= j <= jp and jp<=self.X +self.n )
		# Further prune the statespace to make computations tractable
		if    (   (   i  == ip  or  j == jp   ) and  (  k  ==0 or l ==0 )   )   :
			allow = False
		if (self.m ==0 or self.n ==0 ) and ( ( i == ip   )  or  (j == jp ) ) and  ( ( k < self.X -(2 - int(self.m/3) )  -1) or  (l <  self.X   -(2 - int(self.n/3) )  -1) ):
			allow = False
		if (self.m !=0 and self.n !=0 ) and ( ( i == ip   )  or  (j == jp ) ) and  ( ( k < self.X -1 ) or  (l <  self.X   -1) ):
			allow = False


		if ip < self.X   or jp < self.X:
			allow = False
		if  (i != ip and j != jp ) and  (  abs ( i  -l   )+  abs ( i - k )  + abs(  j - l ) +  abs ( j -  k ) >  4 +   ( 2- self.n/3  )+ ( 2- self.m/3  ) ):
			allow = False

		return allow

	def possible_states( self, state):
		"""Returns the neighbors of state"""
		
		if state in self.PSD:
			return self.PSD[state]
		i, j, ip, jp ,  k , l = state
		states = []
		maxkl  = max ( k, l )
		if (i == ip ) :
			states += [(n , j, n+ 1, jp , k, l ) for n in range ( maxkl , self.X  +  self.m) ]
		else :
			states += [(i - 1, j, ip, jp , k, l ),
					  (i + 1, j, ip, jp , k, l ),
					  (i, j, ip - 1 , jp , k, l ),
					  (i, j, ip  +1 , jp , k, l )  ]
		if (j == jp )  :
			states += [(i , n , ip, n+1 , k, l )  for n in range ( maxkl , self.X  + self.n)]
		else :
			states +=  [(i, j - 1, ip, jp ,k, l ),
					  (i, j, ip, jp - 1 , k, l ),
					  (i, j + 1,ip, jp , k, l ),
					  (i, j, ip, jp  + 1 , k, l ) ]
		if ( i!=ip or j!= jp ) or  ( self.m == 0 or  self.n ==0 )  :
			states += [(i, j, ip, jp , k - 1, l ) ]
			states += [(i, j, ip, jp , k + 1, l ) ]
			states += [  (i, j, ip, jp ,k , l  -1 ) ]
			states += [  (i, j, ip, jp ,k , l  +1 ) ]
		removed = False
		removeList = []
		for s in states :
			if s[0] == s[2] and  0 <= s[0] < self.X+ self.m :
				removeList.append((s[0],s[1], s[2] , s[3] , s[4] , s[5] ))
				removed= True
		for s in removeList:
			states.remove(s )
		if removed == True :
			states.append((self.X+ self.m  , s[1], self.X+ self.m , s[3] , s[4], s[5] ))
		removed = False
		removeList = []
		for s in states :
			if s[1] == s[3] and  0 <= s[1] < self.X+ self.n :
				removeList.append((s[0],s[1], s[2] , s[3] , s[4] , s[5] ))
				removed= True
		for s in removeList:
			states.remove(s )
		if removed == True :
			states.append((s[0], self.X+ self.n , s[2], self.X+ self.n  , s[4], s[5] ))
		filteredStates = filter(self.allowed_state, states)
		self.PSD[state ] = filteredStates
		return filteredStates

	def initial_final_state_config(self ):
		"""sets the initial and final state for four-way strand exchange"""
		initialStateConfig = (self.X  + self.m , self.X + self.n ,  self.X+ self.m , self.X + self.n , self.X , self.X  )
		finalStateConfig= ( 0 ,   0 , self.X  + self.m , self.X + self.n ,  0    ,  0  )
		return [initialStateConfig, finalStateConfig]

def main(myPickles, bi_real_rate, uni_real_rate ,  lenm, lenn, Xstr,toeholdstrm, toeholdstrn   , ms, ns ,  theta, T, concentration, sodium, magnesium  ,dataset_name, docID  , name):
	Xstr= Xstr.encode()
	toeholdstrm= toeholdstrm.encode()
	toeholdstrn= toeholdstrn.encode()
	complex1  = MyStrand(toeholdstrm[len(toeholdstrm)- lenm: ] + Xstr )
	FullToehold = MyStrand(toeholdstrm[len(toeholdstrm)- ms : ] + Xstr )
	reporter1  = FullToehold.complement
	complex2 = MyStrand( MyStrand(Xstr).complement.sequence  + toeholdstrn[ :lenn ])
	FullToehold = MyStrand( MyStrand(Xstr).complement.sequence  + toeholdstrn[ :  ns  ])
	reporter2 = FullToehold.complement
	lenX = len(Xstr)
	fourwaybranch_complex= FourWayBranchComplex (   myPickles,  ms, ns ,  theta,complex1, complex2, reporter1, reporter2,lenX, lenm, lenn, T, concentration, sodium, magnesium , dataset_name, docID )
	bimolTransition = True
	if name == myenums.DatasetName.DABBY.value:
		"""Table 5.1 from  DABBY's thesis report unimolecular rates and Table 5.2 reports bimolecular rates separately. Combining the two rates to obtain an overall reaction rate constant """
		t1 = 1 / (bi_real_rate * concentration)
		t2 = 1 / uni_real_rate
		real_rate = 1 / ((t1 + t2) * concentration)
	return fourwaybranch_complex.find_answers(concentration, real_rate, bimolTransition )

if __name__ == "__main__":
	main( )
