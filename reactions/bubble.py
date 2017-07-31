from __future__ import division
import parent
from parent import *

class HairpinComplex(ParentComplex):
    """Bubble closing reaction! The Altan-Bonnet paper has bubbles in a hairpin! """
    def __init__(self,  myPickles,  theta, strand1, loop, strand2,   T,concentration, sodium, magnesium , flurPosition,  dataset_name, docID ):
        ParentComplex.__init__(self, myPickles, theta ,  T, concentration, sodium, magnesium,  dataset_name, docID )
        self.flurPosition = flurPosition
        assert strand1.complement == strand2
        self.strand1 = strand1
        self.strand2 = strand2
        self.loop = loop
        self.L = len(strand1)
        assert self.L == len(strand2)
        self.i = 0  # the i is always fixed and equal to 0
        self.l = self.L # the l is always fixed  and equal to L
        [self.initialj, self.initialk] = [self.flurPosition, self.flurPosition + 1]
        self.j = self.initialj
        self.k = self.initialk
        
    def dot_paren(self, state):
        #return the dot paranthesis notation
        L = self.L 
        i, j, k , l = state
        strand2 = '.' * (i ) + '(' * (j - i) + '.' * (k - j )+ '(' * (l - k ) + '.' * (L - l ) 
        strand1 =   '.' * (L-l )       + ')' * (l - k ) + '.' * (k-j) +  ')' * (j-i)  + '.' * i 
        loop  = '.' * len(self.loop)
        return strand2 + loop  + strand1         

    def dot_paren_modify(self, state):
        """Insert 	`*' at the start and end of the dot-paren notation, before and after all `+' signs, and also before and after every space"""
        L = self.L
        i, j, k , l = state
        strand2 = '.' * (i ) + '(' * (j - i) + '.' * (k - j )+ '(' * (l - k ) + '.' * (L - l )
        strand1 =   '.' * (L-l )       + ')' * (l - k ) + '.' * (k-j) +  ')' * (j-i)  + '.' * i
        loop  = '.' * len(self.loop)
        return '*' + strand2 + '*' + loop  + '*' + strand1 +'*'

    def sequence(self, state):
        """Returns the sequence of the complex as NUPACK expects. The
           first line is the number of independent strands, and the last
           line determines how many times each strand appears."""
        return ('1\n' + self.strand2.sequence +  self.loop.sequence+  self.strand1.sequence + '\n1 \n')
    
    def num_complex(self) :
        """counts the number of complexes in each state """
        self.n_complex = dict()   
        for state in self.statespace:
            self.n_complex[state]  = 1   
    def calculate_energy(self):
        ParentComplex.calculate_energy(self,  AUXILIARY_NUPACK_FILE +"/bubble")

    def allowed_state(self, state):
        """Check that a state is allowed."""
        i, j, k, l = state
        return (0 <= i <= j <= self.flurPosition <= k <= l <= self.L ) and (0 != j  and  k != self.L) 
    
    
    
    def possible_states(self, state):
        """Returns the neighbors of state"""
    
        i, j, k, l = state
        states = [ (i, j - 1, k, l),
        (i, j + 1, k, l),
        ( i , j, k - 1 , l),
        (i , j, k + 1 , l)]
        removed = False  
        removeList = []
        for s in states : 
            if s[1] == s[2] and s[1] != self.flurPosition : 
                removeList.append((s[0],s[1], s[2], s[3]))  
                removed= True
        for s in removeList: 
            states.remove(s ) 
        if removed == True : 
            states.append((0,self.flurPosition, self.flurPosition, self.L))
        return filter(self.allowed_state, states)
   
    def initial_final_state_config(self ):
        """sets the initial and final state for bubble closing """
        initialStateConfig =(0,self.initialj, self.initialk, self.L )
        finalStateConfig  = (0,  self.flurPosition , self.flurPosition,  self.L)

        return [initialStateConfig, finalStateConfig]

def main(myPickles, real_rate, theta , beta , loop , gamma ,T, concentration, sodium, magnesium  , flurPosition , dataset_name, docID ):
    strand  = MyStrand(beta)
    strandComplement=MyStrand(gamma)
    loop = MyStrand(loop)
    helix_complex = HairpinComplex( myPickles, theta, strand,  loop,strandComplement,  T, concentration, sodium, magnesium, flurPosition, dataset_name, docID )
    bimolTransition  = False 
    return helix_complex.find_answers(concentration , real_rate, bimolTransition  )
   
if __name__ == "__main__":
    main()