from __future__ import division
import warnings
import gc
import numpy as np
import  random, string
from subprocess import Popen, PIPE
import ConfigParser
import math
import cPickle as pickle
from  scipy.sparse.linalg import * 
from scipy.sparse import csr_matrix, coo_matrix
import myenums

configParser = ConfigParser.ConfigParser()
configParser.readfp(open(r'../learndnakinetics/config_file.txt'))
CONFIG_NAME = 'parent'
NUPACK_bin= configParser.get(CONFIG_NAME, 'NUPACK_bin')
AUXILIARY_NUPACK_FILE= 'auxilary'
R = 0.001987 # R is the molar gas constant in kcal/(mol K).
MOLARITY_OF_WATER = 55.14 # mol/L at 37C
NUCLEOTIDES = "ACTG"
TRANSLATION_TABLE = string.maketrans(NUCLEOTIDES, "TGAC")
RETURN_MINUS_INF = None
class MultiStrandState (object):
    def __init__(self,  uniqueID, uniqueDotParanthesis, energy , sequence):
        self.uniqueID = uniqueID 
        self.uniqueDotParanthesis = uniqueDotParanthesis 
        self.sequence = sequence
        self.energy = energy 
        self.neighborsList  =  [] 
class MyStrand(object):
    @classmethod
    def new_random(cls, length):
        """Create a random strand of a certain length"""
        sequence = ''.join(random.choice("ACTG") for i in range(length))
        return cls(sequence)
    
    def __init__(self, sequence, complement=None):
        """Create a strand by providing a sequence from 5' to 3' ends"""
        self.sequence = sequence
        if complement:
            self.complement = complement
        else:
            seq = ''.join(list(reversed(self.sequence))).translate(
                    TRANSLATION_TABLE)
            self.complement = MyStrand(seq, self)
    
    def __len__(self):
        return len(self.sequence)
        
    def __eq__(self, other):
        return self.sequence == other.sequence

class ParentComplex(object):
    """Contains function and variables that different type of reaction have in common"""
    def __init__(self,myPickles, theta,  T, concentration, sodium, magnesium, dataset_name, docID ):
  
        if rate_method == myenums.ModelName.ARRHENIUSMODELNAME.value :
            self.kinetic_parameters = { "stack": (theta[0] , theta[1]) ,
                               "loop": (theta[2] , theta[3]),
                               "end": (theta[4] , theta[5]),
                               "stack+loop": (theta[6] , theta[7]),
                               "stack+end": (theta[8] , theta[9]),
                               "loop+end": (theta[10] , theta[11]),
                               "stack+stack": (theta[12] , theta[13]), 
                               "alpha" : (theta[14]) }
        elif rate_method == myenums.ModelName.METROPOLISMODELNAME.value :
            self.kinetic_parameters ={ "k_uni" : theta[0] , "k_bi" : theta[1] }
        else:
            raise ValueError('Error: Please specify rate_method to be Arrhenius or Metropolis in the configuration file!')
        
        
        self.dataset_name = dataset_name
        self.docID = docID 
        self.T = T
        self.concentration = concentration
        self.sodium = sodium
        self.magnesium = magnesium
        self.fast_access= dict() 
        self.statespace =[] 
        self.multistrand_state_dictionary = dict()   
        self.energies = dict()
        self.PSD = dict()
        self.transition_structure= myPickles [myenums.Permanent_Folder.TRANSITION_STRUCTURE.value]
        self.PSD= myPickles [myenums.Permanent_Folder.PSD.value]
        self.rates ={}
        self.local_context_bi = dict()     
        self.local_context_uni = dict()     
        for i in self.kinetic_parameters: 
            for j in self.kinetic_parameters:
                self.local_context_bi [i , j]  = 0
                self.local_context_uni [i  , j]  = 0
  
    
 
        
    def possible_states(self, state):
        return self.multistrand_state_dictionary[state].neighborsList
            


           
    def local_context(self, state1,  state2) : 
        """ Finds the local context of a base pair forming or breaking in transition from state1 to state2 """
        s1 = self.dot_paren_modify(state1)
        s2= self.dot_paren_modify(state2)
        count = 0 
        
        for i  in  range( len(s1)  ): 
            if s1[i] != s2[i]:
                if count == 0 :
                    found1 = i 
                else : 
                    found2 = i 
                count +=1 
                if count ==2 : 
                    break
        right = self.half_context(s1[found1 + 1] , s2[found2 -1 ]  , found1   + 1 , found2  - 1, s1  )
        left  = self.half_context(s1[found1 - 1] , s2[found2 + 1 ] , found1 - 1, found2 + 1 , s1 )
        return (left, right)

    def half_context( self,  c1, c2 , f1, f2, s1) :
        """ Finds the half context on one side of a base pair forming or breaking """
        if c1 =='*' and c2 =='*' : 
            return "end" 
        if c1 ==   '.' and c2 =='.' : 
            return "loop"
        if c1 == '(' and c2 ==')' : 
            countStack = 0 
            pointer = -1 
            for i in range( f1 , f2 + 1 ) : 
                if s1[i]  == '('   :
                    countStack  = countStack + 1 
                elif s1[i] == ')' : 
                    countStack  = countStack -  1 
                if countStack == 0 : 
                    pointer = i  
                    break 
            if pointer == f2 : 
                return "stack" 
            else: 
                return "stack+stack" 
        if  ( c1 == ')' and c2== '('  ) or   ( c1 == ')' and c2== ')' ) or  (  c1 == '(' and c2== '(' )  : 
            return "stack+stack" 
        if  ( c1 == '(' and c2 == '.' )    or ( c1 == ')' and c2 == '.' )    or ( c1 == '.' and c2 == ')' )    or ( c1 == '.' and c2 == '(' )    : 
            return "stack+loop"  
        if  ( c1 == '(' and c2 == '*' )    or ( c1 == ')' and c2 == '*' )    or ( c1 == '*' and c2 == ')' )    or ( c1 == '*' and c2 == '(' )    : 
            return "stack+end"  
        if  ( c1 == '.' and c2 == '*' )    or ( c1 == '*' and c2 == '.' )       : 
            return "loop+end"     
    
        
    def initial_final_state(self ):
        initialStateConfig, finalStateConfig = self.initial_final_state_config()
        initialState = self.statespace.index(initialStateConfig) 
        finalState = self.statespace.index(finalStateConfig) 
    
        return [initialState, finalState]
    

    def initial_final_state_config(self ) :
        return [initialStateConfig, finalStateConfig]
        
    
     
    def generate_statespace( self   ):
        #Generates the state space
        state = self.initial_final_state_config()[0] 
        self.statespace = [] 
        self.statespace.append(state)
        self.fast_access[state] = len(self.statespace)- 1         
        color= dict()
        color [state] = True
        head =  0  
        tail = 1
        while head < tail    :
            state = self.statespace[head ]
            pstates = self.possible_states( state )
            for s in pstates :
                if s not in color  :
                    color[s] = True
                    self.statespace.append(s ) 
                    self.fast_access[s] = len(self.statespace)- 1 
                    tail +=1 
            head += 1
        return self.statespace
    

    def calculate_energy (self, auxilary_f):
        """ Calculate the energy of all states using NUPACK"""
        try:   
            fs =open( self.dataset_name+"/"+myenums.Permanent_Folder.ENERGY.value+ "/"  + myenums.Permanent_Folder.ENERGY.value+str(self.docID)  , "rb" )
            self.energies  = pickle.load(fs ) 
            fs.close() 
        except:
           
            self.energies = dict() 
            print myenums.Permanent_Folder.ENERGY.value +" doesnt't exist!"
            shell = Popen("bash", stdin=PIPE, stdout=PIPE)   
            for state in self.statespace: 
                filename = auxilary_f +self.docID
                file_contents = self.sequence(state) + self.dot_paren(state)
                energy  = self.nupack_energy(shell, filename,file_contents )
                complex_count =  self.n_complex[state]  
             
                """NUPACK uses mole fraction units when computing free energy. We correct to use molar concentration units:  DeltaG_molar = DeltaG - (N-1)RT*log(molarity(H2O)) where N is the number of interacting complexes"""
                energy -= ((complex_count - 1) * R * (self.T + 273.15) *
                           np.log(MOLARITY_OF_WATER))
                
                self.energies[state] = energy

            shell.stdout.close()
            shell.stdin.close()
            pickle.dump(        self.energies ,  open(self.dataset_name+"/"+myenums.Permanent_Folder.ENERGY.value+ "/"  + myenums.Permanent_Folder.ENERGY.value+str(self.docID) , "wb") )
   
  
       
    def nupack_energy(self, shell, filename, file_contents) :
        """ Calls NUPACK to calculate energy """
        with open(filename + '.in', 'w') as f:
            f.write(file_contents)

        command = ("%senergy -material dna -T %f -sodium %f "
                   "-magnesium %f  -multi -dangles some %s\n" %
                   (NUPACK_bin, self.T,self.sodium, self.magnesium, filename))
        
        shell.stdin.write(command)
        line = '%'
        while line.startswith('%'):
            line = shell.stdout.readline()
      
        energy = float(line)
        return energy   

    def Metropolis_rate(self, state1, state2 ):
        """ Uses the Metropolis kinetic model to calculate the transition  rate from state1 to state and the transition rate from state2 to state1. Only returns the transition rate from state1 to state2 """
        
        transition1 = (state1,state2 )
        transition2 = (state2,state1)
        rate1 = self.rates.get(transition1)
        
        if rate1:
            return rate1
       
        k_uni  = self.kinetic_parameters["k_uni"]
        k_bi  = self.kinetic_parameters["k_bi"]
        DeltaG = (self.energies[state2] -  self.energies[state1])
        DeltaG2   = -DeltaG
        RT = R * (self.T + 273.15)
    
        if   ( self.n_complex[state1] - self.n_complex[state2] ) == 1    :
            rate1 = k_bi * self.concentration
            rate2 = k_bi * np.e ** ( - DeltaG2 / RT)
        
        elif (self.n_complex[state1]  - self.n_complex[state2] ) ==  -1   :
            rate1 = k_bi * np.e ** ( - DeltaG / RT)
            rate2 = k_bi * self.concentration
         
        elif  self.n_complex[state1] == self.n_complex[state2]   :
            
            if DeltaG > 0.0:
                rate1 = k_uni  * np.e **(-DeltaG  / RT)
                rate2 = k_uni
            else:
                rate1 = k_uni
                rate2 = k_uni  * np.e **(-DeltaG2  / RT)
        else :
            raise ValueError('Exception, fix this in Metropolis_rate function.  Check transition rate calculations!')
       
        self.rates[transition1] = rate1
        self.rates[transition2] = rate2
        
        return rate1
   


    
    def Arrhenius_rate(self, state1, state2   ):
        """Uses the Arrhenius kinetic model to calculate the transition  rate from state1 to state and the transition rate from state2 to state1. Only returns the transition rate from state1 to state2 """

        transition1 = (state1,state2 )
        transition2 = (state2,state1)
        rate1 = self.rates.get(transition1)
       
        if rate1:
            return rate1
        try :
            
            left, right  = self.transition_structure[state1, state2  ]
        except:
          
            left, right = self.local_context(state1, state2)
            self.transition_structure[state1 , state2 ]  =  [left, right ]
            
        lnA_left, E_left = self.kinetic_parameters[left]
        lnA_right, E_right = self.kinetic_parameters[right]
        lnA = lnA_left + lnA_right
        E = E_left + E_right
        DeltaG = (self.energies[state2] -  self.energies[state1])
        DeltaG2   = -DeltaG
        RT = R * (self.T + 273.15)
        
   
       
        n_complex1 = self.n_complex[state1]
        n_complex2 = self.n_complex[state2]
        n_complexdiff = n_complex1 - n_complex2
    
        if   n_complexdiff ==0  :
            self.local_context_uni[left, right] += 2
            """Using plus 2 instead of plus 1 since were calculating we're calculating the transition rate from state1 to state2 and from state2 to state1 simultaneously. """
            if left != right :
                self.local_context_uni[right , left ] += 2
            if DeltaG > 0.0:
                rate1 = np.e **(lnA - (DeltaG + E) / RT)
                rate2 = np.e ** (lnA - E / RT)
            else:
                rate1 = np.e ** (lnA - E / RT)
                rate2 = np.e **(lnA - (DeltaG2 + E) / RT)
        elif   n_complexdiff == 1    :
            rate1 = (self.kinetic_parameters["alpha"] * self.concentration) * np.e  ** (lnA - E / RT)
            rate2 = self.kinetic_parameters["alpha"] * np.e ** (lnA - (DeltaG2 + E) / RT)
            self.local_context_bi[left, right] += 2
            if left != right:
                self.local_context_bi[right, left ] += 2
        elif n_complexdiff ==  -1   :
          
            self.local_context_bi[left, right] += 2
            if left != right:
                self.local_context_bi[right, left] += 2
            rate1 = self.kinetic_parameters["alpha"] * np.e ** (lnA - (DeltaG + E) / RT)
            rate2 = (self.kinetic_parameters["alpha"] * self.concentration) * np.e  ** (lnA - E / RT)
        else :
            raise ValueError('Exception, fix this in Arrhenius_rate function.  Check transition rate calculations!')
      
        self.rates[transition1] = rate1
        self.rates[transition2] = rate2
        
        return rate1
        

        
    def find_firstpassagetime(self):
        """finds the first passage time from the initial state to the final state  by solving a system of linear equations """
        
        
        try:
            """save objects to make computations faster in later iterations"""
            
            self.statespace   = pickle.load( open( self.dataset_name+"/"+myenums.Permanent_Folder.STATESPACE.value+"/"+myenums.Permanent_Folder.STATESPACE.value+str(self.docID)  , "rb" ))
            self.fast_access   = pickle.load( open( self.dataset_name+"/" + myenums.Permanent_Folder.FAST_ACCESS.value + "/"  + myenums.Permanent_Folder.FAST_ACCESS.value+str(self.docID)  , "rb" ))
            
        except:
            self.statespace = self.generate_statespace ()
            pickle.dump( self.statespace ,  open( self.dataset_name+"/"+myenums.Permanent_Folder.STATESPACE.value+"/"+myenums.Permanent_Folder.STATESPACE.value+str(self.docID) , "wb") )
            pickle.dump( self.fast_access,  open( self.dataset_name+"/" + myenums.Permanent_Folder.FAST_ACCESS.value + "/"  + myenums.Permanent_Folder.FAST_ACCESS.value+str(self.docID) , "wb") )
            
        try :
            if rate_method == myenums.ModelName.ARRHENIUSMODELNAME.value:
                if len (self.transition_structure ) == 0 :
                    self.transition_structure  = pickle.load( open( self.dataset_name+"/"+ myenums.Permanent_Folder.TRANSITION_STRUCTURE.value+ "/"  +myenums.Permanent_Folder.TRANSITION_STRUCTURE.value+str(self.docID)  , "rb" ))
            if len (self.PSD ) == 0 :
              
                self.PSD = pickle.load( open( self.dataset_name+"/"+myenums.Permanent_Folder.PSD.value+"/"+ myenums.Permanent_Folder.PSD.value +str(self.docID)  , "rb" ))

            savetransition_structure = False
            
           
     
        except :
    
            savetransition_structure  = True
        [initialState,finalState] = self.initial_final_state()
        self.num_complex()
        self.calculate_energy()
        #self.calculate_refinedenergy()
        vals = []
        rows = []
        cols = []
        d = dict()
        d1 =dict()
        diags = [0 for i in range (len(self.statespace) ) ]

        lens = len(self.statespace)
        for s1 in  range(lens)  :
            state = self.statespace[s1]
            ps = self.possible_states(state )
          
            for state2 in ps:
                s2 = self.fast_access[state2]
                if rate_method ==myenums.ModelName.ARRHENIUSMODELNAME.value :
                    myRate = self.Arrhenius_rate(state, state2 )
                elif rate_method == myenums.ModelName.METROPOLISMODELNAME.value:
                    myRate= self.Metropolis_rate(state,state2)
                else:
                    raise ValueError('Error: Please specify rate_method to be Arrhenius or Metropolis!')
          
                sss = (s1, s2, myRate ) 
        
                if (sss[0], sss[1] ) not in d1 : 
                   
                    diags[sss[0] ] = diags[sss[0]]  - sss[2]
                    d1 [sss[0] , sss[1] ] = 1 
                if sss[0] == finalState or sss[1] == finalState: 
                    continue 
                    
                row = sss[0]- (sss[0] > finalState)
                col = sss[1] - (sss[1] > finalState)

                if (row, col ) in d : 
                    continue 
               
                rows.append(row) 
                cols.append(col )
                vals.append(sss[2])
                d [ row, col  ]  =1 
        if savetransition_structure == True and rate_method ==myenums.ModelName.ARRHENIUSMODELNAME.value :
              
            pickle.dump( self.transition_structure,  open( self.dataset_name+"/"+ myenums.Permanent_Folder.TRANSITION_STRUCTURE.value+ "/"  +myenums.Permanent_Folder.TRANSITION_STRUCTURE.value+str(self.docID) , "wb") )
            pickle.dump( self.PSD,  open( self.dataset_name+"/"+myenums.Permanent_Folder.PSD.value+"/"+myenums.Permanent_Folder.PSD.value+str(self.docID) , "wb") )
      
        b= -1 * np.ones(len(self.statespace)-1)
        diags=  np.delete(diags, finalState, 0)
        vals+= list(diags)
        rowtemp = [i for i in range( len(self.statespace) -1 )] 
        rows += rowtemp  
        cols += rowtemp
        warnings.filterwarnings('error')
        try: 
            rate_matrix_coo = coo_matrix((vals, (rows,cols)), shape=(len(self.statespace) -1, len(self.statespace) -1 ) , dtype=np.float64)
            rate_matrix_csr = csr_matrix(rate_matrix_coo)
            firstpassagetimes = spsolve(rate_matrix_csr  , b   )
        
        except RuntimeWarning as  w: 
            s = str(w) 
            if 'overflow' in s : 
                print "Overflow warning :( "
                return RETURN_MINUS_INF
            if 'underflow' in s : 
                print "Underflow warning :( "
                return 2.2250738585072014e-308
            return RETURN_MINUS_INF
        
        except Exception as w :
            print "Singular matrix exception :( "
            return RETURN_MINUS_INF
        except : 
            print  "Exception - Don't know what happend  :( "
            return RETURN_MINUS_INF
        
        if initialState > finalState : 
            firstpassagetime = firstpassagetimes[initialState-1]
        else : 
            firstpassagetime = firstpassagetimes[initialState]
        return firstpassagetime
        
        
    def find_answers(self, concentration, real_rate, bimolTransition):
        """ Computes the error """
        firstpassagetime = self.find_firstpassagetime()
        if firstpassagetime == RETURN_MINUS_INF  or firstpassagetime <= 0 :
            #firstpassagetime should be greater then 0!
            return  [  np.inf, np.inf, np.inf, self.local_context_uni, self.local_context_bi]
        
        #Estimating reaction rate constant from first passage time.
        if bimolTransition == True :
            predicted_rate= 1.0 / (firstpassagetime * concentration)
        else : 
            predicted_rate= 1.0 / firstpassagetime
        warnings.filterwarnings('error')
        try :
            predicted_log_10_rate =np.log(predicted_rate)/np.log(10)
            real_log_10_rate = np.log(real_rate)/np.log(10)
            error  = math.pow( real_log_10_rate - predicted_log_10_rate, 2)
        except :
            print " Exception occurred :( - why? "
            return [   np.inf,   np.inf, np.inf, self.local_context_uni, self.local_context_bi]
        gc.collect()
        return  [   error ,  predicted_log_10_rate , real_log_10_rate , self.local_context_uni, self.local_context_bi]

     
     
