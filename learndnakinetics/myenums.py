from enum import Enum

class ModelName(Enum):
	""" DNA kinetic model"""
	ARRHENIUSMODELNAME = "Arrhenius"
	METROPOLISMODELNAME = "Metropolis"


class Permanent_Folder(Enum):
	"""used for saving computations to prevent repeatable computations (leads to < 6s speed of each iteration after files are produced)"""
	TRANSITION_STRUCTURE = "transition_structure"
	FAST_ACCESS = "fast_access"
	STATESPACE = "statespace"
	PSD = "PSD"
	ENERGY = "energy"

class MethodName(Enum):
	"""  Specifies whether the freamework is using the MCMC ensemble approach or the MAP approach """
	MCMCNAME = "MCMCENSEMBLE"
	MAPNAME = "MAP"
	
class SetName (Enum) :
    TRAIN =  "Train"
    TEST = "Test"

class LearnedParameters(Enum):
	""" Used for plotting options in config_file.txt. See the learndnakinetics.plot_rates function"""
	ARRHENIUSINITIAL= [ "Arrheniusinitial", "Arrhenius with Initial"]
	ARRHENIUSMAP = [ "ArrheniusMAP", "Arrhenius with MAP"]
	ARRHENIUSMCMCMODE =  [ "ArrheniusMCMCmode", "Arrhenius with MCMC Mode"]
	ARRHENIUSMCMC = ["ArrheniusMCMC", "Arrhenius with MCMC Ensemble"]
	METROPOLISINITIAL = ["Metropolisinitial", "Metropolis with Initial"]
	METROPOLISMAP = ["MetropolisMAP", "Metropolis with MAP"]
	METROPOLISMCMCMODE = ["MetropolisMCMCmode", "Metropolis with MCMC Mode"]
	METROPOLISMCMC = ["MetropolisMCMC", "Metropolis with MCMC Ensemble"]

class DatasetName(Enum):
	""" Dataset name"""
	ALTANBONNET= "Altanbonnet"
	MORRISON = "Morrison"
	REYNALDODISSOCIATE = "ReynaldoDissociate"
	REYNALDOSEQUENTIAL= "ReynaldoSequential"
	ZHANG= "Zhang"
	MACHINEK= "Machinek"
	BONNET= "Bonnet"
	GODDARD= "Goddard"
	KIM ="Kim"
	DABBY= "Dabby"