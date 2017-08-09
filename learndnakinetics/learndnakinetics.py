from __future__ import division
import ConfigParser
import numpy as np
from scipy.optimize import minimize
import csv
import cPickle as pickle
import timeit
import os
import multiprocessing
import sys
import math
import shutil
sys.path.insert(0,os.path.realpath('../reactions'))
import parent
import hairpin
import helix
import bubble
import three_waystranddisplacement
import four_waystrandexchange
import myenums

DATASET_PATH = '../dataset'
PATH_AUXILIARY= "simplifiedstatespace"

use_all_data = False
for_plot= False
iter = 0



class ForMultiProcess(object):
	"""This class used for multiprocessing"""
	def __init__(self, function_name, arguments) :
		self.function_name = function_name
		self.arguments = arguments

def open_document(document) :
	"""open a csv file"""
	my_CSV = list(csv.reader(open(document, 'rb')))
	return my_CSV

#Note that for each dataset a separate function is used, started with read_, since the data sets have different fields!
def read_DabbyThesis(ss , counter_cell , document,  theta , done_queue , row, dataset_name , docID, name ) :
	docID = name +docID
	[  error ,  predicted_log_10_rate, real_log_10_rate ,  stuctureCounterUniLocal, half_context_biLocal] = four_waystrandexchange.main(ss,  float( row[counter_cell][8]) ,   float( row[counter_cell][13])  ,  int(row[counter_cell][1]) , int(row[counter_cell][2]) ,  row[counter_cell][3]  ,  row[counter_cell][4]  , row[counter_cell][5] ,  6,6 , theta,  1000/  float (row[counter_cell][6] )- 273.15 , np.max ( ( float (row[counter_cell][16] ), float (row[counter_cell][17]  ) ) ) ,  float (row[counter_cell][11]) ,float (row[counter_cell][12])  ,  dataset_name, docID , name)
	done_queue.put( ( name , error , counter_cell, document,  predicted_log_10_rate, real_log_10_rate , stuctureCounterUniLocal, half_context_biLocal )  )
def read_Machinek  ( ss, counter_cell , document,  theta , done_queue , row, dataset_name , docID, name) :
	docID = name + docID
	real_log_10_rate = float( row[counter_cell][9])
	[  error ,  predicted_log_10_rate, real_log_10_rate ,   stuctureCounterUniLocal, half_context_biLocal] = three_waystranddisplacement.main(ss, True,        row[counter_cell][3] , real_log_10_rate, int(row[counter_cell][1]) ,    "" , "",    theta,  float ( row[counter_cell][7]) , np.max(( float (row[counter_cell][16]) ,  float (row[counter_cell][17])   )),  float(row[counter_cell][13])  , float (row[counter_cell][14]), "" , dataset_name,  docID , name, row[counter_cell][4][16:] ,  row[counter_cell][5],  row[counter_cell][6], int( row[counter_cell][2]))
	done_queue.put( ( name , error ,counter_cell, document,  predicted_log_10_rate, real_log_10_rate, stuctureCounterUniLocal, half_context_biLocal  )  )

def read_Zhang(ss, counter_cell , document,   theta , done_queue , row,dataset_name  , docID, name) :
	docID = name + docID
	real_log_10_rate = math.pow(10, float( row[counter_cell][7])  )
	[  error ,  predicted_log_10_rate, real_log_10_rate ,  stuctureCounterUniLocal, half_context_biLocal] = three_waystranddisplacement.main(ss, True, row[counter_cell][2],  real_log_10_rate,  int ( row[counter_cell][1]   )  ,row[counter_cell][3], row[counter_cell][4], theta,  1000/ float (row[counter_cell][5]) - 273.15 , np.max ( ( float (row[counter_cell][16] ), float (row[counter_cell][17]  ) ) )   ,  float (row[counter_cell][9])  , float (row[counter_cell][10]) , "" ,  dataset_name,  docID , name,  "" , "", "", ""  )
	done_queue.put( ( name , error ,counter_cell, document,  predicted_log_10_rate, real_log_10_rate, stuctureCounterUniLocal, half_context_biLocal  )  )

def read_ReyanldoSequential(ss, counter_cell , document, theta , done_queue , row, dataset_name , docID, name ) :
	docID = name +docID
	real_log_10_rate =float( row[counter_cell][5])
	[  error ,  predicted_log_10_rate, real_log_10_rate ,   stuctureCounterUniLocal, half_context_biLocal] = three_waystranddisplacement.main(ss, False ,"" ,   real_log_10_rate, 0 ,row[counter_cell][2] ,""   , theta,  float( row[counter_cell][3])  , np.max ( ( float (row[counter_cell][16] ), float (row[counter_cell][17]  ) ) )  ,  float (row[counter_cell][7])  ,  float( row[counter_cell][8]) , row[counter_cell][9], dataset_name,  docID , name ,   "" , "", "", ""  )
	done_queue.put( ( name , error , counter_cell, document ,  predicted_log_10_rate, real_log_10_rate   , stuctureCounterUniLocal, half_context_biLocal)  )

def read_AltanBonnet(ss, counter_cell , document,  theta , done_queue,row, dataset_name, docID, name) :
	docID = name + docID
	flurPosition = 17
	real_log_10_rate = 1 / float( row[counter_cell][5])
	[  error ,  predicted_log_10_rate, real_log_10_rate , stuctureCounterUniLocal, half_context_biLocal] = bubble.main(ss,  real_log_10_rate , theta, row[counter_cell][1].rstrip(),row[counter_cell][2].rstrip(),row[counter_cell][3].rstrip(), (1000/ float (row[counter_cell][4] ))-273.15,  float (row[counter_cell][8] ),  float (row[counter_cell][9] ), 0, flurPosition, dataset_name, docID )
	done_queue.put( ( name , error  , counter_cell, document ,  predicted_log_10_rate, real_log_10_rate  , stuctureCounterUniLocal, half_context_biLocal)  )


def read_Morrison(ss, counter_cell , document, theta , done_queue, _zip, row, dataset_name , docID, name ) :
	docID = name + docID
	[  error ,  predicted_log_10_rate, real_log_10_rate,  stuctureCounterUniLocal, half_context_biLocal] = helix.main(ss,  math.pow(10, float (row[counter_cell][5] )) , theta, row[counter_cell][1].rstrip(), _zip, 1000/  float (row[counter_cell][3] ) - 273.15,  np.max ( ( float (row[counter_cell][16] ), float (row[counter_cell][17]  ) ) ) ,  float (row[counter_cell][8] ), 0,  "" , dataset_name, docID  , name  )
	done_queue.put( ( name , error  , counter_cell, document,    predicted_log_10_rate, real_log_10_rate, stuctureCounterUniLocal, half_context_biLocal )  )
	
def read_ReynaldoDissociate(ss, counter_cell , document, theta , done_queue, _zip, row, dataset_name , docID, name) :
	docID = name + docID
	[  error ,  predicted_log_10_rate, real_log_10_rate,stuctureCounterUniLocal, half_context_biLocal] = helix.main(ss,  float( row[counter_cell][5] )  ,  theta, row[counter_cell][2].rstrip(), _zip,  float (row[counter_cell][3] ),  np.max ( ( float (row[counter_cell][16] ), float (row[counter_cell][17]  ) ) ) ,  float (row[counter_cell][7] ),  float (row[counter_cell][8] ),row[counter_cell][9] , dataset_name, docID , name   )
	done_queue.put( ( name , error  , counter_cell, document,    predicted_log_10_rate, real_log_10_rate , stuctureCounterUniLocal, half_context_biLocal)  )

def read_Bonnet(ss, counter_cell , document,theta , done_queue, _zip , row, dataset_name, docID, name ):
	docID = name +docID
	magnesium = 0
	[  error ,  predicted_log_10_rate, real_log_10_rate,  stuctureCounterUniLocal, half_context_biLocal] =   hairpin.main(ss, float (row[counter_cell][5]) ,  theta, row[counter_cell][1].rstrip(),row[counter_cell][2].rstrip(), _zip, 1000/  float( row[counter_cell][3] )- 273.15,  float ( row[counter_cell][7] ) ,  float ( row[counter_cell][8] ) , magnesium , dataset_name, docID )
	done_queue.put( ( name , error  , counter_cell, document,    predicted_log_10_rate, real_log_10_rate, stuctureCounterUniLocal, half_context_biLocal )  )


def read_Kim(ss, counter_cell , document,theta , done_queue, _zip , row, dataset_name, docID  ,name):
	docID = name +docID
	magnesium = 0
	[  error ,  predicted_log_10_rate, real_log_10_rate,  stuctureCounterUniLocal, half_context_biLocal] =   hairpin.main(ss, float (row[counter_cell][5]) ,  theta, row[counter_cell][1].rstrip(),row[counter_cell][2].rstrip(), _zip, 1000/  float( row[counter_cell][3] )- 273.15,  float ( row[counter_cell][7] ) ,  float ( row[counter_cell][8] ) , magnesium , dataset_name, docID )
	done_queue.put( ( name , error  , counter_cell, document,    predicted_log_10_rate, real_log_10_rate, stuctureCounterUniLocal, half_context_biLocal )  )

def read_BonnetThesis(ss, counter_cell , document,  theta , done_queue, _zip,row,  dataset_name, docID , name):
	docID = name +docID
	real_log_10_rate = 1 / float( row[counter_cell][4])
	[  error ,  predicted_log_10_rate, real_log_10_rate,  stuctureCounterUniLocal, half_context_biLocal] =   hairpin.main(ss,  real_log_10_rate,  theta, row[counter_cell][1].rstrip(),row[counter_cell][2].rstrip(), _zip, 1000/  float (row[counter_cell][3] ) - 273.15,   float (row[counter_cell][7] ),  float (row[counter_cell][8] ), 0 , dataset_name, docID  )
	done_queue.put( ( name , error ,counter_cell, document,  predicted_log_10_rate, real_log_10_rate, stuctureCounterUniLocal, half_context_biLocal )  )

def multi_process(done_queue ,  dataset_list ,  iter , countS , local_context_uni, local_context_bi) :
	"""multi processing function """
	global   predicted_logreactionrateconstants , experimental_logreactionrateconstants
	error = 0
	pool = multiprocessing.Pool( processes = n_processors)
	for ds in dataset_list:
		compile_error = pool.apply_async( ds.function_name ,  ds.arguments )
	#print "Errors: " + str(compile_error.get())
	pool.close( )
	pool.join ()

	while not done_queue.empty():

		(name, s  ,  counter_cell, document  , predictedRate ,real_log_10_rate ,  local_context_uni_l, local_context_bi_l )= done_queue.get()
		if for_plot== True :
			predicted_logreactionrateconstants[iter, counter_cell, document ]  = predictedRate
			experimental_logreactionrateconstants [ iter, counter_cell, document ]  = real_log_10_rate
		error += s
		if iter == 0 and  parent.rate_method == myenums.ModelName.ARRHENIUSMODELNAME.value:
			for i in local_context_uni:
				local_context_uni [i] += local_context_uni_l[i]
			for i in local_context_bi:
				local_context_bi[i] += local_context_bi_l[i]
		
		if name in countS :
			countS [name ] += s
		else :
			countS[name ] = s
	return error
def check_directories (directories) :
	for dir in directories:
		if not os.path.exists(dir):
			os.makedirs(dir)
def objective_function(thetap):
	"""For the MCMC approach, receives an  parameter set and returns an approximation of the log  posterior. For the MAP approach it returns an approximation of the negative log posterior"""
	global iter
	start_time = timeit.default_timer()
	theta =[]
	for x in thetap :
		theta.append(x)
	if parent.rate_method == myenums.ModelName.ARRHENIUSMODELNAME.value:
		theta = [thetap[0] , thetap[1] , thetap[2], thetap[3] , thetap[4] , thetap[5], thetap[6] ,thetap[7], thetap[8], thetap[9], thetap[10] , thetap[11],  thetap[12] , thetap[13], thetap[14] ]
		alpha = theta [14]
	elif parent.rate_method == myenums.ModelName.METROPOLISMODELNAME.value :
		theta = [thetap[0] , thetap[1]]
		alpha =1
	else:
		raise ValueError('Error: Please specify rate_method to be Arrhenius or Metropolis!')
	sigma = thetap[len(thetap)-1]

	if alpha  <= 0  or sigma <= 0  or (parent.rate_method ==myenums.ModelName.METROPOLISMODELNAME.value and  ( theta[0] <= 0 or theta[1]  <= 0 )  ) :
		if METHOD == myenums.MethodName.MCMCNAME.value:
			return -np.inf
		elif METHOD ==myenums.MethodName.MAPNAME.value:
			return np.inf
	parameter_file = open(parameter_file_name, 'a')
	parameter_file.write("Iteration " + str(iter) +" "+str(theta) +  " " + str(sigma) + '\n')
	error = 0
	n = 0
	done_queue =  multiprocessing.Manager().Queue()
	dataset_list = []
	directories =[]
 
	
	if use_all_data ==  False :
		set = [myenums.SetName.TRAIN.value]
	elif  use_all_data == True :
		set = [myenums.SetName.TRAIN.value, myenums.SetName.TEST.value]

	# Zhang
	my_name = '/three_waystranddisplacement/Fig3b'
	dataset_name, document, row = initconf(my_name, directories)
	for set_type in set:
		for counter_cell in traintestset[document, set_type]:
			ss = {myenums.Permanent_Folder.PSD.value:dict() , myenums.Permanent_Folder.TRANSITION_STRUCTURE.value:dict( )}

			dataset_list.append(ForMultiProcess(read_Zhang, ( ss, counter_cell, document,theta, done_queue, row, dataset_name , str(counter_cell) ,  myenums.DatasetName.ZHANG.value  )  ))
			n +=1
	
	#Dabby

	my_name= '/four_waystrandexchange/Table5.2'
	dataset_name, document , row = initconf(my_name , directories)

	for set_type in set:
		for counter_cell in traintestset[document, set_type]:
			ss = {myenums.Permanent_Folder.PSD.value:dict() , myenums.Permanent_Folder.TRANSITION_STRUCTURE.value:dict( )}
			dataset_list.append(ForMultiProcess(read_DabbyThesis ,(ss,  counter_cell , document, theta, done_queue, row , dataset_name ,  str(counter_cell) , myenums.DatasetName.DABBY.value )))
			n +=1
	 
	#Reynaldo
	my_name = '/three_waystranddisplacement1/Fig6b'
	dataset_name, document , row = initconf(my_name , directories)

	for set_type in set:
		for counter_cell in traintestset[document, set_type]:
			ss = {myenums.Permanent_Folder.PSD.value:dict() , myenums.Permanent_Folder.TRANSITION_STRUCTURE.value:dict( )}
			dataset_list.append(ForMultiProcess(read_ReyanldoSequential, (ss, counter_cell , document,theta, done_queue, row, dataset_name, str(counter_cell) ,  myenums.DatasetName.REYNALDOSEQUENTIAL.value)  ))
			n +=1
	  

	#ReynaldoDissociate
	my_name = '/helix1/Fig6a'
	dataset_name, document , row = initconf(my_name , directories)
	for _zip in [False]:
		for set_type in set:
			for counter_cell  in traintestset [document, set_type  ] :
				ss = {myenums.Permanent_Folder.PSD.value:dict() , myenums.Permanent_Folder.TRANSITION_STRUCTURE.value:dict( )}
				dataset_list.append( ForMultiProcess( read_ReynaldoDissociate, (ss , counter_cell , document, theta , done_queue, _zip , row , dataset_name ,str(_zip) + str(counter_cell) ,  myenums.DatasetName.REYNALDODISSOCIATE.value)))
				n +=1

	#Morrison
	for _zip in [True , False ]:

		my_name = '/helix/Fig6_' +str(int(_zip))
		dataset_name, document , row = initconf(my_name , directories)
		for set_type in set:
			for counter_cell in traintestset[document, set_type]:
				ss = {myenums.Permanent_Folder.PSD.value:dict() , myenums.Permanent_Folder.TRANSITION_STRUCTURE.value:dict( )}
				dataset_list.append( ForMultiProcess( read_Morrison, (ss, counter_cell , document, theta , done_queue, _zip , row, dataset_name , str(_zip) + str(counter_cell) , myenums.DatasetName.MORRISON.value)))
				n +=1

	#AltanBonnet

	my_name= '/bubble/Fig4'
	dataset_name, document , row = initconf(my_name , directories)

	for set_type in set:
		for counter_cell in traintestset[document, set_type]:
			ss = {myenums.Permanent_Folder.PSD.value:dict() , myenums.Permanent_Folder.TRANSITION_STRUCTURE.value:dict( )}
			dataset_list.append(ForMultiProcess( read_AltanBonnet ,(ss, counter_cell , document, theta, done_queue, row ,dataset_name , str(counter_cell) , myenums.DatasetName.ALTANBONNET.value)))
			n +=1
	#Bonnet

	for j in [4,6] :
		for _zip in [True , False ]:
			my_name = '/hairpin/Fig'+str(j)  + '_' + str(int(_zip))
			dataset_name, document , row = initconf(my_name , directories)

			for set_type in set:
				for counter_cell in traintestset[document, set_type]:
					docID = str(j) + str(_zip) + str(counter_cell )
					ss = {myenums.Permanent_Folder.PSD.value:dict() , myenums.Permanent_Folder.TRANSITION_STRUCTURE.value:dict( )}
					dataset_list.append(ForMultiProcess(read_Bonnet, (ss, counter_cell , document , theta, done_queue, _zip, row,  dataset_name, docID  , myenums.DatasetName.BONNET.value )))
					n +=1

	#Goddard
	for j in ["T" ]:
		for _zip in [True, False  ]:
			my_name = '/hairpin1/Fig3_'+str(j)  + '_' + str(int(_zip))
			dataset_name, document , row = initconf(my_name , directories)

			for set_type in set:
				for counter_cell in traintestset[document, set_type]:
					docID = str(j) + str(_zip) + str(counter_cell )
					ss = {myenums.Permanent_Folder.PSD.value:dict() , myenums.Permanent_Folder.TRANSITION_STRUCTURE.value:dict( )}
					dataset_list.append(ForMultiProcess(read_BonnetThesis, (ss, counter_cell , document , theta, done_queue, _zip, row , dataset_name, docID  , myenums.DatasetName.GODDARD.value )))
					n +=1
	#Kim

	for _zip in [True , False ]:

		my_name = '/hairpin4/Table1_' + str(int(_zip))
		dataset_name, document , row = initconf(my_name , directories)

		for set_type in set:
			for counter_cell in traintestset[document, set_type]:
				docID = str(_zip) + str(counter_cell )
				ss = {myenums.Permanent_Folder.PSD.value:dict() , myenums.Permanent_Folder.TRANSITION_STRUCTURE.value:dict( )}
				dataset_list.append(ForMultiProcess(read_Kim, (ss, counter_cell , document, theta, done_queue, _zip, row,  dataset_name, docID, myenums.DatasetName.KIM.value  )))
				n +=1


	#machinek
	my_name = '/three_waystranddisplacement2/Fig2'
	dataset_name, document , row = initconf(my_name , directories)

	for set_type in set:
		for counter_cell in traintestset[document, set_type]:
			ss = {myenums.Permanent_Folder.PSD.value:dict() , myenums.Permanent_Folder.TRANSITION_STRUCTURE.value:dict( )}
			dataset_list.append(ForMultiProcess(read_Machinek ,( ss, counter_cell, document,theta, done_queue, row, dataset_name, str(counter_cell) ,  myenums.DatasetName.MACHINEK.value)  ))
			n +=1
	
	#Kim
	for _zip in [True , False ]:
		for j in ["a", "b"]:
			my_name = '/hairpin4/Fig5'+ str(j) + "_" + str(int(_zip))
			dataset_name, document , row = initconf(my_name , directories)

			for set_type in set:
				for counter_cell in traintestset[document, set_type]:
					docID = str(_zip) +"fig5"+j+ str(counter_cell )
					ss = {myenums.Permanent_Folder.PSD.value:dict() , myenums.Permanent_Folder.TRANSITION_STRUCTURE.value:dict( )}
					dataset_list.append(ForMultiProcess(read_Kim, (ss, counter_cell , document, theta, done_queue, _zip, row,  dataset_name, docID , myenums.DatasetName.KIM.value )))
					n +=1


	check_directories (directories)
	countS = dict()
	local_context_bi = dict()
	local_context_uni = dict()
	parameter_name = ("stack", "loop", "end", "stack+loop", "stack+end", "loop+end", "stack+stack")
	for i in parameter_name:
		for j in parameter_name:
			local_context_bi[i, j] = 0
			local_context_uni[i, j] = 0
	error += multi_process(done_queue , dataset_list   , iter , countS , local_context_uni, local_context_bi  )
	if iter  ==0  and  parent.rate_method == myenums.ModelName.ARRHENIUSMODELNAME.value:
		output = open('local_context_bi.pkl', 'wb')
		pickle.dump(local_context_bi, output)
		output = open('local_context_uni.pkl', 'wb')
		pickle.dump(local_context_uni, output)


	regularizer = 0
	if parent.rate_method == myenums.ModelName.ARRHENIUSMODELNAME.value :
		for i in range( len(theta) ) :
			regularizer+=  (theta[i] * theta [i] )
	elif parent.rate_method == myenums.ModelName.METROPOLISMODELNAME.value:
		for i in range(0, 2 ) :
			param = math.log (theta[i] )
			regularizer += (  param * param )
		for i in range( 2, len(theta) ):
			regularizer += (theta[i] * theta[i] )

	else:
		raise ValueError('Error: Specify rate_method to be Arrhenius or Metropolis!')
	LAMBDA =50
	regularizer = regularizer/ (2 * LAMBDA)
	lnprob = -(n + 1  )*np.log(sigma) - (error  /( 2 *(sigma ** 2) )) -  regularizer
	negativelnprob = -lnprob

	elapsed = timeit.default_timer() - start_time
	if for_plot== True :
		plot_probs [iter] = sampler.lnprobability[counti][countj]
	parameter_file.write( "Iteration:" + str(iter) +  "		,error:" + str(error) + "       lnprob:" + str(lnprob) + "      negativelnprob:" + str(negativelnprob) +  "		iteration time:" + str(elapsed) +  '\n')
	parameter_file.write(str(countS)  +  '\n\n')
	parameter_file.close()
	print "Iteration:" + str(iter) +  "		,error:" + str(error) + "       lnprob:" + str(lnprob) + "      negativelnprob:" + str(negativelnprob) +  "		iteration time:" + str(elapsed) +  '\n'
	iter += 1   # Do not move this line or you'll get errors later

	#if os.path.exists(parent.AUXILIARY_NUPACK_FILE):
	#    shutil.rmtree(parent.AUXILIARY_NUPACK_FILE)

	if np.isnan(error)  or error == np.inf:
		negativelnprob = np.inf
		lnprob = -np.inf
	
	
	if METHOD == myenums.MethodName.MCMCNAME.value :
		return lnprob
	if METHOD == myenums.MethodName.MAPNAME.value:
		return negativelnprob

def n_csv_rows(csv) :
	# return the number of rows in a csv file
	count = 0
	row = csv[count]
	while row [0] != '' :
		count += 1
		if count >= len (csv)  :
			break
		row= csv [count]
	return count

def initconf(my_name , directories) :
	#creating required directories
	dataset_name =  PATH_AUXILIARY +my_name
	document = DATASET_PATH + my_name + '.csv'
	directories +=  [dataset_name + "/"+  myenums.Permanent_Folder.PSD.value ,dataset_name+"/"+  myenums.Permanent_Folder.TRANSITION_STRUCTURE.value ,dataset_name +"/" + myenums.Permanent_Folder.STATESPACE.value  , dataset_name +"/"+ myenums.Permanent_Folder.FAST_ACCESS.value  ,  dataset_name +"/"+ myenums.Permanent_Folder.ENERGY.value  ]
	if not os.path.exists(dataset_name):
		os.makedirs(dataset_name)
	row =  open_document(document)
	return dataset_name, document , row

def set_traintestset_doc(document , trainortest ) :
	
	global traintestset
	row =  open_document(document)
	traintestset [document, trainortest] = [i for i in range( 1,  n_csv_rows(row))]
	if trainortest == myenums.SetName.TRAIN.value:
		inv = myenums.SetName.TEST.value
	elif trainortest == myenums.SetName.TEST.value:
		inv = myenums.SetName.TRAIN.value
	traintestset[document, inv ]  = [ ]


def set_traintestset():
	
	"""Split the dataset in to a training set and a testing set. Use  myenums.SetName.TRAIN.value for training set and  myenums.SetName.TEST.value for testing set."""
	global traintestset
	traintestset = dict()
	for _zip in [True, False ] :
		set_traintestset_doc(DATASET_PATH + '/helix/Fig6_' + str(int(_zip)) + '.csv', myenums.SetName.TRAIN.value)
		for j in [4,6] :
			set_traintestset_doc(DATASET_PATH + '/hairpin/Fig' + str(j) + '_' + str(int(_zip)) + '.csv', myenums.SetName.TRAIN.value)
	for _zip in [True , False ] :
		for j in ["T" ] :
			set_traintestset_doc(DATASET_PATH + '/hairpin1/Fig3_' + str(j) + '_' + str(int(_zip)) + '.csv', myenums.SetName.TRAIN.value)

	for _zip in [True , False ]:
		set_traintestset_doc(DATASET_PATH + '/hairpin4/Table1_' + str(int(_zip)) + '.csv', myenums.SetName.TRAIN.value)
		set_traintestset_doc(DATASET_PATH + '/hairpin4/Fig5a_' + str(int(_zip)) + '.csv', myenums.SetName.TEST.value)
		set_traintestset_doc(DATASET_PATH + '/hairpin4/Fig5b_' + str(int(_zip)) + '.csv', myenums.SetName.TEST.value)


	set_traintestset_doc(DATASET_PATH + '/three_waystranddisplacement/Fig3b.csv', myenums.SetName.TRAIN.value)
	set_traintestset_doc(DATASET_PATH + '/four_waystrandexchange/Table5.2.csv', myenums.SetName.TRAIN.value)
	set_traintestset_doc(DATASET_PATH + '/three_waystranddisplacement2/Fig2.csv', myenums.SetName.TEST.value)
	set_traintestset_doc(DATASET_PATH + '/three_waystranddisplacement1/Fig6b.csv', myenums.SetName.TRAIN.value)
	set_traintestset_doc(DATASET_PATH + '/helix1/Fig6a.csv', myenums.SetName.TRAIN.value)
	set_traintestset_doc(DATASET_PATH + '/bubble/Fig4.csv', myenums.SetName.TRAIN.value)

	return traintestset

def set_configuration():
	#setting some configurations!
	global parameter_file_name,  parameter_folder ,   n_processors
	set_traintestset()
	configParser = ConfigParser.ConfigParser()
	configParser.readfp(open(r'config_file.txt'))
	CONFIG_NAME = 'learndnakinetics'
	parent.rate_method= configParser.get(CONFIG_NAME, 'rate_method')
	n_processors =    configParser.getint(CONFIG_NAME, 'n_processors')
	parameter_folder = configParser.get(CONFIG_NAME, 'parameter_folder')
	if not os.path.exists(parameter_folder):
		os.makedirs(parameter_folder)
	check_directories([parent.AUXILIARY_NUPACK_FILE])
	parameter_file_name  = parameter_folder + "/parameter_file_name"

def plot_rates(file_name ,   learned_parameterset, use_only_finalstep) :
	#This function is only used for plotting in plot.py. It is not part of the training process!!!!!
	global  METHOD , sampler, for_plot, predicted_logreactionrateconstants, experimental_logreactionrateconstants, plot_probs, counti, countj
	predicted_logreactionrateconstants  = dict()
	experimental_logreactionrateconstants = dict()
	plot_probs = dict()
	for_plot= True
	f_in = open(file_name+".pkl", 'rb')
	ndim = pickle.load(f_in)
	nwalkers = pickle.load(f_in)
	nburn = pickle.load(f_in)
	nsteps = pickle.load(f_in)
	p0 = pickle.load(f_in)
	sampler = pickle.load(f_in)
	f_in.close()
	theta = []
	if learned_parameterset==   myenums.LearnedParameters.ARRHENIUSMCMC.value[0]  or learned_parameterset ==  myenums.LearnedParameters.METROPOLISMCMC.value[0]   :
		METHOD = myenums.MethodName.MCMCNAME.value
		if use_only_finalstep == True:
			steps= [sampler.chain.shape[1]-1]#Only using the last step of each walker!
		elif use_only_finalstep == False :
			steps = [i for i in range(sampler.chain.shape[1])] #use all steps of each walker! c
		for j in steps :
			for i in range(sampler.chain.shape[0]) :
				if sampler.lnprobability [i][j]!= -np.inf :
					th = sampler.chain[i][j]
					counti =  i
					countj = j
					theta.append(th)
				
			
		if learned_parameterset ==  myenums.LearnedParameters.ARRHENIUSMCMC.value[0] :
			title_name   =  myenums.LearnedParameters.ARRHENIUSMCMC.value[1]
		elif learned_parameterset == myenums.LearnedParameters.METROPOLISMCMC.value[0]:
			title_name = myenums.LearnedParameters.METROPOLISMCMC.value[1]
		else:
			raise ValueError('Error: Please specify learned_parameterset to be one of the options in myenums.LearnedParameters!')

	else :
		METHOD = myenums.MethodName.MAPNAME.value
		counti= 0
		countj = 0
	
		if learned_parameterset ==  myenums.LearnedParameters.ARRHENIUSINITIAL.value[0]  :
			#Initial parameter set for the Arrhenius model
			th = [  13.0580, 3,  13.0580, 3,   13.0580, 3,  13.0580 , 3,   13.0580, 3,  13.0580,  3,   13.0580 , 3,    0.0402 ]
			title_name = myenums.LearnedParameters.ARRHENIUSINITIAL.value[1]
		
		elif learned_parameterset == myenums.LearnedParameters.ARRHENIUSMAP.value[0]:
			#Learned parameter set with the MAP approach for the  Arrhenius model
			th = [10.700511073989023, 3.0406751829628016, 14.177641444707664, 3.3210958616087707, 12.960513664971495, 3.420869159636668, 11.88673466110987, 2.9827816021111917, 13.447865151543084, 3.2025632149181731, 14.716257115604998, 3.223036583523915, 13.791307834028169, 3.0974417518433972, 0.043497424113516377]
			title_name = myenums.LearnedParameters.ARRHENIUSMAP.value[1]
		 
		elif learned_parameterset == myenums.LearnedParameters.ARRHENIUSMCMCMODE.value[0] :
			#Parameter set with the highest probability on the training set  with MCMC ensemble method for the Arrhenius model
			th =  [1.41839430e+01,   5.28692038e+00,   1.64236969e+01,          4.46143369e+00,   1.29648159e+01,   3.49798154e+00,          5.81061725e+00,  -1.12763854e+00,   1.75235569e+01,          2.65589869e+00,   2.42237267e+00,   8.49339120e-02,          8.04573830e+00,  -6.27121400e-01,   1.60062641e-02]
			title_name =   myenums.LearnedParameters.ARRHENIUSMCMCMODE.value[1]
			
		
		
		elif learned_parameterset ==  myenums.LearnedParameters.METROPOLISINITIAL.value[0]:
			# Initial parameter set for the Metropolis model
			th = [8.2 *  (10 **6), 3.3  * (10**5) ]
			title_name = myenums.LearnedParameters.METROPOLISINITIAL.value[1]
		
		elif learned_parameterset == myenums.LearnedParameters.METROPOLISMAP.value[0]:
			# Learned parameter set with the MAP approach for the  Metropolis  model
			th =  [2430988.7336683525, 745530.95818480779]
			title_name =myenums.LearnedParameters.METROPOLISMAP.value[1]

		elif learned_parameterset == myenums.LearnedParameters.METROPOLISMCMCMODE.value[0] :
			# Parameter set with the highest probability on the training set  with MCMC ensemble method for the Metropolis model
			th = [ 2.41686715e+06,   8.01171383e+05]
			title_name = myenums.LearnedParameters.METROPOLISMCMCMODE.value[1]
	
		else :
			raise ValueError('Error: Please specify learned_parameterset to be one of the options in myenums.LearnedParameters!')
		theta.append(th)
		
	try:
		predicted_logreactionrateconstants2   = pickle.load( open( learned_parameterset+"/predicted_logreactionrateconstants2"  , "rb" ))
		experimental_logreactionrateconstants2   = pickle.load( open( learned_parameterset+"/experimental_logreactionrateconstants2"  , "rb" ))
		plot_probs2   = pickle.load( open( learned_parameterset+"/plot_probs2"  , "rb" ))
	except:
	
		toKeep = []
		for s  in  range (len(theta)) :
			th = theta [s]
		
	 
			overallerror = objective_function (th )
			if  overallerror !=  -np.inf  :
			 
				toKeep.append(s )
		predicted_logreactionrateconstants2= dict()
		experimental_logreactionrateconstants2 = dict()
		plot_probs2 = dict()
		for s in toKeep:
			th = theta [s]
			for i in predicted_logreactionrateconstants:
				if i[0] == s:
				 
					predicted_logreactionrateconstants2 [str(th)  ,i[1],i[2]] = predicted_logreactionrateconstants[i]
					experimental_logreactionrateconstants2 [str(th), i[1], i[2]] = experimental_logreactionrateconstants[i]
		   
			plot_probs2 [str(th)]  = plot_probs[s]

		#if not os.path.exists(learned_parameterset):
		#    os.makedirs(learned_parameterset)
		#pickle.dump(predicted_logreactionrateconstants2,open( learned_parameterset+"/predicted_logreactionrateconstants2", 'wb'))
		#pickle.dump(experimental_logreactionrateconstants2,open( learned_parameterset+"/experimental_logreactionrateconstants2", 'wb'))
		#pickle.dump(plot_probs2,open( learned_parameterset+"/plot_probs2", 'wb'))
	
	return traintestset,  predicted_logreactionrateconstants2 , experimental_logreactionrateconstants2, plot_probs2, title_name
	
