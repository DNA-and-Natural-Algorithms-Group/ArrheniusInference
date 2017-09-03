## Introduction 

This is a  Python framework  for learning the parameters of the  <a href="https://link.springer.com/chapter/10.1007/978-3-319-66799-7_12">Arrhenius model</a>, an elementary step
model of DNA structure kinetics with locally context-dependent Arrhenius rates [1].   Additonally,  it  can infer parameters for the Metropolis model implemented  in Multistrand [2]. 

Our software employs a reduced state space that is a strict subset of the full Multistrand state space, enabling efficient sparse matrix computations of mean first passage times, from which reaction rate constants are predicted.  

Our software uses two approaches to  train the aforementioned DNA kinetic models: 
- A maximum a prior (MAP) approach that optimizes a single set of pararameters. 

- A  Markov chain Monte Carlo (MCMC) approach that learns an ensemble of parameter sets. 

Fore more information regarding our software and the methods employed, see <a href="https://link.springer.com/chapter/10.1007/978-3-319-66799-7_12">our paper</a> [1]. Also, see our **online appendix** that can be found on this webpage. 

Please contact us at nasimzf@cs.ubc.ca for questions regarding the software. 


Disclaimer: This software is free for use only for research purposes.  **If you make use of this software, please cite <a href="https://link.springer.com/chapter/10.1007/978-3-319-66799-7_12">our paper</a>**:



[1] Zolaktaf S, Dannenberg F, Rudelis X, Condon A, Schaeffer JM, Schmidt M, Thachuk C, Winfree E. Inferring Parameters for an Elementary Step Model of DNA Structure Kinetics with Locally Context-Dependent Arrhenius Rates. InInternational Conference on DNA-Based Computers 2017 Sep 24 (pp. 172-187). Springer, Cham.

## Requirements 


| Dependency       | Notes          
| ------------- |:-------------:|
|Python2 |2.7.12 |
|<a href="http://www.nupack.org/">NUPACK</a>  [3]    | 3.0.4 |
|numpy|1.12.1+| 
|SciPy | 0.19.0+     |
| enum34 | 1.1.6+ | 
|statistics| 1.0.3.5+|
|matplotlib| 2.0.0+|
|seaborn|0.7.1+|
| ConfigParser | 3.5.0+ | 
| <a href="http://dan.iel.fm/emcee/current/user/pt/"> emcee</a> [4]| 2.2.1+      | 

## Setup 
- Download and install the requirements. 
- Set the system environment variables to point NUPACKHOME to the folder where NUPACK is installed.  For example, in openSUSE, you can use the following command: 
`$ export NUPACKHOME=path/to/NUPACK/Folder`
- Clone our software directory into your workspace. 
 
## Software File Description

Our software contains three directories, namely, learndnakinetics, reactions, and dataset. Each directory ontains the following files:
	
- **learndnakinetics**:
  * **map.py**: this file runs the maximum a priori approach with the <a href="https://docs.scipy.org/doc/scipy/reference/optimize.minimize-neldermead.html"> Nelder-Mead optimization technique </a>.
   * **mcmcensemble.py**: this file runs  the MCMC approach with the <a href="http://dan.iel.fm/emcee/current/user/pt/"> emcee software package </a> [4].  
  * **learndnakinetics.py**: this file computes the objective function. 
  * **plot.py**: this file draws experimental plots of the dataset consistent with the literature they were derived from. 
  * **loadmcmc.py**: this file draws boxplots and correlation plots for the samples obtained from the MCMC ensemble method. 
  * **myenums.py**: this file contains enum types used in the framework. 
  * **configfile.txt**: this file contains variables for configularation. 

- **reactions**: 
   * **parent.py**: this file contains code that different types of reactions have in common.
   * **hairpin.py**: this file contains code specific to hairpin closing and opening.
   * **bubble.py**: this file contains code specific to bubble closing.
   * **helix.py**: this file contains code specific to helix association and dissociation.
   * **three_waystranddisplacement**: this file contains code spefic to three-way strand displacement. 
   * **four_waystrandexchange**: this file contains code specific to four-way strand exchange. 

- **dataset**: this directory contains reaction rate constants and timescales that we compiled from published literature
    * **hairpin**: this directory contains hairpin closing and opening experiments from Bonnet, Grégoire, Oleg Krichevsky, and Albert Libchaber. "Kinetics of conformational fluctuations in DNA hairpin-loops." Proceedings of the National Academy of Sciences 95.15 (1998): 8602-8606.
    * **hairpin1**: this directory contains hairpin closing and opening experiments from  Bonnet, Gregoire. "Dynamics of DNA breathing and folding for molecular recognition and computation." (2000) 
    * **hairpin4**: this directory contains hairpin closing and opening experiments from  Kim, Jiho, et al. "The initial step of DNA hairpin folding: a kinetic analysis using fluorescence correlation spectroscopy." Nucleic acids research 34.9 (2006): 2516-2527.
    * **bubble**: this directory contains bubble closing experiments from Altan-Bonnet, G., Libchaber, A., Krichevsky, O.: Bubble dynamics in doublestranded DNA. Physical Review Letters 90, 138101 (2003)
    * **helix**: this directory contains helix association and dissociation experiments from Morrison, L.E., Stols, L.M.: Sensitive fluorescence-based thermodynamic and kinetic measurements of DNA hybridization in solution. Biochemistry 32, 3095–3104 (1993)
    * **helix1**: this directory contains helix association and dissociation experiments from Reynaldo, L.P., Vologodskii, A.V., Neri, B.P., Lyamichev, V.I.: The kinetics of oligonucleotide replacements. Journal of Molecular Biology 297, 511–520 (2000)
    * **three_waystranddisplacement**: this directory contains three-way strand displacement experiments from Zhang, D.Y., Winfree, E.: Control of DNA strand displacement kinetics using toehold exchange. Journal of the American Chemical Society 131, 17303-17314 (2009)
    * **three_waystranddisplacement1**: this directory contains three-way strand displacement experiments from  Reynaldo, L.P., Vologodskii, A.V., Neri, B.P., Lyamichev, V.I.: The kinetics of oligonucleotide replacements. Journal of Molecular Biology 297, 511–520 (2000)
    * **three_waystranddisplacement2**: this directory contains three-way strand displacement experiments from Machinek, R.R., Ouldridge, T.E., Haley, N.E., Bath, J., TurberField, A.J.: Programmable energy landscapes for kinetic control of DNA strand displacement. Nature Communications 5 (2014)
    * **four_waystrandexchange**: this directory contains four-way strand exchange experiments from Dabby, N.L.: Synthetic molecular machines for active self-assembly: prototype algorithms, designs, and experimental study. Ph.D. thesis, California Institute of Technology (2013). 

## Software Usage: 
	
 To execute the  MAP approach: 

  * In **config_file.txt**:
    * Set *rate_method* to be Arrhenius or Metropolis.
    * Set *parameter_folder* to be the path to a  directory to save results. In this folder, the parameter set of each iteration of the optimization will be saved.
    * Set *n_processors* to be the number of processors  for multiprocessing the computation of the objective function.  
    * Set *NUPACK_bin* to  be the path to the bin folder in NUPACK.
  * Run **map.py**. The initial parameter set and other other optimization configurations can be changed in this file.

To execute the MCMC ensemble approach: 
  * In **config_file.txt**:
  
    * Set *rate_method* to be Arrhenius or Metropolis.
    * Set *parameter_folder*  to be the path to a directory to save results.  In this folder,  the sample parameter sets of each iteration of the optimization will be saved. Additionally,  the samples are saved in pickle format (after each walker has taken N_STEPS steps). 
    * Set *n_processors* to be the number of processors  for mutliprocessing the computation of the objective function.  
    * Set *NUPACK_bin* to be the path to the bin folder in NUPACK.
    * Set *N_WALKERS*  to be the number of  walkers. This should be even and at least twice the numbre of parameters. 
    * Set *N_STEPS* to be the number of steps every walker has to take (after emcee has started!) until samples  are saved  in pkl format.  But note that the first few iterations are not saved (see what the N_INITIAL_WALKERS does for explanation)!
    * Set *FINAL_STEP* to be an integer, such that every walker will take *N_STEPS* * *FINAL_STEP* iterations in total.
    * Set *N_INITIAL_WALKERS*  to be an integer.  An initial parameter set might lead to lnprobability of -np.inf ( equivalently probability 0 or error npinf) and the walker might not be able to update to a higher probability during the next iterations. Therefore, only parameter sets which do not have an lnprobability of -np.inf are chosen from a larger pool to run emcee with!  
 
   * Run **mcmcensemble.py**.    The initial parameter sets and other sampling configurations can be changed in this file. 

To reproduce experimental plots  of the literature: 

  * In **config_file.txt**: 
    * Set *rate_method* to be Arrhenius or Metropolis.
    * Set *learned_parameterset* to be one of the options from myenums.LearnedParameters, i.e., Arrheniusinitial, ArrheniusMAP, ArrheniusMCMCmode, ArrheniusMCMC, Metropolisinitial,  MetropolisMAP,  MetropolisMCMCmode,  MetropolisMCMC. 
    * If *learned_parameterset* is equal to ArrheniusMCMC or MetropolisMCMC, set *pkl_file_plot* to be the path to a pkl file  from the MCMC ensemble method. Note that rate_method, learned_parameterset, and the pkl file should  either all correspond to the Arrhenius model or all correspond to the Metropolis model. 
    * Set *n_processors* to be the number of processors  for paralyzing the computation of the objective function.  
    * Set *n_processors_plot* to be the number of processors for multiprocessing the ploting. 
    * Set *NUPACK_bin* to be the path to the bin folder in NUPACK. 
  * Run **plot.py**.  You can add  other parameter sets by editing the plot_rate function in learndnakinetics.py and consequently  adding an option to myenums.LearnedParameters.
	
 To draw boxplots and correlation plots for the samples from the MCMC ensemble method: 
  * In **config_file.txt**:
    * Set *rate_method* to be *Arrhenius* or *Metropolis*.
    *  Set *pkl_file_load* to be the path to a pkl file to read from the MCMC ensemble method. Note that rate_Method should be consistent with the pkl file. 
    * Set *parameter_folder* to be the path to a directory to save results.  
  * Run **loadmcmc.py**. 
		
	
This software has been tested on system with 16 2.93GHz Intel Xeon processors and 64GB RAM, running open-SUSE 42.1. On this system, the first iteration takes approximately 15 minutes and reusable files are saved. After that, each iteration takes less than 6 s.


## References 

[1] Zolaktaf S, Dannenberg F, Rudelis X, Condon A, Schaeffer JM, Schmidt M, Thachuk C, Winfree E. Inferring Parameters for an Elementary Step Model of DNA Structure Kinetics with Locally Context-Dependent Arrhenius Rates. InInternational Conference on DNA-Based Computers 2017 Sep 24 (pp. 172-187). Springer, Cham.

[2] Schaeffer, J.M., Thachuk, C., Winfree, E.: Stochastic simulation of the kinetics of multiple interacting nucleic acid strands. In: Proceedings of the 21st International Conference on DNA Computing and Molecular Programming-Volume 9211 (2015)


[3] Zadeh, J.N., Steenberg, C.D., Bois, J.S., Wolfe, B.R., Pierce, M.B., Khan, A.R., Dirks, R.M., Pierce, N.A.: NUPACK: analysis and design of nucleic acid systems.
Journal of Computational Chemistry 32, 170{173 (2011)

[4]  Foreman-Mackey, D., Hogg, D.W., Lang, D., Goodman, J.: emcee: The MCMC hammer. Publications of the Astronomical Society of the Pacific 125, 306 (2013)

[5] Bonnet, Grégoire, Oleg Krichevsky, and Albert Libchaber. "Kinetics of conformational fluctuations in DNA hairpin-loops." Proceedings of the National Academy of Sciences 95.15 (1998): 8602-8606.
 
 [6]   Bonnet, Gregoire. "Dynamics of DNA breathing and folding for molecular recognition and computation." (2000) 
 
 [7] Kim, Jiho, et al. "The initial step of DNA hairpin folding: a kinetic analysis using fluorescence correlation spectroscopy." Nucleic acids research 34.9 (2006): 2516-2527.
  
  [8] Altan-Bonnet, G., Libchaber, A., Krichevsky, O.: Bubble dynamics in doublestranded DNA. Physical Review Letters 90, 138101 (2003)
  
  [9] Morrison, L.E., Stols, L.M.: Sensitive fluorescence-based thermodynamic and kinetic measurements of DNA hybridization in solution. Biochemistry 32, 3095–3104 (1993)
  
  [10] Reynaldo, L.P., Vologodskii, A.V., Neri, B.P., Lyamichev, V.I.: The kinetics of oligonucleotide replacements. Journal of Molecular Biology 297, 511–520 (2000)
  
  [11] Zhang, D.Y., Winfree, E.: Control of DNA strand displacement kinetics using toehold exchange. Journal of the American Chemical Society 131, 17303-17314 (2009)
  
  [12]  Reynaldo, L.P., Vologodskii, A.V., Neri, B.P., Lyamichev, V.I.: The kinetics of oligonucleotide replacements. Journal of Molecular Biology 297, 511–520 (2000)

[13]  Machinek, R.R., Ouldridge, T.E., Haley, N.E., Bath, J., TurberField, A.J.: Programmable energy landscapes for kinetic control of DNA strand displacement. Nature Communications 5 (2014)

[14] Dabby, N.L.: Synthetic molecular machines for active self-assembly: prototype algorithms, designs, and experimental study. Ph.D. thesis, California Institute of Technology (2013). 
