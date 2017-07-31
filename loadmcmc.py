from __future__ import division
import ConfigParser
import statistics
from scipy import stats
import scipy as sp  
import numpy as np
import math 
import os
import emcee
import cPickle as pickle
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  
from matplotlib.backends.backend_pdf import PdfPages
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns
import pandas as pd
import csv
import matplotlib.cm as cm
import matplotlib.colors as colors
import sys
from math import floor, log10
import learndnakinetics
import parent
from myenums import *
import plot

TEMPERATURE = 25
R = 0.001987
T = 273.15 + TEMPERATURE
ALPHA_DIM =14

def sci_notation(num, decimal_digits=1, precision=None, exponent=None):
    """     Returns a string representation of the scientific     notation of the given number formatted for use with     notation of the given number formatted for use with     LaTeX or Mathtext, with specified number of significant     decimal digits and precision (number of decimal digits     to show). The exponent to be used can also be specified  xplicitly."""
    if not exponent:
        exponent = int(floor(log10(abs(num))))
    coeff = round(num / float(10**exponent), decimal_digits)
    if not precision:
        precision = decimal_digits
    return r"${0:.{2}f}\cdot10^{{{1:d}}}$".format(coeff, exponent, precision)

def delete_unwanted_walkers(keep_walkers_list):
    global sampler_lnprobability, sampler_chain,  sampler_flatlnprobability,  sampler_flatchain , n_walkers
    n_walkers_before_prune, fj , fk  = sampler.chain.shape
    n_walkers_after_prune = len(keep_walkers_list)
    sampler_lnprobability = np.zeros ((len(keep_walkers_list) ,fj))
    sampler_chain = np.zeros ((len(keep_walkers_list) ,fj, fk)) 
    sampler_flatlnprobability = np.zeros ((len(keep_walkers_list)  * fj)) 
    sampler_flatchain = np.zeros ((len(keep_walkers_list)    * fj, fk)) 
    counter = 0 
    for i in keep_walkers_list : 
        sampler_lnprobability[counter] = sampler.lnprobability[i]      
        sampler_chain[counter] = sampler.chain[i] 
        counter +=1
    for ii in range ( n_walkers_after_prune)  :  
        for j in range( fj )   :
            i = keep_walkers_list [ii ]
            sampler_flatlnprobability [ii  * fj  + j ] = sampler.flatlnprobability[ i * fj  +  j  ] 
            sampler_flatchain [ ii *  fj  + j ] = sampler.flatchain[i * fj   + j]
    n_walkers = len (keep_walkers_list)

def use_only_finalstep_func () :
    
    # if use_only_finalstep == 1, only returns the last step of each walker as a sample. Else, returns all steps of each as a sample
    if plot.USE_ONLY_FINALSTEP == True : #When use_only_finalstep  = True , only  the last step of each walker is used!

        
        n_walkers_before_prune, fj , fk  = sampler_chain.shape
        sampler_lnprobabilityTemp = np.zeros ((n_walkers ,1)) 
        sampler_chainTemp = np.zeros ( (n_walkers ,1, fk)) 
        sampler_flatlnprobabilityTemp = np.zeros ((n_walkers * 1)) 
        sampler_flatchainTemp = np.zeros ((n_walkers   * 1, fk)) 
        for i in range(n_walkers) : 
            sampler_lnprobabilityTemp[i] = sampler_lnprobability[i][fj -1 ]
            sampler_chainTemp [i ] = sampler_chain[i ] [fj -1 ]
        for i in range (n_walkers)  :  
            for j in [fj -1 ]   :
                sampler_flatlnprobabilityTemp  [i ] = sampler_flatlnprobability[ i * fj    +  j ] 
                sampler_flatchainTemp [i ] = sampler_flatchain[ i * fj   + j ]
        n_stepsTemp = 1
        return sampler_lnprobabilityTemp, sampler_chainTemp,  sampler_flatlnprobabilityTemp,  sampler_flatchainTemp , n_stepsTemp
    else :
        return sampler_lnprobability, sampler_chain,  sampler_flatlnprobability,  sampler_flatchain , n_steps 


def load_MCMC(filename) :
    global   n_dim, n_walkers, nburn, n_steps, p0, sampler
    f_in = open(filename+".pkl", 'rb')
    n_dim = pickle.load(f_in)
    n_walkers = pickle.load(f_in)
    n_burn= pickle.load(f_in)
    n_steps = pickle.load(f_in)
    p0 = pickle.load(f_in)
    sampler = pickle.load(f_in)
    f_in.close()
    
    #delete walkers which have probability -np.inf in the last step (we only use the last step)
    keep_walkers_list = []
    for i in range(n_walkers) : 
        for j in [ n_steps - 1  ] : 
            if sampler.lnprobability[  i ][j] != -np.inf :  
                keep_walkers_list.append(i)
    delete_unwanted_walkers (keep_walkers_list)

def valuevsiteration( path_name, figsize ):
    """ plot the value of parameters of the kinetic model as a box_plot. """
    if not os.path.exists(path_name):
        os.makedirs(path_name) 
    sampler_lnprobability_finalstep, sampler_chain_finalstep,  sampler_flatlnprobability_finalstep,  sampler_flatchain_finalstep , n_steps_finalstep   = use_only_finalstep_func ()
    x_list = []
    y_list = []
    if parent.rate_method == ModelName.ARRHENIUSMODELNAME.value :
        box_plot_dim = 14
        parameter_name = ( "stack", "stack", "loop", "loop", "end", "end", "stack+loop", "stack+loop", "stack+end", "stack+end",         "loop+end", "loop+end", "stack+stack", "stack+stack", r"$\alpha$")
        csvfile = path_name + "arrheniusvalue.csv"
        file_name_save = "arrheniusvalue.pdf"
        x = "Half context"
        y = "Value" 
        hue = "Parameter Type"
        hue_list  = []
        for j in range(box_plot_dim) :
            if j%2 == 0  : 
                parameter_type = r"$\ln A_l$"
            else : 
                parameter_type = "$E_l$"
            for i in range(n_walkers):
                for k in range  ( n_steps_finalstep ) : 
                    x_list.append (parameter_name[j] )
                    y_list.append(sampler_chain_finalstep[i][k][j])
                    hue_list.append ( parameter_type)
        raw_data = {x: x_list, y: y_list, hue:hue_list }
        df  = pd.DataFrame(raw_data, columns = [x, y, hue ] )
    elif parent.rate_method == ModelName.METROPOLISMODELNAME.value :
        box_plot_dim = 2
        parameter_name = (r"$k_{\rm uni}$", r"$k_{\rm bi}$")
        csvfile = path_name + "metropolisvalue.csv"
        file_name_save = "metropolisvalueorrate.pdf"
        y = r"$\log_{10} k$"
        x = "Transition"
        for j in range( box_plot_dim):
            for i in range(n_walkers):
                for k in range  ( n_steps_finalstep ) : 
                    x_list.append (parameter_name[j] )
                    y_list.append(math.log10 ( sampler_chain_finalstep[i][k][j]) )
        raw_data = {x: x_list, y: y_list }
        df  = pd.DataFrame(raw_data, columns = [x, y ] )
        hue ="nohue"
    else:
        raise ValueError('Error: Specify rate_method to be Arrhenius or Metropolis!')
    df.to_csv(csvfile)  
    box_plot(path_name, x, y, hue, file_name_save, csvfile, figsize )

def ratevsiteration( path_name ,transType , figsize):
    """plot the the rates of the local_contexts (of the Arrhenius model)  as a box_plot"""
    RT = R * T

    if not os.path.exists(path_name):
        os.makedirs(path_name)

    par = 7
    rates =np.zeros ( (n_walkers, n_steps, par, par ) )
    rateChain = dict()
    countDim = -1
    sampler_lnprobability_finalstep, sampler_chain_finalstep,  sampler_flatlnprobability_finalstep,  sampler_flatchain_finalstep , n_steps_finalstep   = use_only_finalstep_func ()
    rate_chain_for_barchart  = np.zeros ((len (sampler_flatchain_finalstep) , 49))
    
    countHash  =dict()             
    for s in range (par) : 
        for j in range(par): 
            countDim += 1 
            countHash [j, s] = countDim
            rateChain[j, s] = [] 
            count = -1
            for i in range(n_walkers):
                for k in range(n_steps) :
                    rates [i, k , j  ,s ] = np.exp(sampler_chain[i,k, 2 * j  ] - (sampler_chain[i,k,2 *j  +1   ] / RT)) *  np.exp(sampler_chain[i,k, 2 * s  ] - (sampler_chain[i,k,2 *s  +1   ] / RT))
                    if transType == "bi" :
                        rates [i, k , j  , s ]  =   sampler_chain[i,k , ALPHA_DIM] * rates[i,k , j  , s ]
                    rateChain[j, s].append( rates [i, k , j, s ] )  
                    if n_steps_finalstep == 1 : 
                        if k == n_steps -1 :     
                            count += 1 
                            rate_chain_for_barchart[count ] [countDim]   = rates   [ i, k , j, s ]
                    else : 
                        count += 1
                        rate_chain_for_barchart [count ] [countDim]   = rates   [ i, k , j, s ]
    if transType == "bi" : 
        local_countext_counter= open('local_context_bi.pkl', 'rb')
    elif transType == "uni" :
        local_countext_counter = open('local_context_uni.pkl', 'rb')
    else:
        raise ValueError('Error: Specify transType to be "bi" or "uni"!')
    local_countext_counter= pickle.load(local_countext_counter)
    parameter_name =( "stack",  "loop"  , "end",  "stack+loop" , "stack+end", "loop+end" ,"stack+stack",r"$\alpha$" )
    csvfile = path_name + "arrheniusrate"
    file_name_save= "arrheniusrate"
    x = "Local context"
    if transType == "bi" :
        y= r"$\log_{10} k_{\rm bi}(l,r)$" +"(M"+r"$^{-1}$"+"s"+r"$^{-1}$"+ ")"
        csvfile = csvfile +"bi.csv"
        file_name_save = file_name_save + "bi.pdf"
    elif transType == "uni" :
        y = r"$\log_{10} k_{\rm uni}(l,r)$"+ "(s" +r"$^{-1}$"+ ")"
        csvfile = csvfile +"uni.csv"
        file_name_save = file_name_save + "uni.pdf"
   
    x_list = []
    y_list = []
    a1 , b1 = rate_chain_for_barchart.shape
    for a in range(a1) :
        for  s in range( par):
            for j in range(s , par) :
                if  (local_countext_counter[(parameter_name[s],parameter_name[j])]  or local_countext_counter[(parameter_name[j],parameter_name[s])] ) > 0 :
                    x_list.append (  parameter_name[j] +  "\n"+parameter_name[s]   )
                    y_list.append(  math.log10(rate_chain_for_barchart [ a , countHash [j, s]]))
    raw_data = {x: x_list, y: y_list }
    df  = pd.DataFrame(raw_data, columns = [x, y ] )
    df.to_csv(csvfile)  
    box_plot(path_name, x, y, "nohue", file_name_save, csvfile, figsize )
    

def box_plot(path_name, x, y, hue,  file_name_save, filenametoread , figsize):
    """ draws the box_plot """
    if not os.path.exists(path_name): 
        os.makedirs(path_name)

    fig, ax = plt.subplots (figsize =  figsize )
    fontsize =23
    ax.set_xlabel(x, fontsize=fontsize )
    ax.set_ylabel(x, fontsize=fontsize  )
    fontsize = 22
    sns.set_style(style='darkgrid')
    data = pd.read_csv(filenametoread)
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(fontsize ) 
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(fontsize ) 
    ax.set_xticklabels(x, rotation=90)    
    if hue =="nohue":
        sns.boxplot(x=x, y=y, data=data, color= "mediumaquamarine"  , saturation =1 ,ax =ax, width = 0.5, fliersize =3 )
        axes = plt.gca()
        axes.set_ylim([1,12])
    else : 
        sns.boxplot(x=x, y=y, hue=hue, data=data, palette="PuBu" , saturation =1 , ax =ax, width = 0.4 , fliersize =3 )
        plt.legend( loc=0, borderaxespad=0., prop={'size':17})
    plt.savefig( path_name + file_name_save, bbox_inches="tight")
    plt.close()

def correlation_plot(path_name):
    """ draws the correlation plot of parameters of a kinetic model"""
    sampler_lnprobability_finalstep, sampler_chain_finalstep, sampler_flatlnprobability_finalstep, sampler_flatchain_finalstep, n_steps_finalstep = use_only_finalstep_func()
    if not os.path.exists(path_name):
        os.makedirs(path_name)
    if parent.rate_method == ModelName.ARRHENIUSMODELNAME.value:
        legscript =( r"$\lnA_{\rm stack}$", r"$E_{\rm stack}$", r"$\lnA_{\rm loop}$",r"$E_{loop}$",  r"$\lnA_{\rm end}$" ,r"$E_{\rm end}$",  r"$\lnA_{\rm stack+loop}$" , r"$E_{\rm stack+loop}$", r"$\lnA_{\rm stack+end}$" , r"$E_{\rm stack+end}$", r"$A_{\rm loop+end}$" , r"$E_{\rm loop+end}$",r"$\lnA_{\rm stack+stack}$" , r"$E_{\rm stack+stack}$",r"$\alpha$")
    elif parent.rate_method ==ModelName.METROPOLISMODELNAME.value:
        legscript = ( r"$k_{\rm uni}$" , r"$k_{\rm bi}$")
    else:
        raise ValueError('Error: Specify rate_method to be Arrhenius or Metropolis!')
    legscript += (r"$\sigma$", )
    coeff = np.corrcoef(sampler_flatchain_finalstep, rowvar=0)
    row_labels = column_labels = legscript
    fig, ax = plt.subplots()
    pc = ax.pcolor(coeff, cmap=plt.cm.PiYG)
    ax.set_xticks(range(len(row_labels)))
    ax.set_yticks(range(len(row_labels)))
    plt.gca().set_xlim((0, len(legscript)))
    plt.gca().set_ylim((0, len(legscript)))
    fontsize=16
    ax.set_yticklabels(column_labels, fontsize=fontsize)
    ax.set_xticklabels(row_labels, rotation='vertical', fontsize=fontsize)
    for axis in [ax.xaxis, ax.yaxis]:
        axis.set(ticks=np.arange(0.5, len(legscript)), ticklabels=legscript)
    plt.colorbar(pc)
    plt.tight_layout()
    plt.savefig(path_name + "heatmap.pdf", bbox_inches="tight")

if __name__ == "__main__":
    learndnakinetics.set_configuration()
    configParser = ConfigParser.ConfigParser()
    configParser.readfp(open(r'config_file.txt'))
    CONFIG_NAME = 'loadmcmc'
    MCMC_pkl_file= configParser.get(CONFIG_NAME, 'pklFile')
    load_MCMC(MCMC_pkl_file)
    print("Mean acceptance fraction: {0:.3f}".format(np.mean(sampler.acceptance_fraction)))
    correlation_plot( MCMC_pkl_file +"/2correlation/")

    if parent.rate_method == ModelName.ARRHENIUSMODELNAME.value:
        valuevsiteration( MCMC_pkl_file +"/5valuevsiteration/" ,   (6,6)  )
    elif  parent.rate_method == ModelName.METROPOLISMODELNAME.value:
        valuevsiteration( MCMC_pkl_file +"/5valuevsiteration/" ,   (2,6) )
        
    if parent.rate_method == ModelName.ARRHENIUSMODELNAME.value :
        ratevsiteration( MCMC_pkl_file +"/8ratesvsiteration/" , "uni" ,  (12.5, 6 ) )
        ratevsiteration( MCMC_pkl_file +"/8ratesvsiteration/", "bi",  (6,6) )