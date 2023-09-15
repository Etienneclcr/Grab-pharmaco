# -*- coding: utf-8 -*-
"""
Created on Mon Feb 20 14:44:31 2023

@author: Etienne.CLAUSS
"""

import numpy as np, pandas as pd, os, matplotlib.pyplot as plt, glob
import scipy.stats as Stat, seaborn as sns
from scipy.optimize import curve_fit



folder_exp = r"C:\Etienne.CLAUSS\Lab\Exp\Imagerie\GrabOTR\Drug\AVP_test"
 
Conditions = ('0,1nM', '1nM', '10nM','50nM','100nM', '250nM', '500nM', 
              '1000nM', 'ACSF', '500nM+Ato')
ind_len = [1,2,3,4,5,6,7]
pd_auc = pd.DataFrame(columns = Conditions)
stock_TC = {}
pd_fit = pd.DataFrame(columns = Conditions)

#Define parameters of the analysis & create empty lists
drug_time = 480

def exp_func(x, a, b, c):
    return a * np.exp(-b * x) + c

#Combine the root and the file to obtain the file location
folder_exp_list = [x for x in os.listdir(folder_exp) if 'pdf' not in x]

for condition in folder_exp_list:
    folder_cond = f'{folder_exp}\{condition}'
    file_glob = glob.glob(rf'{folder_cond}\*.xlsx')
    
    #Folder to stock the final analysis
    Analysis = f'{folder_cond}\Analysis'
    if not os.path.exists(Analysis): os.makedirs(Analysis)
    
    #Create the Output data frame
    Parameters = ('ID n°', 'Data', 'AUC', 'Data raw', 'Max df')
    Output = pd.DataFrame(columns = Parameters)
        
    print('='*30, '\n'*2, condition, '\n'*2, '='*30)
          
    for file in file_glob:
        file_name = file.split('\\')[-1]
        
        data_raw = (pd.read_excel(file, sheet_name=0, header=None, 
                              usecols = 'B')).to_numpy()
        print('='*30, '\n'*2, file_name, '\n'*2, '='*30)
        
        
        '''                      Beginning of the analysis                   '''
        
        
        #Define the coeficient
        data_raw = data_raw.ravel()
        data = data_raw[120:]
        data_tofit = data-min(data)
        xdata = np.arange(len(data_tofit))
        
        # 1st fit
        popt, pcov = curve_fit(exp_func, xdata[np.r_[:drug_time, 1960:len(data_tofit)]], 
                                data_tofit[np.r_[:drug_time, 1960:len(data_tofit)]], 
                                p0=[1,0.0001,1], maxfev = 10000)

        
        exp_inv = exp_func(xdata, *popt)
        data = data_tofit - exp_inv
        fit_diff = np.nanmean(data_tofit[:drug_time] - exp_inv[:drug_time])**2
        
        if fit_diff > 0.01:     
                                    

            
            # 2nd fit
            popt, pcov = curve_fit(exp_func, xdata[:drug_time], data_tofit[:drug_time], 
                                    p0=[1,0.0001,1], maxfev = 10000)
    
            exp_inv = exp_func(xdata, *popt)        
            data = data_tofit - exp_inv        

        plt.figure()
        plt.plot(data)
        plt.savefig(rf'{Analysis}\{file_name}.pdf')
        plt.close() 
        
        auc_cell = data.clip(0)
        auc = (np.trapz(auc_cell[:drug_time]),
                    np.trapz(auc_cell[drug_time:]))


        #Stock section
        #General Stock
        Out = (file_name, data, auc, data_raw)
        Output = Output.append({x:y for x, y in zip(Parameters, Out)}, 
                               ignore_index = True)
        

    #Barplot AUC 
    bl, dr = [],[]
    for auc in Output['AUC']:
        bl.append(auc[0])
        dr.append(auc[1])
    meanz = np.asarray((np.mean(bl), np.mean(dr)))
    val_stat_AUC = np.asarray((bl, dr))
        
    bl_sem, dr_sem = (Stat.sem(bl), Stat.sem(dr))  
    semz = np.asarray((np.mean(bl_sem), np.mean(dr_sem)))
    x_ax = (0,0.2)
    plt.bar(x_ax, meanz, width = 0.1,
            yerr=semz, capsize=10, zorder=0)
    
    plt.savefig(rf'{Analysis}\AUC.pdf')
    plt.close()
    
    if condition in pd_auc:  
        if not len(ind_len) == len(dr):
            dr.extend(['']*(len(ind_len)-len(dr)))
        pd_auc[condition] = dr

    #Time course 
    data_TC = [x for x in Output['Data']]
 
    data_TC = np.asarray(data_TC)
    stock_TC[condition] = data_TC
    mean_bin = np.nanmean(data_TC, axis = 0)
    sem_bin = [Stat.sem(data_TC[:,i]) for i in xdata]
    
    plt.fill_between(xdata, mean_bin-sem_bin, mean_bin+sem_bin)
    plt.axvline(drug_time, c='gold', lw = 2)
    
    plt.plot(xdata, mean_bin, c='r', zorder=2)

    plt.savefig(rf'{Analysis}\TC.pdf')
    plt.close()
    
    #Time course normalised to max value
    data_norm = []
    for data in Output['Data']:
        max_data = max(data)
        data_norm.append([x/max_data for x in data])
                
    data_norm = np.asarray(data_norm)
    mean_bin = np.nanmean(data_norm, axis = 0)
    sem_bin = [Stat.sem(data_norm[:,i]) for i in xdata]
    
    plt.fill_between(xdata, mean_bin-sem_bin, mean_bin+sem_bin)
    plt.axvline(drug_time, c='gold', lw = 2)
    
    plt.plot(xdata, mean_bin, c='r', zorder=2)

    plt.savefig(rf'{Analysis}\TC_norm.pdf')
    plt.close()

    # Output AUC for stat to excel
    df = pd.DataFrame(val_stat_AUC)
    writer_exp = pd.ExcelWriter(f'{Analysis}/AUC_{condition}.xlsx')
    df.to_excel(writer_exp)
    writer_exp.save()

    
    #Heatmap
    def Sorter(data_tot):
        stim = sum(data_tot[drug_time:1000])
        return stim 
    
    data_tot = np.zeros((len(Output['Data']), len(Output['Data'][0])))
    

    for i,x in enumerate(Output['Data']):
        data_tot[i] = x
    
    data_tot = sorted(data_tot, key=Sorter)
    
    
    plt.figure()
    fig = sns.heatmap(data_tot, vmin = -20, vmax = 150, cmap = 'Reds')
    fig.axvline(drug_time, c='Black', lw = 2)
    
    plt.savefig(rf'{Analysis}\HeatMap.pdf')
    plt.close()


# Plot all AUC 
plt.figure()
x_ax = (0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8)
for i, condition in enumerate(Conditions):
    data_indiv = pd_auc[condition].tolist()
    while '' in data_indiv:
        data_indiv.remove('')
    mean = np.mean(data_indiv)        
    sem = Stat.sem(data_indiv)
    plt.bar(x_ax[i], mean, width = 0.1,
            yerr=sem, capsize=10, zorder=0)
    for x in data_indiv:
        plt.scatter(x_ax[i], x, s=20, c='w', marker='o',
                    edgecolor='k', zorder=2)   
plt.savefig(rf'{folder_exp}_auc.pdf')
plt.close()

# All AUC in Excel
writer_exp = pd.ExcelWriter(f'{folder_exp}_AUC.xlsx')
pd_auc.to_excel(writer_exp)
writer_exp.save()       


#Plot all time courses
plt.figure()
for condition in Conditions:
    try:
        data_TC = stock_TC[condition]
    except KeyError: continue
    mean_bin = np.nanmean(data_TC, axis = 0)
    sem_bin = [Stat.sem(data_TC[:,i]) for i in xdata]
    
    plt.fill_between(xdata, mean_bin-sem_bin, mean_bin+sem_bin)    
    plt.plot(xdata, mean_bin, c='r', zorder=2)

plt.axvline(drug_time, c='gold', lw = 2)
plt.savefig(rf'{folder_exp}_TC.pdf')
plt.close()



