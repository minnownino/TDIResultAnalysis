
import pandas as pd
import numpy as np
import csv
import os
import datetime
from TDIfunctions import *

###### prepare filename
print "Prepare filename"
date = datetime.datetime.now().strftime("%Y-%m-%d")
i = 0
while os.path.exists('driverCallPerTumor%s-%s.csv'%(date, i)):
    i += 1
print date, str(i), "finished"
    
###### set path
print "set path"
TDI_FILE_PATH = "data_noCT/TSDtriplet.csv"
# POS_FILE_PATH = "data_noCT/PANCAN.postprobcutoff.perSGA.csv"
POS_FILE_PATH = "data_noCT/PANCAN.postprobthd.p=0.001.perSGA.csv"
SGA_FILE_PATH = "data_noCT/SGAs.csv"

###### read file to get the dataframe above threshold 
df = read_file(TDI_FILE_PATH, POS_FILE_PATH)


###### find driver per tumor and get the driver gene triplet
(driverTriplet, driverPerTumor) = driver(df)
driverPerTumor.to_csv('driverCallPerTumor%s-%s.csv'%(date, i))
driverTriplet.to_csv('')
print 'driverCallPerTumor%s-%s.csv'%(date, i)+" saved"


###### find significant drivers
result = sigDriverCallRate(driverPerTumor, SGA_FILE_PATH)
with open('driverCallRate%s-%s.csv'%(date, i), 'wb') as csvfile:
    wr = csv.writer(csvfile)
    for row in result:
        wr.writerow([row[0], row[1], row[2], round(float(row[3]),2)])
print 'driverCallRate%s-%s.csv'%(date, i)+" saved"
        
sigDrivers = pd.read_csv('driverCallRate%s-%s.csv'%(date, i), header=None)
sigDriverTriplet = driverTriplet.loc[driverTriplet.cause_gene_name.isin(sigDrivers[0].unique())]
sigDriverTriplet.to_csv('sigDriverTriplet%s-%s.csv'%(date, i))
print 'sigDriverTriplet%s-%s.csv'%(date, i)+" saved"

###### generate SGA DEGs
#sigDrivers = pd.read_csv('driverCallPerTumor%s-%s.csv'%(date, i), header=None)
sigSGADEG = getsigSGADEG(driverTriplet, sigDrivers)
with open('sigSGADEG%s-%s.csv'%(date, i), 'wb') as csvfile:
    wr = csv.writer(csvfile)
    for row in sigSGADEG:
        wr.writerow([row[0], row[1], round(float(row[2]),2)])
print 'sigSGADEG%s-%s.csv'%(date, i)+' saved'
        

###### generate significant triplet
sigDriverTriplet = pd.read_csv('sigDriverTriplet%s-%s.csv'%(date, i))
sigSGADEG = pd.read_csv('sigSGADEG%s-%s.csv'%(date, i), header=None)
sigSGADEG.columns = ['cause_gene_name', 'result_gene_name', 'ratio']
# sigSGADEG.head()
sigSGADEGTUMOR = pd.merge(sigDriverTriplet, sigSGADEG, how = 'inner', on=['cause_gene_name', 'result_gene_name'])
del sigSGADEGTUMOR['Unnamed: 0']
del sigSGADEGTUMOR['ratio']
sigSGADEGTUMOR.to_csv('sigSGADEGTUMOR%s-%s.csv'%(date, i))
print 'sigSGADEGTUMOR%s-%s.csv'%(date, i)+' saved'
print "finish"