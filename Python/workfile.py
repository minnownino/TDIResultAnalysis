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
TDI_FILE_PATH = "TSDtriplet.csv"
# POS_FILE_PATH = "data_noCT/PANCAN.postprobcutoff.perSGA.csv"
POS_FILE_PATH = "postprobthd.p=0.05.perSGA.csv"
SGA_FILE_PATH = "SGApair.csv"
print TDI_FILE_PATH
print POS_FILE_PATH
print SGA_FILE_PATH

###### read file to get the dataframe above threshold 
df = read_file(TDI_FILE_PATH, POS_FILE_PATH)


###### find driver per tumor and get the driver gene triplet
(driverTriplet, driverPerTumor) = driver(df)
driverPerTumor.to_csv('driverCallPerTumor%s-%s.csv'%(date, i))
driverTriplet.to_csv('tripletFilter%s-%s.csv'%(date, i))
print 'driverCallPerTumor%s-%s.csv'%(date, i)+" saved"


###### find significant drivers
result = sigDriverCallRate(driverPerTumor, SGA_FILE_PATH)
with open('driverCallRate%s-%s.csv'%(date, i), 'wb') as csvfile:
    wr = csv.writer(csvfile)
    #wr.writerow(['SGA_name', '#Driverevents', '#SGAevents', 'Drivercallrate'])
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

###### get top 20 SGAs for this version of TDI result
df = driverPerTumor
df_tumor = df[['patient_name', 'cause_gene_name']]
df_tumor = df_tumor.drop_duplicates()
grouped = df_tumor.groupby(df_tumor.cause_gene_name).size()
sgatop20 = [x for x in grouped.nlargest(20).index]
with open('top20SGA%s-%s.csv'%(date, i), 'w') as f:
    for sga in sgatop20:
        f.write(sga+'\n')
print grouped.nlargest(20)
print 'top20SGA%s-%s.csv'%(date, i)+' saved'

###########################################33
# calculate the count of each tumor type
df, groups = getCancertypeTotalCount(TDI_FILE_PATH)
with open('totalcountofeachtype%s-%s.csv'%(date, i), 'wb') as csv_file:
    writer = csv.writer(csv_file)
    for key, value in groups.items():
        writer.writerow([key, value])
print 'totalcountofeachtype%s-%s.csv'%(date, i)+' saved'


top20sga = pd.read_csv('top20SGA2017-11-28-0.csv', header = None)
top20sga.columns = ['sganame']
sgas = top20sga['sganame'].tolist()

count = pd.read_csv('totalcountofeachtype2017-11-28-0.csv', header = None)
count.columns = ['cancertype', 'count']
count = count.set_index('cancertype')
typecountmap = count.to_dict()
typecountmap = typecountmap['count']
result = {}
for sga in sgas:
    result[sga] = getSGAcancertypedistribution(sga)

for item in result:
    for i in result[item]:
        for ct in result[item][i]:
            result[item][i][ct] = result[item][i][ct] / float(typecountmap[ct])

returndata = {"list": sgas, "data": result}

print 'printing the top 20 data into the pickle file indexdata%s-%s.plk'%(date, i)
output = open('indexdata%s-%s.plk'%(date, i), 'wb')
#Pickle dictionary using protocol 0.
pickle.dump(returndata, output)
output.close()
print "finish"
