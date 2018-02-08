import os
import time
import pandas as pd
import multiprocessing as mp


class Timer:
  def __init__(self):
    self.start = time.time()

  def restart(self):
    self.start = time.time()

  def get_time_hhmmss(self):
    end = time.time()
    m, s = divmod(end - self.start, 60)
    h, m = divmod(m, 60)
    time_str = "%02d:%02d:%02d" % (h, m, s)
    return time_str

def getTopGtpostprob_parallel(PATH_TDI, FILE):
    df = pd.read_csv(PATH_TDI+FILE,index_col=0)
    tumor = FILE[:-4] #remove '.csv'
    TSDPs = []
    for col in df:
        topIdx = df[col].idxmax()
        SGA = topIdx
        DEG = col
        postprob = df[col].loc[topIdx]
        TSDP = [tumor, SGA, DEG, postprob]
        TSDPs.append(TSDP)
    # df_triplet = pd.DataFrame(SDPs,columns=['SGA','DEG','postprob'])
    #df_triplet.insert(0,"Tumor",tumor)
    return TSDPs

def outputTriplets(triplets,PATHNAME_triplet):
    df_triplets = pd.DataFrame(data=triplets, columns=['Tumor', 'SGA', 'DEG', 'postprob'])
    df_triplets.to_csv(PATHNAME_triplet,index=False)

######################################################################################################
# Main program
#

#triplets = extractTriplets(PATH_TDI,PATHNAME_triplet,fstart,fend)
if __name__ == "__main__":
    PATH_TDI = "C:/Users/XIM33/Documents/NetBeansProjects/TDIC_01_CPU_Tumors_SGAprior/output/"
    PATHNAME_triplet = "C:/Users/XIM33/Documents/NetBeansProjects/TDIC_01_CPU_Tumors_SGAprior/triplets.csv"

    my_timer = Timer()

    num_cores = mp.cpu_count()
    pool = mp.Pool() #defalut is num_cores ==   pool = mp.Pool(num_cores)

#    triplets = [pool.apply_async( getTopGtpostprob_parallel, args=(PATH_TDI,FILE) ) for FILE in os.listdir(PATH_TDI))]
    results = []
    for FILE in os.listdir(PATH_TDI):
        results.append(pool.apply_async( getTopGtpostprob_parallel, args=(PATH_TDI,FILE) ))

    triplets = []
    for result in results:
        triplet = result.get()
        triplets += triplet
    outputTriplets(triplets,PATHNAME_triplet)

    time_hhmmss = my_timer.get_time_hhmmss()
    print("Time elapsed: %s\n" % time_hhmmss)