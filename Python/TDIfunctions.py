
def read_file(TDI_FILE_PATH, POS_FILE_PATH):  
    """
        read the tdi file and posterior file 
        return filtered dataframe of tdi results
    """
    # df = pd.read_csv("data_noCT/TSDtriplet.csv")
    # df_posterior = pd.read_csv("data_noCT/PANCAN.postprobcutoff.perSGA.csv")

    df = pd.read_csv(TDI_FILE_PATH)
    df_posterior = pd.read_csv(POS_FILE_PATH)
    tmp = dict(zip(df_posterior['Unnamed: 0'], df_posterior.PostProbcutoff))
    df['a'] = df.cause_gene_name.map(tmp)
    # cutoff by posterior
    df_updated_cutoff = df.loc[df['posterior'] >= df['a']]
    del df_updated_cutoff['a']
    return df_updated_cutoff

def driver(df):
    """
        input: read_file() output
        output: csv file of driver call per tumor
    """
    grouped = df.groupby('patient_name')
    tmp = []
    for name, group in grouped:
        filtergroup = group.groupby('cause_gene_name').filter(lambda x : len(x) >= 5)
        tmp.append(filtergroup)
    driverTriplet = pd.concat(tmp)
    driverPerTumor = pd.concat(tmp)
    del driverPerTumor['result_gene_name']
    del driverPerTumor['posterior']
    driverPerTumor = driverPerTumor.drop_duplicates()
    return (driverTriplet, driverPerTumor)

def sigDriverCallRate(driverPerTumor, SGA_FILE_PATH):
    """
        SGA which is driver in more than 30 tumors & callrate larger than 0.5
        input: driverPerTumor, sga events
        output: list of sga, drivercallnum, sgaeventsnum, ratio
    """
    df_sga = pd.read_csv(SGA_FILE_PATH)
    result = []
    sgas = df_sga.SGA_name.unique()
    sigDrivers = driverPerTumor.groupby('cause_gene_name').filter(lambda x: len(x) >= 30)
    count_driver = sigDrivers.groupby('cause_gene_name').size()
    count_sga = df_sga.groupby('SGA_name').size()
    for sga, count in count_driver.iteritems():
        if sga in sgas:
            if count_driver[sga] / float(count_sga[sga]) >= 0.5:
                 result.append((sga, count, count_sga[sga], round(float(count)/count_sga[sga],2)))
    return result

def getsigSGADEG(driverTriplet, sigDrivers):
    """
        Generate the significant DEG list of significant SGAs with a cutoff of 0.2
        input: driverTriplet
        output: list of (sga, deg, ratio)
    """
    result = []
    #sigDrivers.info()
    sigDriverTriplet = driverTriplet.loc[driverTriplet.cause_gene_name.isin(sigDrivers[0].unique())]
    # for sga, take the partial, count the total number of tumors, then group by degs, and calc the partition
    sgas = sigDriverTriplet.cause_gene_name.unique()
    for sga in sgas:
        subdf = sigDriverTriplet.loc[sigDriverTriplet.cause_gene_name == sga]
        numOfTumors = subdf.patient_name.unique().size
        deg_count = subdf.groupby('result_gene_name').size()
        for deg, count in deg_count.iteritems():
            if float(count) / numOfTumors >= 0.2:
                result.append((sga, deg, float(count)/numOfTumors))
    return result

def generateSGAcsv(FILE_PATH):
    """
        transform sga-tumor 0/1 matrix into pairs for further use
    """
    #df = pd.read_csv('PanCancer13tts.SGAmatrix.4TCI.csv')
    df = pd.read_csv(FILE_PATH)
    bt = df.apply(lambda x: x > 0)
    cols = df.columns
    bt.apply(lambda x: list(cols[x.values]), axis=1)
    df = df.set_index('Unnamed: 0')
    result = bt.apply(lambda x: list(cols[x.values]), axis=1)
    f = open('SGApair.csv', 'w')
    f.write('SGA_name,tumor_name\n')
    for item in result.iteritems():
        (tumor, sgas) = item
        for sga in sgas:
            f.write(sga + ',' + tumor + '\n')

def getCancertypeTotalCount(TSDTRIPLET_FILE_PATH):
    """
     get total number of tumors in each cancer type for one version of TDI results
     return groups count of each cancer type
    """
    #read the mapping of tumor-cancer type, generated a lookup dictionary
    df_cttumormap = pd.read_csv('TumorID.vs.CancerType.v20160321.csv')
    typedictionary = dict(zip(df_cttumormap.TumorID, df_cttumormap.CancerType))
    
    df = pd.read_csv('TSDtriplet.csv')
    df['tumortype'] = df.patient_name.map(typedictionary)
    df_new = df[['patient_name', 'tumortype']]
    df_new = df_new.drop_duplicates()
    groups = df_new.groupby('tumortype').size()
    return df, groups.to_dict()
