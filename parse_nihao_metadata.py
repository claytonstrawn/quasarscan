import os

#nihao galaxies, ordered by mass
namedict = {'g2.63e10': 'NIHAO_v1_tipsy_01',
 'g2.19e11': 'NIHAO_v1_tipsy_02',
 'g5.02e11': 'NIHAO_v1_tipsy_03',
 'g5.31e11': 'NIHAO_v1_tipsy_04',
 'g5.36e11': 'NIHAO_v1_tipsy_05',
 'g5.38e11': 'NIHAO_v1_tipsy_06',
 'g6.96e11': 'NIHAO_v1_tipsy_07',
 'g7.08e11': 'NIHAO_v1_tipsy_08',
 'g7.44e11': 'NIHAO_v1_tipsy_09',
 'g7.55e11': 'NIHAO_v1_tipsy_10',
 'g7.66e11': 'NIHAO_v1_tipsy_11',
 'g8.06e11': 'NIHAO_v1_tipsy_12',
 'g8.13e11': 'NIHAO_v1_tipsy_13',
 'g8.26e11': 'NIHAO_v1_tipsy_14',
 'g1.12e12': 'NIHAO_v1_tipsy_15',
 'g1.77e12': 'NIHAO_v1_tipsy_16',
 'g1.92e12': 'NIHAO_v1_tipsy_17',
 'g2.79e12': 'NIHAO_v1_tipsy_18'}
inverted_namedict = dict([[v,k] for k,v in namedict.items()])

quantity_dict_Mvir_fbar = {"a":0,"Rvir":1,"gas_Rvir":2,"star_Rvir":3,"dm_Rvir":4,"Mvir":5}
quantity_dict = quantity_dict_Mvir_fbar.copy()


def dict_of_nihao_info(quantity,loud = 0):      #aexp,Rvir,Mvir_gas,Mvir_star,Mvir_dm,Mvir,fb/fbcosmo,fs/fbcosmo,Mgcold,Mgcool,Mgwarm,Mgwhot,Mghot,Mgcold(r>0.1Rvir),cool,warm,whot,hot 
    
    #This is a guide for what info can be found in what column, starting from column 0.
    # NEED SFR
    
    #This calculates the column index as well as the number of data points for one piece of info
    if quantity in quantity_dict_Mvir_fbar.keys():
        index = quantity_dict_Mvir_fbar[quantity]
        numvals = 1
    
    #Setting up the dictionary of the desired property, where the dictionary key is the name of the galaxy folder : the redshift, and the dictionary value is the respective desired property quantity
    ret_dict = {}
    if os.path.isdir("galaxy_catalogs"):
        basepath = "galaxy_catalogs/galaxy_catalogs_nihao/"
    else: 
        basepath = "quasarscan/galaxy_catalogs/galaxy_catalogs_nihao/"
    allfiles = os.listdir(basepath)

    for folder in allfiles:
        if folder.startswith("."):
            continue
        name = namedict[folder]
        ret_dict[name] = {}
        if quantity in quantity_dict_Mvir_fbar.keys():
            pathname = basepath + folder + "/Mvir_fbar_R200.out"
        elif loud > 1: 
            print("where is %s stored? We don't know :("%quantity)
        try:
            f = open(pathname)
        except:
            if loud:
                print("Error reading %s"%pathname)
            continue
        
        #skips the first two metadata lines
        f.readline()
        f.readline()
        
         # the following lines including the while loop retrieve the desired property value by iterating through the textfile line by line 
        line = f.readline()
        a = float(line.split()[0])
        while 1:
            if numvals == 1:
                ret_dict[name][a] = float(line.split()[index])
            else:
                ret_dict[name][a] = [float(line.split()[index[0]]),float(line.split()[index[1]]),float(line.split()[index[2]])]
            line = f.readline()
            try:
                a = float(line.split()[0])
            except:
                break
                f.close()
    return ret_dict

Rdict = dict_of_nihao_info("Rvir")
#Ldict = dict_of_vela_info("L")
adict = dict_of_nihao_info("a")
