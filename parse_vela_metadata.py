import numpy as np
def dict_of_vela_info(quantity, iarr= np.arange(1,35),loud = False):
    quantity_dict_single = {"a":0,"Rvir":1,"Rdisk":2,"Mvir":3,\
                    "gas_Rvir":4,"star_Rvir":5,"dm_Rvir":6,\
                    "gas_.1Rvir":5,"star_.1Rvir":6,"dm_.1Rvir":7,\
                    "gas_10kpc":8,"star_10kpc":9,"dm_10kpc":10,\
                    "gas_Rdisk":11,"star_Rdisk":12,"dm_Rdisk":13}
    quantity_dict_triple = {"L":[8,9,10],"cm":[2,3,4],"vcm":[5,6,7]}
    if quantity in quantity_dict_single.keys():
        index = quantity_dict_single[quantity]
        numvals = 1
    elif quantity in quantity_dict_triple.keys():
        index = quantity_dict_triple[quantity]
        numvals = 3        
    ret_dict = {}
    basepath = "~/quasarscan/galaxy_catalogs/"
    for version in (1,2):
        if version == 1:
            folderstart = "VELA"
        elif version == 1:
            folderstart = "VELA_v2_"
        for i in range(len(iarr)):
            folder = folderstart+"%02i"%iarr[i]
            ret_dict[folder] = {}
            if numvals == 1:
                pathname = basepath + folder + "/galaxy_catalogue/Mstar.txt"
            else:
                pathname = basepath + folder + "/galaxy_catalogue/Nir_disc_cat.txt"
            try:
                f = file(pathname)
            except:
                if loud:
                    print("Error reading %s"%pathname)
                continue
            f.readline()
            line = f.readline()
            a = line.split()[0]
            while 1:
                if numvals == 1:
                    ret_dict[folder][a] = line.split()[index]
                else:
                    ret_dict[folder][a] = [line.split()[index[0]],line.split()[index[1]],line.split()[index[2]]]
                line = f.readline()[1:]
                try:
                    a = line.split()[0]
                except:
                    break
                    f.close()
        return ret_dict
    
Rdict = dict_of_vela_info("Rvir",iarr = np.arange(1,35))
Ldict = dict_of_vela_info("Rvir",iarr = np.arange(1,35))