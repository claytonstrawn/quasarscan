import numpy as np
import os
def dict_of_vela_info(quantity,loud = False):      #aexp,Rvir,Mvir_gas,Mvir_star,Mvir_dm,Mvir,fb/fbcosmo,fs/fbcosmo,Mgcold,Mgcool,Mgwarm,Mgwhot,Mghot,Mgcold(r>0.1Rvir),cool,warm,whot,hot 
    
    #This is a guide for what info can be found in what column, starting from column 0.
    quantity_dict_Mvir_fbar = {"a":0,"Rvir":1,"gas_Rvir":2,"star_Rvir":3,"dm_Rvir":4,"Mvir":5}
    quantity_dict_Nir_disc_cat = {"L":[8,9,10],"cm":[2,3,4],"vcm":[5,6,7],"L_mag":11}
    quantity_dict_Nir_spherical_galaxy_cat = {"SFR":13}
    
    #This calculates the column index as well as the number of data points for one piece of info
    if quantity in quantity_dict_Mvir_fbar.keys():
        index = quantity_dict_Mvir_fbar[quantity]
        numvals = 1
    elif quantity in quantity_dict_Nir_spherical_galaxy_cat.keys():
        print "not implemented for NIHAO"
        return
        index = quantity_dict_Nir_spherical_galaxy_cat[quantity]
        numvals = 1
    elif quantity in quantity_dict_Nir_disc_cat.keys():
        print "not implemented for NIHAO"
        return
        index = quantity_dict_Nir_disc_cat[quantity]
        numvals = 3
        if quantity == "L_mag":
            numvals = 1
    
    #Setting up the dictionary of the desired property, where the dictionary key is the name of the galaxy folder : the redshift, and the dictionary value is the respective desired property quantity
    ret_dict = {}
    if os.path.isdir("galaxy_catalogs_nihao"):
        basepath = "galaxy_catalogs_nihao/"
    else: 
        basepath = "quasarscan/galaxy_catalogs_nihao/"
    allfiles = os.listdir(basepath)

    for folder in allfiles:
        if folder.startswith("."):
            continue
        ret_dict[folder] = {}
        if quantity in quantity_dict_Mvir_fbar.keys():
            pathname = basepath + folder + "/Mvir_fbar_R200.out"
        #elif quantity in quantity_dict_Nir_disc_cat.keys():
            #pathname = basepath + folder + "/galaxy_catalogue/Nir_disc_cat.txt"
        #elif quantity in quantity_dict_Nir_spherical_galaxy_cat.keys():
            #pathname = basepath + folder + "/galaxy_catalogue/Nir_spherical_galaxy_cat.txt"
        else: 
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
        a = line.split()[0]
        while 1:
            if numvals == 1:
                ret_dict[folder][a] = line.split()[index]
            else:
                ret_dict[folder][a] = [line.split()[index[0]],line.split()[index[1]],line.split()[index[2]]]
            line = f.readline()
            try:
                a = line.split()[0]
            except:
                break
                f.close()
    return ret_dict

Rdict = dict_of_vela_info("Rvir")
#Ldict = dict_of_vela_info("L")