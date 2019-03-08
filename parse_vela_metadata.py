import os
def dict_of_vela_info(quantity,loud = False):
    #This is a guide for what info can be found in what column, starting from column 0.
    quantity_dict_Mstar = {"a":0,"Rvir":1,"Rdisk":2,"Mvir":3,\
                    "gas_Rvir":4,"star_Rvir":5,"dm_Rvir":6,\
                    "gas_.1Rvir":5,"star_.1Rvir":6,"dm_.1Rvir":7,\
                    "gas_10kpc":8,"star_10kpc":9,"dm_10kpc":10,\
                    "gas_Rdisk":11,"star_Rdisk":12,"dm_Rdisk":13}
    quantity_dict_Nir_disc_cat = {"L":[8,9,10],"cm":[2,3,4],"vcm":[5,6,7],"L_mag":11}
    quantity_dict_Nir_spherical_galaxy_cat = {"SFR":13}
    
    #This calculates the column index as well as the number of data points for one piece of info
    if quantity in quantity_dict_Mstar.keys():
        index = quantity_dict_Mstar[quantity]
        numvals = 1
    elif quantity in quantity_dict_Nir_spherical_galaxy_cat.keys():
        index = quantity_dict_Nir_spherical_galaxy_cat[quantity]
        numvals = 1
    elif quantity in quantity_dict_Nir_disc_cat.keys():
        index = quantity_dict_Nir_disc_cat[quantity]
        numvals = 3
        if quantity == "L_mag":
            numvals = 1
  
    #Setting up the dictionary of the desired property, where the dictionary key is the name of the galaxy folder : the redshift, and the dictionary value is the respective desired property quantity
    ret_dict = {}
    if os.path.isdir("galaxy_catalogs"):
        basepath = "galaxy_catalogs/galaxy_catalogs_vela/"
    else: 
        basepath = "quasarscan/galaxy_catalogs/galaxy_catalogs_vela/"
    for version in range(1,3):
        if version == 1:
            folderstart = "VELA"
        elif version == 2:
            folderstart = "VELA_v2_"
        for i in range(35):
            folder = folderstart+"%02i"%i
            ret_dict["VELA_v%d_art_%02d"%(version,i)] = {}
            if quantity in quantity_dict_Mstar.keys():
                pathname = basepath + folder + "/galaxy_catalogue/Mstar.txt"
            elif quantity in quantity_dict_Nir_disc_cat.keys():
                pathname = basepath + folder + "/galaxy_catalogue/Nir_disc_cat.txt"
            elif quantity in quantity_dict_Nir_spherical_galaxy_cat.keys():
                pathname = basepath + folder + "/galaxy_catalogue/Nir_spherical_galaxy_cat.txt"
            else: 
                print("where is %s stored? We don't know :("%quantity)
            try:
                f = open(pathname)
            except:
                if loud:
                    print("Error reading %s"%pathname)
                continue
            f.readline() #one-time skips the identification number on the first line
            
            # the following lines including the while loop retrieve the desired property value by iterating through the textfile line by line 
            line = f.readline() 
            a = float(line.split()[0])
            while 1:
                if numvals == 1:
                    ret_dict["VELA_v%d_art_%02d"%(version,i)][a] = float(line.split()[index])
                else:
                    ret_dict["VELA_v%d_art_%02d"%(version,i)][a] = [float(line.split()[index[0]]),float(line.split()[index[1]]),float(line.split()[index[2]])]
                line = f.readline()[:]
                try:
                    a = float(line.split()[0])
                except:
                    break
                    f.close()
    return ret_dict
    
Rdict = dict_of_vela_info("Rvir")
Ldict = dict_of_vela_info("L")
adict = dict_of_vela_info("a")
