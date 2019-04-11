import os

#nihao galaxies, ordered by mass
namedict = {'ENZO': 'enzo',
 'RAMSES': 'ramses',
 'GADGET-3': 'gadget',
 'GEAR': 'gear',
 'GIZMO': 'gizmo',
 'ART': 'art'}

quantity_dict_MsMv = {"a":-1,"z":0,"Rvir":1,"star_Rvir":2,"Mvir":3,"gas_Rvir":4}
quantity_dict_center = {"center":1}

quantity_dict = quantity_dict_MsMv.copy()
quantity_dict.update(quantity_dict_center)

def dict_of_agora_info(quantity,loud = 0):      #aexp,Rvir,Mvir_gas,Mvir_star,Mvir_dm,Mvir,fb/fbcosmo,fs/fbcosmo,Mgcold,Mgcool,Mgwarm,Mgwhot,Mghot,Mgcold(r>0.1Rvir),cool,warm,whot,hot 
    
    #This is a guide for what info can be found in what column, starting from column 0.
    # NEED SFR
    
    #This calculates the column index as well as the number of data points for one piece of info
    if quantity in quantity_dict_MsMv.keys():
        index = quantity_dict[quantity]
        numvals = 1
    elif quantity in quantity_dict_center.keys():
        index = quantity_dict[quantity]
        numvals = 3
        
    #Setting up the dictionary of the desired property, where the dictionary key is the name of the galaxy folder : the redshift, and the dictionary value is the respective desired property quantity
    ret_dict = {'AGORA_v1_enzo_01':{},'AGORA_v1_ramses_01':{},'AGORA_v1_gadget_01':{},'AGORA_v1_gear_01':{},'AGORA_v1_gizmo_01':{},'AGORA_v1_art_01':{}}
    if os.path.isdir("galaxy_catalogs"):
        basepath = "galaxy_catalogs/galaxy_catalogs_agora/"
    else: 
        basepath = "quasarscan/galaxy_catalogs/galaxy_catalogs_agora/"
    allfiles = os.listdir(basepath)

    for z in allfiles:
        if z.startswith("."):
            continue
        if quantity in quantity_dict_MsMv.keys():
            pathname = basepath + z + "/%s_MsMv.txt"%z
        elif quantity in quantity_dict_center.keys():
            pathname = basepath + z + "/%s_center.txt"%z
        elif loud > 1: 
            print("where is %s stored? We don't know :("%quantity)
        try:
            f = open(pathname)
        except:
            if loud:
                print("Error reading %s"%pathname)
            continue
        
        # the following lines including the while loop retrieve the desired property value by iterating through the textfile line by line 
        lines = f.readlines()
        for i,line in enumerate(lines):
            if line.startswith("Code: "):
                code = namedict[line.split(' ')[1]]
                name = "AGORA_v1_%s_01"%code
                start = i
                z = float(line.split("=")[1])
                a = 1./(z+1)
                if index < 0:
                    ret_dict[name][a] = a
                elif numvals == 1:
                    v_str = lines[i+index].split('=')[1].strip(' ').split()[0]
                    if quantity in ["star_Rvir","Mvir","gas_Rvir"]:
                        ret_dict[name][a] = float(v_str)*1e8
                    else:
                        ret_dict[name][a] = float(v_str)
                else:
                    v_str = lines[i+index].split('=')[1].split("code_length")[0].strip(' []').split()
                    ret_dict[name][a] = [float(v_str[0]),float(v_str[1]),float(v_str[2])]
        f.close()
    return ret_dict

Rdict = dict_of_agora_info("Rvir")
adict = dict_of_agora_info("a")
