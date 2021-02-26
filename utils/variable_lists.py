stringcriteria = ["ions","name","simname","version","code","simnum","compaction_stage"]
intensives = ["Z","T","n"]
intensiveslabels = {"Z":"avg metallicity","T":"avg temperature","n":"avg density"}
intensivespositions = {"Z":-1,"T":-2,"n":-3}
sightline_xVars = ["r","rdivR","rMpc","theta","phi"]
param_xVars = ["redshift","a0","Mvir","gas_Rvir","star_Rvir","dm_Rvir","sfr","ssfr","L_mag","Mstar","Mgas","Rvir"]
sightline_unit_labels = {"r":"r (kpc)","r>0":"r (kpc)","rdivR":"r/Rvir","rdivR>0":"r/Rvir",\
           "theta":"viewing angle (rad)","theta_r>0":"viewing angle (rad)","phi" \
           :"azimuthal viewing angle (rad)","rMpc":"r (Mpc)"}
param_unit_labels = {"redshift":"z","a0":"a","Rvir":'Virial radius (kpc)',"Mvir":"Virial Mass (Msun)",\
                    "gas_Rvir":"Gas Mass within Rvir (Msun)","Mgas":"Gas Mass within Rvir (Msun)","star_Rvir":"Stellar Mass within Rvir (Msun)",\
                    "Mstar":"Stellar Mass within Rvir (Msun)","dm_Rvir":"Dark Matter Mass within Rvir (Msun)","sfr":"Star Formation Rate (Msun yr-1)",\
                    "ssfr":"Specific Star Formation Rate (Msun yr-1 Mstar-1)","L_mag":"Magnitude of Angular Momentum"}
all_known_variables = stringcriteria+intensives+sightline_xVars+param_xVars