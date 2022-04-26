stringcriteria = ["ions","name","simname","version","code","simnum","intensives"]

intensives = ["Z","T","nmax","nmean",'n','rho']
intensiveslabels = {"Z":"avg metallicity","T":"avg temperature(K)","nmean":"avg density(1/cm3)",\
                    "nmax":"max density(1/cm3)",'n':'average density(1/cm3)','rho':'average density (g/cm3)'}

sightline_xVars = ["r","rdivR","rMpc","theta","phi","sightline_length"]
sightline_unit_labels = {"r":"r (kpc)","r>0":"r (kpc)","rdivR":"r/Rvir","rdivR>0":"r/Rvir",\
                         "theta":"viewing angle (rad)","theta_r>0":"viewing angle (rad)",\
                         "phi":"azimuthal viewing angle (rad)","rMpc":"r (Mpc)","sightline_length":"sightline length (kpc)"}

projects = ['VELA','MOCK','AGORA']
all_metadata_quantities = {}

all_metadata_quantities['VELA'] = ["a", "center_x", "center_y", "center_z", "L_x", "L_y", "L_z", "bulk_velocity_x",\
                            "bulk_velocity_y", "bulk_velocity_z", "Rvir", "Mvir", "Lmag", "sfr", "Mstar", "Mgas", \
                            "compaction_stage",'redshift']
all_metadata_quantities['MOCK'] = ["a", "Rvir", "center_x", "center_y", "center_z", "L_x", "L_y", "L_z", "model_name",\
                          "box_size", "Mvir", "z", "stream_density_beta", "bulk_density_beta",\
                          "stream_temperature_beta", "bulk_temperature_beta", "stream_metallicity_beta", \
                          "interface_metallicity_beta", "bulk_metallicity_beta", "n", "interface_thickness", \
                          "stream_metallicity", "interface_metallicity", "bulk_metallicity", "stream_rotation", \
                          "endpoint", "dist_method", "n_streams", "stream_size_growth", "stream_width", \
                          "stream_temperature", "bulk_temperature", "stream_density", "bulk_density",'redshift']
all_metadata_quantities['AGORA'] = ["a", "center_x", "center_y", "center_z", "Rvir", "L_x", "L_y", "L_z", \
                                    "bulk_velocity_x", "bulk_velocity_y", "bulk_velocity_z", "unitary_conversion",\
                                    "Mgas", "Mdm", "Lmag", "Mstar", "sfr",'redshift']

all_metadata_quantities_all_projects = []
for project in projects:
    all_metadata_quantities_all_projects+=all_metadata_quantities[project]
all_metadata_quantities_all_projects = list(set(all_metadata_quantities_all_projects))

string_type_metadata_quantities = ["compaction_stage","model_name","dist_method"]

param_xVars = [var for var in all_metadata_quantities_all_projects if var not in string_type_metadata_quantities]

param_unit_labels = {"redshift":"z",\
                     "a":"a",\
                     "Rvir":'Virial radius (kpc)',\
                     "Mvir":"Virial Mass (Msun)",\
                     "Mgas":"Gas Mass within Rvir (Msun)",\
                     "Mstar":"Stellar Mass within Rvir (Msun)",\
                     "sfr":"Star Formation Rate (Msun yr-1)",\
                     "ssfr":"Specific Star Formation Rate (Msun yr-1 Mstar-1)",\
                     "Lmag":"Magnitude of Angular Momentum",\
                     "unitary_conversion":"Unitary conversion"}
def get_param_unit_label(quantity):
    try:
        return param_unit_labels[quantity]
    except KeyError:
        print('param %s is not set up with unique label. Set it up in quasarscan/utils/variable_lists.py'%quantity)
        return quantity
    
all_known_variables = stringcriteria+intensives+sightline_xVars+all_metadata_quantities_all_projects