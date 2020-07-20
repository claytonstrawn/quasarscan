#***import statements***
import sys
import yt
import trident
import numpy as np
import matplotlib.pyplot as plt
# from quasarscan import parse_metadata
import parse_metadata

#define analysis functions
print('starting the script...')

def rahuls_function(galaxy_name,galaxy_file_loc,options):
	print("test")

def sallys_function(simulation, filename, redshift):
    # print("loading file, calculating Rvir and Mvir")
    # loading the file, calculate Rvir and Mvir
    ds = yt.load(filename)
    redshift = float(redshift)
    Rvir = parse_metadata.get_value('Rvir',simulation, redshift = redshift)
    Mvir = '%e'%parse_metadata.get_value('Mvir',simulation, redshift = redshift)
    # print("Rvir" + str(Rvir))
    # print("Mvir" + str(Mvir))
    # add all the oxygen ion fields
    # from quasarscan import code_specific_setup
    import code_specific_setup
    ds,_ = code_specific_setup.load_and_setup(filename,'art',
                                             ['O I', 'O II', 'O III', 'O IV', 'O V', 'O VI', 'O VII', 'O VIII', 'O IX'], 
                                             add_pi_fracs = True)
    # print("adding ion fields")
    trident.add_ion_fields(ds,['O I'])
    trident.add_ion_fields(ds,['O II'])
    trident.add_ion_fields(ds,['O III'])
    trident.add_ion_fields(ds,['O IV'])
    trident.add_ion_fields(ds,['O V'])
    trident.add_ion_fields(ds,['O VI'])
    trident.add_ion_fields(ds,['O VII'])
    trident.add_ion_fields(ds,['O VIII'])
    trident.add_ion_fields(ds,['O IX'])
    
    # adding fields for PI and CI O_mass
    # print("defining and adding new fields for PI CI mass")
    def O_p0_PI_mass(field,data):
        return data['O_p0_mass']*data['PI_OI']
    ds.add_field(('gas','O_p0_PI_mass'),
                   sampling_type="cell",
                   function=O_p0_PI_mass,
                   units='g',force_override=True)

    def O_p1_PI_mass(field,data):
        return data['O_p1_mass']*data['PI_OII']
    ds.add_field(('gas','O_p1_PI_mass'),
                   sampling_type="cell",
                   function=O_p1_PI_mass,
                   units='g',force_override=True)

    def O_p2_PI_mass(field,data):
        return data['O_p2_mass']*data['PI_OIII']
    ds.add_field(('gas','O_p2_PI_mass'),
                   sampling_type="cell",
                   function=O_p2_PI_mass,
                   units='g',force_override=True)

    def O_p3_PI_mass(field,data):
        return data['O_p3_mass']*data['PI_OIV']
    ds.add_field(('gas','O_p3_PI_mass'),
                   sampling_type="cell",
                   function=O_p3_PI_mass,
                   units='g',force_override=True)

    def O_p4_PI_mass(field,data):
        return data['O_p4_mass']*data['PI_OV']
    ds.add_field(('gas','O_p4_PI_mass'),
                   sampling_type="cell",
                   function=O_p4_PI_mass,
                   units='g',force_override=True)

    def O_p5_PI_mass(field,data):
        return data['O_p5_mass']*data['PI_OVI']
    ds.add_field(('gas','O_p5_PI_mass'),
                   sampling_type="cell",
                   function=O_p5_PI_mass,
                   units='g',force_override=True)

    def O_p6_PI_mass(field,data):
        return data['O_p6_mass']*data['PI_OVII']
    ds.add_field(('gas','O_p6_PI_mass'),
                   sampling_type="cell",
                   function=O_p6_PI_mass,
                   units='g',force_override=True)

    def O_p7_PI_mass(field,data):
        return data['O_p7_mass']*data['PI_OVIII']
    ds.add_field(('gas','O_p7_PI_mass'),
                   sampling_type="cell",
                   function=O_p7_PI_mass,
                   units='g',force_override=True)

    def O_p8_PI_mass(field,data):
        return data['O_p8_mass']*data['PI_OIX']
    ds.add_field(('gas','O_p8_PI_mass'),
                   sampling_type="cell",
                   function=O_p8_PI_mass,
                   units='g',force_override=True)
    
    def O_p0_CI_mass(field,data):
        return data['O_p0_mass']*data['CI_OI']
    ds.add_field(('gas','O_p0_CI_mass'),
                   sampling_type="cell",
                   function=O_p0_CI_mass,
                   units='g',force_override=True)

    def O_p1_CI_mass(field,data):
        return data['O_p1_mass']*data['CI_OII']
    ds.add_field(('gas','O_p1_CI_mass'),
                   sampling_type="cell",
                   function=O_p1_CI_mass,
                   units='g',force_override=True)

    def O_p2_CI_mass(field,data):
        return data['O_p2_mass']*data['CI_OIII']
    ds.add_field(('gas','O_p2_CI_mass'),
                   sampling_type="cell",
                   function=O_p2_CI_mass,
                   units='g',force_override=True)

    def O_p3_CI_mass(field,data):
        return data['O_p3_mass']*data['CI_OIV']
    ds.add_field(('gas','O_p3_CI_mass'),
                   sampling_type="cell",
                   function=O_p3_CI_mass,
                   units='g',force_override=True)

    def O_p4_CI_mass(field,data):
        return data['O_p4_mass']*data['CI_OV']
    ds.add_field(('gas','O_p4_CI_mass'),
                   sampling_type="cell",
                   function=O_p4_CI_mass,
                   units='g',force_override=True)

    def O_p5_CI_mass(field,data):
        return data['O_p5_mass']*data['CI_OVI']
    ds.add_field(('gas','O_p5_CI_mass'),
                   sampling_type="cell",
                   function=O_p5_CI_mass,
                   units='g',force_override=True)

    def O_p6_CI_mass(field,data):
        return data['O_p6_mass']*data['CI_OVII']
    ds.add_field(('gas','O_p6_CI_mass'),
                   sampling_type="cell",
                   function=O_p6_CI_mass,
                   units='g',force_override=True)

    def O_p7_CI_mass(field,data):
        return data['O_p7_mass']*data['CI_OVIII']
    ds.add_field(('gas','O_p7_CI_mass'),
                   sampling_type="cell",
                   function=O_p7_CI_mass,
                   units='g',force_override=True)

    def O_p8_CI_mass(field,data):
        return data['O_p8_mass']*data['CI_OIX']
    ds.add_field(('gas','O_p8_CI_mass'),
                   sampling_type="cell",
                   function=O_p8_CI_mass,
                   units='g',force_override=True)
    
    # print("Making CGM and galaxy spheres")
    # making the CGM sphere
    sp = ds.sphere('m', (Rvir, "kpc"))
    
    # making the galaxy sphere
    gal = ds.sphere('m', (0.1*Rvir,"kpc"))
    
    # calculating the PI / CI mass in the CGM
    O_PI_mass = []
    O_CI_mass = []
    O_mass = 0
    # print("Calculating PI and CI mass for total")
    for i in range(9):
        name = "O_p"+str(i)+"_PI_mass"
        oxygen_mass = sp.quantities.total_quantity([name])
        O_PI_mass.append(oxygen_mass)
        name = "O_p"+str(i)+"_CI_mass"
        oxygen_mass = sp.quantities.total_quantity([name])
        O_CI_mass.append(oxygen_mass)
        O_mass = O_mass + O_PI_mass[i] + O_CI_mass[i]
        # print(str(i) + " : PI = " + str(O_PI_mass[i]) + " CI = " + str(O_CI_mass[i]))
    O_PI_mass = 10**(np.log10(O_PI_mass))
    O_CI_mass = 10**(np.log10(O_CI_mass))
    O_mass = 10**(np.log10(O_mass))
    # print(O_PI_mass)
    # print(O_CI_mass)
    # print(O_mass)
    
    # print("Calculating PI and CI mass for galaxy")
    # calculating the PI / CI mass in the galaxy
    O_PI_gal_mass = []
    O_CI_gal_mass = []
    O_gal_mass = 0
    for i in range(9):
        name = "O_p"+str(i)+"_PI_mass"
        oxygen_mass = gal.quantities.total_quantity([name])
        O_PI_gal_mass.append(oxygen_mass)
        name = "O_p"+str(i)+"_CI_mass"
        oxygen_mass = gal.quantities.total_quantity([name])
        O_CI_gal_mass.append(oxygen_mass)
        O_gal_mass = O_gal_mass + O_PI_gal_mass[i] + O_CI_gal_mass[i]
        # print(str(i) + " : PI = " + str(O_PI_gal_mass[i]) + " CI = " + str(O_CI_mass[i]))
    O_PI_gal_mass = 10**(np.log10(O_PI_gal_mass))
    O_CI_gal_mass = 10**(np.log10(O_CI_gal_mass))
    O_gal_mass = 10**(np.log10(O_gal_mass))
    # print(O_PI_gal_mass)
    # print(O_CI_gal_mass)
    # print(O_gal_mass)
    
    # print("calculating the fractions")
    # calculating the fractions 
    O_PI_nogal_fraction = []
    O_CI_nogal_fraction = []
    O_nogal_fraction = []
    for i in range(9):
        O_PI_nogal_fraction.append((O_PI_mass[i] - O_PI_gal_mass[i])/(O_mass-O_gal_mass))
        O_CI_nogal_fraction.append((O_CI_mass[i] - O_CI_gal_mass[i])/(O_mass-O_gal_mass))
        O_nogal_fraction.append(O_PI_nogal_fraction[i]+O_CI_nogal_fraction[i])
    print("PI fractions:")
    print(O_PI_nogal_fraction)
    print("CI fractions:")
    print(O_CI_nogal_fraction)
    print("total fractions:")
    print(O_nogal_fraction)
    
    # print("plotting")
    # plotting the fraction of ion levels linear
    x=[1,2,3,4,5,6,7,8,9]
    
    print("Simulation = " + simulation)
    print("Redshift = " + str(redshift))
    print("Rvir = " + str(Rvir))
    print("Mvir = " + str(Mvir))
    
    # plotting the fraction of ion levels logarithmic
    plt.plot(x, np.log10(O_PI_nogal_fraction), label="PI oxygen mass",color='b',linewidth=1,marker='o',linestyle='-')
    plt.plot(x, np.log10(O_CI_nogal_fraction), label="CI oxygen mass",color='r',linewidth=1,marker='o',linestyle='-')
    plt.plot(x, np.log10(O_nogal_fraction), label="total oxygen mass",color='g',linewidth=1,marker='o',linestyle='-')
    plt.legend(loc=0, fontsize=10)
    # plt.savefig(' ')

if __name__ == '__main__':
	# if you're running this from terminal like "python ionization_state_ytanalysis.py galname file_loc"
	galaxy_name = sys.argv[1]
	galaxy_file_loc = sys.argv[2]
	options = sys.argv[3]
	rahuls_function(galaxy_name,galaxy_file_loc,options)
	sallys_function(galaxy_name,galaxy_file_loc,options)





