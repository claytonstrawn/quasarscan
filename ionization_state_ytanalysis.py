#***import statements***
import sys
import yt
import trident
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
import os
import shutil
# from quasarscan import parse_metadata
import parse_metadata

#define analysis functions
print('starting the script...')

def rahuls_function(simulation,filename,redshift):

    #create folder for images
    dir_1 = 'O_num_density_plots'
    if os.path.exists(dir_1):
      shutil.rmtree(dir_1)
    os.makedirs(dir_1) 
    
    dir_2 = simulation + "_" + str(redshift)
    if os.path.exists(dir_2):
      shutil.rmtree(dir_2)
    os.makedirs(dir_2) 

    dir_3 = 'O_distribution_plots'
    if os.path.exists(dir_3):
      shutil.rmtree(dir_3)
    os.makedirs(dir_3)

    #initialize file
    new_file_name = simulation + "_" + str(redshift) + "/o_number_density_data.txt"
    f = open(new_file_name, "x")
    f.close()
   
    #load file, calculate Rvir and center
    ds = yt.load(filename)
    redshift = float(redshift)
    Rvir = parse_metadata.get_value('Rvir',simulation,redshift)
    center = ds.find_max('density')[1]

    #add fields
    trident.add_ion_fields(ds, ['O I'])
    trident.add_ion_fields(ds, ['O II'])
    trident.add_ion_fields(ds, ['O III'])
    trident.add_ion_fields(ds, ['O IV'])
    trident.add_ion_fields(ds, ['O V'])
    trident.add_ion_fields(ds, ['O VI'])
    trident.add_ion_fields(ds, ['O VII'])
    trident.add_ion_fields(ds, ['O VIII'])
    trident.add_ion_fields(ds, ['O IX'])

    #removing colorbar, tick marks, etc
    def remove_extraneous(plot):
        plot.hide_colorbar()
        plot.hide_axes()
        plot.set_minorticks('all', False)
        plot.save()
        return plot

    #create custom colormap with 15 colors
    yt.make_colormap([('blue',  1), ('green', 1), ('red',   1), ('black', 1), ('gray', 1),
              ('purple', 1), ('orange', 1),('yellow', 1),('dgray', 1),('dblue', 1),
              ('dpurple', 1),('dred', 1),('dorange', 1),('dyellow', 1),('dgreen', 1)], 
               name='15', interpolate=False)

    proj_OI = yt.ProjectionPlot(ds,'x',('gas', 'O_p0_number_density'),center=center, width = (2* Rvir, 'kpc'))
    proj_OI.set_cmap(('gas', 'O_p0_number_density'), ('15'))
    proj_OI.set_zlim(('gas', 'O_p0_number_density'),10**3, 10**18)
    proj_OI = remove_extraneous(proj_OI)
    proj_OI.save('O_num_density_plots/OI.png')

    proj_OII = yt.ProjectionPlot(ds,'x',('gas', 'O_p1_number_density'),center=center, width = (2* Rvir, 'kpc'))
    proj_OII.set_cmap(('gas', 'O_p1_number_density'), ('15'))
    proj_OII.set_zlim(('gas', 'O_p1_number_density'),10**3, 10**18)
    proj_OII = remove_extraneous(proj_OII)
    proj_OII.save('O_num_density_plots/OII.png')

    proj_OIII = yt.ProjectionPlot(ds,'x',('gas', 'O_p2_number_density'),center=center, width = (2* Rvir, 'kpc'))
    proj_OIII.set_cmap(('gas', 'O_p2_number_density'), ('15'))
    proj_OIII.set_zlim(('gas', 'O_p2_number_density'),10**3, 10**18)
    proj_OIII = remove_extraneous(proj_OIII)
    proj_OIII.save('O_num_density_plots/OIII.png')

    proj_OIV = yt.ProjectionPlot(ds,'x',('gas', 'O_p3_number_density'),center=center, width = (2* Rvir, 'kpc'))
    proj_OIV.set_cmap(('gas', 'O_p3_number_density'), ('15'))
    proj_OIV.set_zlim(('gas', 'O_p3_number_density'),10**3, 10**18)
    proj_OIV = remove_extraneous(proj_OIV)
    proj_OIV.save('O_num_density_plots/OIV.png')

    proj_OV = yt.ProjectionPlot(ds,'x',('gas', 'O_p4_number_density'),center=center, width = (2* Rvir, 'kpc'))
    proj_OV.set_cmap(('gas', 'O_p4_number_density'), ('15'))
    proj_OV.set_zlim(('gas', 'O_p4_number_density'),10**3, 10**18)
    proj_OV = remove_extraneous(proj_OV)
    proj_OV.save('O_num_density_plots/OV.png')

    proj_OVI = yt.ProjectionPlot(ds,'x',('gas', 'O_p5_number_density'),center=center, width = (2* Rvir, 'kpc'))
    proj_OVI.set_cmap(('gas', 'O_p5_number_density'), ('15'))
    proj_OVI.set_zlim(('gas', 'O_p5_number_density'),10**3, 10**18)
    proj_OVI = remove_extraneous(proj_OVI)
    proj_OVI.save('O_num_density_plots/OVI.png')

    proj_OVII = yt.ProjectionPlot(ds,'x',('gas', 'O_p6_number_density'),center=center, width = (2* Rvir, 'kpc'))
    proj_OVII.set_cmap(('gas', 'O_p6_number_density'), ('15'))
    proj_OVII.set_zlim(('gas', 'O_p6_number_density'),10**3, 10**18)
    proj_OVII = remove_extraneous(proj_OVII)
    proj_OVII.save('O_num_density_plots/OVII.png')

    proj_OVIII = yt.ProjectionPlot(ds,'x',('gas', 'O_p7_number_density'),center=center, width = (2* Rvir, 'kpc'))
    proj_OVIII.set_cmap(('gas', 'O_p7_number_density'), ('15'))
    proj_OVIII.set_zlim(('gas', 'O_p7_number_density'),10**3, 10**18)
    proj_OVIII = remove_extraneous(proj_OVIII)
    proj_OVIII.save('O_num_density_plots/OVIII.png')

    proj_OIX = yt.ProjectionPlot(ds,'x',('gas', 'O_p8_number_density'),center=center, width = (2* Rvir, 'kpc'))
    proj_OIX.set_cmap(('gas', 'O_p8_number_density'), ('15'))
    proj_OIX.set_zlim(('gas', 'O_p8_number_density'),10**3, 10**18)
    proj_OIX = remove_extraneous(proj_OIX)
    proj_OIX.save('O_num_density_plots/OIX.png')

    extract_colors_to_file('O_num_density_plots/OI.png', new_file_name, OI)
    extract_colors_to_file('O_num_density_plots/OII.png', new_file_name, OII)
    extract_colors_to_file('O_num_density_plots/OIII.png', new_file_name, OIII)
    extract_colors_to_file('O_num_density_plots/OIV.png', new_file_name, OIV)
    extract_colors_to_file('O_num_density_plots/OV.png', new_file_name, OV)
    extract_colors_to_file('O_num_density_plots/OVI.png', new_file_name, OVI)
    extract_colors_to_file('O_num_density_plots/OVII.png', new_file_name, OVII)
    extract_colors_to_file('O_num_density_plots/OVIII.png', new_file_name, OVIII)
    extract_colors_to_file('O_num_density_plots/OIX.png', new_file_name, OIX)

    OI_plot = extract_colors_to_return('O_num_density_plots/OI.png')
    OII_plot = extract_colors_to_return('O_num_density_plots/OII.png')
    OIII_plot = extract_colors_to_return('O_num_density_plots/OIII.png')
    OIV_plot = extract_colors_to_return('O_num_density_plots/OIV.png')
    OV_plot = extract_colors_to_return('O_num_density_plots/OV.png')
    OVI_plot = extract_colors_to_return('O_num_density_plots/OVI.png')
    OVII_plot = extract_colors_to_return('O_num_density_plots/OVII.png')
    OVIII_plot = extract_colors_to_return('O_num_density_plots/OVIII.png')
    OIX_plot = extract_colors_to_return('O_num_density_plots/OIX.png')

    labels = ['3','4','5','6','7','8','9','10','11','12','13','14','15','16','17', 'none']

    OI_graph = plt.plot(labels, OI_plot, label = 'OI')
    plt.legend()
    plt.savefig('O_distribution_plots/OI.png')
    
    plt.plot(labels, OII_plot, label = 'OII')
    plt.legend()
    plt.savefig('O_distribution_plots/OII.png')

    plt.plot(labels, OIII_plot, label = 'OIII')
    plt.legend()
    plt.savefig('O_distribution_plots/OIII.png')
    
    plt.plot(labels, OIV_plot, label = 'OIV')
    plt.legend()
    plt.savefig('O_distribution_plots/OIV.png')
    
    plt.plot(labels, OV_plot, label = 'OV')
    plt.legend()
    plt.savefig('O_distribution_plots/OV.png')
    
    plt.plot(labels, OVI_plot, label = 'OVI')
    plt.legend()
    plt.savefig('O_distribution_plots/OVI.png')
    
    plt.plot(labels, OVII_plot, label = 'OVII')
    plt.legend()
    plt.savefig('O_distribution_plots/OVII.png')
    
    plt.plot(labels, OVIII_plot, label = 'OVIII')
    plt.legend()
    plt.savefig('O_distribution_plots/OVIII.png')
    
    plt.plot(labels, OIX_plot, label = 'OIX')
    plt.legend()
    plt.savefig('O_distribution_plots/OIX.png')
  

#extracting colors from saved images
def extract_colors_to_file(img_path, file, oxygen_state):
    im = Image.open(img_path)   
    width = im.width
    height = im.height
    color_map = {}
    for i in range(width):
        for j in range(height):
            current_color = im.getpixel((i,j))
            if current_color in color_map.keys():
                color_map[current_color]+=1
            else:
                color_map[current_color]=1
    total = 0
    for k in color_map.keys():
        total+=color_map[k]

    color_map = {k: (v/total) * 100 for k, v in color_map.items()}
    
    values = {(0,0,255,255): '10^3-4', (0,255,0,255): '10^4-5',
              (255,0,0,255): '10^5-6', (0,0,0,255): '10^6-7',
              (130,130,130,255): '10^7-8', (100,0,200,255): '10^8-9',
              (255,128,0,255): '10^9-10', (255,255,0,255): '10^10-11',
              (80,80,80,255): '10^11-12', (0,0,160,255): '10^12-13',
              (66,0,133,255): '10^13-14', (160,0,0,255): '10^14-15',
              (200,100,0,255): '10^15-16', (200,200,0,255): '10^16-17',
              (0,160,0,255): '10^17-18', (255, 255, 255, 255): 'none'}
    
    new_values = {}
    for key in color_map:
        new_values[values[key]] = color_map[key]
       
    sorted_order = ['10^3-4', '10^4-5', 
                    '10^5-6', '10^6-7', 
                    '10^7-8', '10^8-9', 
                    '10^9-10', '10^10-11', 
                    '10^11-12', '10^12-13', 
                    '10^13-14', '10^14-15', 
                    '10^15-16', '10^16-17', 
                    '10^17-18', 'none']
    graph_values = []
    output = []
    for key in sorted_order:
        if key in new_values.keys():
            output.append([key, new_values[key]])
            graph_values.append(new_values[key])
        else:
            output.append([key, 0])
            graph_values.append(0)

    f = open(file, 'a')
    f.write(oxygen_state + '\n' + str(output) + '\n \n')
    f.close()

    return(graph_values)

def extract_colors_to_return(img_path):
    im = Image.open(img_path)   
    width = im.width
    height = im.height
    color_map = {}
    for i in range(width):
        for j in range(height):
            current_color = im.getpixel((i,j))
            if current_color in color_map.keys():
                color_map[current_color]+=1
            else:
                color_map[current_color]=1
    total = 0
    for k in color_map.keys():
        total+=color_map[k]

    color_map = {k: (v/total) * 100 for k, v in color_map.items()}
    
    values = {(0,0,255,255): '10^3-4', (0,255,0,255): '10^4-5',
              (255,0,0,255): '10^5-6', (0,0,0,255): '10^6-7',
              (130,130,130,255): '10^7-8', (100,0,200,255): '10^8-9',
              (255,128,0,255): '10^9-10', (255,255,0,255): '10^10-11',
              (80,80,80,255): '10^11-12', (0,0,160,255): '10^12-13',
              (66,0,133,255): '10^13-14', (160,0,0,255): '10^14-15',
              (200,100,0,255): '10^15-16', (200,200,0,255): '10^16-17',
              (0,160,0,255): '10^17-18', (255, 255, 255, 255): 'none'}
    
    new_values = {}
    for key in color_map:
        new_values[values[key]] = color_map[key]
       
    sorted_order = ['10^3-4', '10^4-5', 
                    '10^5-6', '10^6-7', 
                    '10^7-8', '10^8-9', 
                    '10^9-10', '10^10-11', 
                    '10^11-12', '10^12-13', 
                    '10^13-14', '10^14-15', 
                    '10^15-16', '10^16-17', 
                    '10^17-18', 'none']
    graph_values = []
    output = []
    for key in sorted_order:
        if key in new_values.keys():
            output.append([key, new_values[key]])
            graph_values.append(new_values[key])
        else:
            output.append([key, 0])
            graph_values.append(0)

    return(graph_values)


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
    
    O_PI_cgm_mass = []
    O_CI_cgm_mass = []
    O_cgm_mass = []
    for i in range(9):
        O_PI_cgm_mass.append(O_PI_mass[i] - O_PI_gal_mass[i])
        O_CI_cgm_mass.append(O_CI_mass[i] - O_CI_gal_mass[i])
        O_cgm_mass.append(O_PI_cgm_mass[i] + O_CI_cgm_mass[i])
    
    # print("plotting")
    # plotting the fraction of ion levels linear
    x=['O I', 'O II', 'O III', 'O IV', 'O V', 'O VI', 'O VII', 'O VIII', 'O IX']
    
    print("Simulation = " + simulation)
    print("Redshift = " + str(redshift))
    print("Rvir = " + str(Rvir))
    print("Mvir = " + str(Mvir))
    
    # plotting the fraction of ion levels logarithmic
    plt.plot(x, np.log10(O_PI_nogal_fraction), label="PI oxygen mass",color='b',linewidth=1,marker='o',linestyle='-')
    plt.plot(x, np.log10(O_CI_nogal_fraction), label="CI oxygen mass",color='r',linewidth=1,marker='o',linestyle='-')
    plt.plot(x, np.log10(O_nogal_fraction), label="total oxygen mass",color='g',linewidth=1,marker='o',linestyle='-')
    plt.legend(loc=0, fontsize=10)
    plt.xlabel("Oxygen ion")
    plt.ylabel("ion fraction of total oxygen mass in log10")
    plt.title(simulation + " z = " + str(redshift) + " Rvir = " + str(Rvir) + " Mvir = " + str(Mvir))
    #plt.show()
    figname = 'ion_state_ytanalysis/' + simulation + "_" + str(redshift) + "/o_ion_fraction_plot.png"
    plt.savefig(figname)
    plt.clf()
    dataname = 'ion_state_ytanalysis/' + simulation + "_" + str(redshift) + "/o_ion_fraction_data.txt"
    f = open(dataname, "x")
    f.close()
    f = open(dataname, "a")
    f.write("PI fraction:")
    f.write(str(O_PI_nogal_fraction))
    f.write("\n \n CI fraction:")
    f.write(str(O_CI_nogal_fraction))
    f.write("\n \n total fraction:")
    f.write(str(O_nogal_fraction))
    f.write("\n \n Rvir: ")
    f.write(str(Rvir))
    f.write("\n \n Mvir: ")
    f.write(str(Mvir))
    f.close()
    plt.plot(x, O_PI_cgm_mass, label="PI oxygen mass",color='b',linewidth=1,marker='o',linestyle='-')
    plt.plot(x, O_CI_cgm_mass, label="CI oxygen mass",color='r',linewidth=1,marker='o',linestyle='-')
    plt.plot(x, O_cgm_mass, label="total oxygen mass",color='g',linewidth=1,marker='o',linestyle='-')
    plt.legend(loc=0, fontsize=10)
    plt.xlabel("Oxygen ion")
    plt.ylabel("total mass of oxygen ion")
    plt.title(simulation + " z = " + str(redshift) + " Rvir = " + str(Rvir) + " Mvir = " + str(Mvir))
    #plt.show()
    plotname = 'ion_state_ytanalysis/' + simulation + "_" + str(redshift) + "/o_ion_mass_plot.png"
    plt.savefig(plotname)
    textname = 'ion_state_ytanalysis/' + simulation + "_" + str(redshift) + "/o_ion_mass_data.txt"
    f1 = open(textname, "x")
    f1.close()
    f1 = open(textname, "a")
    f1.write("PI mass:")
    f1.write(str(O_PI_cgm_mass))
    f1.write("\n \n CI mass:")
    f1.write(str(O_CI_cgm_mass))
    f1.write("\n \n total mass:")
    f1.write(str(O_cgm_mass))
    f1.write("\n \n Rvir: ")
    f1.write(str(Rvir))
    f1.write("\n \n Mvir: ")
    f1.write(str(Mvir))
    f1.close()
    

if __name__ == '__main__':
	# if you're running this from terminal like "python ionization_state_ytanalysis.py galname file_loc"
	galaxy_name = sys.argv[1]
	galaxy_file_loc = sys.argv[2]
	options = sys.argv[3]
	rahuls_function(galaxy_name,galaxy_file_loc,options)
	sallys_function(galaxy_name,galaxy_file_loc,options)





