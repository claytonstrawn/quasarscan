import os
import sys
try:
    from quasarscan.multi_quasar_sphere_plotter import get_all_textfiles
    from quasarscan.ion_lists import *
    from quasarscan.parse_metadata import get_value
    level = 0
except:
    from multi_quasar_sphere_plotter import get_all_textfiles
    from ion_lists import *
    from parse_metadata import get_value
    level = 1

alltextfiles = get_all_textfiles(0)

#the list of files we want to process at any given time
list_of_files_to_process = [('07','500'),
                            ('07','330'),
                            ('07','250'),
                            ('07','200'),
                            ('08','460'),
                            ('08','330'),
                            ('08','250'),
                            ('08','200'),
                            ('10','500'),
                            ('10','330'),
                            ('10','250'),
                            ('10','200'),
                            ('21','500'),
                            ('21','330'),
                            ('21','250'),
                            ('21','200'),
                            ('22','500'),
                            ('22','330'),
                            ('22','250'),
                            ('22','200'),
                            ('29','500'),
                            ('29','340'),
                            ('29','200'),
                            ('29','250')]

a2z = {"460":"1.1","500":"1.0","400":"1.5","330":"2.0","340":"2.0","250":"3.0","200":"4.0"}

computer = 'nersc'

root_for_data_files = {'nersc':"/global/cscratch1/sd/cstrawn",'nasa':'/nobackupp2/cstrawn/mydir'}[computer]

ionlist = Strawn20
ionlist_name = "Strawn20"

#turn the list into the filenames you expect
filenames_to_make = []
for t in list_of_files_to_process:
    filenames_to_make.append("quasarscan/output/VELA_v2_art_"+t[0]+"coldensinfo/number_of_lines-"+ionlist_name+'_z'+a2z[t[1]]+'.txt')

# if we successfully cleared 250, we can move on and don't need to start over
minimumlines = 250

#this is just to turn off writing to blacklists while testing
final = True

#can do this just to test, by printing the output instead of writing 
#to actual file (in this case the number of processors won't matter)
if len(sys.argv)>2 and sys.argv[2]=='test':
    test = True
else:
    test = False
procs = int(sys.argv[1])

#check if the file we want to make already exists and has enough lines. 
#if it does, return True, if it doesn't, return False. If it does but 
#doesn't have enough lines, return the number to start at
def check_in_allfiles(alltextfiles,name_to_check):
    for fil in alltextfiles:
        afteroutput = fil.split("output/")[1]
        aftercoldensinfo = afteroutput.split("coldensinfo/")[1]
        lines = int(aftercoldensinfo.split("_of_")[0])
        outOf = int(aftercoldensinfo.split("_of_")[1].split("-")[0])
        if name_to_check.replace('number',str(lines)).replace('lines',str(outOf)) == fil:
            if lines >= minimumlines:
                return True
            elif lines < outOf:
                return lines
        else:
            continue
    return False

    

#check if the data file exists that we need, and if it "appears"
#to be there, but it is in the blacklist, or isn't there, return False
def check_validity(tocheck):
    blacklist = open("quasarscan/blacklist.txt")
    blacklist_list = blacklist.readlines()
    for fil in blacklist_list:
        firstpart = fil.split(" ")[0]
        if firstpart == tocheck[0]:
            if float(fil.split(" ")[1]) == tocheck[1]:
                return False
            else:
                continue
        else:
            continue
    my_saved_directories = os.listdir(root_for_data_files)
    foldername = 'VELA_v2_art_%s'%tocheck[0]
    print("looking for %s/%s/%s"%(root_for_data_files,foldername,"10MpcBox_csf512_a0.%s.d"%tocheck[1]))
    if not foldername in my_saved_directories:
        print("didn't see folder")
        return False
    my_saved_data = os.listdir(root_for_data_files+'/'+foldername)
    if not "10MpcBox_csf512_a0.%s.d"%tocheck[1] in my_saved_data:
        print("didn't see file")
        return False
    if not get_value('Rvir','VELA_v2_art_%s'%tocheck[0],a0=float('0.'+tocheck[1]),check_exists = True):
        print("didn't see metadata")
        return False
    return True

def convert_check_to_strings(tocheck):
    simname = "VELA_v2_art_%s"%tocheck[0]
    filename = "%s/%s/%s"%(root_for_data_files,simname,'10MpcBox_csf512_a0.%s.d'%tocheck[1])
    
    #only on NERSC (temporary, need to make consistent)
    #filename.replace('v2.0_art','v2')
    
    redshift = a2z[tocheck[1]]
    return simname, filename, redshift

def add_to_blacklist(dirname,z):
    f = open("quasarscan/blacklist.txt","a+")
    line = dirname + " " + str(z)+'\n'
    f.write(line)
    f.close

def write_files(tocheck,cont = 0):
    if tocheck is None:
        f = open("quasarscan/nextfile.sh","w+")
        f.close()
        return
    print("I think we should work on %s where we've so far gotten to %d"%(tocheck,cont))
    simname,filename,redshift = convert_check_to_strings(tocheck)
    firstline = "#!/bin/bash"
    secondline = "quasarscan/batch_scripts/./run_one_new_snapshot_%s.sh %s %s %s %s %s"%(computer, simname, filename, redshift, cont, procs)
    f = open("quasarscan/nextfile.sh")
    currentfirstline = f.readline()
    currentsecondline = f.readline()
    print(secondline == currentsecondline.strip()),secondline,currentsecondline
    if (secondline == currentsecondline.strip()) and final and cont == 0 and not test:
        print("I already tried that, I guess it didn't work :(")
        add_to_blacklist(tocheck[0],tocheck[1])
        f.close()
        return main_func()
    f.close()
    if test:
        print(firstline)
        print(secondline)
        return
    f = open("quasarscan/nextfile.sh","w+")
    f.write(firstline+'\n')
    f.write(secondline+'\n')
    f.close()
    
def main_func():
    #figure out what is next necessary file to scan
    #write a bash script to go get it, and to delete it after
    for i,tocheck in enumerate(list_of_files_to_process):
        print(tocheck)
        isValid = check_in_allfiles(alltextfiles, filenames_to_make[i]) 
        if type(isValid) is int and check_validity(tocheck):
            write_files(tocheck, cont = isValid)
            return
        elif type(isValid) is bool and isValid:
            print("already done")
            continue
        if check_validity(tocheck):
            write_files(tocheck)
            return
        print("can't do it")
    write_files(None)

main_func()

