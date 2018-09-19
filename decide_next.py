import os
import quasar_scan
from multi_quasar_sphere_plotter import get_all_textfiles
from parse_vela_metadata import Rdict
import sys
from ion_lists import *

alltextfiles = get_all_textfiles(False)

#NOTE: missing [7,8,9,10,11,12,13,14,15,18]
numsnersc = [1,2,3,4,5,6,16,17,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35]

nerscsimnames = []
for num in numsnersc:
    nerscsimnames.append("VELA_v2_%02d"%num)

knownredshifts = [1.0,1.5,2.0,3.0,4.0]

joeions = ['C III', 'H I', 'Mg II', 'Mg X', 'N II', 'N III', \
           'N IV', 'N V', 'Ne VIII', 'O II', 'O III', 'O IV', \
           'O V', 'O VI', 'S II', 'S III', 'S IV', 'S V', 'S VI']
allions = ['Al II', 'Al III', 'Ar I', 'Ar II', 'Ar VII', 'C I', 'C II', \
           'C III', 'C IV', 'Ca X', 'Fe II', 'H I', 'Mg II', 'Mg X', 'N I', \
           'N II', 'N III', 'N IV', 'N V', 'Na IX', 'Ne V', 'Ne VI', 'Ne VII', \
           'Ne VIII', 'O I', 'O II', 'O III', 'O IV', 'O V', 'O VI', 'P IV', \
           'P V', 'S II', 'S III', 'S IV', 'S V', 'S VI', 'S XIV', 'Si II', \
           'Si III', 'Si IV', 'Si XII']

minimumlines = 250
final = True
if len(sys.argv)>1:
    test = True
else:
    test = False

def check_in_allfiles(tocheck,alltextfiles,ionlist):
    startAt = 0
    for fil in alltextfiles:
        afteroutput = fil.split("output2.0/")[1]
        aftercoldensinfo = afteroutput.split("coldensinfo/")[1]
        lines = int(aftercoldensinfo.split("_of_")[0])
        outOf = int(aftercoldensinfo.split("_of_")[1].split("-")[0])
        velaname = afteroutput.split("coldensinfo")[0]
        if velaname == tocheck[0]:
            afterz = fil.split("z")[1]
            file_redshift = float(afterz.split(".t")[0])
            if abs(file_redshift - tocheck[1]) <= 0.04:
                splitunderscore = afteroutput.split("_")
                ion_ok = True
                for ion in ionlist:
                    if velaname == "VELA_v2_24":
                        print ion
                        print splitunderscore
                    if not ion in allions:
                        print("that ion %s cannot be made in TRIDENT (yet)"%ion)
                        Error
                    if "allions" in splitunderscore[-2]:
                        continue
                    if "joeions" in splitunderscore[-2] and ion in joeions:
                        continue
                    ion = ion.replace(" ","")
                    if ion in splitunderscore:
                        continue
                    ion_ok = False
                if ion_ok and lines >= minimumlines and lines >= outOf:
                    return True
                elif lines < minimumlines and lines >= outOf:
                    continue
                elif lines < outOf and lines > 0:
                    startAt = max(lines,startAt)
    if startAt > 0:
        return startAt
    else:
        return False
    
z2a = {"1.0":"0.500","1.5":"0.400","2.0":"0.330","3.0":"0.250","4.0":"0.200"}
a2z = {"0.500":"1.0","0.400":"1.5","0.330":"2.0","0.250":"3.0","0.200":"4.0"}

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
    my_saved_directories = os.listdir("/global/cscratch1/sd/cstrawn")
    if not tocheck[0] in my_saved_directories:
        return False
    my_saved_data = os.listdir("/global/cscratch1/sd/cstrawn/%s"%tocheck[0])
    if not "10MpcBox_csf512_a%s.d"%z2a[str(tocheck[1])] in my_saved_data:
        return False
    try:
        table_metadata = Rdict[tocheck[0]][z2a[str(tocheck[1])]]
        return True
    except:
        return False
    
def convert_check_to_strings(tocheck):
    v = "2"
    num = tocheck[0][-2:]
    a = z2a[str(tocheck[1])]
    z = str(tocheck[1])
    return v,num,a,z

def add_to_blacklist(dirname,z):
    f = open("quasarscan/blacklist.txt","a+")
    line = dirname + " " + str(z)+'\n'
    f.write(line)
    f.close

def write_files(tocheck,cont = 0):
    print("I think we should work on %s where we've so far gotten to %d"%(tocheck,cont))
    v,num,a,z = convert_check_to_strings(tocheck)
    firstline = "#!/bin/bash"
    secondline = "quasarscan/./run_one_new_snapshot_nersc.sh %s %s %s %s %s"%(v,num,a[2:], z, cont)
    f = open("quasarscan/nextfile.sh")
    currentfirstline = f.readline()
    currentsecondline = f.readline()
    print secondline
    print currentsecondline.strip()
    print secondline == currentsecondline.strip()
    if (secondline == currentsecondline.strip()) and final and cont == 0 and not test:
        print("I already tried that, I guess it didn't work :(")
        add_to_blacklist(tocheck[0],tocheck[1])
        f.close()
        return main_func()
    f.close()
    if test:
        print firstline
        print secondline
        return
    f = open("quasarscan/nextfile.sh","w+")
    f.write(firstline+'\n')
    f.write(secondline+'\n')
    f.close()
    
def main_func():
    #figure out what is next necessary file to scan
    #write a bash script to go get it, and to delete it after
    ionlists = [alloxygens]
    for i in range(len(ionlists)):
        ionlist = ionlists[i]
        for file in nerscsimnames:
            for redshift in knownredshifts:
                tocheck = (file,redshift)
                print tocheck
                isValid = check_in_allfiles(tocheck,alltextfiles,ionlist)
                if type(isValid) is int and check_validity(tocheck):
                    write_files(tocheck, cont = isValid)
                    return
                elif type(isValid) is bool and isValid:
                    print "already done"
                    continue
                if check_validity(tocheck):
                    write_files(tocheck)
                    return
                print "can't do it"
        if i<len(ionlists)-1:
            print("Moving on to larger ion list")
    print("Finished, all available files satisfy strictest criteria!")

main_func()
