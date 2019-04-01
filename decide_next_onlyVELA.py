import os
from multi_quasar_sphere_plotter import get_all_textfiles
from parse_metadata import get_value
import sys
from ion_lists import *

alltextfiles = get_all_textfiles(False)

#NOTE: missing [7,8,9,10,11,12,13,14,15,18]
numsnersc = [1,2,3,4,5,6,16,17,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35]

nerscsimnames = []
for num in numsnersc:
    nerscsimnames.append("VELA_v2_art_%02d"%num)

knownredshifts = [1.0,1.5,2.0,3.0,4.0]

minimumlines = 250
final = True
if len(sys.argv)>1:
    test = True
else:
    test = False

def check_in_allfiles(tocheck,alltextfiles,ionlist):
    startAt = 0
    for fil in alltextfiles:
        afteroutput = fil.split("output/")[1]
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
    foldername = tocheck[0].replace("_art","")
    if not foldername in my_saved_directories:
        print "didn't see folder"
        return False
    my_saved_data = os.listdir("/global/cscratch1/sd/cstrawn/%s"%foldername)
    if not "10MpcBox_csf512_a%s.d"%z2a[str(tocheck[1])] in my_saved_data:
        print "didn't see file"
        return False
    if not get_value('Rvir',tocheck[0],tocheck[1],check_exists = True):
        print "didn't see metadata"
        return False
    return True
    
def convert_check_to_strings(tocheck):
    simname = tocheck[0]
    foldername = tocheck[0].replace("_art","")
    filename = "/global/cscratch1/sd/cstrawn/%s/%s"%(foldername,'10MpcBox_csf512_a%s.d'%z2a[str(tocheck[1])])
    redshift = tocheck[1]
    return simname, filename, redshift

def add_to_blacklist(dirname,z):
    f = open("quasarscan/blacklist.txt","a+")
    line = dirname + " " + str(z)+'\n'
    f.write(line)
    f.close

def write_files(tocheck,cont = 0):
    print("I think we should work on %s where we've so far gotten to %d"%(tocheck,cont))
    simname,filename,redshift = convert_check_to_strings(tocheck)
    firstline = "#!/bin/bash"
    secondline = "quasarscan/./run_one_new_snapshot_nersc.sh %s %s %s %s"%( simname, filename, redshift, cont)
    f = open("quasarscan/nextfile.sh")
    currentfirstline = f.readline()
    currentsecondline = f.readline()
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
    f = open("quasarscan/nextfile.sh","w+")
    f.write("first attempt\n")
    f.write("first attempt\n")
    f.close()
    #figure out what is next necessary file to scan
    #write a bash script to go get it, and to delete it after
    ionlists = [agoraions]
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