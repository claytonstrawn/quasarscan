import os
import quasar_scan
from multi_quasar_sphere_plotter import get_all_textfiles
from parse_vela_metadata import Rdict
import sys

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
final = False

def check_in_allfiles(tocheck,alltextfiles,ionlist):
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
                for ion in ionlist:
                    if not ion in allions:
                        print("that ion %s cannot be made in TRIDENT"%ion)
                        return
                    if "allions" in splitunderscore:
                        continue
                    if "joeions" in splitunderscore and ion in joeions:
                        continue
                    ion = ion.replace(" ","")
                    if ion in splitunderscore:
                        continue
                if lines >= minimumlines and lines >= outOf:
                    return True
                elif lines < minimumlines and lines >= outOf:
                    return False
                elif lines < outOf:
                    return lines
                else:
                    return False
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
    print("do I find %s in %s?"%(tocheck[0],my_saved_directories))
    if not tocheck[0] in my_saved_directories:
        print("no")
        return False
    print("yes")
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
    line = dirname + " " + z
    f.write(line)
    f.close

def write_files(tocheck,cont = 0):
    print("I think we should work on %s where we've so far gotten to %d"%(tocheck,cont))
    v,num,a,z = convert_check_to_strings(tocheck)
    firstline = "#!/bin/bash\n"
    secondline = "quasarscan/./run_one_new_snapshot_nersc.sh %s %s %s %s %s \n"%(v,num,a[2:], z, cont)
    f = open("quasarscan/nextfile.sh")
    currentfirstline = f.read()
    currentsecondline = f.read()
    print currentfirstline
    print secondline
    print secondline in currentfirstline
    if (secondline in currentfirstline) and final:
        print("I already tried that, I guess it didn't work :(")
        add_to_blacklist(tocheck[0],tocheck[1])
        return main_func()
    f.close()
    f = open("quasarscan/nextfile.sh","w+")
    f.write(firstline)
    f.write(secondline)
    f.close()
    
def main_func():
    #figure out what is next necessary file to scan
    #write a bash script to go get it, and to delete it after
    for ionlist in [joeions,allions]:
        for file in nerscsimnames:
            for redshift in knownredshifts:
                tocheck = (file,redshift)
                print tocheck
                isValid = check_in_allfiles(tocheck,alltextfiles,ionlist)
                if type(isValid) is int and check_validity(tocheck):
                    write_files(tocheck, cont = isValid)
                    return
                elif isValid:
                    print "already done"
                    continue
                if check_validity(tocheck):
                    write_files(tocheck)
                    return
                print "can't do it"
        print("Moving on to larger ion list")
    print("Finished, all available files satisfy strictest criteria!")

main_func()
