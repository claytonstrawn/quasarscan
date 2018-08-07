import os
import quasar_scan_testMPI as quasar_scan
from multi_quasar_sphere_plotter import get_all_textfiles
from parse_vela_metadata import Rdict
import sys

alltextfiles = get_all_textfiles()

#NOTE: missing [7,8,9,10,11,12,13,14,15,18]
numsnersc = [1,2,3,4,5,6,16,17,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35]

nerscsimnames = []
for num in numsnersc:
    nerscsimnames.append("VELA_v2_%2d"%num)

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

def check_in_allfiles(tocheck,alltextfiles,ionlist):
    for fil in alltextfiles:
        afteroutput = fil.split("output/")[1]
        aftercoldensinfo = afteroutput.split("coldensinfo/")[1]
        lines = int(aftercoldensinfo.split("_of")[0])
        if lines >= minimumlines:
            velaname = afteroutput.split("coldensinfo")[0]
            if velaname == tocheck[0]:
                afterz = fil.split("z")[1]
                file_redshift = float(afterz.split(".t")[0])
                if abs(file_redshift - tocheck[1]) <= 0.04:
                    splitunderscore = afteroutput.split("_")
                    for ion in ionlist:
                        ion = ion.replace(" ","")
                        if ion in splitunderscore:
                            continue
                        else:
                            return False
                    return True
    return False
    
z2a = {"1.0":"0.500","1.5":"0.400","2.0":"0.330","3.0":"0.250","4.0":"0.200"}
a2z = {"0.500":"1.0","0.400":"1.5","0.330":"2.0","0.250":"3.0","0.200":"4.0"}

def check_validity(tocheck):
    blacklist = open(blacklist.txt)
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
    try:
        table_metadata = Rdict[tocheck[0]][z2a[str(tocheck[1])]]
        return True
    except:
        return False
    
def convert_check_to_strings(tocheck):
    print("convert_check_to_strings not implemented")
    #return str(dirname), str(a)

def add_to_blacklist(dirname,z):
    f = open("blacklist.txt","w+")
    line = dirname + " " + z
    f.write(line)
    f.close

def write_files(tocheck):
    v,num,a,z = convert_check_to_strings(tocheck)
    f = open("nextfile.sh","w+")
    firstline = "#!/bin/bash"
    secondline = "quasarscan/./run_one_new_snapshot_nersc.sh %s %s %s %s"%(v,num,a z)
    f.write(firstline)
    f.write(secondline)
    f.close()

    f = open("deletefile.sh","w+")
    firstline = "#!/bin/bash"
    secondline = "rm -rf %s"%tocheck[0]
    f.write(firstline)
    f.write(secondline)
    f.close()

def main_func():
    if len(sys.argv)>1 and sys.argv[1] == "check":
        #check we successfully made the folder and populated
        #if it doesn't exist, it probably is missing from hsi
        #put into list of failures, to skip in future
        f = open("nextfile.sh")
        f.readline()
        f.readline()
        f.readline()
        linetoparse = f.readline()
        dirname = linetoparse.split(" ")[2][:-1]
        a = linetoparse.split("*a0.")[1][:-1]
        dirs = os.listdir("~")
        if dirname in dirs:
            folderDirs = os.listdir("~/%s"%dirname)
            if len(folderDirs)==0:
                add_to_blacklist(dirname,a)
            elif folderDirs[0].split("a0.")[1][:3] != a:
                add_to_blacklist(dirname,a)
        else:
            print("folder does not exist")

    else:
        #figure out what is next necessary file to scan
        #write a bash script to go get it, and to delete it after
        for ionlist in [joeions,allions]:
            for file in nerscsimnames:
                for redshift in knownredshifts:
                    tocheck = (file,redshift)
                    if check_in_allfiles(tocheck,alltextfiles,ionlist):
                        continue
                    if check_validity(tocheck):
                        write_files(tocheck)
                        return
            print("Moving on to larger ion list")
        print("Finished, all available files satisfy criteria!")

main_func()
