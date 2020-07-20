#***import statements***
import sys

#define analysis functions
print('starting the script...')

def rahuls_function(galaxy_name,galaxy_file_loc,options):
	print("test")

def sallys_function(galaxy_name,galaxy_file_loc,options):
	#do other things
	#save some images and 
	# save some raw data values

if __name__ == '__main__':
	# if you're running this from terminal like "python ionization_state_ytanalysis.py galname file_loc"
	galaxy_name = sys.argv[1]
	galaxy_file_loc = sys.argv[2]
	rahuls_function(galaxy_name,galaxy_file_loc)
	sallys_function(galaxy_name,galaxy_file_loc)





