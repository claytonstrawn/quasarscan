#nihao_script.py
"""
def main():
	print "attempting to import yt, trident, numpy" 
	import yt
	import trident
	import numpy as np

	filename = "~/quasarscan/sample_file/g1.12e12.00672"

	print "attempting to load file %s"%filename
	ds = yt.load(filename)
	c = ds.find_max(('gas','density'))[1].in_units("kpc")

	print "attempting to add derived gas mass field"
	def gas_mass(field,data):
		return data['deposit','Gas_mass']
	ds.add_field(('gas','mass'),units = 'g', function = gas_mass,sampling_type = 'cell')

	print "attempting to make very small sightline"
	ray_start = c
	ray_end = c+yt.YTArray([1,1,1],'kpc')
	ray = trident.make_simple_ray(ds,start_position = ray_start,\
							end_position = ray_end, fields = [('gas','mass'),('gas','metallicity')])

	print "attempting to add ion fields to sightline"
	ions = ['H I','O VI']
	trident.add_ion_fields(ray,ions)
	ad = ray.all_data()

	print "attempting to import from quasar_scan"
	try:
		from quasar_scan import ion_to_field_name
		print "getting from quasar_scan"
	except:
		from quasarscan.quasar_scan import ion_to_field_name
		print "getting from quasarscan.quasar_scan"

	print "calculating column density"
	output = np.zeros(2)
	for i,ion in enumerate(ions):
		cdens = np.sum(ad[('gas',ion_to_field_name(ion))]*ad['dl'])
		output[i] = cdens
	print output

if __name__ == "__main__":
	main()


print "attempting to import yt, trident, numpy" 
import yt
import trident
import numpy as np

filename = "files_to_process/g2.19e11.00672"

print "attempting to load file %s"%filename
ds = yt.load(filename)
c = ds.find_max(('gas','density'))[1].in_units("kpc")

print "attempting to add derived gas mass field"
def gas_mass(field,data):
	return data['deposit','Gas_mass']
ds.add_field(('gas','mass'),units = 'g', function = gas_mass,sampling_type = 'cell')

print "attempting to make very small sightline, from center"
ray_start,ray_end = c + yt.YTArray([0,1,0],'kpc'), c+yt.YTArray([1,0,0],'kpc')
ray = trident.make_simple_ray(ds,start_position = ray_start,end_position = ray_end, fields = [('gas','mass'),('gas','metallicity')])

print "attempting to make very small sightline, through center"
ray_start,ray_end = c + yt.YTArray([0,1,0],'kpc'), c+yt.YTArray([0,-1,0],'kpc')
ray = trident.make_simple_ray(ds,start_position = ray_start,end_position = ray_end, fields = [('gas','mass'),('gas','metallicity')])

print "attempting to make very small sightline, not through center"
ray_start,ray_end = c + yt.YTArray([0,1,0],'kpc'), c+yt.YTArray([1,0,0],'kpc')
ray = trident.make_simple_ray(ds,start_position = ray_start,end_position = ray_end, fields = [('gas','mass'),('gas','metallicity')])

print "attempting to add ion fields to sightline"
ions = ['H I','O VI']
trident.add_ion_fields(ray,ions)
ad = ray.all_data()

print "attempting to import from quasar_sphere"
try:
	from quasar_sphere import ion_to_field_name
	print "getting from quasar_sphere"
except:
	from quasarscan.quasar_sphere import ion_to_field_name
	print "getting from quasarscan.quasar_sphere"

print "calculating column density"
output = np.zeros(2)
for i,ion in enumerate(ions):
	cdens = np.sum(ad[('gas',ion_to_field_name(ion))]*ad['dl'])
	output[i] = cdens
print output
"""
import yt
import numpy as np

myrange = range(16,1025,16)
mytable = np.zeros((3,len(myrange)))
print "will open %d files"%len(myrange)
for i in myrange:
    mytable[0,i//16] = i
    str_num = str(i).zfill(5)
    ds = yt.load("g2.63e10.%s"%str_num)
    mytable[1,i//16] = ds.current_redshift
    mytable[2,i//16] = ds.current_time
    
    print mytable
    