#nihao_script.py
def main():
	print "attempting to import yt, trident, numpy" 
	import yt
	import trident
	import numpy as np

	filename = "quasarscan/sample_file/g1.12e12.00672"

	print "attempting to load file %s"%filename
	ds = yt.load(filename)
	_,c = ds.find_max('gas','density').in_units("kpc")

	print "attempting to add derived gas mass field"
	def gas_mass(field,data):
		return data['deposit','Gas_mass']
	ds.add_field(('gas','mass'),units = 'g' function = gas_mass,sampling_type = 'cell')

	print "attempting to make very small sightline"
	ray_start = c
	ray_end = c+yt.YTArray([1,1,1],'kpc')
	ray = trident.make_simple_ray(start = ray_start,\
							end = ray_end)
	
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
		cdens = np.sum(field_data[('gas',ion_to_field_name(ion))]*field_data['dl'])
		output[i] = cdens
	print output

if __name__ == "__main__":
	main()
