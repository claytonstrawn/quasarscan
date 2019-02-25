"""
-ART-I (Still running in this folder)

/project/projectdirs/agora/ART_I/Cosmo_v8/

-RAMSES (Still running in this folder)

/project/projectdirs/agora/RAMSES/Cosmo_run_v2/

-ENZO

/project/projectdirs/agora/paper_CGM/Cal_step2/ENZO/DD0008/

-GADGET-3

/project/projectdirs/agora/paper_CGM/Cal_step2/GADGET-3/v2/zXX

-GEAR

/project/projectdirs/agora/paper_CGM/Cal_step2/GEAR/v3/own/

-GIZMO

/project/projectdirs/agora/paper_CGM/Cal_step2/GIZMO/cos-mech/output/
"""

filenames = {"art-I":"$SCRATCH/AGORAfiles/art-I/10MpcBox_csf512_a0.050.d",\
"ramses":"$SCRATCH/AGORAfiles/ramses/output_00002/info_00002.txt",\
"enzo":"$SCRATCH/AGORAfiles/enzo/DD0008/DD0008",\
"gadget3":"$SCRATCH/AGORAfiles/gadget3/snapshot_021.0.hdf5",\
"gear":"$SCRATCH/AGORAfiles/gear/snapshot_0047.hdf5",\
"gizmo":"$SCRATCH/AGORAfiles/gizmo/snapshot_010.hdf5"}

#agora_script.py
def main(filetype):
	print "attempting to import yt, trident, numpy" 
	import yt
	import trident
	import numpy as np

	filename = filenames[filetype]

	print "attempting to load file %s of type %s"%(filename,filetype)
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
	for filetype in filenames.keys():
		try:
			main(filetype)
		except Exception as e:
			print e
