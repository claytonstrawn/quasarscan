#this file uses both agora_analysis and quasarscan. Important for forcing all codes to use the same parameters as in AGORA
import agora_analysis
from agora_analysis.field_setup.main import load_necessary_fields

def agora_load_and_setup(path,code,redshift,ions=None,add_pi_fracs=True):
    snap = agora_analysis.AgoraSnapshot(code,redshift)
    assert snap.lookup_snap_path() == path,f'path found for redshift {redshift} was {snap.lookup_snap_path()} (expected {path})'
    snap.load_snapshot()
    list_of_fields = ['radial_distance','radial_velocity','agora_cell_volume',
                      'agora_particle_volume']
    if ions is not None:
        for ion in ions:
            if ion not in ['H I','H II']:
                list_of_fields+=['PI_'+ion.replace(' ','')]
    load_necessary_fields(snap,list_of_fields)
    if snap.sampling_type == 'particle':
        fields_to_keep = [('gas','smoothing_length')]
    elif snap.sampling_type == 'cell':
        fields_to_keep = [('gas','cell_volume')]
    fields_to_keep += [('gas','agora_metallicity'),\
                          ('gas','agora_density'),\
                          ('gas','density'),\
                          ('gas','radial_distance'),\
                          ('gas','radial_velocity')]
    if ions is not None:
        for ion in ions:
            if ion not in ['H I','H II']:
                PI_field_name = ('gas','PI_%s'%(ion.replace(' ','')))
                CI_field_name = ('gas','CI_%s'%(ion.replace(' ','')))
                fields_to_keep += [PI_field_name,CI_field_name]    
    return snap.ds,fields_to_keep
