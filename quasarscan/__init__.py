__version__ = "3.0.2"

print('Checking for yt...',end = '')
try:
    import yt
    from yt import load
    hasyt = True
    print('done.')
except ImportError:
    hasyt = False
    print('not found.')
    print('Running quasarscan in readonly mode. Cannot load simulations, create metadata, or process new sightlines. Quasarscan can still plot existing data.')

if hasyt:
    print('Checking for trident and yt_astro_analysis...',end = '')
    try:
        import trident
        from trident import verify
        import yt_astro_analysis
        from yt_astro_analysis import halo_analysis
        hastrident = True
        print('done.')
    except ImportError:
        hastrident = False
        print('not found.')
        print('Running quasarscan in loadonly mode. Can load simulations, but not create new metadata or process new sightlines. Quasarscan can still plot existing data.')

print('Checking for saved output data...',end = '')
path = "quasarscan_data"
import os
if os.path.exists(path):
    print("folder 'quasarscan_data' found.")
else:
    print("no data found.")
    if not hasyt:
        print('Quasarscan has no functionality in readonly mode without data.')

if hasyt and hastrident:
    mode = 'readloadwrite'
elif hasyt and not hastrident:
    mode = 'readload'
else:
    mode = 'read'

if 'write' in mode:
    from quasarscan.preprocessing.main import create_metadata_table
    from quasarscan.processing.main import run_sightlines
if 'load' in mode:
    from quasarscan.preprocessing.main import create_qso_endpoints,get_value
if 'read' in mode:
    from quasarscan.plotting.main import create_mq