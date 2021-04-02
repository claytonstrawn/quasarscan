try:
    import yt
    hasyt = True
except ImportError:
    hasyt = False
try:
    import trident
    hastrident = True
except ImportError:
    hastrident = False

#depends on Sean's code (which does not exist yet, may depend on yt)
if hasyt:
    from quasarscan.preprocessing.main import create_metadata_table
else:
    print('yt not found. Cannot create new metadata')

from quasarscan.preprocessing.main import create_qso_endpoints,get_value

#depends on yt, trident. Needs existing unfilled file
#in "quasarscan_data/output." Will run in parallel if given option
if hasyt and hastrident:
    from quasarscan.processing.main import run_sightlines
else:
    print('yt or trident not found. Cannot run new sightlines')

#depends on matplotlib. Needs (several) existing filled files
#in quasarscan_data/output, will load all of them
from quasarscan.plotting.main import create_mq