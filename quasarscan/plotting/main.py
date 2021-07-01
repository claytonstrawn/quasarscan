#skeleton code for quasarscan plotting

from quasarscan.plotting.multi_quasar_sphere_plotter import MultiQuasarSpherePlotter

def create_mq(loadsim='all',**kwargs):
    return MultiQuasarSpherePlotter(loadsim=loadsim,**kwargs)
