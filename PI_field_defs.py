import numpy as np
import yt 

def rhoz_func(z):
    k="""0<z<0.05 -- 0/120
    0.05<z<0.15 -- 16/120
    0.15<z<0.25 -- 30/120
    0.25<z<0.35 -- 44/120
    0.35<z<0.45 -- 54/120
    0.45<z<0.55 -- 63/120
    0.55<z<0.65 -- 72/120
    0.65<z<0.75 -- 80/120
    0.75<z<0.85 -- 88/120
    0.85<z<0.95 -- 95/120
    0.95<z<1.05 -- 102/120
    1.05<z<1.15 -- 108/120
    1.15<z<1.25 -- 113/120
    1.25<z<1.35 -- 118/120
    1.35<z<1.45 -- 122/120
    1.45<z<1.55 -- 125/120
    1.55<z<1.65 -- 128/120
    1.65<z<1.75 -- 130/120
    1.75<z<1.85 -- 130/120
    1.85<z<1.95 -- 130/120
    1.95<z<2.05 -- 130/120
    2.05<z<2.15 -- 128/120
    2.15<z<2.25 -- 125/120
    2.25<z<2.35 -- 122/120
    2.35<z<2.45 -- 120/120
    2.45<z<2.55 -- 118/120
    2.55<z<2.65 -- 116/120
    2.65<z<2.75 -- 113/120
    2.75<z<2.85 -- 110/120
    2.85<z<2.95 -- 106/120
    2.95<z<3.05 -- 102/120
    3.05<z<3.15 -- 99/120
    3.15<z<3.25 -- 96/120
    3.25<z<3.35 -- 93/120
    3.35<z<3.45 -- 90/120
    3.45<z<3.55 -- 87/120
    3.55<z<3.65 -- 84/120
    3.65<z<3.75 -- 81/120
    3.75<z<3.85 -- 78/120
    3.85<z<3.95 -- 75/120
    3.95<z<np.inf -- 72/120"""
    zs  = []
    rhos  = []
    for line in k.split("\n"):
        key = float(line.split("<")[0])
        num = float((line.split("- ")[1]).split('/')[0])
        denom = float((line.split("- ")[1]).split('/')[1])
        value = num/denom
        zs.append(key)
        rhos.append(value)
    zs = np.array(zs)
    rhos = np.array(rhos)
    pos = len(zs[z>=zs])-1
    return rhos[pos]


def make_funcs(ds=None,z=None,add_fields=False):
    if not z and not ds:
        print('assuming z=1')
        z=1.0
    elif z is None:
        z = ds.current_redshift
    def PI_OIV(field, data):
        #0 if CI, 1 if PI
        rhoz = np.log10(yt.YTQuantity(rhoz_func(z),'1/cm**3'))
        temps = np.log10(data['gas','temperature'])
        """logT(K)<4.6, for all ρ
        – 4.6< logT(K)<4.7, logρ(cm−3)> -2.75+∆ρz
        – 4.7< logT(K)<4.8, logρ(cm−3)> -3.0+∆ρz
        – 4.8< logT(K)<4.9, logρ(cm−3)> -3.25+∆ρz
        – 4.9< logT(K)<5.0, logρ(cm−3)> -3.5+∆ρz"""
        mask1 = temps < 4.6
        mask2 = np.logical_and(4.6 <= temps,temps< 4.7)
        mask3 = np.logical_and(4.7 <= temps,temps< 4.8)
        mask4 = np.logical_and(4.8 <= temps,temps< 4.9)
        mask5 = np.logical_and(4.9 <= temps,temps< 5.0)
        mask6 = temps > 5.0
        comp = np.zeros(temps.shape)
        comp[mask1] = -np.inf
        comp[mask2] = -2.75+rhoz
        comp[mask3] = -3.0+rhoz
        comp[mask4] = -3.25+rhoz
        comp[mask5] = -3.5+rhoz
        comp[mask6] = np.inf
        tr = yt.YTArray((np.log10(data['gas','density'])>comp).astype(float))
        return tr
    def CI_OIV(field,data):
        tr = data['gas','PI_OIV']
        return 1.-tr

    def PI_OV(field, data):
        #0 if CI, 1 if PI
        rhoz = np.log10(yt.YTQuantity(rhoz_func(z),'1/cm**3'))
        temps = np.log10(data['gas','temperature'])
        """– logT(K)<4.85, all ρ
        – 4.85< logT(K)<4.95, logρ(cm−3)> -3.0+∆ρz
        – 4.95< logT(K)<5.05, logρ(cm−3)> -3.25+∆ρz
        – 5.15< logT(K)<5.25, logρ(cm−3)> -3.5+∆ρz
        – 5.25< logT(K)<5.35, logρ(cm−3)> -3.75+∆ρz"""
        mask1 = temps < 4.85
        mask2 = np.logical_and(4.85 <= temps,temps< 4.95)
        mask3 = np.logical_and(4.95 <= temps,temps< 5.05)
        mask4 = np.logical_and(5.15 <= temps,temps< 5.25)
        mask5 = np.logical_and(5.25 <= temps,temps< 5.35)
        mask6 = temps > 5.35
        comp = np.zeros(temps.shape)
        comp[mask1] = -np.inf
        comp[mask2] = -3.0+rhoz
        comp[mask3] = -3.25+rhoz
        comp[mask4] = -3.5+rhoz
        comp[mask5] = -3.75+rhoz
        comp[mask6] = np.inf
        tr = yt.YTArray((np.log10(data['gas','density'])>comp).astype(float))
        return tr
    def CI_OV(field,data):
        tr = data['gas','PI_OV']
        return 1.-tr

    def PI_OVI(field, data):
        #0 if CI, 1 if PI
        """– logT(K)<5.05, all ρ
        – 5.05< logT(K)<5.15, logρ(cm−3)> -3.75+∆ρz
        – 5.15< logT(K)<5.25, logρ(cm−3)> -4.0+∆ρz
        – 5.25< logT(K)<5.35, logρ(cm−3)> -4.15+∆ρz
        – 5.35< logT(K)<5.45, logρ(cm−3)> -4.2+∆ρz"""
        rhoz = np.log10(yt.YTQuantity(rhoz_func(z),'1/cm**3'))
        temps = np.log10(data['gas','temperature'])
        mask1 = temps < 5.05
        mask2 = np.logical_and(5.05 <= temps,temps< 5.15)
        mask3 = np.logical_and(5.15 <= temps,temps< 5.25)
        mask4 = np.logical_and(5.25 <= temps,temps< 5.35)
        mask5 = np.logical_and(5.35 <= temps,temps< 5.45)
        mask6 = temps > 5.45
        comp = np.zeros(temps.shape)
        comp[mask1] = -np.inf
        comp[mask2] = -3.75+rhoz
        comp[mask3] = -4.0+rhoz
        comp[mask4] = -4.15+rhoz
        comp[mask5] = -4.2+rhoz
        comp[mask6] = np.inf
        tr = yt.YTArray((np.log10(data['gas','density'])>comp).astype(float))
        return tr
    def CI_OVI(field,data):
        tr = data['gas','PI_OVI']
        return 1.-tr

    def PI_OVII(field, data):
        #0 if CI, 1 if PI
        """– logT(K)<5.2, all ρ
        – 5.2< logT(K)<5.3, logρ(cm−3)> -3.8+∆ρz
        – 5.3< logT(K)<5.4, logρ(cm−3)> -4.0+∆ρz
        – 5.4< logT(K)<5.5, logρ(cm−3)> -4.2+∆ρz
        – 5.5< logT(K)<5.6, logρ(cm−3)> -4.4+∆ρz"""
        rhoz = np.log10(yt.YTQuantity(rhoz_func(z),'1/cm**3'))
        temps = np.log10(data['gas','temperature'])
        mask1 = temps < 5.2
        mask2 = np.logical_and(5.2 <= temps,temps< 5.3)
        mask3 = np.logical_and(5.3 <= temps,temps< 5.4)
        mask4 = np.logical_and(5.4 <= temps,temps< 5.5)
        mask5 = np.logical_and(5.5 <= temps,temps< 5.6)
        mask6 = temps > 5.6
        comp = np.zeros(temps.shape)
        comp[mask1] = -np.inf
        comp[mask2] = -3.8+rhoz
        comp[mask3] = -4.0+rhoz
        comp[mask4] = -4.2+rhoz
        comp[mask5] = -4.4+rhoz
        comp[mask6] = np.inf
        tr = yt.YTArray((np.log10(data['gas','density'])>comp).astype(float))
        return tr
    def CI_OVII(field,data):
        tr = data['gas','PI_OVII']
        return 1.-tr

    def PI_OVIII(field, data):
        #0 if CI, 1 if PI
        """– logT(K)<5.85, all ρ
        – 5.85< logT(K)<5.95, logρ(cm−3)> -4.0+∆ρz
        – 5.95< logT(K)<6.05, logρ(cm−3)> -4.25+∆ρz
        – 6.05< logT(K)<6.15, logρ(cm−3)> -4.5+∆ρz
        – 6.15< logT(K)<6.25, logρ(cm−3)> -4.75+∆ρz"""
        rhoz = np.log10(yt.YTQuantity(rhoz_func(z),'1/cm**3'))
        temps = np.log10(data['gas','temperature'])
        mask1 = temps < 5.85
        mask2 = np.logical_and(5.85 <= temps,temps< 5.95)
        mask3 = np.logical_and(5.95 <= temps,temps< 6.05)
        mask4 = np.logical_and(6.05 <= temps,temps< 6.15)
        mask5 = np.logical_and(6.15 <= temps,temps< 6.25)
        mask6 = temps > 6.25
        comp = np.zeros(temps.shape)
        comp[mask1] = -np.inf
        comp[mask2] = -4.0+rhoz
        comp[mask3] = -4.25+rhoz
        comp[mask4] = -4.5+rhoz
        comp[mask5] = -4.75+rhoz
        comp[mask6] = np.inf
        tr = yt.YTArray((np.log10(data['gas','density'])>comp).astype(float))
        return tr
    def CI_OVIII(field,data):
        tr = data['gas','PI_OVIII']
        return 1.-tr
    if ds and add_fields:
        ds.add_field(('gas','PI_OIV'),
               sampling_type="cell",
               function=PI_OIV,
               units='')
        ds.add_field(('gas','CI_OIV'),
               sampling_type="cell",
               function=CI_OIV,
               units='')
        ds.add_field(('gas','PI_OV'),
               sampling_type="cell",
               function=PI_OV,
               units='')
        ds.add_field(('gas','CI_OV'),
               sampling_type="cell",
               function=CI_OV,
               units='')
        ds.add_field(('gas','PI_OVI'),
               sampling_type="cell",
               function=PI_OVI,
               units='')
        ds.add_field(('gas','CI_OVI'),
               sampling_type="cell",
               function=CI_OVI,
               units='')
        ds.add_field(('gas','PI_OVII'),
               sampling_type="cell",
               function=PI_OVII,
               units='')
        ds.add_field(('gas','CI_OVII'),
               sampling_type="cell",
               function=CI_OVII,
               units='')
        ds.add_field(('gas','PI_OVIII'),
               sampling_type="cell",
               function=PI_OVIII,
               units='')
        ds.add_field(('gas','CI_OVIII'),
               sampling_type="cell",
               function=CI_OVIII,
               units='')
    
    return PI_OIV,CI_OIV,PI_OV,CI_OV,PI_OVI,CI_OVI,PI_OVII,CI_OVII,PI_OVIII,CI_OVIII
