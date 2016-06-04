import numpy

def calc_total_energy_density(p):

    return p['density']*(0.5*p['velocity']**2+p['energy'])

def primitive2conserved(primitive):

    return numpy.array([(primitive['density'],
                         primitive['density']*primitive['velocity'],
                         calc_total_energy_density(primitive))],
                         dtype=[('mass','f4'),('momentum','f4'),('energy','f4')])
            
def primitive2flux(primitive):

    return numpy.array([primitive['density']*primitive['velocity'],
                        primitive['pressure']+primitive['density']*primitive['velocity']**2,
                        (calc_total_energy_density(primitive)+primitive['pressure'])*primitive['velocity']],
                        dtype=[('mass','f4'),('momentum','f4'),('energy','f4')])
            
def conserved2primitive(c, eos):

    res = {}
    res['density'] = c['mass']
    res['velocity'] = c['momentum']/c['mass']
    res['energy'] = c['energy']/c['mass'] - 0.5*res['velocity']**2
    res['pressure'] = eos.de2p(res['density'], res['energy'])
    res['sound_speed'] = eos.dp2c(res['density'], res['pressure'])
    return res    