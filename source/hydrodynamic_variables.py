def calc_total_energy_density(p):

    return p['density']*(0.5*p['velocity']**2+p['energy'])

def primitive2conserved(primitive):

    return {'mass':primitive['density'],
            'momentum':primitive['density']*primitive['velocity'],
            'energy':calc_total_energy_density(primitive)}
            
def primitive2flux(primitive):

    return {'mass':primitive['density']*primitive['velocity'],
            'momentum':primitive['pressure']+primitive['density']*primitive['velocity']**2,
            'energy':(calc_total_energy_density(primitive)+primitive['pressure'])*primitive['velocity']}
            
def conserved2primitive(c, eos):

    res = {}
    res['density'] = c['mass']
    res['velocity'] = c['momentum']/c['mass']
    res['energy'] = c['energy']/c['mass'] - 0.5*res['velocity']**2
    res['pressure'] = eos.de2p(res['density'], res['energy'], res)
    res['sound_speed'] = eos.dp2c(res['density'], res['pressure'], res)
    return res    