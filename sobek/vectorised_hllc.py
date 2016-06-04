import numpy

def primitives2conserveds(primitives):

    densities = primitives['density']
    velocities = primitives['velocity']
    energies = primitives['energy']
    return numpy.array(zip(densities,
                           densities*velocities,
                           densities*(velocities**2+energies)),
                       dtype=[('mass','d'),('momentum','d'),('energy','d')])
                       
def calc_center_wave_speeds(left, right, sl, sr):

    dl = left['density']
    pl = left['pressure']
    vl = left['velocity']
    cl = left['sound_speed']
    dr = right['density']
    pr = right['pressure']
    vr = right['velocity']
    cr = right['sound_speed']
    
    return (pr-pl+dl*vl*(sl-vl)-dr*vr*(sr-vr))/(dl*(sl-vl)-dr*(sr-vr))
    
def calc_starred_states(states, sk, ss):

    dk = states['density']
    pk = states['pressure']
    vk = states['velocity']
    ds = dk*(sk-vk)/(sk-ss)
    ek = dk*(states['energy']+0.5*vk**2)
    return numpy.array(zip(ds,ds*ss,ek*ds/dk+ds*(ss-vk)*(ss+pk/dk/(sk-vk))),
                        dtype=[('mass','d'),('momentum','d'),('energy','d')])
                        
def primitives2fluxes(primitives):

    densities = primitives['density']
    velocities = primitives['velocity']
    pressures = primitives['pressure']
    energies = primitives['energy']
    return numpy.array(zip(densities*velocities,
                           densities*velocities**2+pressures,
                           densities*velocities*(velocities**2+energies+pressures/densities)),
                       dtype=[('mass','d'),('momentum','d'),('energy','d')])

def calc_vectorised_hllc(left_states, right_states, velocities):

    local_left_states = left_states.copy()
    local_left_states['velocity'] -= velocities
    local_right_states = right_states.copy()
    local_right_states['velocity'] -= velocities
    
    ul_list = primitives2conserveds(local_left_states)
    ur_list = primitives2conserveds(local_right_states)
    
    lws_list = numpy.minimum(local_left_states['velocity']-local_left_states['sound_speed'],
                             local_right_states['velocity']-local_right_states['sound_speed'])
    rws_list = numpy.maximum(local_left_states['velocity']+local_left_states['sound_speed'],
                             local_right_states['velocity']+local_right_states['sound_speed'])
    cws_list = calc_center_wave_speeds(local_left_states, local_right_states, lws_list, rws_list)
    
    usl_list = calc_starred_states(local_left_states, lws_list, cws_list)
    usr_list = calc_starred_states(local_right_states, rws_list, cws_list)
    
    fl_list = primitives2fluxes(local_left_states)
    fr_list = primitives2fluxes(local_right_states)
    
    ql_list = numpy.array(zip(*(fl_list[field]+lws_list*(usl_list[field]-ul_list[field]) for field in fl_list.dtype.names)),
                          dtype=fl_list.dtype)
    qr_list = numpy.array(zip(*(fr_list[field]+rws_list*(usr_list[field]-ur_list[field]) for field in fr_list.dtype.names)),
                          dtype=fr_list.dtype)
    res = numpy.where(lws_list>0, fl_list,
                       numpy.where(cws_list>=0, ql_list,
                                   numpy.where(cws_list<0,qr_list,fr_list)))
    res['energy'] += res['momentum']*velocities+0.5*res['mass']*velocities**2
    res['momentum'] += velocities*res['mass']
    return res
    
def test():

    n = 10
    hydro_variables = ['density','pressure','velocity','energy','sound_speed']
    dtype = [(field,'d') for field in hydro_variables]
    left_states = numpy.array(zip(*(10**(4*(numpy.random.rand(n)-0.5)) for c in range(5))),
                              dtype=dtype)
    right_states = numpy.array(zip(*(10**(4*(numpy.random.rand(n)-0.5)) for c in range(5))),
                               dtype=dtype)
    velocities = 10**(4*(numpy.random.rand(n)-0.5))
    
    print calc_vectorised_hllc(left_states, right_states, velocities)
    
def static_test():

    n = 10
    hydro_variables = ['density','pressure','velocity','energy','sound_speed']
    dtype = [(field,'d') for field in hydro_variables]
    left_states = numpy.array(zip(*(10**(4*(numpy.random.rand(n)-0.5)) for c in range(5))),
                              dtype=dtype)
    left_states['velocity'] = numpy.zeros(n)
    right_states = left_states
    velocities = numpy.zeros(n)
    
    res_1 = calc_vectorised_hllc(left_states, right_states, velocities)
    
    print res_1
    
def advection_test():

    n = 10
    hydro_variables = ['density','pressure','velocity','energy','sound_speed']
    dtype = [(field,'d') for field in hydro_variables]
    left_states = numpy.array(zip(*(10**(4*(numpy.random.rand(n)-0.5)) for c in range(5))),
                              dtype=dtype)
    left_states['velocity'] = numpy.ones(n)
    left_states['density'] = numpy.array(range(1,n+1))
    right_states = left_states
    velocities = numpy.zeros(n)
    
    res_1 = calc_vectorised_hllc(left_states, right_states, velocities)
    
    print res_1
    
def symmetric_collision_test():

    n = 1
    
    hydro_variables = ['density','pressure','velocity','energy','sound_speed']
    dtype = [(field,'d') for field in hydro_variables]
    left_states = numpy.array(zip(*(10**(4*(numpy.random.rand(n)-0.5)) for c in range(5))),
                              dtype=dtype)
    right_states = left_states.copy()
    right_states['velocity'] = -left_states['velocity']
    velocities = numpy.zeros(n)
    
    res_1 = calc_vectorised_hllc(left_states, right_states, velocities)
    
    print 'final flux',res_1
    
def compare_to_old_hllc():

    from hllc import HLLC

    n = 3
    hydro_variables = ['density','pressure','velocity','energy','sound_speed']
    dtype = [(field,'d') for field in hydro_variables]
    left_states = numpy.array(zip(*(10**(4*(numpy.random.rand(n)-0.5)) for c in range(5))),
                              dtype=dtype)
    right_states = numpy.array(zip(*(10**(4*(numpy.random.rand(n)-0.5)) for c in range(5))),
                               dtype=dtype)
    velocities = 10**(4*(numpy.random.rand(n)-0.5))
    
    rs = HLLC()
    
    res_1 = calc_vectorised_hllc(left_states, right_states, velocities)
    res_2 = [rs(left, right, vel) for left, right, vel in zip(left_states, right_states, velocities)]
    
    print res_1
    print
    print res_2
    
if __name__ == '__main__':

    symmetric_collision_test()
    