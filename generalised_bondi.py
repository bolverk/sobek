import time

def mid_array(a):

    import numpy

    ar = numpy.array(a)
    
    return 0.5*(ar[1:]+ar[:-1])
    
def calc_gravity(cells, grid, extensives, dt):

    GM = 1
    r_list = mid_array(grid)
    for r, ext in zip(r_list, extensives):
        ext['energy'] -= dt*ext['momentum']*GM/r**2
        ext['momentum'] -= dt*ext['mass']*GM/r**2
        
def calc_spherical_complementary(cells, grid, extensives, dt):

    import numpy

    volume_list = numpy.diff([(4.0/3.0)*numpy.pi*r**3 for r in grid])
    r_list = mid_array(grid)
    for r, ext, cell, volume in zip(r_list, extensives, cells, volume_list):
        p = cell['pressure']
        ext['momentum'] += 2*volume*p*dt/r
        
def calc_mass_injection(cells, grid, extensives, dt):

    import numpy
    
    vw = 1.0
    r_list = mid_array(grid)
    volume_list = numpy.diff([(4.0/3.0)*numpy.pi*r**3 for r in grid])
    for r, volume, ext in zip(r_list, volume_list, extensives):
        mass_addition = r**(-2.5)
        ext['mass'] += mass_addition*volume
        ext['energy'] += 0.5*mass_addition*vw**2*volume
        
def calc_net_source(cells, grid, extensives, dt):

    calc_gravity(cells, grid, extensives, dt)
    calc_spherical_complementary(cells, grid, extensives, dt)
    calc_mass_injection(cells, grid, extensives, dt)

def custom_extensive_updater(grid, cells, extensive_list, flux_list, geometry, dt):

    from sobek.simple_extensive_updater import simple_extensive_updater
    import numpy

    res = simple_extensive_updater(grid, cells, extensive_list, flux_list, geometry, dt)

    # Gravity
    r_list = mid_array(grid)
    GM = 1e-3
    res['energy'] -= dt*extensive_list['momentum']*GM/r_list**2
    res['momentum'] -= dt*extensive_list['mass']*GM/r_list**2

    # Spherical complementary
    volume_list = numpy.diff((4.0/3.0)*numpy.pi*grid**3)
    res['momentum'] += 2*volume_list*cells['pressure']*dt/r_list

    # Mass source
    vw = 1.0
    res['mass'] += dt*r_list**(-2.5)*volume_list
    res['energy'] += dt*0.5*vw**2*r_list**(-2.5)*volume_list

    return res

def main():

    import numpy
    from sobek.simulation import Simulation
    from sobek.ideal_gas import IdealGas
    from sobek.simple_cfl import SimpleCFL
    from sobek.eulerian import Eulerian
    from sobek.flux_condition_action import FluxConditionAction
    from sobek.hllc import HLLC
    from sobek.vectorised_hllc import calc_vectorised_hllc
    from sobek.simple_cell_updater import simple_cell_updater
    import matplotlib.pyplot as plt
    plt.ion()
    
    data = {}
    data['grid'] = numpy.logspace(-1, 1, 30)
    #data['snapshot']['cells'] = [{'density':1e-9, 'pressure':1e-9, 'velocity': 0} for r in mid_array(data['snapshot']['grid'])]
    data['cells'] = numpy.array([(1e-9,1e-9,0.0) for r in mid_array(data['grid'])],
                                dtype=[('density','d'),('pressure','d'),('velocity','d')])
    data['equation_of_state'] = IdealGas(5./3.)
    data['physical_geometry'] = {'area': lambda r: 4.0*numpy.pi*r**2, 'volume': lambda r: 4.0*numpy.pi*r**3/3.0}
    data['time_step_function'] = SimpleCFL(0.3)
    data['grid_motion'] = Eulerian()
    data['riemann_solver'] = HLLC()

    def flux_calculator(grid, cells, velocity_list):

        leftmost = {'density':cells['density'][0]/10,
                    'pressure':cells['pressure'][0]/10,
                    'energy':cells['energy'][0],
                    'sound_speed':cells['sound_speed'][0],
                    'velocity':0}
        rightmost = {'density':cells['density'][-1]/10,
                     'pressure':cells['pressure'][-1]/10,
                     'energy':cells['energy'][-1],
                     'sound_speed':cells['sound_speed'][-1],
                     'velocity':0}
        left_states = numpy.array(zip(*(numpy.concatenate(([leftmost[field]],
                                                           cells[field]))
                                        for field in cells.dtype.names)),
                                  dtype = cells.dtype)
        right_states = numpy.array(zip(*(numpy.concatenate((cells[field],
                                                            [rightmost[field]]))
                                         for field in cells.dtype.names)),
                                   dtype=cells.dtype)
        return calc_vectorised_hllc(left_states, right_states, velocity_list)

    data['flux_calculator'] = flux_calculator
    data['extensive_updater'] = custom_extensive_updater
    data['cell_updater'] = simple_cell_updater
    #data['source_term'] = calc_net_source
    sim = Simulation(data)

    # Visualise
    r_list = mid_array(sim.data['grid'])
    fig = plt.figure()
    axes_list = [fig.add_subplot(3,1,n+1) for n in range(3)]
    plots = {}
    for n, field in enumerate(['density','sound_speed','velocity']):
        #pylab.subplot(3,1,n+1)
        y_list = map(lambda x:x[field], sim.data['cells'])
        plots[field], = axes_list[n].semilogx(r_list, map(lambda x:x[field], sim.data['cells']))
        axes_list[n].set_ylabel(field)

    def update_figure():
        for n, field in enumerate(['density','sound_speed','velocity']):
            y_list = sim.data['cells'][field]
            plots[field].set_ydata(y_list)
            if field == 'velocity':
                slack = 0.1*(numpy.max(y_list)-numpy.min(y_list))
                axes_list[n].set_ylim((numpy.min(y_list)-slack, numpy.max(y_list)+slack))
            else:
                axes_list[n].set_ylim((numpy.min(y_list), numpy.max(y_list)))
        plt.suptitle('t = '+str(sim.data['time'])+', c = '+str(sim.data['cycle']))
        plt.pause(0.01)

    start = time.time()
    for i in range(10000):
        sim.timeAdvance()
        if i%100==0:
            update_figure()
    print time.time() - start
    
def profile_main():

    import cProfile
    
    cProfile.run('main()','restats')
    
if __name__ == '__main__':

    main()
