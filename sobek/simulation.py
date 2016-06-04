class Simulation:

    def __init__(self, 
                 data):
        import numpy
        from hydrodynamic_variables import primitive2conserved
    
        self.data = data
        eos = data['equation_of_state']
        cells = data['cells']
        self.data['cells'] = numpy.array(zip(cells['density'],
                                             cells['pressure'],
                                             cells['velocity'],
                                             numpy.vectorize(eos.dp2e)(cells['density'],
                                                                       cells['pressure']),
                                             numpy.vectorize(eos.dp2c)(cells['density'],
                                                                       cells['pressure'])),
                                         dtype=[('density','d'),
                                                ('pressure','d'),
                                                ('velocity','d'),
                                                ('energy','d'),
                                                ('sound_speed','d')])
        volume_list = numpy.diff(numpy.vectorize(data['physical_geometry']['volume'])(data['grid']))
        
        cells = self.data['cells']
        intensive_list = numpy.array(zip(cells['density'], 
                                         cells['density']*cells['velocity'],
                                         0.5*cells['density']*(cells['velocity']**2+cells['energy'])),
                                     dtype=[('mass','d'),('momentum','d'),('energy','d')])
        self.data['extensive'] = numpy.array(zip(intensive_list['mass']*volume_list,
                                                 intensive_list['momentum']*volume_list,
                                                 intensive_list['energy']*volume_list),
                                             dtype=[('mass','d'),('momentum','d'),('energy','d')])
        self.data['time'] = 0
        self.data['cycle'] = 0
        
    def timeAdvance(self):
    
        dt = self.data['time_step_function'](self.data['grid'], self.data['cells'])
        
        grid_velocity = self.data['grid_motion'](self.data['grid'], self.data['cells'])
        
        fluxes = self.data['flux_calculator'](self.data['grid'], self.data['cells'], grid_velocity)
                
        self.data['extensive'] = self.data['extensive_updater'](self.data['grid'], 
                                                                self.data['extensive'],
                                                                fluxes,
                                                                self.data['physical_geometry'],
                                                                dt)
                                                          
        self.data['cells'] = self.data['cell_updater'](self.data['grid'],
                                                       self.data['extensive'],
                                                       self.data['equation_of_state'],
                                                       self.data['physical_geometry'],
                                                       self.data['cells'])
                                                                   
        self.data['time'] += dt
        self.data['cycle'] += 1
    
def flip_velocity(p):

    return {field:-p[field] if field=='velocity' else p[field]}
    
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
        
def test():

    import matplotlib.pyplot as plt
    plt.ion()
    import numpy
    from simple_cfl import SimpleCFL
    from ideal_gas import IdealGas
    from flux_condition_action import FluxConditionAction
    from hllc import HLLC
    from eulerian import Eulerian
    from vectorised_hllc import calc_vectorised_hllc
    from simple_extensive_updater import simple_extensive_updater
    from simple_cell_updater import simple_cell_updater
    from mid_array import mid_array
    
    data = {}
    #data['physical_geometry'] = {'area': lambda r: 1, 'volume': lambda r: r}
    data['physical_geometry'] = {'area': lambda r: 4.0*numpy.pi*r**2, 'volume': lambda r: 4.0*numpy.pi*r**3/3.0}
    data['time_step_function'] = SimpleCFL(0.3)
    data['grid'] = numpy.logspace(-1, 1, 200)
    data['cells'] = numpy.array(zip(1e-9*numpy.ones(len(data['grid'])-1),
                                    1e-9*numpy.ones(len(data['grid'])-1),
                                    0*numpy.ones(len(data['grid'])-1)),
                                 dtype=[('density','d'),('pressure','d'),('velocity','d')])
    data['equation_of_state'] = IdealGas(5./3.)
    data['grid_motion'] = Eulerian()

    def flux_calculator(grid, cells, velocity_list):
    
        left_states = numpy.array(zip(*(numpy.concatenate(([cells[field][0]],cells[field])) for field in cells.dtype.names)),
                                  dtype=cells.dtype)
        right_states = numpy.array(zip(*(numpy.concatenate((cells[field],[cells[field][-1]])) for field in cells.dtype.names)),
                                   dtype=cells.dtype)
        return calc_vectorised_hllc(left_states, right_states, velocity_list)
    data['flux_calculator'] = flux_calculator
    data['extensive_updater'] = simple_extensive_updater
    data['cell_updater'] = simple_cell_updater
    data['source_term'] = calc_net_source
    sim = Simulation(data)
    
    # Visualise
    r_list = mid_array(sim.data['grid'])
    fig = plt.figure()
    axes_list = [fig.add_subplot(3,1,n+1) for n in range(3)]
    plots = {}
    for n, field in enumerate(['density','pressure','velocity']):
        #pylab.subplot(3,1,n+1)
        y_list = map(lambda x:x[field], sim.data['cells'])
        plots[field], = axes_list[n].plot(r_list, map(lambda x:x[field], sim.data['cells']))
        axes_list[n].set_ylabel(field)
    
    def update_figure():
        for n, field in enumerate(['density','pressure','velocity']):
            y_list = sim.data['cells'][field]
            plots[field].set_ydata(y_list)
            if field == 'velocity':
                slack = 0.1*(numpy.max(y_list)-numpy.min(y_list))
                axes_list[n].set_ylim((numpy.min(y_list)-slack, numpy.max(y_list)+slack))
            else:
                axes_list[n].set_ylim((numpy.min(y_list), numpy.max(y_list)))
        plt.suptitle('t = '+str(sim.data['time'])+', c = '+str(sim.data['cycle']))
        plt.pause(0.01)
    
    while sim.data['time']<0.1:
        sim.timeAdvance()
        update_figure()
        
if __name__ == '__main__':

    test()