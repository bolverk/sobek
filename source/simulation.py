class Simulation:

    def __init__(self, 
                 data):
        import numpy
        from hydrodynamic_variables import primitive2conserved
    
        self.data = data
        for cell in self.data['snapshot']['cells']:
            cell['sound_speed'] = self.data['equation_of_state'].dp2c(cell['density'], cell['pressure'], cell)
            cell['energy'] = self.data['equation_of_state'].dp2e(cell['density'], cell['pressure'], cell)
        volume_list = numpy.diff([self.data['physical_geometry']['volume'](r) for r in self.data['snapshot']['grid']])
        self.data['intensive'] = [primitive2conserved(cell) for cell in self.data['snapshot']['cells']]
        self.data['extensive'] = [{field:volume*intensive[field] for field in intensive} 
                                    for volume,intensive in zip(volume_list,self.data['intensive'])]
        self.data['time'] = 0
        self.data['cycle'] = 0
        
    def timeAdvance(self):
    
        dt = self.data['time_step_function'](self.data['snapshot'])
        
        grid_velocity = self.data['grid_motion'](self.data['snapshot'])
        
        fluxes = self.data['flux_calculator'](self.data['snapshot'], grid_velocity)
        
        self.data['extensive'] = self.data['extensive_updater'](self.data['snapshot']['grid'], 
                                                                self.data['extensive'],
                                                                fluxes,
                                                                self.data['physical_geometry'],
                                                                dt)
        self.data['source_term'](self.data['snapshot']['cells'], 
                                 self.data['snapshot']['grid'],
                                 self.data['extensive'], 
                                 dt)
                                                          
        self.data['snapshot']['cells'] = self.data['cell_updater'](self.data['snapshot']['grid'],
                                                                   self.data['extensive'],
                                                                   self.data['equation_of_state'],
                                                                   self.data['physical_geometry'])
                                                                   
        self.data['time'] += dt
        self.data['cycle'] += 1
        
def mid_array(a):

    import numpy

    ar = numpy.array(a)
    
    return 0.5*(ar[1:]+ar[:-1])
    
def flip_velocity(p):

    import copy

    res = copy.deepcopy(p)
    res['velocity'] = -res['velocity']
    return res
    
def simple_extensive_updater(grid, extensive_list, flux_list, geometry, dt):

    import copy

    res = copy.deepcopy(extensive_list)
    area_list = [geometry['area'](r) for r in grid]
    for i in range(len(flux_list)):
        for field in extensive_list[0]:
            diff = flux_list[i][field]*area_list[i]*dt
            if i>0:
                res[i-1][field] -= diff
            if i<len(flux_list)-1:
                res[i][field] += diff
    return res
    
def simple_cell_updater(grid, extensive_list, eos, pg):

    import numpy
    import copy
    from hydrodynamic_variables import conserved2primitive
    
    volume_list = numpy.diff([pg['volume'](r) for r in grid])
    intensive_list = copy.deepcopy(extensive_list)
    for intensive, volume in zip(intensive_list,volume_list):
        for field in intensive:
            intensive[field] /= volume
    return [conserved2primitive(intensive, eos) for intensive in intensive_list]
    
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
    
    data = {}
    #data['physical_geometry'] = {'area': lambda r: 1, 'volume': lambda r: r}
    data['physical_geometry'] = {'area': lambda r: 4.0*numpy.pi*r**2, 'volume': lambda r: 4.0*numpy.pi*r**3/3.0}
    data['time_step_function'] = SimpleCFL(0.3)
    data['snapshot'] = {'grid':numpy.logspace(-1, 1, 200)}
    data['snapshot']['cells'] = [{'density':1e-9, 'pressure':1e-9, 'velocity': 0} for r in mid_array(data['snapshot']['grid'])]
    data['equation_of_state'] = IdealGas(5./3.)
    data['grid_motion'] = Eulerian()
    data['riemann_solver'] = HLLC()
    data['flux_calculator'] = FluxConditionAction([{'condition': lambda snapshot, velocity_list, i: len(velocity_list)-1>i>0,
                                                    'action': lambda snapshot, velocity_list, i: data['riemann_solver'](snapshot['cells'][i-1], snapshot['cells'][i], velocity_list[i])},
                                                   {'condition': lambda snapshot, velocity_list, i: i==0,
                                                    'action': lambda snapshot, velocity_list, i: data['riemann_solver'](snapshot['cells'][0],snapshot['cells'][0], velocity_list[0])},
                                                   {'condition': lambda snapshot, velocity_list, i: i==len(velocity_list)-1,
                                                    'action': lambda snapshot, velocity_list, i: data['riemann_solver'](snapshot['cells'][i-1], snapshot['cells'][i-1], velocity_list[0])}])
    data['extensive_updater'] = simple_extensive_updater
    data['cell_updater'] = simple_cell_updater
    data['source_term'] = calc_net_source
    sim = Simulation(data)
    
    # Visualise
    r_list = mid_array(sim.data['snapshot']['grid'])
    fig = plt.figure()
    axes_list = [fig.add_subplot(3,1,n+1) for n in range(3)]
    plots = {}
    for n, field in enumerate(['density','pressure','velocity']):
        #pylab.subplot(3,1,n+1)
        y_list = map(lambda x:x[field], sim.data['snapshot']['cells'])
        if field == 'velocity':
            plots[field], = axes_list[n].semilogx(r_list, map(lambda x:x[field], sim.data['snapshot']['cells']))
        else:
            plots[field], = axes_list[n].loglog(r_list, map(lambda x:x[field], sim.data['snapshot']['cells']))
        axes_list[n].set_ylabel(field)
    
    def update_figure():
        for n, field in enumerate(['density','pressure','velocity']):
            y_list = map(lambda x:x[field], sim.data['snapshot']['cells'])
            plots[field].set_ydata(y_list)
            if field == 'velocity':
                slack = 0.1*(numpy.max(y_list)-numpy.min(y_list))
                axes_list[n].set_ylim((numpy.min(y_list)-slack, numpy.max(y_list)+slack))
            else:
                axes_list[n].set_ylim((numpy.min(y_list), numpy.max(y_list)))
        plt.suptitle('t = '+str(sim.data['time'])+', c = '+str(sim.data['cycle']))
        plt.pause(0.01)
        
    
    while sim.data['time']<1:
        sim.timeAdvance()
        update_figure()
            
        
    raw_input()
        
if __name__ == '__main__':

    test()