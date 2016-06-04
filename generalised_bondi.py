import time

def mid_array(a):

    import numpy

    ar = numpy.array(a)
    
    return 0.5*(ar[1:]+ar[:-1])
    
def simple_extensive_updater(grid, extensive_list, flux_list, geometry, dt):

    area_list = [geometry['area'](r) for r in grid]
    for i in range(len(flux_list)):
        for field in extensive_list[0]:
            diff = flux_list[i][field]*area_list[i]*dt
            if i>0:
                extensive_list[i-1][field] -= diff
            if i<len(flux_list)-1:
                extensive_list[i][field] += diff
    
def simple_cell_updater(grid, extensive_list, eos, pg, cells):

    import numpy
    from sobek.hydrodynamic_variables import conserved2primitive
    
    volume_list = numpy.diff([pg['volume'](r) for r in grid])
    intensive_list = extensive_list/volume
    
    return numpy.vectorize(conserved2primitive)(intensive_list, eos)
    
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

def main():

    import numpy
    from sobek.simulation import Simulation
    from sobek.ideal_gas import IdealGas
    from sobek.simple_cfl import SimpleCFL
    from sobek.eulerian import Eulerian
    from sobek.flux_condition_action import FluxConditionAction
    from sobek.hllc import HLLC
    
    data = {}
    data['snapshot'] = {'grid':numpy.logspace(-1, 1, 1000)}
    #data['snapshot']['cells'] = [{'density':1e-9, 'pressure':1e-9, 'velocity': 0} for r in mid_array(data['snapshot']['grid'])]
    data['snapshot']['cells'] = numpy.array([(1.0,1.0,0.0) for r in mid_array(data['snapshot']['grid'])],
                                            dtype=[('density','f4'),('pressure','f4'),('velocity','f4')])
    data['equation_of_state'] = IdealGas(5./3.)
    data['physical_geometry'] = {'area': lambda r: 4.0*numpy.pi*r**2, 'volume': lambda r: 4.0*numpy.pi*r**3/3.0}
    data['time_step_function'] = SimpleCFL(0.3)
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
    
    start = time.time()
    for i in range(100):
        sim.timeAdvance()
    print time.time() - start
    
def profile_main():

    import cProfile
    
    cProfile.run('main()','restats')
    
if __name__ == '__main__':

    main()