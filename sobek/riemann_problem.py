def main():

    from simple_extensive_updater import simple_extensive_updater
    from simple_cell_updater import simple_cell_updater

    import matplotlib.pyplot as plt
    plt.ion()
    import numpy
    from simple_cfl import SimpleCFL
    from ideal_gas import IdealGas
    from flux_condition_action import FluxConditionAction
    from hllc import HLLC
    from eulerian import Eulerian
    from vectorised_hllc import calc_vectorised_hllc
    from simulation import Simulation
    from mid_array import mid_array
    from simple_extensive_updater import simple_extensive_updater
    from physical_geometry import planar_geometry
    
    data = {}
    data['physical_geometry'] = planar_geometry
    data['time_step_function'] = SimpleCFL(0.3)
    data['grid'] = numpy.linspace(0, 1, 100)
    r_list = mid_array(data['grid'])
    data['cells'] = numpy.array(zip(1*numpy.ones(len(data['grid'])-1),
                                    numpy.where(r_list<0.5,2,1),
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

    main()