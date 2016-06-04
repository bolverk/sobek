def simple_cell_updater(grid, extensive_list, eos, pg, cells):

    import numpy
    from hydrodynamic_variables import conserved2primitive
    
    volume_list = numpy.diff([pg['volume'](r) for r in grid])
    intensive_list = numpy.array(zip(*(extensive_list[field]/volume_list for field in extensive_list.dtype.names)),
                                 dtype=extensive_list.dtype)
    density_list = intensive_list['mass']
    velocity_list = intensive_list['momentum']/intensive_list['mass']
    thermal_energy_list = intensive_list['energy']/density_list - 0.5*velocity_list**2
    pressure_list = numpy.vectorize(eos.de2p)(density_list, thermal_energy_list)
    sound_speed_list = numpy.vectorize(eos.dp2c)(density_list, pressure_list)
    return numpy.array(zip(density_list,pressure_list, velocity_list, thermal_energy_list, sound_speed_list),
                       dtype=[('density','d'),('pressure','d'),('velocity','d'),('energy','d'),('sound_speed','d')])