def simple_extensive_updater(grid, cells, extensive_list, flux_list, geometry, dt):

    import numpy

    area_list = numpy.vectorize(geometry['area'])(grid)
    current_list = numpy.array(zip(*(flux_list[field]*area_list*dt for field in flux_list.dtype.names)),
                            dtype=flux_list.dtype)
    diff_list = numpy.array(zip(*(numpy.diff(current_list[field]) for field in current_list.dtype.names)),
                            dtype=current_list.dtype)
    return numpy.array(zip(*(extensive_list[field] - diff_list[field] for field in extensive_list.dtype.names)),
                       dtype=extensive_list.dtype)
