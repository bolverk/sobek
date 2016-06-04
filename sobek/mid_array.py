def mid_array(a):

    import numpy

    ar = numpy.array(a)
    
    return 0.5*(ar[1:]+ar[:-1])