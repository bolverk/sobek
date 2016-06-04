def main():

    import numpy
    
    density = numpy.linspace(0,1,100)
    velocity = numpy.zeros_like(density)
    
    mass_list = density
    momentum_list = density*velocity
    #temp = [(mass,momentum) for mass,momentum ]
    #temp = numpy.array(temp,dtype=[('mass','d'),('momentum','d')])

    print numpy.array(zip(mass_list, momentum_list),dtype=[('mass','d'),('momentum','d')])
    
if __name__ == '__main__':

    main()