class FluxConditionAction:

    def __init__(self, ca_list):
    
        self.ca_list = ca_list
        
    def calcSingleFlux(self, grid, cells, velocity_list, i):
    
        for case in self.ca_list:
            if case['condition'](grid, cells, velocity_list, i):
                return case['action'](grid, cells, velocity_list, i)
        
    def __call__(self, grid, cells, velocity_list):
    
        import numpy
    
        def wrapper(i):
            return self.calcSingleFlux(grid,cells,velocity_list, i)
        #return [self.calcSingleFlux(grid, cells, velocity_list, i) for i in range(len(velocity_list))]
        return numpy.vectorize(wrapper)(range(len(velocity_list)))