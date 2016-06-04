class SimpleCFL:

    def __init__(self, cfl):
    
        self.cfl = cfl
        
    def __call__(self, grid, cells):
    
        import numpy
    
        cell_widths = numpy.diff(grid)
        inverse_time_steps = (cells['sound_speed']+numpy.absolute(cells['velocity']))/cell_widths
        return self.cfl/numpy.max(inverse_time_steps)