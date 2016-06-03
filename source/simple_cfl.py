class SimpleCFL:

    def __init__(self, cfl):
    
        self.cfl = cfl
        
    def __call__(self, snapshot):
    
        import numpy
    
        res = 0
        for dx, cell in zip(numpy.diff(snapshot['grid']),snapshot['cells']):
            idt = (cell['sound_speed']+abs(cell['velocity']))/dx
            res = max(res, idt)
        return self.cfl/res