class Eulerian:

    def __init__(self):
    
        pass
        
    def __call__(self, snapshot):
    
        import numpy
    
        return numpy.zeros_like(snapshot['grid'])