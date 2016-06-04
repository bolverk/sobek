import math

def planar_area(r):

    return 1
    
def planar_volume(r):

    return r
    
def cylindrical_area(r):

    return 2.0*math.pi*r
    
def cylindrical_volume(r):

    return math.pi*r**2
    
def spherical_area(r):

    return 4.0*math.pi*r**2
    
def spherical_volume(r):

    return 4.0*math.pi*r**3/3.0
    
planar_geometry = {'area':planar_area,'volume':planar_volume}
cylindrical_geometry = {'area':cylindrical_area,'volume':cylindrical_volume}
spherical_geometry = {'area':spherical_area,'volume':spherical_volume}