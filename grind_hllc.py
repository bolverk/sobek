def main():

    from sobek.hllc import HLLC
    from sobek.ideal_gas import IdealGas
    import numpy
    import random
    
    eos = IdealGas(5./3.)
    rs = HLLC()
    
    left = {'density':10**random.randrange(-2,2),
            'pressure':10**random.randrange(-2,2),
            'velocity':10**random.randrange(-2,2)}
    right = {'density':10**random.randrange(-2,2),
             'pressure':10**random.randrange(-2,2),
             'velocity':10**random.randrange(-2,2)}
    for p in [left, right]:
        p['energy'] = eos.dp2e(p['density'], p['pressure'], p)
        p['sound_speed'] = eos.dp2c(p['density'], p['pressure'], p)
    velocity = 10**random.randrange(-2,2)
    
    for i in range(10000):
        rs(left, right, velocity)
        
if __name__ == '__main__':

    main()