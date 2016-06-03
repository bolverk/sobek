import copy
import unittest
from hydrodynamic_variables import primitive2conserved, calc_total_energy_density, primitive2flux

def calc_wave_speeds(left, right):

    dl = left['density']
    pl = left['pressure']
    vl = left['velocity']
    cl = left['sound_speed']
    dr = right['density']
    pr = right['pressure']
    vr = right['velocity']
    cr = right['sound_speed']
    
    sl = min(vl-cl,vr-cr)
    sr = max(vl+cl,vr+cr)
    ss = (pr-pl+dl*vl*(sl-vl)-dr*vr*(sr-vr))/(dl*(sl-vl)-dr*(sr-vr))
    return {'left':sl,
            'right':sr,
            'center':ss}
            
def calc_starred_state(state, sk, ss):

    dk = state['density']
    pk = state['pressure']
    vk = state['velocity']
    ds = dk*(sk-vk)/(sk-ss)
    ek = calc_total_energy_density(state)
    return {'mass':ds,
            'momentum':ds*ss,
            'energy':ek*ds/dk+ds*(ss-vk)*(ss+pk/dk/(sk-vk))}

class HLLC:

    def __init__(self):
    
        pass
        
    def __call__(self, left, right, velocity):
        
        local_left = copy.deepcopy(left)
        local_right = copy.deepcopy(right)
        
        local_left['velocity'] -= velocity
        local_right['velocity'] -= velocity
        
        ul = primitive2conserved(local_left)
        ur = primitive2conserved(local_right)
        
        ws = calc_wave_speeds(local_left, local_right)
        
        usl = calc_starred_state(local_left, ws['left'], ws['center'])
        usr = calc_starred_state(local_right, ws['right'], ws['center'])
        
        fl = primitive2flux(local_left)
        fr = primitive2flux(local_right)
        
        if ws['left']>0:
            res = fl
        elif ws['left']<=0 and ws['center']>=0:
            #res = fl+ws['center']*(usl-ul)
            #res = {field: fl[field]+ws['center']*(usl[field]-ul[field] for field in fl}
            res = {}
            for field in fl:
                res[field] = fl[field] + ws['left']*(usl[field]-ul[field])
        elif ws['center']<0 and ws['right']>=0:
            #res = fr + ws['right']*(usr-ur)
            res = {}
            for field in fl:
                res[field] = fr[field] + ws['right']*(usr[field]-ur[field])
        else:
            res = fr
            
        res['energy'] += res['momentum']*velocity + 0.5*res['mass']*velocity**2
        res['momentum'] += velocity*res['mass']
        return res
        
        
class TestHLLC(unittest.TestCase):

    def test_static_uniform_wavespeeds(self):
    
        from ideal_gas import IdealGas
    
        eos = IdealGas(5./3.)        
        d = 1.0
        primitive = {'density':d,'pressure':d/eos.g,'velocity':0}
        primitive['energy'] = eos.dp2e(primitive['density'], primitive['pressure'], None)
        primitive['sound_speed'] = eos.dp2c(primitive['density'], primitive['pressure'], None)
        temp = calc_wave_speeds(primitive, primitive)
        self.assertEqual(temp['center'],0)
        self.assertEqual(temp['left'],-1)
        self.assertEqual(temp['right'],1)
        
    def test_calc_total_energy_density(self):
    
        p = {'density':1.0,'velocity':2.0,'energy':3.0}
        e = p['density']*(0.5*p['velocity']**2+p['energy'])
        self.assertEqual(calc_total_energy_density(p),e)
        
    def test_hllc_static_uniform(self):
    
        from ideal_gas import IdealGas
    
        eos = IdealGas(5./3.)        
        d = 1.0
        primitive = {'density':d,'pressure':d/eos.g,'velocity':0}
        primitive['energy'] = eos.dp2e(primitive['density'], primitive['pressure'], None)
        primitive['sound_speed'] = eos.dp2c(primitive['density'], primitive['pressure'], None)
        rs = HLLC()
        temp = rs(primitive, primitive, 0)
        self.assertEqual(temp['mass'],0)
        
if __name__ == '__main__':

    unittest.main()