import unittest

class IdealGas:

    def __init__(self, g):
    
        self.g = g
        
    def de2p(self, d, e, extra):
    
        return (self.g-1)*d*e
        
    def dp2e(self, d, p, extra):
    
        return p/d/(self.g-1)
        
    def dp2c(self, d, p, extra):
    
        import math
        
        if p<0 or d<0:
            import pdb
            pdb.set_trace()
    
        return math.sqrt(self.g*p/d)
        

class TestIdealGas(unittest.TestCase):

    def test_de2p(self):
    
        g = 5./3.
        d = 5.0
        e = 7.0
        eos = IdealGas(g)
        self.assertEqual(eos.de2p(d,e,None),(g-1)*d*e)
        
    def test_dp2e(self):
    
        g = 5./3.
        d = 4.0
        p = 11.0
        eos = IdealGas(g)
        self.assertEqual(eos.dp2e(d,p,None),p/d/(g-1))
        
    def test_dp2c(self):
    
        import math
    
        g = 5./3.
        d = 9.0
        p = 2.5
        eos = IdealGas(g)
        self.assertEqual(eos.dp2c(d,p,None),math.sqrt(g*p/d))
        
if __name__ == '__main__':

    unittest.main()