class FluxConditionAction:

    def __init__(self, ca_list):
    
        self.ca_list = ca_list
        
    def calcSingleFlux(self, snapshot, velocity_list, i):
    
        for case in self.ca_list:
            if case['condition'](snapshot, velocity_list, i):
                return case['action'](snapshot, velocity_list, i)
        
    def __call__(self, snapshot, velocity_list):
        
        return [self.calcSingleFlux(snapshot, velocity_list, i) for i in range(len(velocity_list))]