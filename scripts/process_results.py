import numpy as np


def process_results(v, bus):
    for element in bus:
        i = element.Bus-1
        print ('Bus:', element.Bus, ' voltage =', np.sqrt(v[i*2]+v[i*2+1]), ' angle =', 180/np.pi*np.arctan(v[i*2+1]/v[i*2]) )

    
    pass