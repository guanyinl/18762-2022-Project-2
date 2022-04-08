import numpy as np


def process_results(v, bus):
    avg_voltage = 0
    avg_angle = 0
    bus_number =0
    for element in bus:
        bus_number = bus_number+1
        i = element.Bus-1
        voltage = np.sqrt(v[i*2]+v[i*2+1])
        angle = 180/np.pi*np.arctan(v[i*2+1]/v[i*2])
        print ('Bus:', element.Bus, ' voltage =', voltage, ' angle =', angle )
        avg_voltage = avg_voltage + voltage
        avg_angle = avg_angle + angle
    
    avg_voltage = avg_voltage/bus_number 
    avg_angle = avg_angle/bus_number

    print ('Avg_voltage =', avg_voltage, ' Avg_angle =', avg_angle )
    
    
    pass