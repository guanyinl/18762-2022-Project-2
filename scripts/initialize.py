
from pickle import NONE
import numpy as numpy


def initialize(parsed_data, case_name):
    if case_name == 'testcases/GS-4_prior_solution.RAW':
        #The v vector is defined as vr1, vi1, vr2, vi2, vr3, vi3, vr4, vi4, qg, ir1, ii1
        #The definition is based on how the circuits of GS-4 look like.
        v_init = [0]*11
        i = 0 
        #Get initial values from the slack bus and put them into the V vector
        for slack in parsed_data['slack']:
            v_init[i] = slack.Vset*numpy.cos(slack.ang)
            i = i+1       
            v_init[i] = slack.Vset*numpy.sin(slack.ang)
            i = i+1             
            
        #Get initial values from loads and put them into the V vector
        for load in parsed_data['loads']:
            for bus in parsed_data['buses']:
                if (load.Bus == bus.Bus) and (bus.Type == 1):
                    v_init[i] = bus.Vm_init*numpy.cos(bus.Va_init) #vr2, 3 = pq on bus 2, 3. vm =1, va =0;
                    i = i+1 
                    v_init[i] = bus.Vm_init*numpy.sin(bus.Va_init) #vi2, 3 = pq on bus 2, 3. vm =1, va =0;
                    i = i+1 

        #Get initial values from PVs and put them into the V vector
        for gen in parsed_data['generators']:
            for bus in parsed_data['buses']:
                if (gen.Bus == bus.Bus) and (bus.Type == 2):
                    v_init[i] = bus.Vm_init*numpy.cos(bus.Va_init) #vr4 = pv on bus 4. vm =1, va =0;
                    i = i+1 
                    v_init[i] = bus.Vm_init*numpy.sin(bus.Va_init) #vi4 = pv on bus 4. vm =1, va =0;
                    i = i+1

        #Set initial values of Qg into the V vector
        for gen in parsed_data['generators']:
            v_init[i] = gen.Qinit #qg = pv's reactive power. 
            i = i+1

        #set the initial values of current of the slack bus into the V vector
        for slack in parsed_data['slack']:
            #slack.Pinit = 50/100
            #slack.Qinit = 30.99/100
            v_init[i] = -1*slack.Pinit/slack.Vset #ir1 = the real current flowing into the sb. -1 means is actually flowing out.
            i = i+1
            v_init[i] = -1*slack.Qinit/slack.Vset #ii1 = the imaginary current flowing into the sb. -1 means is actually flowing out.

        print ('v_init =', v_init)

    return (v_init)