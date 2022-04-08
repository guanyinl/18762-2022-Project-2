
from pickle import NONE
import numpy as numpy


def initialize(parsed_data, case_name, size_Y):
    if case_name != None:
        #The v vector is defined as vr1, vi1, vr2, vi2, vr3, vi3, vr4, vi4, qg, ir1, ii1
        #The definition is based on how the circuits of GS-4 look like.

        for slack in parsed_data['slack']:
            slack.Pinit = slack.Pinit/100
            slack.Qinit = slack.Qinit/100

        for gen in parsed_data['generators']:
            gen.P = gen.P/100
            gen.Qinit = gen.Qinit/100

        for load in parsed_data['loads']:
            load.P = load.P/100 
            load.Q = load.Q/100

        for shunt in parsed_data['shunts']:
            shunt.G_MW = shunt.G_MW/100 
            shunt.B_MVAR = shunt.B_MVAR/100

        bus_number = 0
        for bus in parsed_data['buses']:
            bus_number = bus_number+1
        
        gen_number = 0
        for gen in parsed_data['generators']:
            gen_number = gen_number+1


        v_init = [0]*size_Y
        #print (size_Y)
        i = 0 
        #Get initial values from the slack bus and put them into the V vector
        #for slack in parsed_data['slack']:
        #    for bus in parsed_data['buses']:
        #        if (slack.Bus == bus.Bus):
        #            v_init[(bus.Bus-1)*2] = slack.Vset*numpy.cos(slack.ang*numpy.pi/180)      
        #            v_init[(bus.Bus-1)*2+1] = slack.Vset*numpy.sin(slack.ang*numpy.pi/180)     
                 
            
        #Get initial values from loads and put them into the V vector
        #for load in parsed_data['loads']:
        #    for bus in parsed_data['buses']:
        #        if (load.Bus == bus.Bus) and (bus.Type == 1):
        #            v_init[(bus.Bus-1)*2] = bus.Vm_init*numpy.cos(bus.Va_init*numpy.pi/180) #vr2, 3 = pq on bus 2, 3. vm =1, va =0;
        #            v_init[(bus.Bus-1)*2+1] = bus.Vm_init*numpy.sin(bus.Va_init*numpy.pi/180) #vi2, 3 = pq on bus 2, 3. vm =1, va =0;
        #            #print ((bus.Bus-1)*2, v_init[(bus.Bus-1)*2],  (bus.Bus-1)*2+1 ,v_init[(bus.Bus-1)*2+1])

        for bus in parsed_data['buses']:
            v_init[(bus.Bus-1)*2] = bus.Vm_init*numpy.cos(bus.Va_init*numpy.pi/180)
            v_init[(bus.Bus-1)*2+1] = bus.Vm_init*numpy.sin(bus.Va_init*numpy.pi/180)
            #print((bus.Bus-1)*2, (bus.Bus-1)*2+1)


        #Get initial values from PVs and put them into the V vector
        #for gen in parsed_data['generators']:
        #    for bus in parsed_data['buses']:
        #        if (gen.Bus == bus.Bus) and (bus.Type == 2):
        #            v_init[(bus.Bus-1)*2] = bus.Vm_init*numpy.cos(bus.Va_init*numpy.pi/180) #vr4 = pv on bus 4. vm =1, va =0;
        #            v_init[(bus.Bus-1)*2+1] = bus.Vm_init*numpy.sin(bus.Va_init*numpy.pi/180) #vi4 = pv on bus 4. vm =1, va =0;
                    #print ((bus.Bus-1)*2, v_init[(bus.Bus-1)*2],  (bus.Bus-1)*2+1 ,v_init[(bus.Bus-1)*2+1])

        #Set initial values of Qg into the V vector
        gen_index = 0
        for gen in parsed_data['generators']:
            v_init[bus_number*2+gen_index] = gen.Qinit #qg = pv's reactive power.  
            #print (bus_number*2+gen_index, gen.Qinit)
            gen_index = gen_index+1
            

        #set the initial values of current of the slack bus into the V vector
        slack_index = 0
        for slack in parsed_data['slack']:
            #slack.Pinit = 50/100?
            #slack.Qinit = 30.99/100?
            v_init[bus_number*2+gen_number+slack_index] = -1*slack.Pinit/slack.Vset #ir1 = the real current flowing into the sb. -1 means is actually flowing out.
            v_init[bus_number*2+gen_number+slack_index+1] = -1*slack.Qinit/slack.Vset #ii1 = the imaginary current flowing into the sb. -1 means is actually flowing out.
            slack_index = slack_index+1


        print ('v_init =', v_init, '\n')

    return (v_init)