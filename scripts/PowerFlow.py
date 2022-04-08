from re import X
from sqlite3 import Row
import numpy as np
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import spsolve


class PowerFlow:

    def __init__(self,
                 case_name,
                 tol,
                 max_iters,
                 enable_limiting):
        """Initialize the PowerFlow instance.

        Args:
            case_name (str): A string with the path to the test case.
            tol (float): The chosen NR tolerance.
            max_iters (int): The maximum number of NR iterations.
            enable_limiting (bool): A flag that indicates if we use voltage limiting or not in our solver.
        """
        # Clean up the case name string
        case_name = case_name.replace('.RAW', '')
        case_name = case_name.replace('testcases/', '')

        self.case_name = case_name
        self.tol = tol
        self.max_iters = max_iters
        self.enable_limiting = enable_limiting
    

    def solver(self, matrix_Y, matrix_v, matrix_J, parsed_data):
        matrix_v = spsolve(matrix_Y, matrix_J)
        #print ('matrix_v = ', matrix_v)

        bus_number =  0
        slack_number =  0
        gen_number = 0
        for bus in parsed_data['buses']:
            bus_number = bus_number+1 #each bus contributes 2 variables in V (Vr, Vi)
        for slack in parsed_data['slack']:
            slack_number = slack_number+1 #each slack contributes 2 variables in V. (Ir, Ii)
        for gen in parsed_data['generators']:
            gen_number = gen_number+1 #each generators contributes 1 variable in V. (Qg)  

        #put results back to each class
        i = 0 
        for bus in parsed_data['buses']:
            i = i+1
        for gen in parsed_data['generators']:
            #only Qinit needs to be updated. 
            gen.Qinit = matrix_v[i*2]
        
        #Update for J
        #for slack bus
        i_slack = 0
        v_index = 0
        for slack in parsed_data['slack']:
            for bus in parsed_data['buses']:
                if (slack.Bus == bus.Bus):
                    angle = slack.ang*np.pi/180  
                    Vset_r = slack.Vset*np.cos(angle)
                    Vset_l = slack.Vset*np.sin(angle)
            
                    matrix_J[bus_number*2+gen_number*1+(i_slack)*2] = Vset_r
                    matrix_J[bus_number*2+gen_number*1+(i_slack)*2+1] = Vset_l
                    
                    i_slack = i_slack+1
        
        #slack bus also has P and Q connected to them, and P and Q need to be added into the matrix.
        for slack in parsed_data['slack']:
            for load in parsed_data['loads']:
                if (load.Bus == slack.Bus):
                    for bus in parsed_data['buses']:
                        if (load.Bus == bus.Bus): 
                            Vrl = matrix_v[v_index*2]
                            Vil = matrix_v[v_index*2+1]
                            Pl = load.P
                            Ql = load.Q  

                            #Irl part of PQ
                            dIrl_dVrl = (Pl*(np.square(Vrl)+np.square(Vil))-(Pl*Vrl+1*Ql*Vil)*2*Vrl)/np.square(np.square(Vrl)+np.square(Vil))
                            dIrl_dVil = (Ql*(np.square(Vrl)+np.square(Vil))-(Pl*Vrl+1*Ql*Vil)*2*Vil)/np.square(np.square(Vrl)+np.square(Vil))
                            Irl = (Pl*Vrl+Ql*Vil)/(np.square(Vrl)+np.square(Vil))

                            Vrl_hist = -1*(Irl-dIrl_dVrl*Vrl-dIrl_dVil*Vil)

                            #Iil part of PQ
                            dIil_dVrl = (-1*Ql*(np.square(Vrl)+np.square(Vil))-(Pl*Vil-1*Ql*Vrl)*2*Vrl)/np.square(np.square(Vrl)+np.square(Vil))
                            dIil_dVil = (Pl*(np.square(Vrl)+np.square(Vil))-(Pl*Vil-1*Ql*Vrl)*2*Vil)/np.square(np.square(Vrl)+np.square(Vil))
                            Iil = (Pl*Vil-Ql*Vrl)/(np.square(Vrl)+np.square(Vil))

                            Vil_hist = -1*(Iil-dIil_dVrl*Vrl-dIil_dVil*Vil)
                            matrix_J[(bus.Bus-1)*2] = Vrl_hist
                            matrix_J[(bus.Bus-1)*2+1] = Vil_hist
                            v_index = v_index+1

        #for PQ bus
        for load in parsed_data['loads']:
            for bus in parsed_data['buses']:
                if (load.Bus == bus.Bus) and (bus.Type == 1):
                    #abandoned the angle calculation because it creates "sign flips".
                    #angle = bus.Va_init*np.pi/180 #Va is not just init, it will be updated. 
                    #Vrl = bus.Vm_init*np.cos(angle) #Vm is not just init, it will be updated.
                    #Vil = bus.Vm_init*np.sin(angle)

                    #Used the last v matrx to update J
                    Vrl = matrix_v[v_index*2]
                    Vil = matrix_v[v_index*2+1]
                    Pl = load.P
                    Ql = load.Q  

                    #Irl part of PQ
                    dIrl_dVrl = (Pl*(np.square(Vrl)+np.square(Vil))-(Pl*Vrl+1*Ql*Vil)*2*Vrl)/np.square(np.square(Vrl)+np.square(Vil))
                    dIrl_dVil = (Ql*(np.square(Vrl)+np.square(Vil))-(Pl*Vrl+1*Ql*Vil)*2*Vil)/np.square(np.square(Vrl)+np.square(Vil))
                    Irl = (Pl*Vrl+Ql*Vil)/(np.square(Vrl)+np.square(Vil))

                    Vrl_hist = -1*(Irl-dIrl_dVrl*Vrl-dIrl_dVil*Vil)

                    #Iil part of PQ
                    dIil_dVrl = (-1*Ql*(np.square(Vrl)+np.square(Vil))-(Pl*Vil-1*Ql*Vrl)*2*Vrl)/np.square(np.square(Vrl)+np.square(Vil))
                    dIil_dVil = (Pl*(np.square(Vrl)+np.square(Vil))-(Pl*Vil-1*Ql*Vrl)*2*Vil)/np.square(np.square(Vrl)+np.square(Vil))
                    Iil = (Pl*Vil-Ql*Vrl)/(np.square(Vrl)+np.square(Vil))

                    Vil_hist = -1*(Iil-dIil_dVrl*Vrl-dIil_dVil*Vil)
                    
                    matrix_J[(bus.Bus-1)*2] = Vrl_hist
                    matrix_J[(bus.Bus-1)*2+1] = Vil_hist
                    #print (v_index*2,v_index*2+1, (bus.Bus-1)*2, (bus.Bus-1)*2+1)
                    v_index = v_index+1
        gen_number = 0
        for gen in parsed_data['generators']:
            for bus in parsed_data['buses']:
                if bus.Bus == gen.Bus:
                    #abandoned the angle calculation because it creates "sign flips".
                    #angle = bus.Va_init*np.pi/180 #it's not just init, it will be updated. 
                    #Vrg = bus.Vm_init*np.cos(angle)
                    #Vig = bus.Vm_init*np.sin(angle)

                    #Used the last v matrx to update J
                    Vrg = matrix_v[v_index*2]
                    Vig = matrix_v[v_index*2+1]
                    Pg = gen.P
                    Qg = gen.Qinit #it's not just init, it will be updated. 
                    VSet = gen.Vset
                    #print (v_index*2,v_index*2+1)
                    #gen_angle = np.arctan(Qg/Pg)
                    
                    #Irg part of PV
                    dIrg_dVrg = (-1*Pg*(np.square(Vrg)+np.square(Vig))-(-1*Pg*Vrg-1*Qg*Vig)*2*Vrg)/np.square(np.square(Vrg)+np.square(Vig))
                    dIrg_dVig = (-1*Qg*(np.square(Vrg)+np.square(Vig))-(-1*Pg*Vrg-1*Qg*Vig)*2*Vig)/np.square(np.square(Vrg)+np.square(Vig))
                    dIrg_dQg = -1*Vig/(np.square(Vrg)+np.square(Vig))
                    Irg = (-1*Pg*Vrg-Qg*Vig)/(np.square(Vrg)+np.square(Vig))

                    Vrg_hist = -1*(Irg-dIrg_dQg*Qg-dIrg_dVrg*Vrg-dIrg_dVig*Vig)

                    #Iirg part of PV
                    dIig_dVig = (-1*Pg*(np.square(Vrg)+np.square(Vig))-(-1*Pg*Vig+Qg*Vrg)*2*Vig)/np.square(np.square(Vrg)+np.square(Vig))
                    dIig_dVrg = (1*Qg*(np.square(Vrg)+np.square(Vig))-(-1*Pg*Vig+Qg*Vrg)*2*Vrg)/np.square(np.square(Vrg)+np.square(Vig))
                    dIig_dQg = 1*Vrg/(np.square(Vrg)+np.square(Vig))
                    Iig = (-1*Pg*Vig+Qg*Vrg)/(np.square(Vrg)+np.square(Vig))

                    Vig_hist = -1*(Iig-dIig_dQg*Qg-dIig_dVrg*Vrg-dIig_dVig*Vig)

                    matrix_J[(bus.Bus-1)*2] = Vrg_hist
                    matrix_J[(bus.Bus-1)*2+1] = Vig_hist

                    #gen bus also has P and Q connected to them, and P and Q need to be added into the matrix.
                    for load in parsed_data['loads']:
                        if (load.Bus == gen.Bus):
                            Vrl = matrix_v[v_index*2]
                            Vil = matrix_v[v_index*2+1]
                            Pl = load.P
                            Ql = load.Q   

                            #Irl part of PQ
                            dIrl_dVrl = (Pl*(np.square(Vrl)+np.square(Vil))-(Pl*Vrl+1*Ql*Vil)*2*Vrl)/np.square(np.square(Vrl)+np.square(Vil))
                            dIrl_dVil = (Ql*(np.square(Vrl)+np.square(Vil))-(Pl*Vrl+1*Ql*Vil)*2*Vil)/np.square(np.square(Vrl)+np.square(Vil))
                            Irl = (Pl*Vrl+Ql*Vil)/(np.square(Vrl)+np.square(Vil))

                            Vrl_hist = -1*(Irl-dIrl_dVrl*Vrl-dIrl_dVil*Vil)

                            #Iil part of PQ
                            dIil_dVrl = (-1*Ql*(np.square(Vrl)+np.square(Vil))-(Pl*Vil-1*Ql*Vrl)*2*Vrl)/np.square(np.square(Vrl)+np.square(Vil))
                            dIil_dVil = (Pl*(np.square(Vrl)+np.square(Vil))-(Pl*Vil-1*Ql*Vrl)*2*Vil)/np.square(np.square(Vrl)+np.square(Vil))
                            Iil = (Pl*Vil-Ql*Vrl)/(np.square(Vrl)+np.square(Vil))

                            Vil_hist = -1*(Iil-dIil_dVrl*Vrl-dIil_dVil*Vil)
                    
                            matrix_J[(bus.Bus-1)*2] += Vrl_hist
                            matrix_J[(bus.Bus-1)*2+1] += Vil_hist
                            #print (v_index*2,v_index*2+1)
                    
                    #Veq_hist, Vrg and Vig
                    Veg_hist = np.square(VSet)+np.square(Vrg)+np.square(Vig)
                    matrix_J[bus_number*2+gen_number] = Veg_hist
                    #print (bus_number*2+gen_number)
                    v_index = v_index+1
            gen_number = gen_number+1
        return (matrix_v, matrix_J)

    def apply_limiting(self, v_last, matrix_v, v_delta_limit, bus_amount, size_Y):
        temp_v = np.copy(v_last)
        voltage_amount = bus_amount*2
        for index in range (0, size_Y):     
            #handle voltage bues
            if index < voltage_amount:
                if (matrix_v[index]-v_last[index])> v_delta_limit:
                    v_last[index] = v_last[index]+v_delta_limit
                    print ('v[',index,']', 'step too big', matrix_v[index]-v_last[index], 'before limiting:', matrix_v[index], 'after limiting:', v_last[index], 'previous v:', temp_v[index])
                if (v_last[index]-matrix_v[index])> v_delta_limit:
                    v_last[index] = v_last[index]-v_delta_limit
                    print ('v[',index,']', 'step too big', -1*(v_last[index]-matrix_v[index]), 'before limiting:', matrix_v[index], 'after limiting:', v_last[index], 'previous v:', temp_v[index])
                else:
                    v_last[index] = np.copy(matrix_v[index])
            #ignore Q and I of slacks
            else: 
                #print (index)
                v_last[index] = np.copy(matrix_v[index])
            pass

        return(v_last)

    def bus_number(self, parsed_data):
        bus_amount = 0
        for ele in parsed_data['buses']:
            bus_amount = bus_amount+1
        return (bus_amount)

    def check_error(self, matrix_v, v_last):
        error = np.amax(np.abs(v_last - matrix_v))
        return (error)

    def stamp_linear(self, v_init, row_Y, columun_Y, value_Y, parsed_data, size_Y, matrix_J):
        #update the Y matrix and the first J matrix
        bus_number =  0
        slack_number =  0
        gen_number = 0
        for bus in parsed_data['buses']:
            bus_number = bus_number+1 #each bus contributes 2 variables in V (Vr, Vi)
        for slack in parsed_data['slack']:
            slack_number = slack_number+1 #each slack contributes 2 variables in V. (Ir, Ii)
        for gen in parsed_data['generators']:
            gen_number = gen_number+1 #each generators contributes 1 variable in V. (Qg)   
        
        normal_Y_size = size_Y

        value_temp = 0
        i = 0
        
        #update Y for real variables in branches. （The main diagonal)
        #(0,0)(2,2)(4,4)(6,6)
        for bus in parsed_data['buses']:
            for branch in parsed_data['branches']:
                if bus.Bus == branch.from_bus:
                    value_temp = (branch.r/((branch.r)*(branch.r)+(branch.x)*(branch.x))) #G part of the branch
                    row_Y[i]= (bus.Bus-1)*2
                    columun_Y[i] = (bus.Bus-1)*2
                    value_Y[i] = value_Y[i]+value_temp 
                    #print (row_Y[i], columun_Y[i], i, branch.b)
                if bus.Bus == branch.to_bus:
                    value_temp = (branch.r/((branch.r)*(branch.r)+(branch.x)*(branch.x))) #G part of the branch
                    row_Y[i]= (bus.Bus-1)*2
                    columun_Y[i] = (bus.Bus-1)*2
                    value_Y[i] = value_Y[i]+value_temp 
                    #print (row_Y[i], columun_Y[i], i, branch.b)
            i = i+1

        #update Y for imaginary variables in branches. （The main diagonal)
        #(1,1)(3,3)(5,5)(7,7)
        for bus in parsed_data['buses']:
            for branch in parsed_data['branches']:
                if bus.Bus == branch.from_bus:
                    value_temp = (branch.r/((branch.r)*(branch.r)+(branch.x)*(branch.x))) #G part of the branch
                    row_Y[i]= bus.Bus*2-1
                    columun_Y[i] = bus.Bus*2-1
                    value_Y[i] = value_Y[i]+value_temp 
                if bus.Bus == branch.to_bus:
                    value_temp = (branch.r/((branch.r)*(branch.r)+(branch.x)*(branch.x))) #G part of the branch
                    row_Y[i]= bus.Bus*2-1
                    columun_Y[i] = bus.Bus*2-1
                    value_Y[i] = value_Y[i]+value_temp 

            i = i+1

        #update Y for real variables's imaginary partners in branches. （diagonal)
        #(0,1)(2,3)(4,5)(6,7)
        for bus in parsed_data['buses']:
            for branch in parsed_data['branches']:
                if bus.Bus == branch.from_bus:
                    value_temp = (branch.x/((branch.r)*(branch.r)+(branch.x)*(branch.x)))-(branch.b/2)#-B part of the branch
                    row_Y[i]= (bus.Bus-1)*2
                    columun_Y[i] = (bus.Bus-1)*2+1
                    value_Y[i] = value_Y[i]+value_temp 
                    #print (row_Y[i], columun_Y[i], i, branch.b)
                if bus.Bus == branch.to_bus:
                    value_temp = (branch.x/((branch.r)*(branch.r)+(branch.x)*(branch.x)))-(branch.b/2) #-B part of the branch
                    row_Y[i]= (bus.Bus-1)*2
                    columun_Y[i] = (bus.Bus-1)*2+1
                    value_Y[i] = value_Y[i]+value_temp 
                    #print (row_Y[i], columun_Y[i], i, branch.b)

            i = i+1

        #update Y for imaginary variables's real partners in branches. （diagonal)
        #(1,0)(3,2)(5,4)(7,6)
        for bus in parsed_data['buses']:
            for branch in parsed_data['branches']:
                if bus.Bus == branch.from_bus:
                    value_temp = -1*(branch.x/((branch.r)*(branch.r)+(branch.x)*(branch.x)))+(branch.b/2) #B part of the branch
                    row_Y[i]= bus.Bus*2-1
                    columun_Y[i] = bus.Bus*2-1-1
                    value_Y[i] = value_Y[i]+value_temp 
                    #print (row_Y[i], columun_Y[i], i, branch.b)
                if bus.Bus == branch.to_bus:
                    value_temp = -1*(branch.x/((branch.r)*(branch.r)+(branch.x)*(branch.x)))+(branch.b/2) #B part of the branch
                    row_Y[i]= bus.Bus*2-1
                    columun_Y[i] = bus.Bus*2-1-1
                    value_Y[i] = value_Y[i]+value_temp 
                    #print (row_Y[i], columun_Y[i], i, branch.b)

            i = i+1

        #update Y for real variables in branches. （Non-diagonal)
        for bus in parsed_data['buses']:
            for branch in parsed_data['branches']:
                if bus.Bus == branch.from_bus:
                    value_temp = -1*(branch.r/((branch.r)*(branch.r)+(branch.x)*(branch.x))) #G part of the branch
                    row_Y[i]= (bus.Bus-1)*2
                    columun_Y[i] = (branch.to_bus-1)*2
                    value_Y[i] = value_Y[i]+value_temp 
                    i = i+1
                    

        for bus in parsed_data['buses']:
            for branch in parsed_data['branches']:
                if bus.Bus == branch.to_bus:
                    value_temp = -1*(branch.r/((branch.r)*(branch.r)+(branch.x)*(branch.x))) #G part of the branch
                    row_Y[i]= (bus.Bus-1)*2
                    columun_Y[i] = (branch.from_bus-1)*2
                    value_Y[i] = value_Y[i]+value_temp 
                    i = i+1
                    

        for bus in parsed_data['buses']:
            for branch in parsed_data['branches']:
                if bus.Bus == branch.from_bus:
                    value_temp = -1*(branch.x/((branch.r)*(branch.r)+(branch.x)*(branch.x))) #B part of the branch
                    row_Y[i]= (bus.Bus-1)*2
                    columun_Y[i] = (branch.to_bus-1)*2+1
                    value_Y[i] = value_Y[i]+value_temp 
                    i = i+1

        for bus in parsed_data['buses']:
            for branch in parsed_data['branches']:
                if bus.Bus == branch.to_bus:
                    value_temp = -1*(branch.x/((branch.r)*(branch.r)+(branch.x)*(branch.x))) #B part of the branch
                    row_Y[i]= (bus.Bus-1)*2
                    columun_Y[i] = (branch.from_bus-1)*2+1
                    value_Y[i] = value_Y[i]+value_temp 
                    i = i+1        

        #update Y for imaginary variables in branches. （Non-diagonal)
        #(1, 2)(1, 4)(3, 6)(5, 6)
        for bus in parsed_data['buses']:#Fixed. G-> -B
            for branch in parsed_data['branches']:
                if bus.Bus == branch.from_bus:
                    value_temp = 1*(branch.x/((branch.r)*(branch.r)+(branch.x)*(branch.x))) #G part of the branch
                    row_Y[i]= (bus.Bus)*2-1
                    columun_Y[i] = (branch.to_bus-1)*2
                    value_Y[i] = value_Y[i]+value_temp 
                    i = i+1

        #(3, 0)(5, 0)(7, 2)(7, 4)
        for bus in parsed_data['buses']:#Fixe -G-> -B
            for branch in parsed_data['branches']:
                if bus.Bus == branch.to_bus:
                    value_temp = 1*(branch.x/((branch.r)*(branch.r)+(branch.x)*(branch.x))) #G part of the branch
                    row_Y[i]= (bus.Bus-1)*2+1
                    columun_Y[i] = (branch.from_bus-1)*2
                    value_Y[i] = value_Y[i]+value_temp 
                    i = i+1        

        #(1, 3)(1, 5)(3, 7)(5, 7) 
        for bus in parsed_data['buses']:#Fixe -B-> -G
            for branch in parsed_data['branches']:
                if bus.Bus == branch.from_bus:
                    value_temp = -1*(branch.r/((branch.r)*(branch.r)+(branch.x)*(branch.x))) #B part of the branch
                    row_Y[i]= (bus.Bus-1)*2+1
                    columun_Y[i] = (branch.to_bus-1)*2+1
                    value_Y[i] = value_Y[i]+value_temp 
                    i = i+1

        #(3, 1)(5, 1)(7, 3)(7, 5)
        for bus in parsed_data['buses']:#Fixe -B-> -G
            for branch in parsed_data['branches']:
                if bus.Bus == branch.to_bus:
                    value_temp = -1*(branch.r/((branch.r)*(branch.r)+(branch.x)*(branch.x))) #B part of the branch
                    row_Y[i]= (bus.Bus-1)*2+1
                    columun_Y[i] = (branch.from_bus-1)*2+1
                    value_Y[i] = value_Y[i]+value_temp 
                    i = i+1   
        
        bus_number =  0
        slack_number =  0
        gen_number = 0
        trans_number = 0
        for bus in parsed_data['buses']:
            bus_number = bus_number+1 #each bus contributes 2 variables in V (Vr, Vi)
        for slack in parsed_data['slack']:
            slack_number = slack_number+1 #each slack contributes 2 variables in V. (Ir, Ii)
        for gen in parsed_data['generators']:
            gen_number = gen_number+1 #each generators contributes 1 variable in V. (Qg)  
        for trans in parsed_data['xfmrs']:
            trans_number = trans_number+1

        #update Y for Ir1, Ii1 of the slack
        for slack in parsed_data['slack']:#fixed for multiple slack buses
            # for Ir
            row_Y[i] = (slack.Bus-1)*2
            columun_Y[i] = bus_number*2+gen_number+(slack.Bus-1)*2
            value_Y[i] = 1
            i = i+1
            
            # for Ii
            row_Y[i] = (slack.Bus-1)*2+1
            columun_Y[i] = bus_number*2+gen_number+(slack.Bus-1)*2+1
            value_Y[i] = 1
            i = i+1     
            
            # for Vset_r
            row_Y[i] = bus_number*2+gen_number+(slack.Bus-1)*2 
            columun_Y[i] = (slack.Bus-1)*2
            value_Y[i] = 1
            i = i+1

            # for Vset_i
            row_Y[i] = bus_number*2+gen_number+(slack.Bus-1)*2+1 
            columun_Y[i] = (slack.Bus-1)*2+1
            value_Y[i] = 1
            i = i+1 

        #add empty sapce for generators
        gen_index = 0
        for gen in parsed_data['generators']:
            for bus in parsed_data['buses']:
                #add empty sapce for dIrg_dQ and dIig_dQ
                if bus.Bus == gen.Bus: #(2, 28)(4, 29)(10, 30)(14, 31)
                    row_Y[i] = (bus.Bus-1)*2
                    columun_Y[i] = (bus_number)*2+gen_index
                    value_Y[i] = 0
                    #print(row_Y[i], columun_Y[i])
                    i = i+1 
                
                if bus.Bus == gen.Bus:#(3 28)(5 29)(11 30)(15 31)
                    row_Y[i] = (bus.Bus-1)*2+1
                    columun_Y[i] = (bus_number)*2+gen_index
                    value_Y[i] = 0
                    i = i+1 
        
                #add empty sapce for 2Vrg and 2Vig
                if bus.Bus == gen.Bus:
                    row_Y[i] = (bus_number)*2+gen_index
                    columun_Y[i] = (bus.Bus-1)*2
                    value_Y[i] = 0
                    #print(row_Y[i], columun_Y[i])
                    i = i+1 

                if bus.Bus == gen.Bus: #(28 3)(29 5)(30 11)(31 15)
                    row_Y[i] = (bus_number)*2+gen_index
                    columun_Y[i] = (bus.Bus-1)*2+1
                    value_Y[i] = 0
                    #print(, row_Y[i], columun_Y[i])
                    i = i+1 
            gen_index = gen_index+1

        #add empty sapce for transformer
        trans_index = 0
        for trans in parsed_data['xfmrs']:
            for bus in parsed_data['buses']:
                if (trans.to_bus) == (bus.Bus): 
                    
                    #(12, 34)(16, 36)(10. 38)
                    row_Y[i] = (bus.Bus-1)*2
                    columun_Y[i] = (bus_number*2+gen_number+slack_number*2)+(trans_index)*2
                    value_Y[i] = 0
                    i = i+1
                    
                    #(12, 35)(16, 37)(10. 39)
                    row_Y[i] = (bus.Bus-1)*2
                    columun_Y[i] = (bus_number*2+gen_number+slack_number*2)+(trans_index)*2+1
                    value_Y[i] = 0
                    i = i+1

                    #(13, 34)(17, 36)(11, 38)
                    row_Y[i] = (bus.Bus-1)*2+1
                    columun_Y[i] = (bus_number*2+gen_number+slack_number*2)+(trans_index)*2
                    value_Y[i] = 0
                    i = i+1

                    #(13, 35)(17, 37)(11, 39)
                    row_Y[i] = (bus.Bus-1)*2+1
                    columun_Y[i] = (bus_number*2+gen_number+slack_number*2)+(trans_index)*2+1
                    value_Y[i] = 0
                    i = i+1

                    #(34, 34)(36, 36)(38. 38)
                    row_Y[i] = (bus_number*2+gen_number+slack_number*2)+(trans_index)*2
                    columun_Y[i] = (bus_number*2+gen_number+slack_number*2)+(trans_index)*2
                    value_Y[i] = 0
                    i = i+1

                    #(34, 35)(36, 37)(38. 39)
                    row_Y[i] = (bus_number*2+gen_number+slack_number*2)+(trans_index)*2
                    columun_Y[i] = (bus_number*2+gen_number+slack_number*2)+(trans_index)*2+1
                    value_Y[i] = 0
                    i = i+1

                    #(35, 34)(37, 36)(39. 38)
                    row_Y[i] = (bus_number*2+gen_number+slack_number*2)+(trans_index)*2+1
                    columun_Y[i] = (bus_number*2+gen_number+slack_number*2)+(trans_index)*2
                    value_Y[i] = 0
                    i = i+1

                    #(35, 35)(37, 37)(39. 39)
                    row_Y[i] = (bus_number*2+gen_number+slack_number*2)+(trans_index)*2+1
                    columun_Y[i] = (bus_number*2+gen_number+slack_number*2)+(trans_index)*2+1
                    value_Y[i] = 0
                    i = i+1

                    #(34, 12)(36, 16)(38. 10)
                    row_Y[i] = (bus_number*2+gen_number+slack_number*2)+(trans_index)*2
                    columun_Y[i] = (bus.Bus-1)*2
                    value_Y[i] = 0
                    i = i+1

                    #(34, 13)(36, 17)(38. 11)
                    row_Y[i] = (bus_number*2+gen_number+slack_number*2)+(trans_index)*2
                    columun_Y[i] = (bus.Bus-1)*2+1
                    value_Y[i] = 0
                    i = i+1

                    #(35, 12)(37, 16)(39. 10)
                    row_Y[i] = (bus_number*2+gen_number+slack_number*2)+(trans_index)*2+1
                    columun_Y[i] = (bus.Bus-1)*2
                    value_Y[i] = 0
                    i = i+1

                    #(35, 13)(37, 17)(39. 11)
                    row_Y[i] = (bus_number*2+gen_number+slack_number*2)+(trans_index)*2+1
                    columun_Y[i] = (bus.Bus-1)*2+1
                    value_Y[i] = 0
                    i = i+1

                    #(34, 40)(36, 42)(38, 44)
                    row_Y[i] = (bus_number*2+gen_number+slack_number*2)+(trans_index)*2
                    columun_Y[i] = (bus_number*2+gen_number+slack_number*2)+trans_number*2+(trans_index)*2
                    value_Y[i] = 0
                    i = i+1

                    #(34, 41)(36, 43)(38, 45)
                    row_Y[i] = (bus_number*2+gen_number+slack_number*2)+(trans_index)*2
                    columun_Y[i] = (bus_number*2+gen_number+slack_number*2)+trans_number*2+(trans_index)*2+1
                    value_Y[i] = 0
                    i = i+1

                    #(35, 40)(37, 42)(39, 44)
                    row_Y[i] = (bus_number*2+gen_number+slack_number*2)+(trans_index)*2+1
                    columun_Y[i] = (bus_number*2+gen_number+slack_number*2)+trans_number*2+(trans_index)*2
                    value_Y[i] = 0
                    i = i+1

                    #(35, 41)(37, 43)(39, 45)
                    row_Y[i] = (bus_number*2+gen_number+slack_number*2)+(trans_index)*2+1
                    columun_Y[i] = (bus_number*2+gen_number+slack_number*2)+trans_number*2+(trans_index)*2+1
                    value_Y[i] = 0
                    i = i+1
                    
                    #(40, 12)(42, 16)(44, 10)
                    row_Y[i] = (bus_number*2+gen_number+slack_number*2)+trans_number*2+(trans_index)*2
                    columun_Y[i] = (bus.Bus-1)*2
                    value_Y[i] = 0
                    i = i+1

                    #(40, 13)(42, 17)(44, 11)
                    row_Y[i] = (bus_number*2+gen_number+slack_number*2)+trans_number*2+(trans_index)*2
                    columun_Y[i] = (bus.Bus-1)*2+1
                    value_Y[i] = 0
                    i = i+1

                    #(41, 12)(43, 16)(45, 10)
                    row_Y[i] = (bus_number*2+gen_number+slack_number*2)+trans_number*2+(trans_index)*2+1
                    columun_Y[i] = (bus.Bus-1)*2
                    value_Y[i] = 0
                    i = i+1

                    #(41, 13)(43, 17)(45, 11)
                    row_Y[i] = (bus_number*2+gen_number+slack_number*2)+trans_number*2+(trans_index)*2+1
                    columun_Y[i] = (bus.Bus-1)*2+1
                    value_Y[i] = 0
                    i = i+1

                    trans_index= trans_index+1
        


        trans_index = 0
        for trans in parsed_data['xfmrs']:
            for bus in parsed_data['buses']:
                if (trans.from_bus) == (bus.Bus):    
                    #(40, 6) (42, 6)(44, 8)
                    row_Y[i] = (bus_number*2+gen_number+slack_number*2)+trans_number*2+(trans_index)*2 
                    columun_Y[i] = (bus.Bus-1)*2
                    value_Y[i] = 1  
                    i = i+1 

                    #(41, 7) (43, 7)(45, 9)
                    row_Y[i] = (bus_number*2+gen_number+slack_number*2)+trans_number*2+(trans_index)*2+1 
                    columun_Y[i] = (bus.Bus-1)*2+1
                    value_Y[i] = 1  
                    i = i+1

                    #(6, 40) (6, 42)(8, 44)
                    row_Y[i] = (bus.Bus-1)*2 
                    columun_Y[i] = (bus_number*2+gen_number+slack_number*2)+trans_number*2+(trans_index)*2
                    value_Y[i] = 1  
                    i = i+1

                    #(7, 41) (7, 43)(9, 45)
                    row_Y[i] = (bus.Bus-1)*2+1 
                    columun_Y[i] = (bus_number*2+gen_number+slack_number*2)+trans_number*2+(trans_index)*2+1
                    value_Y[i] = 1  
                    #print (i, row_Y[i], columun_Y[i], value_Y[i])
                    i = i+1

                    trans_index= trans_index+1 
        
        trans_index = 0
        for trans in parsed_data['xfmrs']:
            for bus in parsed_data['buses']:
                if (trans.to_bus) == (bus.Bus):
                    for i1 in range (0, i): 
                        if row_Y[i1] == (bus.Bus-1)*2: #(10, 10)(12, 12)(16, 16)
                            if columun_Y[i1] == (bus.Bus-1)*2:
                                value_temp = trans.r/(np.square(trans.r)+np.square(trans.x))
                                value_Y[i1] = value_Y[i1]+ value_temp
                                
                        if row_Y[i1] == (bus.Bus-1)*2: #(12, 13)(16, 17)(10, 11)
                            if columun_Y[i1] == (bus.Bus-1)*2+1:
                                value_temp = trans.x/(np.square(trans.r)+np.square(trans.x))
                                value_Y[i1] = value_Y[i1]+ value_temp
                                
                        if row_Y[i1] == (bus.Bus-1)*2+1: #(13, 12)(17, 16)(11, 10)
                            if columun_Y[i1] == (bus.Bus-1)*2:
                                value_temp = -1*trans.x/(np.square(trans.r)+np.square(trans.x))
                                value_Y[i1] = value_Y[i1]+ value_temp
                                
                        if row_Y[i1] == (bus.Bus-1)*2+1: #(13, 13)(17, 17)(11, 11)
                            if columun_Y[i1] == (bus.Bus-1)*2+1:
                                value_temp = trans.r/(np.square(trans.r)+np.square(trans.x))
                                value_Y[i1] = value_Y[i1]+ value_temp
                                
                        if row_Y[i1] == (bus.Bus-1)*2: #(12, 34)(16, 36)(10, 38)
                            if columun_Y[i1] == (bus_number*2+gen_number+slack_number*2)+(trans_index)*2:
                                value_temp = -1*trans.r/(np.square(trans.r)+np.square(trans.x))
                                value_Y[i1] = value_Y[i1]+ value_temp

                        if row_Y[i1] == (bus.Bus-1)*2: #(12, 35)(16, 37)(10, 39)
                            if columun_Y[i1] == (bus_number*2+gen_number+slack_number*2)+(trans_index)*2+1:
                                value_temp = -1*trans.x/(np.square(trans.r)+np.square(trans.x))
                                value_Y[i1] = value_Y[i1]+ value_temp
                        
                        if row_Y[i1] == (bus.Bus-1)*2+1: #(13, 34)(17, 36)(11, 38)
                            if columun_Y[i1] == (bus_number*2+gen_number+slack_number*2)+(trans_index)*2:
                                value_temp = trans.x/(np.square(trans.r)+np.square(trans.x))
                                value_Y[i1] = value_Y[i1]+ value_temp

                        if row_Y[i1] == (bus.Bus-1)*2+1: #(13, 35)(17, 37)(11, 39)
                            if columun_Y[i1] == (bus_number*2+gen_number+slack_number*2)+(trans_index)*2+1:
                                value_temp = -1*trans.r/(np.square(trans.r)+np.square(trans.x))
                                value_Y[i1] = value_Y[i1]+ value_temp

                        if row_Y[i1] == (bus_number*2+gen_number+slack_number*2)+(trans_index)*2: #(34, 34)(36, 36)(38, 38)
                            if columun_Y[i1] == (bus_number*2+gen_number+slack_number*2)+(trans_index)*2:
                                value_temp = trans.r/(np.square(trans.r)+np.square(trans.x))
                                value_Y[i1] = value_Y[i1]+ value_temp

                        if row_Y[i1] == (bus_number*2+gen_number+slack_number*2)+(trans_index)*2: #(34, 35)(36, 37)(38, 39)
                            if columun_Y[i1] == (bus_number*2+gen_number+slack_number*2)+(trans_index)*2+1:
                                value_temp = trans.x/(np.square(trans.r)+np.square(trans.x))
                                value_Y[i1] = value_Y[i1]+ value_temp

                        if row_Y[i1] == (bus_number*2+gen_number+slack_number*2)+(trans_index)*2+1: #(35, 34)(37, 36)(39, 38)
                            if columun_Y[i1] == (bus_number*2+gen_number+slack_number*2)+(trans_index)*2:
                                value_temp = -1*trans.x/(np.square(trans.r)+np.square(trans.x))
                                value_Y[i1] = value_Y[i1]+ value_temp

                        if row_Y[i1] == (bus_number*2+gen_number+slack_number*2)+(trans_index)*2+1: #(35, 34)(37, 36)(39, 38)
                            if columun_Y[i1] == (bus_number*2+gen_number+slack_number*2)+(trans_index)*2+1:
                                value_temp = trans.r/(np.square(trans.r)+np.square(trans.x))
                                value_Y[i1] = value_Y[i1]+ value_temp

                        if row_Y[i1] == (bus_number*2+gen_number+slack_number*2)+(trans_index)*2: #(34, 12)(36, 16)(38, 10)
                            if columun_Y[i1] == (bus.Bus-1)*2:
                                value_temp = -1*trans.r/(np.square(trans.r)+np.square(trans.x))
                                value_Y[i1] = value_Y[i1]+ value_temp

                        if row_Y[i1] == (bus_number*2+gen_number+slack_number*2)+(trans_index)*2: #(34, 13)(36, 17)(38, 11)
                            if columun_Y[i1] == (bus.Bus-1)*2+1:
                                value_temp = -1*trans.x/(np.square(trans.r)+np.square(trans.x))
                                value_Y[i1] = value_Y[i1]+ value_temp

                        if row_Y[i1] == (bus_number*2+gen_number+slack_number*2)+(trans_index)*2+1: #(35, 12)(37, 16)(39, 10)
                            if columun_Y[i1] == (bus.Bus-1)*2:
                                value_temp = trans.x/(np.square(trans.r)+np.square(trans.x))
                                value_Y[i1] = value_Y[i1]+ value_temp

                        if row_Y[i1] == (bus_number*2+gen_number+slack_number*2)+(trans_index)*2+1: #(35, 13)(37, 17)(39, 11)
                            if columun_Y[i1] == (bus.Bus-1)*2+1:
                                value_temp = -1*trans.r/(np.square(trans.r)+np.square(trans.x))
                                value_Y[i1] = value_Y[i1]+ value_temp
                                
                        if row_Y[i1] == (bus_number*2+gen_number+slack_number*2)+(trans_index)*2: #(34, 40)(36, 42)(38, 44)
                            if columun_Y[i1] == (bus_number*2+gen_number+slack_number*2)+trans_number*2+(trans_index)*2:
                                angle = trans.ang*np.pi/180
                                value_temp = -1*trans.tr*np.cos(angle)
                                value_Y[i1] = value_Y[i1]+ value_temp

                        if row_Y[i1] == (bus_number*2+gen_number+slack_number*2)+(trans_index)*2: #(34, 41)(36, 43)(38, 45)
                            if columun_Y[i1] == (bus_number*2+gen_number+slack_number*2)+trans_number*2+(trans_index)*2+1:
                                angle = trans.ang*np.pi/180
                                value_temp = -1*trans.tr*np.sin(angle)
                                value_Y[i1] = value_Y[i1]+ value_temp

                        if row_Y[i1] == (bus_number*2+gen_number+slack_number*2)+(trans_index)*2+1: #(35, 40)(37, 42)(39, 44)
                            if columun_Y[i1] == (bus_number*2+gen_number+slack_number*2)+trans_number*2+(trans_index)*2:
                                angle = trans.ang*np.pi/180
                                value_temp = trans.tr*np.sin(angle)
                                value_Y[i1] = value_Y[i1]+ value_temp
                        
                        if row_Y[i1] == (bus_number*2+gen_number+slack_number*2)+(trans_index)*2+1: #(35, 41)(37, 43)(39, 45)
                            if columun_Y[i1] == (bus_number*2+gen_number+slack_number*2)+trans_number*2+(trans_index)*2+1:
                                angle = trans.ang*np.pi/180
                                value_temp = -1*trans.tr*np.cos(angle)
                                value_Y[i1] = value_Y[i1]+ value_temp

                        if row_Y[i1] == (bus_number*2+gen_number+slack_number*2)+trans_number*2+(trans_index)*2: #(40, 12)(42, 16)(44, 10)
                            if columun_Y[i1] == (bus.Bus-1)*2:
                                angle = trans.ang*np.pi/180
                                value_temp = trans.tr*np.cos(angle)*(-1)
                                value_Y[i1] = value_Y[i1]+ value_temp

                        if row_Y[i1] == (bus_number*2+gen_number+slack_number*2)+trans_number*2+(trans_index)*2: #(40, 13)(42, 17)(44, 11)
                            if columun_Y[i1] == (bus.Bus-1)*2+1:
                                angle = trans.ang*np.pi/180
                                value_temp = -1*trans.tr*np.sin(angle)*(-1)
                                value_Y[i1] = value_Y[i1]+ value_temp

                        if row_Y[i1] == (bus_number*2+gen_number+slack_number*2)+trans_number*2+(trans_index)*2+1: #(41, 12)(43, 16)(45, 10)
                            if columun_Y[i1] == (bus.Bus-1)*2:
                                angle = trans.ang*np.pi/180
                                value_temp = trans.tr*np.sin(angle)*(-1)
                                value_Y[i1] = value_Y[i1]+ value_temp

                        if row_Y[i1] == (bus_number*2+gen_number+slack_number*2)+trans_number*2+(trans_index)*2+1: #(41, 13)(43, 17)(45, 11)
                            if columun_Y[i1] == (bus.Bus-1)*2+1:
                                angle = trans.ang*np.pi/180
                                value_temp = trans.tr*np.cos(angle)*(-1)
                                value_Y[i1] = value_Y[i1]+ value_temp
                                #print (row_Y[i1], columun_Y[i1], value_Y[i1], angle)


                    trans_index = trans_index+1

        for shunt in parsed_data['shunts']:
            for bus in parsed_data['buses']:
                if shunt.Bus == bus.Bus:    
                    for i1 in range (0, i):
                        if row_Y[i1] == (bus.Bus-1)*2:#(16 16)
                            if columun_Y[i1] == (bus.Bus-1)*2:
                                value_Y[i1] = value_Y[i1]+ shunt.G_MW #changed to p.u already.
                        if row_Y[i1] == (bus.Bus-1)*2:#(16 17)
                            if columun_Y[i1] == ((bus.Bus-1)*2+1):
                                value_Y[i1] = value_Y[i1]+shunt.B_MVAR*(-1) #changed to p.u already.
                        if row_Y[i1] == (bus.Bus-1)*2+1:#(17 16)
                            if columun_Y[i1] == ((bus.Bus-1)*2):
                                value_Y[i1] = value_Y[i1]+shunt.B_MVAR #changed to p.u already.
                        if row_Y[i1] == (bus.Bus-1)*2+1:#(17 16)
                            if columun_Y[i1] == ((bus.Bus-1)*2)+1:
                                value_Y[i1] = value_Y[i1]+shunt.G_MW #changed to p.u already.
                                #print ('!!!!!', row_Y[i1], columun_Y[i1], shunt.G_MW)

        #Create the first J
        #for slack bus
        i_slack = 0
        for slack in parsed_data['slack']:
            angle = slack.ang*np.pi/180 #Va is a constant for the slack bus.
            Vset_r = slack.Vset*np.cos(angle)
            Vset_l = slack.Vset*np.sin(angle)
            
            matrix_J[bus_number*2+gen_number*1+(i_slack)*2] = Vset_r
            matrix_J[bus_number*2+gen_number*1+(i_slack)*2+1] = Vset_l

            i_slack = i_slack+1

        #for PQ bus
        i_index = 0
        #slack bus also has P and Q connected to them, and P and Q need to be added into the matrix.
        for slack in parsed_data['slack']:
            for load in parsed_data['loads']:
                if (load.Bus == slack.Bus):
                    for bus in parsed_data['buses']:
                        if (load.Bus == bus.Bus):
                            angle = bus.Va_init*np.pi/180 #Va is not just init, it will be updated. 
                            Vrl = bus.Vm_init*np.cos(angle) #Vm is not just init, it will be updated.
                            Vil = bus.Vm_init*np.sin(angle)
                            Pl = load.P
                            Ql = load.Q
                            #print (Vrl, Vil, Pl, Ql)

                            #Irl part of PQ
                            dIrl_dVrl = (Pl*(np.square(Vrl)+np.square(Vil))-(Pl*Vrl+1*Ql*Vil)*2*Vrl)/np.square(np.square(Vrl)+np.square(Vil))
                            dIrl_dVil = (Ql*(np.square(Vrl)+np.square(Vil))-(Pl*Vrl+1*Ql*Vil)*2*Vil)/np.square(np.square(Vrl)+np.square(Vil))
                            Irl = (Pl*Vrl+Ql*Vil)/(np.square(Vrl)+np.square(Vil))

                            Vrl_hist = -1*(Irl-dIrl_dVrl*Vrl-dIrl_dVil*Vil)

                            #Iil part of PQ
                            dIil_dVrl = (-1*Ql*(np.square(Vrl)+np.square(Vil))-(Pl*Vil-1*Ql*Vrl)*2*Vrl)/np.square(np.square(Vrl)+np.square(Vil))
                            dIil_dVil = (Pl*(np.square(Vrl)+np.square(Vil))-(Pl*Vil-1*Ql*Vrl)*2*Vil)/np.square(np.square(Vrl)+np.square(Vil))
                            Iil = (Pl*Vil-Ql*Vrl)/(np.square(Vrl)+np.square(Vil))

                            Vil_hist = -1*(Iil-dIil_dVrl*Vrl-dIil_dVil*Vil)
                            matrix_J[(bus.Bus-1)*2] = Vrl_hist
                            matrix_J[(bus.Bus-1)*2+1] = Vil_hist
                            #print ((bus.Bus-1)*2, (bus.Bus-1)*2+1)

            i_index = i_index+1

        for load in parsed_data['loads']:
            for bus in parsed_data['buses']:
                if (load.Bus == bus.Bus) and (bus.Type == 1):
                    angle = bus.Va_init*np.pi/180 #Va is not just init, it will be updated. 
                    Vrl = bus.Vm_init*np.cos(angle) #Vm is not just init, it will be updated.
                    Vil = bus.Vm_init*np.sin(angle)
                    Pl = load.P
                    Ql = load.Q  

                    #Irl part of PQ
                    dIrl_dVrl = (Pl*(np.square(Vrl)+np.square(Vil))-(Pl*Vrl+1*Ql*Vil)*2*Vrl)/np.square(np.square(Vrl)+np.square(Vil))
                    dIrl_dVil = (Ql*(np.square(Vrl)+np.square(Vil))-(Pl*Vrl+1*Ql*Vil)*2*Vil)/np.square(np.square(Vrl)+np.square(Vil))
                    Irl = (Pl*Vrl+Ql*Vil)/(np.square(Vrl)+np.square(Vil))

                    Vrl_hist = -1*(Irl-dIrl_dVrl*Vrl-dIrl_dVil*Vil)

                    #Iil part of PQ
                    dIil_dVrl = (-1*Ql*(np.square(Vrl)+np.square(Vil))-(Pl*Vil-1*Ql*Vrl)*2*Vrl)/np.square(np.square(Vrl)+np.square(Vil))
                    dIil_dVil = (Pl*(np.square(Vrl)+np.square(Vil))-(Pl*Vil-1*Ql*Vrl)*2*Vil)/np.square(np.square(Vrl)+np.square(Vil))
                    Iil = (Pl*Vil-Ql*Vrl)/(np.square(Vrl)+np.square(Vil))

                    Vil_hist = -1*(Iil-dIil_dVrl*Vrl-dIil_dVil*Vil)
                    
                    matrix_J[(bus.Bus-1)*2] = Vrl_hist
                    matrix_J[(bus.Bus-1)*2+1] = Vil_hist
                    #print ((bus.Bus-1)*2, (bus.Bus-1)*2+1)
                    i_index = i_index+1
        
        gen_index = 0
        for gen in parsed_data['generators']:
            for bus in parsed_data['buses']:
                if bus.Bus == gen.Bus:
                    angle = bus.Va_init*np.pi/180 #it's not just init, it will be updated. 
                    Vrg = bus.Vm_init*np.cos(angle)
                    Vig = bus.Vm_init*np.sin(angle)
                    Pg = gen.P
                    Qg = gen.Qinit #it's not just init, it will be updated. 
                    VSet = gen.Vset

                    #Irg part of PV
                    dIrg_dVrg = (-1*Pg*(np.square(Vrg)+np.square(Vig))-(-1*Pg*Vrg-1*Qg*Vig)*2*Vrg)/np.square(np.square(Vrg)+np.square(Vig))
                    dIrg_dVig = (-1*Qg*(np.square(Vrg)+np.square(Vig))-(-1*Pg*Vrg-1*Qg*Vig)*2*Vig)/np.square(np.square(Vrg)+np.square(Vig))
                    dIrg_dQg = -1*Vig/(np.square(Vrg)+np.square(Vig))
                    Irg = (-1*Pg*Vrg-Qg*Vig)/(np.square(Vrg)+np.square(Vig))

                    Vrg_hist = -1*(Irg-dIrg_dQg*Qg-dIrg_dVrg*Vrg-dIrg_dVig*Vig)

                    #Iig part of PV
                    dIig_dVig = (-1*Pg*(np.square(Vrg)+np.square(Vig))-(-1*Pg*Vig+Qg*Vrg)*2*Vig)/np.square(np.square(Vrg)+np.square(Vig))
                    dIig_dVrg = (1*Qg*(np.square(Vrg)+np.square(Vig))-(-1*Pg*Vig+Qg*Vrg)*2*Vrg)/np.square(np.square(Vrg)+np.square(Vig))
                    dIig_dQg = 1*Vrg/(np.square(Vrg)+np.square(Vig))
                    Iig = (-1*Pg*Vig+Qg*Vrg)/(np.square(Vrg)+np.square(Vig))

                    Vig_hist = -1*(Iig-dIig_dQg*Qg-dIig_dVrg*Vrg-dIig_dVig*Vig)

                    matrix_J[(bus.Bus-1)*2] = Vrg_hist
                    matrix_J[(bus.Bus-1)*2+1] = Vig_hist

                    #print ((bus.Bus-1)*2, (bus.Bus-1)*2+1)

                    #Veq_hist, Vrg and Vig
                    Veg_hist = np.square(VSet)+np.square(Vrg)+np.square(Vig)
                    matrix_J[bus_number*2+gen_index] = Veg_hist
                    #print (bus_number*2+gen_index)

                    i_index = i_index+1

                    #gen bus also has P and Q connected to them, and P and Q need to be added into the matrix.
                    for load in parsed_data['loads']:
                        if (load.Bus == gen.Bus):
                            angle = bus.Va_init*np.pi/180 #Va is not just init, it will be updated. 
                            Vrl = bus.Vm_init*np.cos(angle) #Vm is not just init, it will be updated.
                            Vil = bus.Vm_init*np.sin(angle)
                            Pl = load.P
                            Ql = load.Q  

                            #Irl part of PQ
                            dIrl_dVrl = (Pl*(np.square(Vrl)+np.square(Vil))-(Pl*Vrl+1*Ql*Vil)*2*Vrl)/np.square(np.square(Vrl)+np.square(Vil))
                            dIrl_dVil = (Ql*(np.square(Vrl)+np.square(Vil))-(Pl*Vrl+1*Ql*Vil)*2*Vil)/np.square(np.square(Vrl)+np.square(Vil))
                            Irl = (Pl*Vrl+Ql*Vil)/(np.square(Vrl)+np.square(Vil))

                            Vrl_hist = -1*(Irl-dIrl_dVrl*Vrl-dIrl_dVil*Vil)

                            #Iil part of PQ
                            dIil_dVrl = (-1*Ql*(np.square(Vrl)+np.square(Vil))-(Pl*Vil-1*Ql*Vrl)*2*Vrl)/np.square(np.square(Vrl)+np.square(Vil))
                            dIil_dVil = (Pl*(np.square(Vrl)+np.square(Vil))-(Pl*Vil-1*Ql*Vrl)*2*Vil)/np.square(np.square(Vrl)+np.square(Vil))
                            Iil = (Pl*Vil-Ql*Vrl)/(np.square(Vrl)+np.square(Vil))

                            Vil_hist = -1*(Iil-dIil_dVrl*Vrl-dIil_dVil*Vil)
                    
                            matrix_J[(bus.Bus-1)*2] += Vrl_hist
                            matrix_J[(bus.Bus-1)*2+1] += Vil_hist
                            #print ((bus.Bus-1)*2, (bus.Bus-1)*2+1)
                    
            gen_index = gen_index+1

            

        return (row_Y, columun_Y, value_Y, matrix_J)

    def stamp_nonlinear(self, err_max, v_init, row_Y, columun_Y, value_Y, parsed_data, v_last, size_Y_csc):
        bus_number =  0
        slack_number =  0
        gen_number = 0
        trans_number = 0
        for bus in parsed_data['buses']:
            bus_number = bus_number+1 #each bus contributes 2 variables in V (Vr, Vi)
        for slack in parsed_data['slack']:
            slack_number = slack_number+1 #each slack contributes 2 variables in V. (Ir, Ii)
        for gen in parsed_data['generators']:
            gen_number = gen_number+1 #each generators contributes 1 variable in V. (Qg)  
        for trans in parsed_data['xfmrs']:
            trans_number = trans_number+1
        
        v_index = 0
        for slack in parsed_data['slack']:
            for bus in parsed_data['buses']:
                if (slack.Bus == bus.Bus):
                    for load in parsed_data['loads']:
                        if (load.Bus == bus.Bus): 
                            Vrl = v_last[v_index*2]
                            Vil = v_last[v_index*2+1]
                            Pl = load.P
                            Ql = load.Q  
                            #print ("here!!!", Vrl, Vil, Pl, Ql)

                            #Irl part of PQ
                            dIrl_dVrl = (Pl*(np.square(Vrl)+np.square(Vil))-(Pl*Vrl+1*Ql*Vil)*2*Vrl)/np.square(np.square(Vrl)+np.square(Vil))
                            dIrl_dVil = (Ql*(np.square(Vrl)+np.square(Vil))-(Pl*Vrl+1*Ql*Vil)*2*Vil)/np.square(np.square(Vrl)+np.square(Vil))
                            Irl = (Pl*Vrl+Ql*Vil)/(np.square(Vrl)+np.square(Vil))

                            #Iil part of PQ
                            dIil_dVrl = (-1*Ql*(np.square(Vrl)+np.square(Vil))-(Pl*Vil-1*Ql*Vrl)*2*Vrl)/np.square(np.square(Vrl)+np.square(Vil))
                            dIil_dVil = (Pl*(np.square(Vrl)+np.square(Vil))-(Pl*Vil-1*Ql*Vrl)*2*Vil)/np.square(np.square(Vrl)+np.square(Vil))
                            Iil = (Pl*Vil-Ql*Vrl)/(np.square(Vrl)+np.square(Vil))

                            
                            for i1 in range (0, size_Y_csc):
                                #add dIrl_dVrl to the Y matrix (0, 0)
                                if row_Y[i1] == (bus.Bus-1)*2:
                                    if columun_Y[i1] == (bus.Bus-1)*2:
                                        value_Y[i1]=value_Y[i1]+dIrl_dVrl
                                        #print("here!!!!",row_Y[i1], columun_Y[i1])
                                #add dIrl_dVil to the Y matrix (0, 1)
                                if row_Y[i1] == (bus.Bus-1)*2:
                                    if columun_Y[i1] == (bus.Bus-1)*2+1:
                                        value_Y[i1]=value_Y[i1]+dIrl_dVil
                                        #print("here!!!!", row_Y[i1], columun_Y[i1])
                                #add dIil_dVrl to the Y matrix (1, 0)
                                if row_Y[i1] == (bus.Bus-1)*2+1:
                                    if columun_Y[i1] == (bus.Bus-1)*2:
                                        value_Y[i1]=value_Y[i1]+dIil_dVrl
                                        #print("here!!!!", row_Y[i1], columun_Y[i1])
                                #add dIil_dVil to the Y matrix(1, 1)
                                if row_Y[i1] == (bus.Bus-1)*2+1:
                                    if columun_Y[i1] == (bus.Bus-1)*2+1:
                                        value_Y[i1]=value_Y[i1]+dIil_dVil
                                        #print("here!!!!", row_Y[i1], columun_Y[i1])

                            #print ("here!!!", (bus.Bus-1)*2, (bus.Bus-1)*2+1)
                            v_index = v_index+1

        for load in parsed_data['loads']:
            for bus in parsed_data['buses']:
                if (load.Bus == bus.Bus) and (bus.Type == 1):
                    #abandoned the angle calculation because it creates "sign flips".
                    #angle = bus.Va_init*np.pi/180 #Va is not just init, it will be updated. 
                    #Vrl = bus.Vm_init*np.cos(angle) #Vm is not just init, it will be updated.
                    #Vil = bus.Vm_init*np.sin(angle)

                    #Used the last v matrx to update Y
                    Vrl = v_last[v_index*2]
                    Vil = v_last[v_index*2+1]
                    Pl = load.P
                    Ql = load.Q  
                    
                    #Irl part of PQ
                    dIrl_dVrl = (Pl*(np.square(Vrl)+np.square(Vil))-(Pl*Vrl+1*Ql*Vil)*2*Vrl)/np.square(np.square(Vrl)+np.square(Vil))
                    dIrl_dVil = (Ql*(np.square(Vrl)+np.square(Vil))-(Pl*Vrl+1*Ql*Vil)*2*Vil)/np.square(np.square(Vrl)+np.square(Vil))

                    #Iil part of PQ
                    dIil_dVrl = (-1*Ql*(np.square(Vrl)+np.square(Vil))-(Pl*Vil-1*Ql*Vrl)*2*Vrl)/np.square(np.square(Vrl)+np.square(Vil))
                    dIil_dVil = (Pl*(np.square(Vrl)+np.square(Vil))-(Pl*Vil-1*Ql*Vrl)*2*Vil)/np.square(np.square(Vrl)+np.square(Vil))

                    #add PQ bus parameters
                    for i1 in range (0, size_Y_csc):
                    #add dIrl_dVrl to the Y matrix
                        if row_Y[i1] == (bus.Bus-1)*2:
                            if columun_Y[i1] == (bus.Bus-1)*2:
                                value_Y[i1]=value_Y[i1]+dIrl_dVrl
                                #print(row_Y[i1], columun_Y[i1])
                    #add dIrl_dVil to the Y matrix
                        if row_Y[i1] == (bus.Bus-1)*2:
                            if columun_Y[i1] == (bus.Bus-1)*2+1:
                                value_Y[i1]=value_Y[i1]+dIrl_dVil
                                #print(row_Y[i1], columun_Y[i1])
                    #add dIil_dVrl to the Y matrix
                        if row_Y[i1] == (bus.Bus-1)*2+1:
                            if columun_Y[i1] == (bus.Bus-1)*2:
                                value_Y[i1]=value_Y[i1]+dIil_dVrl
                                #print(row_Y[i1], columun_Y[i1])
                    #add dIil_dVil to the Y matrix
                        if row_Y[i1] == (bus.Bus-1)*2+1:
                            if columun_Y[i1] == (bus.Bus-1)*2+1:
                                value_Y[i1]=value_Y[i1]+dIil_dVil
                                #print(row_Y[i1], columun_Y[i1])
                    v_index = v_index+1
        
        gen_index = 0
        for gen in parsed_data['generators']:
            for bus in parsed_data['buses']:
                if bus.Bus == gen.Bus:
                    #abandoned the angle calculation because it creates "sign flips".
                    #angle = bus.Va_init*np.pi/180 #it's not just init, it will be updated. 
                    #Vrg = bus.Vm_init*np.cos(angle)
                    #Vig = bus.Vm_init*np.sin(angle)

                    #Used the last v matrx to update Y
                    Vrg = v_last[v_index*2]
                    Vig = v_last[v_index*2+1]
                    Pg = gen.P
                    Qg = gen.Qinit #it's not just init, it will be updated. 
                    #print (Pg, Qg)
                    #Irg part of PV
                    dIrg_dVrg = (-1*Pg*(np.square(Vrg)+np.square(Vig))-(-1*Pg*Vrg-1*Qg*Vig)*2*Vrg)/np.square(np.square(Vrg)+np.square(Vig))
                    dIrg_dVig = (-1*Qg*(np.square(Vrg)+np.square(Vig))-(-1*Pg*Vrg-1*Qg*Vig)*2*Vig)/np.square(np.square(Vrg)+np.square(Vig))
                    dIrg_dQg = -1*Vig/(np.square(Vrg)+np.square(Vig))

                    #Iirg part of PV
                    dIig_dVig = (-1*Pg*(np.square(Vrg)+np.square(Vig))-(-1*Pg*Vig+Qg*Vrg)*2*Vig)/np.square(np.square(Vrg)+np.square(Vig))
                    dIig_dVrg = (1*Qg*(np.square(Vrg)+np.square(Vig))-(-1*Pg*Vig+Qg*Vrg)*2*Vrg)/np.square(np.square(Vrg)+np.square(Vig))
                    dIig_dQg = 1*Vrg/(np.square(Vrg)+np.square(Vig))
                    
                    #add PV bus parameters
                    for i1 in range (0, size_Y_csc):
                        #add dIrg_dVrg to the Y matrix
                        if row_Y[i1] == (bus.Bus-1)*2:
                            if columun_Y[i1] == (bus.Bus-1)*2:
                                value_Y[i1]=value_Y[i1]+dIrg_dVrg
                                #print (row_Y[i1], columun_Y[i1], value_Y[i1])

                        #add dIrg_dVig to the Y matrix
                        if row_Y[i1] == (bus.Bus-1)*2:
                            if columun_Y[i1] == (bus.Bus-1)*2+1:
                                value_Y[i1]=value_Y[i1]+dIrg_dVig
                                #print (row_Y[i1], columun_Y[i1], value_Y[i1])
                    
                        #add dIrg_dQ to the Y matrix
                        if row_Y[i1] == (bus.Bus-1)*2: #(2 28)(4 29)(10 30)(14 31)
                            if columun_Y[i1] == (bus_number)*2+gen_index:
                                value_Y[i1]=value_Y[i1]+dIrg_dQg
                                #print (row_Y[i1], columun_Y[i1], value_Y[i1])

                        #add dIig_dVrg to the Y matrix
                        if row_Y[i1] == (bus.Bus-1)*2+1:
                            if columun_Y[i1] == (bus.Bus-1)*2:
                                value_Y[i1]=value_Y[i1]+dIig_dVrg
                                #print (row_Y[i1], columun_Y[i1], value_Y[i1])

                        #add dIig_dVig to the Y matrix
                        if row_Y[i1] == (bus.Bus-1)*2+1:
                            if columun_Y[i1] == (bus.Bus-1)*2+1:
                                value_Y[i1]=value_Y[i1]+dIig_dVig
                                #print (row_Y[i1], columun_Y[i1], value_Y[i1])

                        #add dIig_dQ to the Y matrix
                        if row_Y[i1] == (bus.Bus-1)*2+1: #(3 28)(5 29)(11 30)(15 31)
                            if columun_Y[i1] == (bus_number)*2+gen_index:
                                value_Y[i1]=value_Y[i1]+dIig_dQg
                                #print (row_Y[i1], columun_Y[i1], value_Y[i1])
                    
                        #add 2Vrg and 2Vig into the Y matrix
                        if row_Y[i1] == (bus_number)*2+gen_index:
                            if columun_Y[i1] == (bus.Bus-1)*2:
                                value_Y[i1]=value_Y[i1]+2*Vrg
                                #print (row_Y[i1], columun_Y[i1], value_Y[i1])

                        if row_Y[i1] == (bus_number)*2+gen_index:
                            if columun_Y[i1] == (bus.Bus-1)*2+1:
                                value_Y[i1]=value_Y[i1]+2*Vig
                                #print (row_Y[i1], columun_Y[i1], value_Y[i1])
                    
                    for load in parsed_data['loads']:
                        if (load.Bus == bus.Bus): 
                            Vrl = v_last[v_index*2]
                            Vil = v_last[v_index*2+1]
                            Pl = load.P
                            Ql = load.Q  
                            #print ("here!!!", Vrl, Vil, Pl, Ql)

                            #Irl part of PQ
                            dIrl_dVrl = (Pl*(np.square(Vrl)+np.square(Vil))-(Pl*Vrl+1*Ql*Vil)*2*Vrl)/np.square(np.square(Vrl)+np.square(Vil))
                            dIrl_dVil = (Ql*(np.square(Vrl)+np.square(Vil))-(Pl*Vrl+1*Ql*Vil)*2*Vil)/np.square(np.square(Vrl)+np.square(Vil))
                            Irl = (Pl*Vrl+Ql*Vil)/(np.square(Vrl)+np.square(Vil))

                            #Iil part of PQ
                            dIil_dVrl = (-1*Ql*(np.square(Vrl)+np.square(Vil))-(Pl*Vil-1*Ql*Vrl)*2*Vrl)/np.square(np.square(Vrl)+np.square(Vil))
                            dIil_dVil = (Pl*(np.square(Vrl)+np.square(Vil))-(Pl*Vil-1*Ql*Vrl)*2*Vil)/np.square(np.square(Vrl)+np.square(Vil))
                            Iil = (Pl*Vil-Ql*Vrl)/(np.square(Vrl)+np.square(Vil))

                            #gen bus also has P and Q connected to them, and P and Q need to be added into the matrix.
                            for i1 in range (0, size_Y_csc):
                                #add dIrg_dVrg to the Y matrix (6, 6)
                                if row_Y[i1] == (bus.Bus-1)*2:
                                    if columun_Y[i1] == (bus.Bus-1)*2:
                                        value_Y[i1]=value_Y[i1]+dIrl_dVrl
                                        #print ("here!!!", row_Y[i1], columun_Y[i1], value_Y[i1])
                                #add dIrg_dVig to the Y matrix (6, 7)
                                if row_Y[i1] == (bus.Bus-1)*2:
                                    if columun_Y[i1] == (bus.Bus-1)*2+1:
                                        value_Y[i1]=value_Y[i1]+dIrl_dVil
                                        #print ("here!!!", row_Y[i1], columun_Y[i1], value_Y[i1])
                                #add dIig_dVrg to the Y matrix (7, 6)
                                if row_Y[i1] == (bus.Bus-1)*2+1:
                                    if columun_Y[i1] == (bus.Bus-1)*2:
                                        value_Y[i1]=value_Y[i1]+dIil_dVrl
                                        #print ("here!!!", row_Y[i1], columun_Y[i1], value_Y[i1])
                                #add dIig_dVig to the Y matrix (7, 7)
                                if row_Y[i1] == (bus.Bus-1)*2+1:
                                    if columun_Y[i1] == (bus.Bus-1)*2+1:
                                        value_Y[i1]=value_Y[i1]+dIil_dVil
                                        #print ("here!!!", row_Y[i1], columun_Y[i1], value_Y[i1])
                    gen_index = gen_index +1
                                
        return (err_max, row_Y, columun_Y, value_Y)

    def run_powerflow(self,
                      v_init,
                      bus,
                      slack,
                      generator,
                      transformer,
                      branch,
                      shunt,
                      load,
                      row_Y, 
                      columun_Y, 
                      value_Y,
                      parsed_data,
                      tol_setting,
                      size_Y,
                      matrix_v, 
                      matrix_J,
                      size_Y_csc):
        """Runs a positive sequence power flow using the Equivalent Circuit Formulation.

        Args:
            v_init (np.array): The initial solution vector which has the same number of rows as the Y matrix.
            bus (list): Contains all the buses in the network as instances of the Buses class.
            slack (list): Contains all the slack generators in the network as instances of the Slack class.
            generator (list): Contains all the generators in the network as instances of the Generators class.
            transformer (list): Contains all the transformers in the network as instance of the Transformers class.
            branch (list): Contains all the branches in the network as instances of the Branches class.
            shunt (list): Contains all the shunts in the network as instances of the Shunts class.
            load (list): Contains all the loads in the network as instances of the Load class.

        Returns:
            v(np.array): The final solution vector.

        """

        # # # Copy v_init into the Solution Vectors used during NR, v, and the final solution vector v_sol # # #
        v = np.copy(v_init)
        v_sol = np.copy(v)
        v_last = np.copy(v_init)
       
        # # # Stamp Linear Power Grid Elements into Y matrix # # #
        # TODO: PART 1, STEP 2.1 - Complete the stamp_linear function which stamps all linear power grid elements.
        #  This function should call the stamp_linear function of each linear element and return an updated Y matrix.
        #  You need to decide the input arguments and return values.
        (row_Y, columun_Y, value_Y, matrix_J) = self.stamp_linear(v_init, row_Y, columun_Y, value_Y, parsed_data, size_Y, matrix_J)

        row_Y_linear = np.copy(row_Y)
        columun_Y_linear = np.copy(columun_Y)
        value_Y_linear = np.copy(value_Y)

        # # # Initialize While Loop (NR) Variables # # #
        # TODO: PART 1, STEP 2.2 - Initialize the NR variables
        err_max = 10  # maximum error at the current NR iteration. Guan-Ying: Set to 10 to enable the while loop.
        tol = tol_setting  # chosen NR tolerance
        NR_count = 1  # current NR iteration

        bus_amount = self.bus_number(parsed_data)

        # # # Begin Solving Via NR # # #
        # TODO: PART 1, STEP 2.3 - Complete the NR While Loop
        while err_max > tol:
            
            row_Y = np.copy(row_Y_linear)
            columun_Y = np.copy(columun_Y_linear)
            value_Y = np.copy(value_Y_linear)
            # # # Stamp Nonlinear Power Grid Elements into Y matrix # # #
            # TODO: PART 1, STEP 2.4 - Complete the stamp_nonlinear function which stamps all nonlinear power grid
            #  elements. This function should call the stamp_nonlinear function of each nonlinear element and return
            #  an updated Y matrix. You need to decide the input arguments and return values.
            (err_max, row_Y, columun_Y, value_Y) = self.stamp_nonlinear(err_max, v_init, row_Y, columun_Y, value_Y, parsed_data, v_last, size_Y_csc)
            
            matrix_Y = csc_matrix((value_Y, (row_Y, columun_Y)), shape=(size_Y, size_Y))
            matrix_Y_dense = matrix_Y.todense()
            #print('Y =', matrix_Y_dense)

            #for i1 in range (0, size_Y):
            #    for i in range(0,size_Y_csc-1):
            #        if columun_Y[i]==i1:
            #            print (columun_Y[i], row_Y[i], value_Y[i])

            #print ('J = ', matrix_J)


            # # # Solve The System # # #
            # TODO: PART 1, STEP 2.5 - Complete the solve function which solves system of equations Yv = J. The
            #  function should return a new v_sol.
            #  You need to decide the input arguments and return values.
            (matrix_v, matrix_J) = self.solver(matrix_Y, matrix_v, matrix_J, parsed_data)
            print ('NR iteration =', NR_count)
            print ('v = ', matrix_v)

            # # # Compute The Error at the current NR iteration # # #
            # TODO: PART 1, STEP 2.6 - Finish the check_error function which calculates the maximum error, err_max
            #  You need to decide the input arguments and return values.
            err_max = self.check_error(matrix_v, v_last)
            print ('err = ', err_max, '\n')
            
            #err_max = 0

            # # # Compute The Error at the current NR iteration # # #
            # TODO: PART 2, STEP 1 - Develop the apply_limiting function which implements voltage and reactive power
            #  limiting. Also, complete the else condition. Do not complete this step until you've finished Part 1.
            #  You need to decide the input arguments and return values.
            
            #This number defines the max. delta voltage (p.u.) in the v matrix. 
            v_delta_limit = 0.005
            if self.enable_limiting and err_max > tol:
                v_last = self.apply_limiting(v_last, matrix_v, v_delta_limit, bus_amount, size_Y)
            else:
                v_last = matrix_v
            

            v = np.copy(v_last)
            NR_count = NR_count+1 

        return v
