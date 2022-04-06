from curses import raw
from multiprocessing.sharedctypes import Value
from pickle import NONE
import time as time
from numpy import matrix
from parsers.parser import parse_raw
from scripts.PowerFlow import PowerFlow
from scripts.process_results import process_results
from scripts.initialize import initialize
from scripts.csc_size import csc_size
from models.Buses import Buses


def solve(TESTCASE, SETTINGS):
    print ("time.time(): %f " %  time.time())
    print ("Simulation starts at = ", time.asctime( time.localtime(time.time())))
    """Run the power flow solver.

    Args:
        TESTCASE (str): A string with the path to the test case.
        SETTINGS (dict): Contains all the solver settings in a dictionary.

    Returns:
        None
    """
    # TODO: PART 1, STEP 0 - Initialize all the model classes in the models directory (models/) and familiarize
    #  yourself with the parameters of each model. Use the docs/DataFormats.pdf for assistance.

    # # # Parse the Test Case Data # # #
    case_name = TESTCASE
    parsed_data = parse_raw(case_name)

    # # # Assign Parsed Data to Variables # # #
    bus = parsed_data['buses']
    slack = parsed_data['slack']
    generator = parsed_data['generators']
    transformer = parsed_data['xfmrs']
    branch = parsed_data['branches']
    shunt = parsed_data['shunts']
    load = parsed_data['loads']

    # # # Solver Settings # # #
    tol_setting = SETTINGS['Tolerance']  # NR solver tolerance
    max_iters = SETTINGS['Max Iters']  # maximum NR iterations
    enable_limiting = SETTINGS['Limiting']  # enable/disable voltage and reactive power limiting

    # # # Assign System Nodes Bus by Bus # # #
    # We can use these nodes to have predetermined node number for every node in our Y matrix and J vector.
    for ele in bus:
        ele.assign_nodes()

    # Assign any slack nodes
    for ele in slack:
        ele.assign_nodes()

    # # # Initialize Solution Vector - V and Q values # # #

    # determine the size of the Y matrix by looking at the total number of nodes in the system
    size_Y = Buses._node_index.__next__()
    #print (size_Y)

    # TODO: PART 1, STEP 1 - Complete the function to initialize your solution vector v_init.
    v_init = None  # create a solution vector filled with zeros of size_Y
    v_init = initialize(parsed_data, case_name, size_Y)
    
    # # # Run Power Flow # # #
    powerflow = PowerFlow(case_name, tol_setting, max_iters, enable_limiting)

    # TODO: PART 1, STEP 2 - Complete the PowerFlow class and build your run_powerflow function to solve Equivalent
    #  Circuit Formulation powerflow. The function will return a final solution vector v. Remove run_pf and the if
    #  condition once you've finished building your solver.
    
    size_Y_csc = csc_size(parsed_data)
    #print (size_Y_csc)

    row_Y = [0]*size_Y_csc
    columun_Y = [0]*size_Y_csc
    value_Y = [0]*size_Y_csc

    matrix_J = [0]*size_Y
    matrix_v = [0]*size_Y

    run_pf = True
    if run_pf:
        v = powerflow.run_powerflow(v_init, bus, slack, generator, transformer, branch, shunt, load, 
                                    row_Y, columun_Y, value_Y, parsed_data, tol_setting, size_Y, 
                                    matrix_v, matrix_J)

    # # # Process Results # # #
    # TODO: PART 1, STEP 3 - Write a process_results function to compute the relevant results (voltages, powers,
    #  and anything else of interest) and find the voltage profile (maximum and minimum voltages in the case).
    #  You can decide which arguments to pass to this function yourself.
    process_results(v, bus)
    print ("time.time(): %f " %  time.time())
    print ("Simulation stops at = ", time.asctime( time.localtime(time.time())))
