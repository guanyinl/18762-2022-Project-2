from pickle import NONE
import numpy as np

def matrix_size (parsed_data, size_Y):
    size_Y_csc = 0
    for bus in parsed_data['buses']:
        size_Y_csc += 12
    for slack in parsed_data['slack']:
        size_Y_csc += 4
    for gen in parsed_data['generators']:
        size_Y_csc += 4
    for trans in parsed_data['xfmrs']:
        size_Y_csc += 32
        size_Y = size_Y+4 

    return (size_Y_csc, size_Y)