from pickle import NONE
import numpy as np

def csc_size (parsed_data):
    size_Y_csc = 0
    for bus in parsed_data['buses']:
        size_Y_csc += 12
    for slack in parsed_data['slack']:
        size_Y_csc += 4
    for gen in parsed_data['generators']:
        size_Y_csc += 4
    for trans in parsed_data['xfmrs']:
        size_Y_csc += 10 #this is my guess. Don't know why.
    return (size_Y_csc)