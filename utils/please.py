'''
Bunch of helper functions to make my life easier. It's always 
nice to say please when asking for help.
'''
import os

import numpy as np
import pandas as pd
from collections import defaultdict

import matplotlib.patches as patches

def digits_to_int(num_list: list[int]) -> int:
    """
    Convert a list of numbers into a single integer.
    
    Args:
    List of numbers to be concatenated.
    
    Returns:
    The concatenated integer.
    """
    # Convert each number to a string and concatenate them
    concatenated_str = ''.join(map(str, num_list))
    
    # Convert the concatenated string back to an integer
    return int(concatenated_str)

def int_to_digits(num: int) -> list[int]:
    """
    Convert an integer into a list of numbers.
    
    Args:
    The integer to be converted.
    
    Returns:
    The list of numbers.
    """
    # Convert the integer to a string
    num_str = str(num)
    
    # Convert each character to an integer and store in a list
    return [int(char) for char in num_str]

def moving_average(data, window_size):
    """
    Compute the moving average along the first axis of an N by K input list or np.ndarray.

    Parameters:
    data (list or np.ndarray): Input data, an N by K list or ndarray.
    window_size (int): Size of the moving window.

    Returns:
    np.ndarray: Array of moving averages with shape (N - window_size + 1, K).
    """
    if isinstance(data, list):
        data = np.array(data)
    
    if window_size < 1:
        raise ValueError("Window size must be at least 1.")
    if window_size > data.shape[0]:
        raise ValueError("Window size must be less than or equal to the length of the data along the first axis.")
    
    # Create a window of ones, normalized by the window size
    window = np.ones(window_size) / window_size
    
    # Apply the moving average along the first axis for each column
    moving_avg = np.apply_along_axis(lambda m: np.convolve(m, window, mode='valid'), axis=0, arr=data)
    
    # Pad the result to keep the same shape as the input
    pad_width = (window_size - 1) // 2
    if window_size % 2 == 0:
        pad_width_end = pad_width + 1
    else:
        pad_width_end = pad_width
    
    moving_avg_padded = np.pad(moving_avg, ((pad_width, pad_width_end), (0, 0)), mode='edge')
    
    return moving_avg_padded

def find_closest_index(vector, value):
    """
    Find the index of the closest value in a vector to a given float value.

    Parameters:
    vector (np.ndarray): Input vector.
    value (float): The value to find the closest index to.

    Returns:
    int: The index of the closest value in the vector.
    """
    vector = np.asarray(vector)  # Ensure the input is a NumPy array
    index = np.argmin(np.abs(vector - value))
    return index

def draw_board(ax, radius=5):
    """
        Draw the experiment display with circle targets on the given axes.
    """
    x_offset = radius+2
    y_offset = 0
    
    target_board_r = [{'radius':0, 'angle':0, 'sz':0.1, 'color':[0,0,0], 'x_offset':x_offset, 'y_offset':y_offset},
                  {'radius':radius, 'angle':0, 'sz':0.15, 'color':[0.5,0.5,0.5], 'x_offset':x_offset, 'y_offset':y_offset},
                  {'radius':radius, 'angle':60, 'sz':0.15, 'color':[0.5,0.5,0.5], 'x_offset':x_offset, 'y_offset':y_offset},
                  {'radius':radius, 'angle':120, 'sz':0.15, 'color':[0.5,0.5,0.5], 'x_offset':x_offset, 'y_offset':y_offset},
                  {'radius':radius, 'angle':180, 'sz':0.15, 'color':[0.5,0.5,0.5], 'x_offset':x_offset, 'y_offset':y_offset},
                  {'radius':radius, 'angle':240, 'sz':0.15, 'color':[0.5,0.5,0.5], 'x_offset':x_offset, 'y_offset':y_offset},
                  {'radius':radius, 'angle':300, 'sz':0.15, 'color':[0.5,0.5,0.5], 'x_offset':x_offset, 'y_offset':y_offset}]

    target_board_l = [{'radius':0, 'angle':0, 'sz':0.1, 'color':[0,0,0], 'x_offset':-x_offset, 'y_offset':y_offset},
                    {'radius':radius, 'angle':0, 'sz':0.15, 'color':[0.5,0.5,0.5], 'x_offset':-x_offset, 'y_offset':y_offset},
                    {'radius':radius, 'angle':60, 'sz':0.15, 'color':[0.5,0.5,0.5], 'x_offset':-x_offset, 'y_offset':y_offset},
                    {'radius':radius, 'angle':120, 'sz':0.15, 'color':[0.5,0.5,0.5], 'x_offset':-x_offset, 'y_offset':y_offset},
                    {'radius':radius, 'angle':180, 'sz':0.15, 'color':[0.5,0.5,0.5], 'x_offset':-x_offset, 'y_offset':y_offset},
                    {'radius':radius, 'angle':240, 'sz':0.15, 'color':[0.5,0.5,0.5], 'x_offset':-x_offset, 'y_offset':y_offset},
                    {'radius':radius, 'angle':300, 'sz':0.15, 'color':[0.5,0.5,0.5], 'x_offset':-x_offset, 'y_offset':y_offset}]

    # plot circles at target locs:
    for target in target_board_r:
        circle = patches.Circle((target['radius']*np.cos(np.deg2rad(target['angle'])) + target['x_offset'],
                                target['radius']*np.sin(np.deg2rad(target['angle'])) + target['y_offset']),
                                radius=target['sz'], color=target['color'])
        ax.add_patch(circle)
    # plot circles at target locs:
    for target in target_board_l:
        circle = patches.Circle((target['radius']*np.cos(np.deg2rad(target['angle'])) + target['x_offset'],
                                target['radius']*np.sin(np.deg2rad(target['angle'])) + target['y_offset']),
                                radius=target['sz'], color=target['color'])
        ax.add_patch(circle)

def make_it_pretty(ax):
    # Make it pretty:
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(1)
    ax.spines['bottom'].set_linewidth(1)

    ax.spines["left"].set_bounds(ax.get_ylim()[0], ax.get_ylim()[-1])
    ax.spines["bottom"].set_bounds(ax.get_xticks()[0], ax.get_xticks()[-1])
