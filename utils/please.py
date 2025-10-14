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
    center_size = 0.3
    target_size = 0.3

    target_board_r = [{'radius':0, 'angle':0, 'sz':center_size, 'color':[0,0,0], 'x_offset':x_offset, 'y_offset':y_offset},
                  {'radius':radius, 'angle':0, 'sz':target_size, 'color':[0.5,0.5,0.5], 'x_offset':x_offset, 'y_offset':y_offset},
                  {'radius':radius, 'angle':60, 'sz':target_size, 'color':[0.5,0.5,0.5], 'x_offset':x_offset, 'y_offset':y_offset},
                  {'radius':radius, 'angle':120, 'sz':target_size, 'color':[0.5,0.5,0.5], 'x_offset':x_offset, 'y_offset':y_offset},
                  {'radius':radius, 'angle':180, 'sz':target_size, 'color':[0.5,0.5,0.5], 'x_offset':x_offset, 'y_offset':y_offset},
                  {'radius':radius, 'angle':240, 'sz':target_size, 'color':[0.5,0.5,0.5], 'x_offset':x_offset, 'y_offset':y_offset},
                  {'radius':radius, 'angle':300, 'sz':target_size, 'color':[0.5,0.5,0.5], 'x_offset':x_offset, 'y_offset':y_offset}]

    target_board_l = [{'radius':0, 'angle':0, 'sz':center_size, 'color':[0,0,0], 'x_offset':-x_offset, 'y_offset':y_offset},
                    {'radius':radius, 'angle':0, 'sz':target_size, 'color':[0.5,0.5,0.5], 'x_offset':-x_offset, 'y_offset':y_offset},
                    {'radius':radius, 'angle':60, 'sz':target_size, 'color':[0.5,0.5,0.5], 'x_offset':-x_offset, 'y_offset':y_offset},
                    {'radius':radius, 'angle':120, 'sz':target_size, 'color':[0.5,0.5,0.5], 'x_offset':-x_offset, 'y_offset':y_offset},
                    {'radius':radius, 'angle':180, 'sz':target_size, 'color':[0.5,0.5,0.5], 'x_offset':-x_offset, 'y_offset':y_offset},
                    {'radius':radius, 'angle':240, 'sz':target_size, 'color':[0.5,0.5,0.5], 'x_offset':-x_offset, 'y_offset':y_offset},
                    {'radius':radius, 'angle':300, 'sz':target_size, 'color':[0.5,0.5,0.5], 'x_offset':-x_offset, 'y_offset':y_offset}]

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

def matrix_to_vector(mat, include_diagonal=True):
    """Convert upper triangle of symmetric matrix to 1D vector"""
    triu_indices = np.triu_indices_from(mat, k=0 if include_diagonal else 1)
    return mat[triu_indices], triu_indices

def vector_to_matrix(vec, indices, N):
    """Reconstruct full symmetric matrix from vector and upper-triangle indices."""
    mat = np.zeros((N, N))
    r_idx, c_idx = indices
    mat[r_idx, c_idx] = vec
    mat[c_idx, r_idx] = vec  # ensure symmetry
    return mat

def double_center(mat):
    """Double center a square matrix (subtract row/col means, add grand mean)."""
    row_mean = np.mean(mat, axis=1, keepdims=True)
    col_mean = np.mean(mat, axis=0, keepdims=True)
    grand_mean = np.mean(mat)
    return mat - row_mean - col_mean + grand_mean

def indicator(conds, cols):
    '''
    Create an indicator matrix Z where each row corresponds to a condition in `conds`,
    and each column corresponds to a condition name in `cols`. An element Z[i, j] = 1
    if conds[i] == cols[j], otherwise 0.

    Parameters
    ----------
    conds : list or 1D numpy array of str
        A sequence of condition names (length n).

    cols : list or 1D array-like of str
        A list or array of unique condition names to be used as columns in the indicator matrix.

    Returns
    -------
    Z : numpy array of shape (n, len(cols))
        Indicator matrix with 1 where `conds[i] == cols[j]`, else 0.
    '''

    # Ensure conds and cols are numpy arrays
    conds = np.asarray(conds)
    cols = list(cols)  # convert cols to list if it's an ndarray or other iterable

    Z = np.zeros((len(conds), len(cols)), dtype=int)

    for i, cond in enumerate(conds):
        if cond in cols:
            j = cols.index(cond)
            Z[i, j] = 1
    
    return Z

def flatten_matrix(mat, colnames, prefix):
    """Returns a dict mapping 'prefix_colname' to matrix row"""
    flat = {}
    for i, col in enumerate(colnames):
        flat[f"{prefix}_{col}"] = mat[i]
    return flat

def retrieve_matrix(row, prefix, colnames):
    """
    Reconstructs an MxM matrix from a row of a DataFrame using a prefix and colnames.

    Parameters
    ----------
    row : pd.Series
        A single row (from df.loc[...] or df.iloc[...]).

    prefix : str
        Prefix used in column names (e.g., 'G', 'D').

    colnames : list of str
        The names of the matrix columns/rows.

    Returns
    -------
    np.ndarray
        The reconstructed matrix (M x M).
    """
    mat = np.vstack([row[f"{prefix}_{col}"].values[0].flatten() for col in colnames])
    return np.array(mat)

def cka(G1, G2):
    '''
    Compute the Centered Kernel Alignment (CKA) between two kernel (similarity) matrices G1 and G2.
    '''
    # Center the matrices
    n = G1.shape[0]
    H = np.eye(n) - np.ones((n, n)) / n
    G1_centered = H @ G1 @ H
    G2_centered = H @ G2 @ H

    # Compute the CKA
    cka = np.sum(G1_centered * G2_centered) / (np.sqrt(np.sum(G1_centered ** 2)) * np.sqrt(np.sum(G2_centered ** 2)))
    return cka