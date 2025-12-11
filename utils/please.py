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


def draw_sig_lines(ax, df, x_col, y_col, pairs, plot_type='boxplot', line_height_gap=2, line_height_increase=3, height_tick=2):
    """
    Draws significance lines between categories on a plot.
    Works for boxplot and barplot.

    Args:
        ax: The matplotlib axes object.
        df: The pandas DataFrame used for plotting.
        x_col: The name of the column for the x-axis categories.
        y_col: The name of the column for the y-axis values.
        pairs: A list of pairs of category indices to connect.
        plot_type: The type of plot ('boxplot' or 'barplot').
        line_height_gap: Gap above the highest point in a pair.
        line_height_increase: Vertical distance between stacked lines.
    """
    # Get category names in the order they are plotted
    cats = df[x_col].unique()
    
    # Keep track of the y-level for lines to avoid overlaps
    line_levels = []

    for p in pairs:
        cat1_idx, cat2_idx = p
        cat1_name, cat2_name = cats[cat1_idx], cats[cat2_idx]

        # Get x-coordinates for the centers of the categories
        x1, x2 = cat1_idx, cat2_idx

        # Determine the y-level for the line
        if plot_type == 'boxplot':
            # For boxplots, find the max value in the data for the two categories
            max_val = df[df[x_col].isin([cat1_name, cat2_name])][y_col].max()
        elif plot_type == 'barplot':
            # For barplots, get the height of the bars
            bar_heights = [b.get_height() for b in ax.patches]
            max_val = max(bar_heights[cat1_idx], bar_heights[cat2_idx])
        else:
            raise ValueError("plot_type must be 'boxplot' or 'barplot'")

        y = max_val + line_height_gap
        while any(abs(y - level) < line_height_increase for level in line_levels):
            y += line_height_increase
        line_levels.append(y)

        # Draw the horizontal line and vertical ticks
        ax.plot([x1, x1, x2, x2], [y - (height_tick), y, y, y - (height_tick)], color='k', linewidth=1)

def analyze_r(D, cond_order, specify_conditions=None, within_corr_conditions=None):
    """
    gets the pearson r between the first and second half of the conditions. 
    
    Params:
    D (list): List of PCM Datasets
    cond_order (list): The order of conditions such as [0,1,2,3,4,5,6,7,8,9,10,11] must be passed.
    specify_conditions (np.array): array of condition indices to only include in the analysis. If None, all conditions are used. 
                                   the indices are based on the first half of the conditions (e.g., [0,1,2,3,4,5])

    Returns:
    r (np.array): list of pearson r values for each subject
    """
    N = len(D) # number of subjects
    nParts = len(np.unique(D[0].obs_descriptors['part_vec']))
    r = np.zeros(N)
    r_within_contra = np.zeros(N)
    r_within_ipsi = np.zeros(N)
    for i in range(N):
        measurements = D[i].measurements
        cond_vec = np.array(cond_order * nParts)
        num_voxels = measurements.shape[1]
        ncond = 12

        # Create an array to store the averaged patterns for each condition
        averaged_patterns = np.zeros((ncond, num_voxels))

        if specify_conditions is not None:
            mask = np.array(specify_conditions)
        else:
            mask = np.arange(6)

        # Loop through each condition and calculate the average pattern
        for c in range(ncond):
            # Find rows corresponding to the current condition
            condition_indices = cond_vec == c
            # Calculate the mean pattern for the current condition
            averaged_patterns[c, :] = measurements[condition_indices, :].mean(axis=0)
        
        # get the contra conditions: 
        y_contra_avg = averaged_patterns[0:6, :]
        y_contra_avg = y_contra_avg[mask, :]
        # remove mean across conditions:
        y_contra_avg = y_contra_avg - y_contra_avg.mean(axis=0)
        # flatten the averaged patterns:
        y_contra_vec = y_contra_avg.flatten()
        
        # get the ipsi conditions:
        y_ipsi_avg = averaged_patterns[6:12, :]
        y_ipsi_avg = y_ipsi_avg[mask, :]
        # remove mean across conditions:
        y_ipsi_avg = y_ipsi_avg - y_ipsi_avg.mean(axis=0)
        # flatten the averaged patterns:
        y_ipsi_vec = y_ipsi_avg.flatten()

        if within_corr_conditions is not None:
            y1 = y_contra_avg[np.array([]), :]
            y2 = y_contra_avg[within_corr_conditions, :]
            y1_vec = y1.flatten()
            y2_vec = y2.flatten()
            r_within_contra[i] = np.corrcoef(y1_vec, y2_vec)[0,1]

            y1 = y_ipsi_avg[within_corr_conditions, :]
            y2 = y_ipsi_avg[within_corr_conditions, :]
            y1_vec = y1.flatten()
            y2_vec = y2.flatten()
            r_within_ipsi[i] = np.corrcoef(y1_vec, y2_vec)[0,1]
            
        # concatenate the contra and ipsi patterns:
        r[i] = np.corrcoef(y_contra_vec, y_ipsi_vec)[0,1]

    if within_corr_conditions is not None:
        return r, r_within_contra, r_within_ipsi
    return r

