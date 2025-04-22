import numpy as np

def get_RT(mov: np.ndarray, baseline_threshold: float, fGain: list[float], global_gain: float):
    '''
    RT is defined as the time from the beginning of the execution period (when the targets signal the go-cue)
    to the time that any finger exits the baseline zone (baseline_threshold) for the first time.

    params:
        mov: the mov data for a single trial
        baseline_threshold: the threshold for the baseline zone
        fGain: the gain for force of each finger
        global_gain: the global gain for all finger forces
    '''

    WAIT_EXEC = 3

    # find the beginning of the execution period:
    start_idx = np.where(mov[:, 0] == WAIT_EXEC)[0][0]
    end_idx = np.where(mov[:, 0] == WAIT_EXEC)[0][-1]

    # get the differential forces - five columns:
    force = mov[:, 13:18]
    
    # apply the gains:
    force = force * fGain * global_gain

    # find the first time any finger exits the baseline zone:
    for i in range(start_idx, end_idx+1):
        if np.any(np.abs(force[i, :]) > baseline_threshold):
            RT_idx = i
            break
    
    return mov[RT_idx, 2] - mov[start_idx, 2]

def get_ET(mov: np.ndarray):
    '''
    ET is defined as the time from the the beginning of the execution period (when the targets signal the go-cue)
    to the beginning of the holding period (The first time when the fingers are in the correct position).

    params:
        mov: 2D the mov data for a single trial
    '''

    WAIT_EXEC = 3

    # find the beginning of the execution period:
    start_idx = np.where(mov[:, 0] == WAIT_EXEC)[0][0]
    end_idx = np.where(mov[:, 0] == WAIT_EXEC)[0][-1]

    return mov[end_idx, 2] - mov[start_idx, 2] - 600

def get_MD(f, straight_trajectory):
    '''
    Mean deviation captures the simultaneity of the forces of the fingers in trial. It is defined as the average of 
    the norm of the deviation of the forces from the ideal trajectory.

    sum from t=1 to T of: norm(F_t - (C' * F_t)/norm(C)^2 . C) / T
    where F_t is the force of the 5 fingers at each time point t and C 
    is the ideal trajectory that you can take to reach the target
    position (Here it is assumed that the ideal trajetory is a straight 
    line from the starting position to the ending position).

    Look Waters-Metenier et al. 2014 for more details.

    params:
        mov: the mov data for a single trial
        baseline_threshold: the threshold for the baseline zone
        fGain: the gain for force of each finger
        global_gain: the global gain for all finger forces
        fs: the sampling frequency of the force data
        hold_time: the time to hold the target in ms
    '''
    c = straight_trajectory
    
    # calculate mean deviation:
    deviation = []
    for i in range(1,f.shape[0]):
        # force vector:
        tmp = f[i, :] - f[0, :]

        # projection:
        projection = np.dot(tmp, c) / np.dot(c, c) * c

        # deviation:
        deviation.append(np.linalg.norm(tmp - projection))

    return np.mean(deviation)

