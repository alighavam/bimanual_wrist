import os
import numpy as np
import pandas as pd
from tqdm import tqdm
from .please import *
from .measures import *

def movload(fname):
    # loads .mov files given the path of the file. The .mov files have a specific custom hence the need for a custom function
    A = []
    fid = open(fname, 'rt')
    if fid == -1:
        raise Exception('Could not open ' + fname)

    trial = 0
    for line in fid:
        if line[0] == 'T':
            #print('Trial: ', line.split()[1])
            a = int(line.split()[1])
            trial += 1
            if a != trial:
                print('Trials out of sequence')
                trial = a
            A.append(None)  # placeholder for now
        else:
            lineData = line.strip().split('\t')
            a = np.array([float(x) for x in lineData], ndmin=2)
            if A[trial - 1] is None:
                A[trial - 1] = np.empty((0, a.shape[1]))  # infer shape dynamically

            A[trial - 1] = np.vstack((A[trial - 1], a))

    fid.close()
    return A

def trial_routine(row):
    '''
    This function performs all the necessary preprocessing for a single trial
    '''
    C = pd.DataFrame()

    # add row to the dataframe:
    C = pd.concat([C, row], ignore_index=True)
    # C = C.rename(columns={'trialCorr': 'trial_correct'})

    return C
    
def subject_routine_train(sn, path, smooth_win_sz=0, fs=200):
    """
    This function is used to preprocess the behavioural data of a subject
    
    params:
        sn: int, the subject number
        path: dict, the path to directories to access data
        smooth_win_sz: int, the window size for the moving average filter for forces/kinematics
        fs: int, the sampling frequency of the force/kinematics data
    """

    # ===========================================================================
    #  HANDLING THE TRAINING BEHAVIOURAL DATA:
    # ===========================================================================
    # empty dataframe to store the data:
    D = pd.DataFrame()
    D_mov = None

    # states on mov file:
    # WAIT_TRIAL,				///< 0 Trial is not started yet, is default mode
	# START_TRIAL,			///< 1 Start trial	
	# WAIT_TR,				///< 2 Wait for TR pulse to begin the trial	
	# WAIT_CENTER_HOLD,		///< 3 Wait for the hand to be start position 	
	# SHOW_TARGET,			///< 4 Show target but don't move yet
	# GO_CUE,					///< 5 Show go signal
	# MOVING,					///< 6 During the movement 
	# GIVE_FEEDBACK,			///< 7 Show feedback 
	# ACQUIRE_HRF,			///< 8 Non mandatory state, for use in scanner to hold recording
	# END_TRIAL,				///< 9 Absorbant State: Trial is finished
    
    tmp_left = {
        'cond': ['left_0', 'left_60', 'left_120', 'left_180', 'left_240', 'left_300'],
        'reach_type': ['unimanual', 'unimanual', 'unimanual', 'unimanual', 'unimanual', 'unimanual'],
        'Uni_or_Bi': [0, 0, 0, 0, 0, 0],
        'Hand': [0, 0, 0, 0, 0, 0],
        'targetAngle_L': [0, 60, 120, 180, 240, 300],
        'targetAngle_R': [-1, -1, -1, -1, -1, -1],
    }

    tmp_right = {
        'cond': ['right_0', 'right_60', 'right_120', 'right_180', 'right_240', 'right_300'],
        'reach_type': ['unimanual', 'unimanual', 'unimanual', 'unimanual', 'unimanual', 'unimanual'],
        'Uni_or_Bi': [0, 0, 0, 0, 0, 0],
        'Hand': [1, 1, 1, 1, 1, 1],
        'targetAngle_L': [-1, -1, -1, -1, -1, -1],
        'targetAngle_R': [0, 60, 120, 180, 240, 300],
    }

    tmp_bimanual = {
        'cond': ['bimanual_0_0', 'bimanual_0_60', 'bimanual_0_120', 'bimanual_0_180', 'bimanual_0_240', 'bimanual_0_300',
                'bimanual_60_0', 'bimanual_60_60', 'bimanual_60_120', 'bimanual_60_180', 'bimanual_60_240', 'bimanual_60_300',
                'bimanual_120_0', 'bimanual_120_60', 'bimanual_120_120', 'bimanual_120_180', 'bimanual_120_240', 'bimanual_120_300',
                'bimanual_180_0', 'bimanual_180_60', 'bimanual_180_120', 'bimanual_180_180', 'bimanual_180_240', 'bimanual_180_300',
                'bimanual_240_0', 'bimanual_240_60', 'bimanual_240_120', 'bimanual_240_180', 'bimanual_240_240', 'bimanual_240_300',
                'bimanual_300_0', 'bimanual_300_60', 'bimanual_300_120', 'bimanual_300_180', 'bimanual_300_240', 'bimanual_300_300'],
        'reach_type': ['matched', 'unrelated', 'unrelated', 'mirror', 'unrelated', 'unrelated',
                    'unrelated', 'matched', 'mirror', 'unrelated', 'mirror-diagonal', 'unrelated',
                    'unrelated', 'mirror', 'matched', 'unrelated', 'unrelated', 'mirror-diagonal',
                    'mirror', 'unrelated', 'unrelated', 'matched', 'unrelated', 'unrelated',
                    'unrelated', 'mirror-diagonal', 'unrelated', 'unrelated', 'matched', 'mirror',
                    'unrelated', 'unrelated', 'mirror-diagonal', 'unrelated', 'mirror', 'matched'],
        'Uni_or_Bi': [1, 1, 1, 1, 1, 1,
                    1, 1, 1, 1, 1, 1,
                    1, 1, 1, 1, 1, 1,
                    1, 1, 1, 1, 1, 1,
                    1, 1, 1, 1, 1, 1,
                    1, 1, 1, 1, 1, 1],
        'Hand': [-1, -1, -1, -1, -1, -1,
                -1, -1, -1, -1, -1, -1,
                -1, -1, -1, -1, -1, -1,
                -1, -1, -1, -1, -1, -1,
                -1, -1, -1, -1, -1, -1,
                -1, -1, -1, -1, -1, -1],
        'targetAngle_L': [0, 0, 0, 0, 0, 0,
                        60, 60, 60, 60, 60, 60,
                        120, 120, 120, 120, 120, 120,
                        180, 180, 180, 180, 180, 180,
                        240, 240, 240, 240, 240, 240,
                        300, 300, 300, 300, 300, 300],
        'targetAngle_R': [0, 60, 120, 180, 240, 300,
                        0, 60, 120, 180, 240, 300,
                        0, 60, 120, 180, 240, 300,
                        0, 60, 120, 180, 240, 300,
                        0, 60, 120, 180, 240, 300,
                        0, 60, 120, 180, 240, 300],
    }
    
    # combine the dictionaries:
    combined = {k: tmp_left[k] + tmp_right[k] + tmp_bimanual[k] for k in tmp_left}
    conds = pd.DataFrame(combined)
    
    # Load the training .dat file:
    dat_file_name = os.path.join(path['train_behavDir'], f's{sn}', f's{sn}_train.dat')
    dat = pd.read_table(dat_file_name)

    fMRI_sess = []
    idx_preRT = []
    idx_postRT = []
    idx_gocue = []
    idx_endReach = []
    MD_left = []
    MD_right = []

    oldblock = -1
    # loop on trials:
    for i in tqdm(range(dat.shape[0])):
        if dat['BN'][i] != oldblock:
            # print(f'Processing block {dat["BN"][i]}')
            # load the .mov file:
            ext = int(dat['BN'][i])
            mov = movload(os.path.join(path['train_behavDir'], f's{sn}', f's{sn}_train_{ext:02d}.mov'))
            oldblock = dat['BN'][i]
        # print(f'Processing trial {dat["TN"][i]}')
        # trial routine:
        C = trial_routine(dat.iloc[[i]])

        # append the trial to the dataframe:
        D = pd.concat([D, C], ignore_index=True)
        
        # Whether data is from fMRI or training:
        fMRI_sess.append(int(0))         

        # movemean the kinematics:
        trial_mov = mov[dat['TN'][i]-1]
        if smooth_win_sz > 0:
            # smooth the radius and angle signals:
            trial_mov[:, 5:9] = moving_average(trial_mov[:, 5:9], smooth_win_sz)
        # add the mov trial in the move dataframe:
        tmp = pd.DataFrame({'fMRI_sess': [0], 
                            'sn': [sn], 
                            'BN': [dat['BN'][i]], 
                            'TN': [dat['TN'][i]], 
                            'GoodMovement': [dat['GoodMovement'][i]], 
                            'mov': pd.Series([trial_mov], dtype='object')})
        
        # find the index of the go cue and end of reaching:
        if C.GoodMovement[0]:
            state = trial_mov[:,2]
            idx_gocue.append(np.where(state==5)[0][0])
            
            dEndRadius_L = C.dEndRadius_L.values[0]
            dEndRadius_R = C.dEndRadius_R.values[0]

            if dEndRadius_L > 0:
                tmp_idx = find_closest_index(trial_mov[:, 5], dEndRadius_L)
            
            if dEndRadius_R > 0:
                tmp_idx = find_closest_index(trial_mov[:, 6], dEndRadius_R)
            
            idx_endReach.append(tmp_idx)
        else:
            idx_gocue.append(0)
            idx_endReach.append(0)

        # Find pre_RT and post_RT indices:
        # MAY DEVELOP LATER.

        # Calculate mean deviation for the trial:
        if C.GoodMovement[0]:
            radius_l = trial_mov[:,5].flatten()[idx_gocue[i]:idx_endReach[i]]
            radius_r = trial_mov[:,6].flatten()[idx_gocue[i]:idx_endReach[i]]
            angle_l = trial_mov[:,7].flatten()[idx_gocue[i]:idx_endReach[i]]
            angle_r = trial_mov[:,8].flatten()[idx_gocue[i]:idx_endReach[i]]

            x_l = radius_l * np.cos(np.deg2rad(angle_l))
            y_l = radius_l * np.sin(np.deg2rad(angle_l))
            x_r = radius_r * np.cos(np.deg2rad(angle_r))
            y_r = radius_r * np.sin(np.deg2rad(angle_r))

            # left hand coords:
            f_left = np.vstack((x_l, y_l)).T    # position of left hand
            c_left = f_left[-1] - f_left[0]     # straight trajectory defined as position of go-cue and end of reach
            MD_left.append(get_MD(f_left, c_left))
            
            # right hand coords:
            f_right = np.vstack((x_r, y_r)).T
            c_right = f_right[-1] - f_right[0]
            MD_right.append(get_MD(f_right, c_right))
        else:
            MD_left.append(-1)
            MD_right.append(-1)

        # flatten the mov dataframe:
        flat_rows = []
        for _, row in tmp.iterrows():
            mat = row['mov']  # shape (n, m)
            meta = row.drop('mov').to_dict()  # get all metadata except 'mov'

            for r in mat:
                combined_row = meta.copy()
                for j, val in enumerate(r):
                    combined_row[f'mov_{j}'] = val
                flat_rows.append(combined_row)

        # Turn into flat DataFrame
        if D_mov is None:
            D_mov = pd.DataFrame(flat_rows)
        else:
            D_mov = pd.concat([D_mov, pd.DataFrame(flat_rows)], ignore_index=True)

    D['fMRI_sess'] = fMRI_sess
    D['idx_gocue'] = idx_gocue
    D['idx_endReach'] = idx_endReach
    D['MD_left'] = MD_left
    D['MD_right'] = MD_right

    cond_name = np.empty(len(D), dtype=object)
    reach_type = np.zeros(len(D), dtype=object)
    for index, row in conds.iterrows():
        if row.Uni_or_Bi == 0:  # unimanual
            if row.Hand == 0:   # left hand
                rows = ((D.Uni_or_Bi == 0) & (D.Hand == 0) & (D.targetAngle_L == row.targetAngle_L)).values.flatten()
                cond_name[rows] = row.cond
                reach_type[rows] = row.reach_type
            else:   # right hand
                rows = ((D.Uni_or_Bi == 0) & (D.Hand == 1) & (D.targetAngle_R == row.targetAngle_R)).values.flatten()
                cond_name[rows] = row.cond
                reach_type[rows] = row.reach_type
        else:  # bimanual
            rows = ((D.Uni_or_Bi == 1) & (D.targetAngle_L == row.targetAngle_L) & (D.targetAngle_R == row.targetAngle_R)).values.flatten()
            cond_name[rows] = row.cond
            reach_type[rows] = row.reach_type
    D['cond_name'] = cond_name
    D['reach_type'] = reach_type
    D['cond_name'] = D['cond_name'].astype(str)
    D['reach_type'] = D['reach_type'].astype(str)
    
    # save the data frames:
    D.to_csv(os.path.join(path['anaDir'], f's{sn}_train.csv'), index=False)    
    D_mov.to_csv(os.path.join(path['anaDir'], f's{sn}_train_mov.csv'), index=False)

def subject_routine_fMRI(sn, path, smooth_win_sz=0, fs=200):
    # ===========================================================================
    #  HANDLING THE fMRI BEHAVIOURAL DATA:
    # ===========================================================================
    # empty dataframe to store the data:
    D = pd.DataFrame()
    D_mov = None

    # states on mov file:
    # WAIT_TRIAL,				///< 0 Trial is not started yet, is default mode
	# START_TRIAL,			///< 1 Start trial	
	# WAIT_TR,				///< 2 Wait for TR pulse to begin the trial	
	# WAIT_CENTER_HOLD,		///< 3 Wait for the hand to be start position 	
	# SHOW_TARGET,			///< 4 Show target but don't move yet
	# GO_CUE,					///< 5 Show go signal
	# MOVING,					///< 6 During the movement 
	# GIVE_FEEDBACK,			///< 7 Show feedback 
	# ACQUIRE_HRF,			///< 8 Non mandatory state, for use in scanner to hold recording
	# END_TRIAL,				///< 9 Absorbant State: Trial is finished
    
    tmp_left = {
        'cond': ['left_0', 'left_60', 'left_120', 'left_180', 'left_240', 'left_300'],
        'reach_type': ['unimanual', 'unimanual', 'unimanual', 'unimanual', 'unimanual', 'unimanual'],
        'Uni_or_Bi': [0, 0, 0, 0, 0, 0],
        'Hand': [0, 0, 0, 0, 0, 0],
        'targetAngle_L': [0, 60, 120, 180, 240, 300],
        'targetAngle_R': [-1, -1, -1, -1, -1, -1],
    }

    tmp_right = {
        'cond': ['right_0', 'right_60', 'right_120', 'right_180', 'right_240', 'right_300'],
        'reach_type': ['unimanual', 'unimanual', 'unimanual', 'unimanual', 'unimanual', 'unimanual'],
        'Uni_or_Bi': [0, 0, 0, 0, 0, 0],
        'Hand': [1, 1, 1, 1, 1, 1],
        'targetAngle_L': [-1, -1, -1, -1, -1, -1],
        'targetAngle_R': [0, 60, 120, 180, 240, 300],
    }

    tmp_bimanual = {
        'cond': ['bimanual_0_0', 'bimanual_0_60', 'bimanual_0_120', 'bimanual_0_180', 'bimanual_0_240', 'bimanual_0_300',
                'bimanual_60_0', 'bimanual_60_60', 'bimanual_60_120', 'bimanual_60_180', 'bimanual_60_240', 'bimanual_60_300',
                'bimanual_120_0', 'bimanual_120_60', 'bimanual_120_120', 'bimanual_120_180', 'bimanual_120_240', 'bimanual_120_300',
                'bimanual_180_0', 'bimanual_180_60', 'bimanual_180_120', 'bimanual_180_180', 'bimanual_180_240', 'bimanual_180_300',
                'bimanual_240_0', 'bimanual_240_60', 'bimanual_240_120', 'bimanual_240_180', 'bimanual_240_240', 'bimanual_240_300',
                'bimanual_300_0', 'bimanual_300_60', 'bimanual_300_120', 'bimanual_300_180', 'bimanual_300_240', 'bimanual_300_300'],
        'reach_type': ['matched', 'unrelated', 'unrelated', 'mirror', 'unrelated', 'unrelated',
                    'unrelated', 'matched', 'mirror', 'unrelated', 'mirror-diagonal', 'unrelated',
                    'unrelated', 'mirror', 'matched', 'unrelated', 'unrelated', 'mirror-diagonal',
                    'mirror', 'unrelated', 'unrelated', 'matched', 'unrelated', 'unrelated',
                    'unrelated', 'mirror-diagonal', 'unrelated', 'unrelated', 'matched', 'mirror',
                    'unrelated', 'unrelated', 'mirror-diagonal', 'unrelated', 'mirror', 'matched'],
        'Uni_or_Bi': [1, 1, 1, 1, 1, 1,
                    1, 1, 1, 1, 1, 1,
                    1, 1, 1, 1, 1, 1,
                    1, 1, 1, 1, 1, 1,
                    1, 1, 1, 1, 1, 1,
                    1, 1, 1, 1, 1, 1],
        'Hand': [-1, -1, -1, -1, -1, -1,
                -1, -1, -1, -1, -1, -1,
                -1, -1, -1, -1, -1, -1,
                -1, -1, -1, -1, -1, -1,
                -1, -1, -1, -1, -1, -1,
                -1, -1, -1, -1, -1, -1],
        'targetAngle_L': [0, 0, 0, 0, 0, 0,
                        60, 60, 60, 60, 60, 60,
                        120, 120, 120, 120, 120, 120,
                        180, 180, 180, 180, 180, 180,
                        240, 240, 240, 240, 240, 240,
                        300, 300, 300, 300, 300, 300],
        'targetAngle_R': [0, 60, 120, 180, 240, 300,
                        0, 60, 120, 180, 240, 300,
                        0, 60, 120, 180, 240, 300,
                        0, 60, 120, 180, 240, 300,
                        0, 60, 120, 180, 240, 300,
                        0, 60, 120, 180, 240, 300],
    }
    
    # combine the dictionaries:
    combined = {k: tmp_left[k] + tmp_right[k] + tmp_bimanual[k] for k in tmp_left}
    conds = pd.DataFrame(combined)

    # Load the training .dat file:
    dat_file_name = os.path.join(path['fMRI_behavDir'], f's{sn}', f's{sn}_scan.dat')
    dat = pd.read_table(dat_file_name)
    
    fMRI_sess = []
    idx_preRT = []
    idx_postRT = []
    idx_gocue = []
    idx_endReach = []
    MD_left = []
    MD_right = []

    oldblock = -1
    # loop on trials:
    for i in tqdm(range(dat.shape[0])):
        if dat['BN'][i] != oldblock:
            # print(f'Processing block {dat["BN"][i]}')
            # load the .mov file:
            ext = int(dat['BN'][i])
            mov = movload(os.path.join(path['fMRI_behavDir'], f's{sn}', f's{sn}_scan_{ext:02d}.mov'))
            oldblock = dat['BN'][i]
        # print(f'Processing trial {dat["TN"][i]}')
        # trial routine:
        C = trial_routine(dat.iloc[[i]])

        # append the trial to the dataframe:
        D = pd.concat([D, C], ignore_index=True)
        
        # Whether data is from fMRI or training:
        fMRI_sess.append(int(1))

        # movemean the kinematics:
        trial_mov = mov[dat['TN'][i]-1]
        if smooth_win_sz > 0:
            # smooth the radius and angle signals:
            trial_mov[:, 5:9] = moving_average(trial_mov[:, 5:9], smooth_win_sz)
        # add the mov trial in the move dataframe:
        tmp = pd.DataFrame({'fMRI_sess': [1], 
                            'sn': [sn], 
                            'BN': [dat['BN'][i]], 
                            'TN': [dat['TN'][i]], 
                            'GoodMovement': [dat['GoodMovement'][i]], 
                            'mov': pd.Series([trial_mov], dtype='object')})
        
        # find the index of the go cue and end of reaching:
        if C.GoodMovement[0]:
            state = trial_mov[:,2]
            idx_gocue.append(np.where(state==5)[0][0])
            
            dEndRadius_L = C.dEndRadius_L.values[0]
            dEndRadius_R = C.dEndRadius_R.values[0]

            if dEndRadius_L > 0:
                tmp_idx = find_closest_index(trial_mov[:, 5], dEndRadius_L)
            
            if dEndRadius_R > 0:
                tmp_idx = find_closest_index(trial_mov[:, 6], dEndRadius_R)
            
            idx_endReach.append(tmp_idx)
        else:
            idx_gocue.append(0)
            idx_endReach.append(0)
        
        # Calculate mean deviation for the trial:
        if C.GoodMovement[0]:
            radius_l = trial_mov[:,5].flatten()[idx_gocue[i]:idx_endReach[i]]
            radius_r = trial_mov[:,6].flatten()[idx_gocue[i]:idx_endReach[i]]
            angle_l = trial_mov[:,7].flatten()[idx_gocue[i]:idx_endReach[i]]
            angle_r = trial_mov[:,8].flatten()[idx_gocue[i]:idx_endReach[i]]

            x_l = radius_l * np.cos(np.deg2rad(angle_l))
            y_l = radius_l * np.sin(np.deg2rad(angle_l))
            x_r = radius_r * np.cos(np.deg2rad(angle_r))
            y_r = radius_r * np.sin(np.deg2rad(angle_r))

            # left hand coords:
            f_left = np.vstack((x_l, y_l)).T    # position of left hand
            c_left = f_left[-1] - f_left[0]     # straight trajectory defined as position of go-cue and end of reach
            MD_left.append(get_MD(f_left, c_left))
            
            # right hand coords:
            f_right = np.vstack((x_r, y_r)).T
            c_right = f_right[-1] - f_right[0]
            MD_right.append(get_MD(f_right, c_right))
        else:
            MD_left.append(-1)
            MD_right.append(-1)

        # flatten the mov dataframe:
        flat_rows = []
        for _, row in tmp.iterrows():
            mat = row['mov']  # shape (n, m)
            meta = row.drop('mov').to_dict()  # get all metadata except 'mov'

            for r in mat:
                combined_row = meta.copy()
                for j, val in enumerate(r):
                    combined_row[f'mov_{j}'] = val
                flat_rows.append(combined_row)

        # Turn into flat DataFrame
        if D_mov is None:
            D_mov = pd.DataFrame(flat_rows)
        else:
            D_mov = pd.concat([D_mov, pd.DataFrame(flat_rows)], ignore_index=True)

    D['fMRI_sess'] = fMRI_sess
    D['idx_gocue'] = idx_gocue
    D['idx_endReach'] = idx_endReach
    D['MD_left'] = MD_left
    D['MD_right'] = MD_right
    
    cond_name = np.empty(len(D), dtype=object)
    reach_type = np.zeros(len(D), dtype=object)
    for index, row in conds.iterrows():
        if row.Uni_or_Bi == 0:  # unimanual
            if row.Hand == 0:   # left hand
                rows = ((D.Uni_or_Bi == 0) & (D.Hand == 0) & (D.targetAngle_L == row.targetAngle_L)).values.flatten()
                cond_name[rows] = row.cond
                reach_type[rows] = row.reach_type
            else:   # right hand
                rows = ((D.Uni_or_Bi == 0) & (D.Hand == 1) & (D.targetAngle_R == row.targetAngle_R)).values.flatten()
                cond_name[rows] = row.cond
                reach_type[rows] = row.reach_type
        else:  # bimanual
            rows = ((D.Uni_or_Bi == 1) & (D.targetAngle_L == row.targetAngle_L) & (D.targetAngle_R == row.targetAngle_R)).values.flatten()
            cond_name[rows] = row.cond
            reach_type[rows] = row.reach_type
    D['cond_name'] = cond_name
    D['reach_type'] = reach_type
    D['cond_name'] = D['cond_name'].astype(str)
    D['reach_type'] = D['reach_type'].astype(str)

    # save the data frames:
    D.to_csv(os.path.join(path['anaDir'], f's{sn}_scan.csv'), index=False)    
    D_mov.to_csv(os.path.join(path['anaDir'], f's{sn}_scan_mov.csv'), index=False)

def make_all_dataframe(df_list, sn_list, path):
    """
    This function takes a list of dataframes and a list of subject numbers and concatenates them into a single dataframe.
    
    params:
        df_list: list of dataframes
        sn_list: list of subject numbers
        path: dict, the path to directories to access data
    """
    # Concatenate all dataframes in the list
    bmw_all = pd.concat(df_list, ignore_index=True)
    
    # make sn column from sn_list:
    bmw_all['sn'] = np.concatenate([[sn] * len(df) for sn, df in zip(sn_list, df_list)])

    # save the dataframe:
    bmw_all.to_csv(os.path.join(path['anaDir'], 'bmw_all.csv'), index=False)
    
    return bmw_all

def make_summary_dataframe(path):
    """
        makes a summary dataframe from the bmw_all dataframe
    """
    # Concatenate all dataframes in the list
    bmw_all = pd.read_csv(os.path.join(path['anaDir'], 'bmw_all.csv'))
    bmw_all = bmw_all[bmw_all['GoodMovement'] == 1]
    
    # group by:
    bmw = bmw_all.groupby(['sn', 'cond_name']).agg(
        Uni_or_Bi=('Uni_or_Bi', 'first'),
        Hand=('Hand', 'first'),
        targetAngle_L=('targetAngle_L', 'first'),
        targetAngle_R=('targetAngle_R', 'first'),
        time2plan=('time2plan', 'first'),
        dEndRadius_L=('dEndRadius_L', 'first'),
        dEndRadius_R=('dEndRadius_R', 'first'),
        dEndAngle_L=('dEndAngle_L', 'first'),
        dEndAngle_R=('dEndAngle_R', 'first'),
        fMRI_sess=('fMRI_sess', 'first'),
        reach_type=('reach_type', 'first'),
        RT=('RT', lambda x: np.median(x)),
        MT=('MT', lambda x: np.median(x)),
        MD_left=('MD_left', lambda x: np.mean(x)),
        MD_right=('MD_right', lambda x: np.mean(x))
    ).reset_index()
    
    # save the dataframe:
    bmw.to_csv(os.path.join(path['anaDir'], 'bmw.csv'), index=False)
    
    return bmw