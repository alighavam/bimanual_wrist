import os
import numpy as np
import pandas as pd
from .please import *

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
    
def subject_routine(sn, path, smooth_win_sz=0, fs=200):
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
    dat_file_name = os.path.join(path['train_behavDir'], f's{sn:02d}', f'BimanualWrist_MR_{sn}.dat')
    dat = pd.read_table(dat_file_name)

    fMRI_sess = []
    oldblock = -1
    # loop on trials:
    for i in range(dat.shape[0]):
        if dat['BN'][i] != oldblock:
            print(f'Processing block {dat["BN"][i]}')
            # load the .mov file:
            ext = int(dat['BN'][i])
            mov = movload(os.path.join(path['train_behavDir'], f's{sn:02d}', f'BimanualWrist_MR_{sn}_{ext:02d}.mov'))
            oldblock = dat['BN'][i]
        print(f'Processing trial {dat["TN"][i]}')
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

    cond_name = np.zeros(len(D))
    reach_type = np.zeros(len(D))
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
    D.to_csv(os.path.join(path['anaDir'], f'bmw_train_{sn}.csv'), index=False)    
    D_mov.to_csv(os.path.join(path['anaDir'], f'bmw_train_{sn}_mov.csv'), index=False)

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

    # Load the training .dat file:
    dat_file_name = os.path.join(path['fMRI_behavDir'], f's{sn:02d}', f'BimanualWrist_MR_{sn}.dat')
    dat = pd.read_table(dat_file_name)
    
    fMRI_sess = []

    oldblock = -1
    # loop on trials:
    for i in range(dat.shape[0]):
        if dat['BN'][i] != oldblock:
            print(f'Processing block {dat["BN"][i]}')
            # load the .mov file:
            ext = int(dat['BN'][i])
            mov = movload(os.path.join(path['train_behavDir'], f's{sn:02d}', f'BimanualWrist_MR_{sn}_{ext:02d}.mov'))
            oldblock = dat['BN'][i]
        print(f'Processing trial {dat["TN"][i]}')
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
    
    # sort the dataframes by :
    # D = D.sort_values(by='day', kind='mergesort')
    # df_mov = df_mov.sort_values(by='day', kind='mergesort')

    # save the data frames:
    D.to_csv(os.path.join(path['anaDir'], f'bmw_fMRI_{sn}.csv'), index=False)    
    D_mov.to_csv(os.path.join(path['anaDir'], f'bmw_fMRI_{sn}_mov.csv'), index=False)

    
                


