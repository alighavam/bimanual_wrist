import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import PcmPy as pcm
import getpass

baseDir = os.path.join('/Users', getpass.getuser(), 'Desktop', 'Projects', 'bimanual_wrist', 'data', 'fMRI')
bidsDir = 'BIDS'
anatomicalDir = 'anatomicals'
freesurferDir = 'surfaceFreesurfer'
surfacewbDir = 'surfaceWB' 
behavDir = 'behavioural'
regDir = 'ROI'
atlasDir = '/Volumes/diedrichsen_data$/data/Atlas_templates/fs_LR_32'
analysisDir = os.path.join(os.path.dirname(os.path.dirname(baseDir)), 'analysis')

# non-structured bimanual component model:
def non_structured_bimanual_components(verbose=True):
    '''
        we will make a non-structured contra, ipsi, and interaction components.
    '''
    ncond = 36
    G_model = {}
    # Contra component:
    cov = np.kron(np.eye(6), np.ones((6, 6)))
    cov = pcm.centering(36) @ cov @ pcm.centering(36)
    cov = cov/np.abs(np.trace(cov))/ncond
    G_model['contra'] = cov

    # Ipsilateral component:
    cov = np.kron(np.ones((6, 6)), np.eye(6))
    cov = pcm.centering(36) @ cov @ pcm.centering(36)
    cov = cov/np.abs(np.trace(cov))/ncond
    G_model['ipsi'] = cov  

    # interaction term:
    cov = np.eye(ncond)
    cov = pcm.centering(36) @ cov @ pcm.centering(36)
    cov = cov/np.abs(np.trace(cov))/ncond
    G_model['interaction'] = cov

    if verbose:
        fig, axes = plt.subplots(1, 3, figsize=(6,2))
        # plot the components:
        for i in range(3):
            G_tmp = G_model[list(G_model.keys())[i]]
            vmax = np.max(np.abs(G_tmp))
            vmin = -vmax
            axes[i].imshow(G_tmp, cmap='RdBu_r', vmin=vmin, vmax=vmax)
            axes[i].set_xticklabels('')
            axes[i].set_yticklabels('')
        plt.suptitle('Non-structured bimanual components')
        plt.tight_layout()
        plt.show()

    return G_model

# structured bimanual component model:
def structured_bimanual_components(region, verbose=True):
    '''
        we will make a structured contra, ipsi, and interaction components.
        the structured components are:
        - contra: 
        - ipsi: 
        - interaction: 
    '''
    labels = ['flx_flx',    'flx_flxup',   'flx_extup',   'flx_ext',   'flx_extdn',   'flx_flxdn',
          'flxup_flx',  'flxup_flxup', 'flxup_extup', 'flxup_ext', 'flxup_extdn', 'flxup_flxdn',
          'extup_flx',  'extup_flxup', 'extup_extup', 'extup_ext', 'extup_extdn', 'extup_flxdn',
          'ext_flx',    'ext_flxup',   'ext_extup',   'ext_ext',   'ext_extdn',   'ext_flxdn',
          'extdn_flx',  'extdn_flxup', 'extdn_extup', 'extdn_ext', 'extdn_extdn', 'extdn_flxdn',
          'flxdn_flx',  'flxdn_flxup', 'flxdn_extup', 'flxdn_ext', 'flxdn_extdn', 'flxdn_flxdn']
    ncond = 36

    # ================= Load unimanual 6by6 =================
    file_path_unimanual = os.path.join(analysisDir, f'pcm_dataset_6by6_{region}.npz')
    U = np.load(file_path_unimanual, allow_pickle=True)
    U = U['Y']
    Guni = np.zeros((6, 6))
    for i in range(len(U)):
        tmp,_ = pcm.est_G_crossval(U[i].measurements,
                                    U[i].obs_descriptors['cond_vec'],
                                    U[i].obs_descriptors['part_vec'],
                                    X=pcm.matrix.indicator(U[i].obs_descriptors['part_vec']))
        Guni += tmp
    Guni /= len(U)
    Guni = Guni / np.trace(np.abs(Guni))
    
    G_model = {}
    
    # ================= MAKE CONTRA COMPONENT =================
    F_contra = np.zeros((36,6))
    cnd2idx = {'flx':0, 'flxup':1, 'extup':2, 'ext':3, 'extdn':4, 'flxdn':5}
    for i in range(36):
        cond_pair = labels[i].split('_')
        cnd_contra = cond_pair[0]
        cnd_ipsi = cond_pair[1]
        idx_contra = cnd2idx[cnd_contra]
        idx_ipsi = cnd2idx[cnd_ipsi]
        F_contra[i, idx_contra] = 1
    # center the features:
    F_contra -= np.mean(F_contra, axis=0)
    G_model['contra'] = F_contra @ Guni @ F_contra.T
    G_model['contra'] = G_model['contra'] / np.trace(np.abs(G_model['contra']))
    
    # ================= MAKE IPSILATERAL COMPONENT ==================
    F_ipsi = np.zeros((36,6))
    for i in range(36):
        cond_pair = labels[i].split('_')
        cnd_contra = cond_pair[0]
        cnd_ipsi = cond_pair[1]
        idx_contra = cnd2idx[cnd_contra]
        idx_ipsi = cnd2idx[cnd_ipsi]
        F_ipsi[i, idx_ipsi] = 1
    # center the features:
    F_ipsi -= np.mean(F_ipsi, axis=0)
    G_model['ipsi'] = F_ipsi @ Guni @ F_ipsi.T
    G_model['ipsi'] = G_model['ipsi'] / np.trace(np.abs(G_model['ipsi']))

    # interaction term:
    cov = np.eye(ncond)
    cov = pcm.centering(36) @ cov @ pcm.centering(36)
    cov = cov/np.abs(np.trace(cov))/ncond
    G_model['interaction'] = cov
    G_model['interaction'] = G_model['interaction'] / np.trace(np.abs(G_model['interaction']))

    if verbose:
        fig, axes = plt.subplots(1, 3, figsize=(6,2))
        # plot the components:
        for i in range(3):
            G_tmp = G_model[list(G_model.keys())[i]]
            vmax = np.max(np.abs(G_tmp))
            vmin = -vmax
            axes[i].imshow(G_tmp, cmap='RdBu_r', vmin=vmin, vmax=vmax)
            axes[i].set_xticklabels('')
            axes[i].set_yticklabels('')
        plt.suptitle('Structured bimanual components')
        plt.tight_layout()
        plt.show()

    return G_model