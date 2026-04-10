import os, sys
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'src'))
import numpy as np
import PcmPy as pcm
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as stats
import pandas as pd
from .please import *
from mpl_toolkits.axes_grid1 import make_axes_locatable
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

def make_data_component_eye(nsim=50, sig=1, theta_true=np.array([-10e10,1,1,1,-10e10]), verbose=True):
    labels = ['flx', 'flxup', 'extup', 'ext', 'extdn', 'flxdn', 'flx', 'flxup', 'extup', 'ext', 'extdn', 'flxdn']
    signal = [sig]*nsim
    ncond = 12
    num_items = 6
    # theta_true = np.array([-10e10,1,1,-10e10,1]) # extrinsic
    # theta_true = np.array([-10e10,1,1,10,-10e10]) # intrinsic
    # theta_true = np.array([-10e10,1,1,1,-0.2]) # intrinsic + extrinsic
    # theta_true = np.array([-10e10,1,0.8,-10e10,-10e10]) # no covariance

    G_model = {}
    # hand component:
    cov = np.zeros((ncond, ncond))
    cov[0:6,0:6] = 1
    cov[6:12,6:12] = 1
    cov = pcm.centering(ncond) @ cov @ pcm.centering(ncond)
    cov = cov / np.trace(cov)
    G_model['hand'] = cov

    # contra hand:
    cov_eye = np.eye(num_items)
    # cov_eye = np.eye(num_items)
    cov = np.zeros((ncond, ncond))
    cov[0:6,0:6] = cov_eye
    G_model['contra'] = cov

    # ipsi hand:
    cov = np.zeros((ncond, ncond))
    cov[6:12,6:12] = cov_eye
    G_model['ipsi'] = cov

    # intrinsic:
    cov = np.zeros((ncond, ncond))
    cov[0:6,6:12] = cov_eye
    cov[6:12,0:6] = cov_eye.T
    G_model['intrinsic'] = cov

    # extrinsic:
    cov = np.zeros((ncond, ncond))
    cov[0:6,6:12] = cov_eye[:, [3,2,1,0,5,4]]
    cov[6:12,0:6] = cov_eye[:, [3,2,1,0,5,4]].T
    G_model['extrinsic'] = cov

    Mtrue = pcm.ComponentModel('trueG', [G_model['hand'], 
                                        G_model['contra'], G_model['ipsi'], 
                                        G_model['intrinsic'], G_model['extrinsic']])

    # ============= Visualize The Generative Model =============
    if verbose:
        H = Mtrue.n_param
        for i in range(H):
            plt.subplot(1, H, i+1)
            vmin = -np.max(np.abs(Mtrue.Gc[i,:,:]))
            vmax = np.max(np.abs(Mtrue.Gc[i,:,:]))
            plt.imshow(Mtrue.Gc[i,:,:], vmin=vmin, vmax=vmax, cmap='RdBu_r')
        plt.tight_layout()
        plt.show()

    # ============= Generate The Data =============
    cond_vec, part_vec = pcm.sim.make_design(n_cond=12, n_part=10)
    D = pcm.sim.make_dataset(model=Mtrue, \
        theta=theta_true,
        cond_vec=cond_vec,
        part_vec=part_vec,
        n_sim=nsim,
        n_channel=300,
        signal=signal)

    # ============= Estimate The Data G Matrix =============
    N = len(D)
    G_hat = np.zeros((N, ncond, ncond))
    Dist = np.zeros((N,ncond,ncond))
    for i in range(N):
        G_hat[i, :, :], _ = pcm.est_G_crossval(D[i].measurements,
                                                D[i].obs_descriptors['cond_vec'],
                                                D[i].obs_descriptors['part_vec'],
                                                X=pcm.matrix.indicator(D[i].obs_descriptors['part_vec']))
        Dist[i,:,:] = pcm.G_to_dist(G_hat[i,:,:])

    # ============= Visualize The Estimated G Matrix =============
    if verbose:
        plt.rcParams.update({'font.size': 7})
        G_mean = np.mean(G_hat, axis=0)
        plt.figure(figsize=(2.5,2.5))
        vmin = -np.max(np.abs(G_mean))
        vmax = np.max(np.abs(G_mean))
        plt.imshow(G_mean, cmap='RdBu_r', vmin=vmin, vmax=vmax)
        plt.colorbar()
        plt.title(f'Estimated G, SNR={signal[0]}')
        plt.xlabel('Conditions')
        plt.ylabel('Conditions')
        plt.xticks(ticks=np.arange(ncond), labels=labels, rotation=45)
        plt.yticks(ticks=np.arange(ncond), labels=labels)
        plt.axhline(5.5, color='k', linestyle='--')
        plt.axvline(5.5, color='k', linestyle='--')
        plt.tight_layout()
        plt.show()
    
    return D, G_hat

def make_data_correlation_eye(nsim=50, sig=1, r=1, theta_true=np.array([0,1,-1.5,-1]), verbose=True):
    labels = ['flx', 'flxup', 'extup', 'ext', 'extdn', 'flxdn', 'flx', 'flxup', 'extup', 'ext', 'extdn', 'flxdn']
    signal = [sig]*nsim
    ncond = 12
    num_items = 6
    C = pcm.centering(num_items) @ np.eye(num_items) @ pcm.centering(num_items)
    C = np.eye(num_items)
    within_cov = C.reshape(1,num_items,num_items)
    Mtrue = pcm.CorrelationModel('corr', num_items=num_items, corr=r, cond_effect=True,
                                within_cov=within_cov)

    # ============= Visualize The Generative Model =============
    if verbose:
        plt.figure(figsize=(7,2))
        H = Mtrue.n_param
        for i in range(H):
            plt.subplot(1, H, i+1)
            vmin = -np.max(np.abs(Mtrue.Gc[i,:,:]))
            vmax = np.max(np.abs(Mtrue.Gc[i,:,:]))
            plt.imshow(Mtrue.Gc[i,:,:], vmin=vmin, vmax=vmax, cmap='RdBu_r')
            plt.colorbar()
        plt.suptitle('True model components')
        plt.tight_layout()
        plt.show()

    # ============= Generate The Data =============
    cond_vec, part_vec = pcm.sim.make_design(n_cond=12, n_part=10)
    D = pcm.sim.make_dataset(model=Mtrue, \
        theta=theta_true,
        cond_vec=cond_vec,
        part_vec=part_vec,
        n_sim=nsim,
        n_channel=300,
        signal=signal)

    # ============= Estimate The Data G Matrix =============
    N = len(D)
    G_hat = np.zeros((N, ncond, ncond))
    Dist = np.zeros((N,ncond,ncond))
    for i in range(N):
        G_hat[i, :, :], _ = pcm.est_G_crossval(D[i].measurements,
                                                D[i].obs_descriptors['cond_vec'],
                                                D[i].obs_descriptors['part_vec'],
                                                X=pcm.matrix.indicator(D[i].obs_descriptors['part_vec']))

    # ============= Visualize The Estimated G Matrix =============
    if verbose:
        plt.rcParams.update({'font.size': 7})
        G_mean = np.mean(G_hat, axis=0)
        plt.figure(figsize=(3,3))
        vmin = -np.max(np.abs(G_mean))
        vmax = np.max(np.abs(G_mean))
        plt.imshow(G_mean, cmap='RdBu_r', vmin=vmin, vmax=vmax)
        plt.colorbar()
        plt.suptitle('Estimated G_hat mean')
        plt.xlabel('Conditions')
        plt.ylabel('Conditions')
        plt.xticks(ticks=np.arange(ncond), labels=labels, rotation=45)
        plt.yticks(ticks=np.arange(ncond), labels=labels)
        plt.axhline(5.5, color='k', linestyle='--')
        plt.axvline(5.5, color='k', linestyle='--')
        plt.tight_layout()
        plt.show()

    return D, G_hat

def make_data_component_reach_return(nsim=50, sig=1, theta_true=np.array([-10e10,1,1,1,0.5]), verbose=True):
    labels = ['flx', 'flxup', 'extup', 'ext', 'extdn', 'flxdn', 'flx', 'flxup', 'extup', 'ext', 'extdn', 'flxdn']
    signal = [sig]*nsim
    ncond = 12
    num_items = 6
    # theta_true = np.array([-10e10,1,1,-10e10,1]) # extrinsic
    # theta_true = np.array([-10e10,1,1,10,-10e10]) # intrinsic
    # theta_true = np.array([-10e10,1,1,1,-0.2]) # intrinsic + extrinsic
    # theta_true = np.array([-10e10,1,0.8,-10e10,-10e10]) # no covariance

    G_model = {}
    # hand component:
    cov = np.zeros((ncond, ncond))
    cov[0:6,0:6] = 1
    cov[6:12,6:12] = 1
    cov = pcm.centering(ncond) @ cov @ pcm.centering(ncond)
    cov = cov / np.trace(cov)
    G_model['hand'] = cov

    # contra hand:
    cov_rr = np.eye(num_items)
    cov_rr[0,3] = 1
    cov_rr[1,4] = 1
    cov_rr[2,5] = 1
    cov_rr[3,0] = 1
    cov_rr[4,1] = 1
    cov_rr[5,2] = 1
    cov_rr = pcm.centering(num_items) @ cov_rr @ pcm.centering(num_items)
    cov = np.zeros((ncond, ncond))
    cov[0:6,0:6] = cov_rr
    G_model['contra'] = cov

    # ipsi hand:
    cov = np.zeros((ncond, ncond))
    cov[6:12,6:12] = cov_rr
    G_model['ipsi'] = cov

    # intrinsic:
    cov_eye = np.eye(num_items)
    cov = np.zeros((ncond, ncond))
    cov[0:6,6:12] = cov_eye
    cov[6:12,0:6] = cov_eye.T
    G_model['intrinsic'] = cov

    # extrinsic:
    cov = np.zeros((ncond, ncond))
    cov[0:6,6:12] = cov_eye[:, [3,2,1,0,5,4]]
    cov[6:12,0:6] = cov_eye[:, [3,2,1,0,5,4]].T
    G_model['extrinsic'] = cov

    Mtrue = pcm.ComponentModel('trueG', [G_model['hand'], 
                                        G_model['contra'], G_model['ipsi'], 
                                        G_model['intrinsic'], G_model['extrinsic']])

    # ============= Visualize The Generative Model =============
    if verbose:
        H = Mtrue.n_param
        for i in range(H):
            plt.subplot(1, H, i+1)
            vmin = -np.max(np.abs(Mtrue.Gc[i,:,:]))
            vmax = np.max(np.abs(Mtrue.Gc[i,:,:]))
            plt.imshow(Mtrue.Gc[i,:,:], vmin=vmin, vmax=vmax, cmap='RdBu_r')
        plt.tight_layout()
        plt.show()

    # ============= Generate The Data =============
    cond_vec, part_vec = pcm.sim.make_design(n_cond=12, n_part=10)
    D = pcm.sim.make_dataset(model=Mtrue, \
        theta=theta_true,
        cond_vec=cond_vec,
        part_vec=part_vec,
        n_sim=nsim,
        n_channel=300,
        signal=signal)

    # ============= Estimate The Data G Matrix =============
    N = len(D)
    G_hat = np.zeros((N, ncond, ncond))
    Dist = np.zeros((N,ncond,ncond))
    for i in range(N):
        G_hat[i, :, :], _ = pcm.est_G_crossval(D[i].measurements,
                                                D[i].obs_descriptors['cond_vec'],
                                                D[i].obs_descriptors['part_vec'],
                                                X=pcm.matrix.indicator(D[i].obs_descriptors['part_vec']))
        Dist[i,:,:] = pcm.G_to_dist(G_hat[i,:,:])

    # ============= Visualize The Estimated G Matrix =============
    if verbose:
        plt.rcParams.update({'font.size': 7})
        G_mean = np.mean(G_hat, axis=0)
        plt.figure(figsize=(2.5,2.5))
        vmin = -np.max(np.abs(G_mean))
        vmax = np.max(np.abs(G_mean))
        plt.imshow(G_mean, cmap='RdBu_r', vmin=vmin, vmax=vmax)
        plt.colorbar()
        plt.title(f'Estimated G, SNR={signal[0]}')
        plt.xlabel('Conditions')
        plt.ylabel('Conditions')
        plt.xticks(ticks=np.arange(ncond), labels=labels, rotation=45)
        plt.yticks(ticks=np.arange(ncond), labels=labels)
        plt.axhline(5.5, color='k', linestyle='--')
        plt.axvline(5.5, color='k', linestyle='--')
        plt.tight_layout()
        plt.show()
    
    return D, G_hat

def sort_as_extrinsic(D, verbose=True):
    D_extinrisic = []
    for d in D:
        obs_des = d.obs_descriptors
        new_obs_des = obs_des.copy()
        new_obs_des['cond_vec'] = np.tile(np.array([0,1,2,3,4,5, 9,8,7,6,11,10]),10)
        D_extinrisic.append(pcm.dataset.Dataset(d.measurements, obs_descriptors=new_obs_des))
    
    if verbose:
        labels = ['flx', 'flxup', 'extup', 'ext', 'extdn', 'flxdn', 'flx', 'flxup', 'extup', 'ext', 'extdn', 'flxdn']
        # ============= Visualize The Estimated G Matrix =============
        N = len(D_extinrisic)
        G_hat = np.zeros((N, 12, 12))
        for i in range(N):
            G_hat[i, :, :], _ = pcm.est_G_crossval(D_extinrisic[i].measurements,
                                                    D_extinrisic[i].obs_descriptors['cond_vec'],
                                                    D_extinrisic[i].obs_descriptors['part_vec'],
                                                    X=pcm.matrix.indicator(D_extinrisic[i].obs_descriptors['part_vec']))
        plt.rcParams.update({'font.size': 7})
        G_mean = np.mean(G_hat, axis=0)
        plt.figure(figsize=(3,3))
        vmin = -np.max(np.abs(G_mean))
        vmax = np.max(np.abs(G_mean))
        plt.imshow(G_mean, cmap='RdBu_r', vmin=vmin, vmax=vmax)
        plt.colorbar()
        plt.suptitle('G_hat mean sorted as EXTRINSIC')
        plt.xlabel('Conditions')
        plt.ylabel('Conditions')
        plt.xticks(ticks=np.arange(12), labels=labels, rotation=45)
        plt.yticks(ticks=np.arange(12), labels=labels)
        plt.axhline(5.5, color='k', linestyle='--')
        plt.axvline(5.5, color='k', linestyle='--')
        plt.tight_layout()
        plt.show()

    return D_extinrisic

def pearson_r(D, verbose=True):
    # ============= Find Intrinsic and Extrinsic r =============
    # specify_conditions = np.array([1,2,4,5])
    intrinsic = [0,1,2,3,4,5, 6,7,8,9,10,11]
    r_intrinsic = analyze_r(D, intrinsic)
    extrinsic = [0,1,2,3,4,5, 9,8,7,6,11,10]
    r_extrinsic = analyze_r(D, extrinsic)
    print(f'Intrinsic: mean={np.mean(r_intrinsic):.2f} +/- {stats.sem(r_intrinsic):.2f}')
    print(f'Extrinsic: mean={np.mean(r_extrinsic):.2f} +/- {stats.sem(r_extrinsic):.2f}')

    # ============= Visualize The Pearson r Distribution =============
    if verbose:
        # fig, ax = plt.subplots(1,2, figsize=(4,2))
        # sns.histplot(r_intrinsic, kde=True, color='gray', ax=ax[0], bins=10)
        # ax[0].axvline(np.mean(r_intrinsic), color='red', linestyle=':')
        # ax[0].set_xlabel('Pearson r')
        # ax[0].set_ylabel('Count')
        # ax[0].set_title(f'Intrinsic Pearson r, SNR={signal[0]}')
        # ax[0].set_xlim(-1.1,1.1)
        # sns.histplot(r_extrinsic, kde=True, color='gray', ax=ax[1], bins=10)
        # ax[1].axvline(np.mean(r_extrinsic), color='red', linestyle=':')
        # ax[1].set_xlabel('Pearson r')
        # ax[1].set_ylabel('Count')
        # ax[1].set_title(f'Extrinsic Pearson r, SNR={signal[0]}')
        # ax[1].set_xlim(-1.1,1.1)
        # plt.tight_layout()
        # plt.show()

        df = pd.DataFrame({'r': np.concatenate([r_intrinsic, r_extrinsic]),
                        'region': ['Intrinsic']*len(r_intrinsic) + ['Extrinsic']*len(r_extrinsic)})
        # boxplot:
        plt.figure(figsize=(2.5,2.5))
        plt.axhline(0, color='gray', linestyle=':')
        # sns.swarmplot(x='region', y='r', data=df, color='gray', size=1, alpha=0.5)
        sns.stripplot(x='region', y='r', data=df, jitter=True, size=2, alpha=0.6)
        sns.boxplot(x='region', y='r', data=df, width=0.5, color='k', fill=False, showcaps=False, linewidth=1, fliersize=0)
        plt.xlabel('Region')
        plt.ylabel('Pearson r')
        plt.title(f'SNR=')
        plt.ylim(-1.1,1.1)
        # rotate x-tick labels:
        plt.xticks(rotation=45)
        plt.tight_layout()
        plt.show()

    return r_intrinsic, r_extrinsic

def pcm_correlation_model(D, M:list, verbose=True):
    ## ================ Visualize Model Components ================
    if verbose:
        # visualize the model components:
        H = M[0].n_param-1
        plt.figure(figsize=(5,2))
        for i in range(H):
            plt.subplot(1, H, i+1)
            vmin = -np.max(np.abs(M[0].Gc[i,:,:]))
            vmax = np.max(np.abs(M[0].Gc[i,:,:]))
            plt.imshow(M[0].Gc[i,:,:], vmin=vmin, vmax=vmax, cmap='RdBu_r')
            plt.colorbar()
        plt.suptitle('fitting model components')
        plt.tight_layout()
        plt.show()

    ## ================ MODEL FITTING ================
    # If we want to remove the condition effect as a fixed effect, we need an indicator for that
    cond_vec = D[0].obs_descriptors['cond_vec']
    c_vec,item_vec= pcm.cond_to_item(cond_vec)
    X = pcm.matrix.indicator(c_vec)
    T, theta = pcm.fit_model_individ(D, M, fixed_effect=X, fit_scale=False, verbose=False)
    # fit group for visualization:
    T_group, theta_group = pcm.fit_model_group(D, M, fixed_effect=X, fit_scale=True, verbose=False)

    ## ================ PLOT PREDICTED G MATRIX ================
    if verbose:
        G_pred,_ = M[0].predict(theta_group[0][:M[0].n_param])
        # visualize the predicted G:
        plt.figure(figsize=(6,2))
        vmin = -np.max(np.abs(G_pred))
        vmax = np.max(np.abs(G_pred))
        plt.imshow(G_pred, vmin=vmin, vmax=vmax, cmap='RdBu_r')
        plt.axhline(5.5, color='k', linestyle='--')
        plt.axvline(5.5, color='k', linestyle='--')
        plt.colorbar()
        plt.title('Predicted G')
        plt.show()

    # Get the maximum likelihood estimate of each parameter
    if verbose:
        maxr = M[-1].get_correlation(theta[-1])
        fsnr = M[0].get_fSNR(theta[0], 10, separate=False)
        
        print(f'maxr = {np.mean(maxr):.3f} +/- {stats.sem(maxr):.3f}')

        # load results:
        fig, ax = plt.subplots(1,1, figsize=(4,4))
        ax.axhline(0, color='gray', linestyle=':', linewidth=1)
        ax.axhline(1, color='gray', linestyle=':', linewidth=1)
        ax.axhline(-1, color='gray', linestyle=':', linewidth=1)
        ax.scatter(fsnr, maxr, color='k', s=5)

        # --- vertical marginal histogram of correlation MLE (maxr) ---
        divider = make_axes_locatable(ax)
        ax_marg = divider.append_axes("right", size="20%", pad=0.02, sharey=ax)
        ax_marg.hist(maxr, bins=5, orientation='horizontal', color='gray', alpha=0.5, edgecolor='none')
        # clean up marginal axis
        ax_marg.set_xticks([])
        ax_marg.set_xlabel('')
        plt.setp(ax_marg.get_yticklabels(), visible=False)
        for spine in ax_marg.spines.values():
            spine.set_visible(False)
        # -------------------------------------------------------------
        ax.set_xlabel('fsnr')
        ax.set_title(f'correlation MLE')
        plt.show()

    return maxr, fsnr

def cka_model_fit_eye(D, verbose=True):
    ncond = 12
    num_items = 6

    # ================= SETUP MODEL ================
    G_model = {}
    # hand component:
    cov = np.zeros((ncond, ncond))
    cov[0:6,0:6] = 1
    cov[6:12,6:12] = 1
    cov = pcm.centering(ncond) @ cov @ pcm.centering(ncond)
    cov = cov / np.trace(cov)
    G_model['hand'] = cov

    # contra hand:
    cov_eye = np.eye(num_items)
    cov = np.zeros((ncond, ncond))
    cov[0:6,0:6] = cov_eye
    G_model['contra'] = cov

    # ipsi hand:
    cov = np.zeros((ncond, ncond))
    cov[6:12,6:12] = cov_eye
    G_model['ipsi'] = cov

    # intrinsic:
    cov = np.zeros((ncond, ncond))
    cov[0:6,6:12] = cov_eye
    cov[6:12,0:6] = cov_eye.T
    G_model['intrinsic'] = cov

    # extrinsic:
    cov = np.zeros((ncond, ncond))
    cov[0:6,6:12] = cov_eye[:, [3,2,1,0,5,4]]
    cov[6:12,0:6] = cov_eye[:, [3,2,1,0,5,4]].T
    G_model['extrinsic'] = cov

    Mfull = pcm.ComponentModel('full', [G_model['hand'], 
                                        G_model['contra'], G_model['ipsi'], 
                                        G_model['intrinsic'], G_model['extrinsic']])
    M_within = pcm.ComponentModel('within', [G_model['hand'], 
                                        G_model['contra'], G_model['ipsi']])
    M_intrinsic = pcm.ComponentModel('intrinsic', [G_model['hand'], 
                                        G_model['contra'], G_model['ipsi'], 
                                        G_model['intrinsic']])
    M_extrinsic = pcm.ComponentModel('extrinsic', [G_model['hand'], 
                                        G_model['contra'], G_model['ipsi'], 
                                        G_model['extrinsic']])
    
    ## ================ MODEL FITTING ================
    M = []
    M.append(Mfull)
    M.append(M_within)
    M.append(M_intrinsic)
    M.append(M_extrinsic)
    T, theta, ceil = pcm.fit_CKA_group_crossval(D, M, theta0=None, verbose=False, ceil=True)

    # Do model comparison stats on CKA:
    CKA = T['CKA']
    # ================ PLOT THE CKA OF MODELS ================
    if verbose:
        CKA_long = CKA.melt(var_name='model', value_name='CKA')
        plt.figure(figsize=(3,2))
        sns.boxplot(data=CKA_long, x='model', y='CKA', color='k', 
                    fill=False, order=['within','intrinsic','extrinsic','full'])
        plt.xlabel('Model')
        plt.ylabel('CKA')
        plt.title('CKA of Models')
        plt.xticks(rotation=45)
        plt.title('cross-validated CKA of Models')
        plt.tight_layout()
        plt.show()
    
    # ttest intrinsic > within:
    t, p = stats.ttest_rel(CKA['intrinsic'].values, CKA['within'].values, alternative='greater')
    print(f"t-test intrinsic > within: t({len(CKA['intrinsic'])-1})={t:.3f}, p={p:.4f}")
    # ttest extrinsic > within:
    t, p = stats.ttest_rel(CKA['extrinsic'], CKA['within'], alternative='greater')
    print(f"t-test extrinsic > within: t({len(CKA['extrinsic'])-1})={t:.3f}, p={p:.4f}")
    # ttest full > intrinsic:
    t, p = stats.ttest_rel(CKA['full'], CKA['intrinsic'], alternative='greater')
    print(f"t-test full > intrinsic: t({len(CKA['full'])-1})={t:.3f}, p={p:.4f}")
    # ttest full > extrinsic:
    t, p = stats.ttest_rel(CKA['full'], CKA['extrinsic'], alternative='greater')
    print(f"t-test full > extrinsic: t({len(CKA['full'])-1})={t:.3f}, p={p:.4f}")

    # ================ PLOT THE WEIGHT OF EACH COMPONENT ================
    if verbose:
        # plot the weights of the full model:
        w = np.asarray(theta[0])
        n_comp, n_cv = w.shape
        comp_names = ('hand', 'contra', 'ipsi', 'intrinsic', 'extrinsic')
        df_weights = pd.DataFrame({
            'component': np.repeat(comp_names, n_cv),
            'weight': w.ravel(),
        })
        df_weights['weight'] = np.exp(df_weights['weight'])
        plt.figure(figsize=(3, 2))
        sns.boxplot(data=df_weights, x='component', y='weight', order=comp_names, fliersize=0)
        plt.xlabel('Component')
        plt.ylabel('Weight')
        plt.title('Component weights')
        plt.xticks(rotation=15, ha='right')
        plt.tight_layout()
        plt.show()
    
    _, theta_group = pcm.fit_CKA_group(D, M, theta0=None, verbose=False, X=None)
    ## ================ PLOT PREDICTED G MATRIX ================
    if verbose:
        G_pred, _ = M[0].predict(theta_group[0][:M[0].n_param])
        # visualize the predicted G:
        plt.figure(figsize=(2.5,2.5))
        vmin = -np.max(np.abs(G_pred))
        vmax = np.max(np.abs(G_pred))
        plt.imshow(G_pred, vmin=vmin, vmax=vmax, cmap='RdBu_r')
        plt.axhline(5.5, color='k', linestyle='--')
        plt.axvline(5.5, color='k', linestyle='--')
        plt.colorbar()
        plt.title('Predicted G CKA')
        plt.tight_layout()
        plt.show()

    return T, theta, ceil

def cka_model_fit_rr(D, verbose=True):
    ncond = 12
    num_items = 6
    
    # ================= SETUP MODEL ================
    G_model = {}
    # hand component:
    cov = np.zeros((ncond, ncond))
    cov[0:6,0:6] = 1
    cov[6:12,6:12] = 1
    cov = pcm.centering(ncond) @ cov @ pcm.centering(ncond)
    cov = cov / np.trace(cov)
    G_model['hand'] = cov

    # contra hand:
    cov_rr = np.eye(num_items)
    cov_rr[0,3] = 1
    cov_rr[1,4] = 1
    cov_rr[2,5] = 1
    cov_rr[3,0] = 1
    cov_rr[4,1] = 1
    cov_rr[5,2] = 1
    cov = np.zeros((ncond, ncond))
    cov[0:6,0:6] = cov_rr
    G_model['contra'] = cov

    # ipsi hand:
    cov = np.zeros((ncond, ncond))
    cov[6:12,6:12] = cov_rr
    G_model['ipsi'] = cov

    # intrinsic:
    cov_eye = np.eye(num_items)
    cov = np.zeros((ncond, ncond))
    cov[0:6,6:12] = cov_eye
    cov[6:12,0:6] = cov_eye.T
    G_model['intrinsic'] = cov

    # extrinsic:
    cov = np.zeros((ncond, ncond))
    cov[0:6,6:12] = cov_eye[:, [3,2,1,0,5,4]]
    cov[6:12,0:6] = cov_eye[:, [3,2,1,0,5,4]].T
    G_model['extrinsic'] = cov

    Mfull = pcm.ComponentModel('full', [G_model['hand'], 
                                        G_model['contra'], G_model['ipsi'], 
                                        G_model['intrinsic'], G_model['extrinsic']])
    M_within = pcm.ComponentModel('within', [G_model['hand'], 
                                        G_model['contra'], G_model['ipsi']])
    M_intrinsic = pcm.ComponentModel('intrinsic', [G_model['hand'], 
                                        G_model['contra'], G_model['ipsi'], 
                                        G_model['intrinsic']])
    M_extrinsic = pcm.ComponentModel('extrinsic', [G_model['hand'], 
                                        G_model['contra'], G_model['ipsi'], 
                                        G_model['extrinsic']])
    
    ## ================ MODEL FITTING ================
    M = []
    M.append(Mfull)
    M.append(M_within)
    M.append(M_intrinsic)
    M.append(M_extrinsic)
    T, theta, ceil = pcm.fit_CKA_group_crossval(D, M, theta0=None, verbose=False, ceil=True)

    # ================ PLOT EACH MODEL G MATRIX ================
    if verbose:
        theta_example = [np.array([-10e10,1,1,1,0.7]),np.array([-10e10,1,1]),
                         np.array([-10e10,1,1,1]),np.array([-10e10,1,1,1])]
        fig, axes = plt.subplots(1, len(M), figsize=(10,2))
        for i in range(len(M)):
            G_pred, _ = M[i].predict(theta_example[i])
            # visualize the predicted G:
            vmin = -np.max(np.abs(G_pred))
            vmax = np.max(np.abs(G_pred))
            axes[i].imshow(G_pred, vmin=vmin, vmax=vmax, cmap='RdBu_r')
            axes[i].axhline(5.5, color='k', linestyle='--')
            axes[i].axvline(5.5, color='k', linestyle='--')
            axes[i].set_title(f'Model {i}')
            # axes[i].colorbar()
        fig.tight_layout()
        plt.show()
    
    # Do model comparison stats on CKA:
    CKA = T['CKA']
    # ================ PLOT THE CKA OF MODELS ================
    if verbose:
        CKA_long = CKA.melt(var_name='model', value_name='CKA')
        plt.figure(figsize=(3,2))
        sns.boxplot(data=CKA_long, x='model', y='CKA', color='k', 
                    fill=False, order=['within','intrinsic','extrinsic','full'])
        plt.xlabel('Model')
        plt.ylabel('CKA')
        plt.title('CKA of Models')
        plt.xticks(rotation=45)
        plt.title('cross-validated CKA of Models')
        plt.tight_layout()
        plt.show()
    
    # ttest intrinsic > within:
    t, p = stats.ttest_rel(CKA['intrinsic'].values, CKA['within'].values, alternative='greater')
    print(f"t-test intrinsic > within: t({len(CKA['intrinsic'])-1})={t:.3f}, p={p:.4f}")
    # ttest extrinsic > within:
    t, p = stats.ttest_rel(CKA['extrinsic'], CKA['within'], alternative='greater')
    print(f"t-test extrinsic > within: t({len(CKA['extrinsic'])-1})={t:.3f}, p={p:.4f}")
    # ttest full > intrinsic:
    t, p = stats.ttest_rel(CKA['full'], CKA['intrinsic'], alternative='greater')
    print(f"t-test full > intrinsic: t({len(CKA['full'])-1})={t:.3f}, p={p:.4f}")
    # ttest full > extrinsic:
    t, p = stats.ttest_rel(CKA['full'], CKA['extrinsic'], alternative='greater')
    print(f"t-test full > extrinsic: t({len(CKA['full'])-1})={t:.3f}, p={p:.4f}")

    # ================ PLOT THE WEIGHT OF EACH COMPONENT ================
    if verbose:
        # plot the weights of the full model:
        w = np.asarray(theta[0])
        n_comp, n_cv = w.shape
        comp_names = ('hand', 'contra', 'ipsi', 'intrinsic', 'extrinsic')
        df_weights = pd.DataFrame({
            'component': np.repeat(comp_names, n_cv),
            'weight': w.ravel(),
        })
        df_weights['weight'] = np.exp(df_weights['weight'])
        plt.figure(figsize=(3, 2))
        sns.boxplot(data=df_weights, x='component', y='weight', order=comp_names, fliersize=0)
        plt.xlabel('Component')
        plt.ylabel('Weight')
        plt.title('Component weights')
        plt.xticks(rotation=15, ha='right')
        plt.tight_layout()
        plt.show()
    
    _, theta_group = pcm.fit_CKA_group(D, M, theta0=None, verbose=False, X=None)
    ## ================ PLOT PREDICTED and DATA G MATRIX ================
    if verbose:
        N = len(D)
        G_hat = np.zeros((N, ncond, ncond))
        for i in range(N):
            G_hat[i, :, :], _ = pcm.est_G_crossval(D[i].measurements,
                                                    D[i].obs_descriptors['cond_vec'],
                                                    D[i].obs_descriptors['part_vec'],
                                                    X=pcm.matrix.indicator(D[i].obs_descriptors['part_vec']))
        G_mean = np.mean(G_hat, axis=0)
        # visualize the data G:
        fig, ax = plt.subplots(1,2, figsize=(5,2.5))
        vmin = -np.max(np.abs(G_mean))
        vmax = np.max(np.abs(G_mean))
        im0 = ax[0].imshow(G_mean, vmin=vmin, vmax=vmax, cmap='RdBu_r')
        ax[0].axhline(5.5, color='k', linestyle='--')
        ax[0].axvline(5.5, color='k', linestyle='--')
        fig.colorbar(im0, ax=ax[0], fraction=0.046, pad=0.04)
        ax[0].set_title('Estimated Data G')

        G_pred, _ = M[0].predict(theta_group[0][:M[0].n_param])
        # visualize the predicted G:
        vmin = -np.max(np.abs(G_pred))
        vmax = np.max(np.abs(G_pred))
        im1 = ax[1].imshow(G_pred, vmin=vmin, vmax=vmax, cmap='RdBu_r')
        ax[1].axhline(5.5, color='k', linestyle='--')
        ax[1].axvline(5.5, color='k', linestyle='--')
        fig.colorbar(im1, ax=ax[1], fraction=0.046, pad=0.04)
        ax[1].set_title('Predicted G CKA')
        fig.tight_layout()
        plt.show()

    return T, theta, ceil

def cka_model_fit_structure(D, region, verbose=True):
    ncond = 12
    num_items = 6

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
    
    # ================= SETUP MODEL ================
    G_model = {}
    # hand component:
    cov = np.zeros((ncond, ncond))
    cov[0:6,0:6] = 1
    cov[6:12,6:12] = 1
    cov = pcm.centering(ncond) @ cov @ pcm.centering(ncond)
    cov = cov / np.trace(cov)
    G_model['hand'] = cov

    # contra hand:
    cov = np.zeros((ncond, ncond))
    cov[0:6,0:6] = Guni
    G_model['contra'] = cov

    # ipsi hand:
    cov = np.zeros((ncond, ncond))
    cov[6:12,6:12] = Guni
    G_model['ipsi'] = cov

    # intrinsic:
    cov_eye = np.eye(num_items)/6
    cov = np.zeros((ncond, ncond))
    cov[0:6,6:12] = cov_eye
    cov[6:12,0:6] = cov_eye.T
    G_model['intrinsic'] = cov
     
    # extrinsic:
    cov = np.zeros((ncond, ncond))
    cov[0:6,6:12] = cov_eye[:, [3,2,1,0,5,4]]
    cov[6:12,0:6] = cov_eye[:, [3,2,1,0,5,4]].T
    G_model['extrinsic'] = cov

    Mfull = pcm.ComponentModel('full', [G_model['hand'], 
                                        G_model['contra'], G_model['ipsi'], 
                                        G_model['intrinsic'], G_model['extrinsic']])
    M_within = pcm.ComponentModel('within', [G_model['hand'], 
                                        G_model['contra'], G_model['ipsi']])
    M_intrinsic = pcm.ComponentModel('intrinsic', [G_model['hand'], 
                                        G_model['contra'], G_model['ipsi'], 
                                        G_model['intrinsic']])
    M_extrinsic = pcm.ComponentModel('extrinsic', [G_model['hand'], 
                                        G_model['contra'], G_model['ipsi'], 
                                        G_model['extrinsic']])
    
    ## ================ MODEL FITTING ================
    M = []
    M.append(Mfull)
    M.append(M_within)
    M.append(M_intrinsic)
    M.append(M_extrinsic)
    T, theta, ceil = pcm.fit_CKA_group_crossval(D, M, theta0=None, verbose=False, ceil=True)

    # ================ PLOT EACH MODEL G MATRIX ================
    if verbose:
        theta_example = [np.array([-10e10,1,1,1,0.7]),np.array([-10e10,1,1]),
                         np.array([-10e10,1,1,1]),np.array([-10e10,1,1,1])]
        fig, axes = plt.subplots(1, len(M), figsize=(10,2))
        for i in range(len(M)):
            G_pred, _ = M[i].predict(theta_example[i])
            # visualize the predicted G:
            vmin = -np.max(np.abs(G_pred))
            vmax = np.max(np.abs(G_pred))
            axes[i].imshow(G_pred, vmin=vmin, vmax=vmax, cmap='RdBu_r')
            axes[i].axhline(5.5, color='k', linestyle='--')
            axes[i].axvline(5.5, color='k', linestyle='--')
            axes[i].set_title(f'Model {i}')
            # axes[i].colorbar()
        fig.tight_layout()
        plt.show()
    
    # Do model comparison stats on CKA:
    CKA = T['CKA']
    # ================ PLOT THE CKA OF MODELS ================
    if verbose:
        CKA_long = CKA.melt(var_name='model', value_name='CKA')
        plt.figure(figsize=(3,2))
        sns.boxplot(data=CKA_long, x='model', y='CKA', color='k', 
                    fill=False, order=['within','intrinsic','extrinsic','full'])
        plt.xlabel('Model')
        plt.ylabel('CKA')
        plt.title('CKA of Models')
        plt.xticks(rotation=45)
        plt.title('cross-validated CKA of Models')
        plt.tight_layout()
        plt.show()
    
    # ttest intrinsic > within:
    t, p = stats.ttest_rel(CKA['intrinsic'].values, CKA['within'].values, alternative='greater')
    print(f"t-test intrinsic > within: t({len(CKA['intrinsic'])-1})={t:.3f}, p={p:.4f}")
    # ttest extrinsic > within:
    t, p = stats.ttest_rel(CKA['extrinsic'], CKA['within'], alternative='greater')
    print(f"t-test extrinsic > within: t({len(CKA['extrinsic'])-1})={t:.3f}, p={p:.4f}")
    # ttest full > intrinsic:
    t, p = stats.ttest_rel(CKA['full'], CKA['intrinsic'], alternative='greater')
    print(f"t-test full > intrinsic: t({len(CKA['full'])-1})={t:.3f}, p={p:.4f}")
    # ttest full > extrinsic:
    t, p = stats.ttest_rel(CKA['full'], CKA['extrinsic'], alternative='greater')
    print(f"t-test full > extrinsic: t({len(CKA['full'])-1})={t:.3f}, p={p:.4f}")

    # ================ PLOT THE WEIGHT OF EACH COMPONENT ================
    if verbose:
        # plot the weights of the full model:
        w = np.asarray(theta[0])
        n_comp, n_cv = w.shape
        comp_names = ('hand', 'contra', 'ipsi', 'intrinsic', 'extrinsic')
        df_weights = pd.DataFrame({
            'component': np.repeat(comp_names, n_cv),
            'weight': w.ravel(),
        })
        df_weights['weight'] = np.exp(df_weights['weight'])
        plt.figure(figsize=(3, 2))
        sns.boxplot(data=df_weights, x='component', y='weight', order=comp_names, fliersize=0)
        plt.xlabel('Component')
        plt.ylabel('Weight')
        plt.title('Component weights')
        plt.xticks(rotation=15, ha='right')
        plt.tight_layout()
        plt.show()
    
    _, theta_group = pcm.fit_CKA_group(D, M, theta0=None, verbose=False, X=None)
    ## ================ PLOT PREDICTED and DATA G MATRIX ================
    if verbose:
        N = len(D)
        G_hat = np.zeros((N, ncond, ncond))
        for i in range(N):
            G_hat[i, :, :], _ = pcm.est_G_crossval(D[i].measurements,
                                                    D[i].obs_descriptors['cond_vec'],
                                                    D[i].obs_descriptors['part_vec'],
                                                    X=pcm.matrix.indicator(D[i].obs_descriptors['part_vec']))
        G_mean = np.mean(G_hat, axis=0)
        # visualize the data G:
        fig, ax = plt.subplots(1,2, figsize=(5,2.5))
        vmin = -np.max(np.abs(G_mean))
        vmax = np.max(np.abs(G_mean))
        im0 = ax[0].imshow(G_mean, vmin=vmin, vmax=vmax, cmap='RdBu_r')
        ax[0].axhline(5.5, color='k', linestyle='--')
        ax[0].axvline(5.5, color='k', linestyle='--')
        fig.colorbar(im0, ax=ax[0], fraction=0.046, pad=0.04)
        ax[0].set_title('Estimated Data G')

        G_pred, _ = M[0].predict(theta_group[0][:M[0].n_param])
        # visualize the predicted G:
        vmin = -np.max(np.abs(G_pred))
        vmax = np.max(np.abs(G_pred))
        im1 = ax[1].imshow(G_pred, vmin=vmin, vmax=vmax, cmap='RdBu_r')
        ax[1].axhline(5.5, color='k', linestyle='--')
        ax[1].axvline(5.5, color='k', linestyle='--')
        fig.colorbar(im1, ax=ax[1], fraction=0.046, pad=0.04)
        ax[1].set_title('Predicted G CKA')
        fig.tight_layout()
        plt.show()

    return T, theta, ceil