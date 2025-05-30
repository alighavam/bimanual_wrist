%% Fit hrf params - redo GLM:
usr_path = userpath;
usr_path = usr_path(1:end-17);
% baseDir = fullfile(usr_path,'Desktop','Projects','bimanual_wrist','data','fMRI');

% UCL:
baseDir = fullfile(usr_path,'Desktop','Projects','bimanual_wrist','data','UCL');
sn = 109;
glm = 1;

pinfo = dload(fullfile(baseDir,'participants.tsv'));

% get participant row from participant.tsv
participant_row = getrow(pinfo, pinfo.sn==sn);

% get subj_id
participant_id = participant_row.participant_id{1};

% spm_file = load(['/Users/alighavampour/Desktop/Projects/bimanual_wrist/data/fMRI/glm' num2str(glm) '/s' num2str(sn,'%.2d') '/SPM.mat']);
spm_file = load(fullfile(baseDir,['glm' num2str(glm)], participant_id, 'SPM.mat'));
SPM = spm_file;
% load ROI definition (R)
R = load(fullfile(baseDir, 'ROI', participant_id, sprintf('%s_ROI_glm%d_region.mat', participant_id, glm))); 
R=R.R;

region_data = region_getdata(SPM.xY.VY,R);
% region_data = load('/Users/ali/Desktop/Projects/bimanual_wrist/data/region_data_s04.mat'); region_data=region_data.region_data;

%%
r = 2;
hrf_params = [6 16 1 1 6 0 32];
pre = 8;
post = 22;
run = 1;

Yraw = region_data{r}; % Pick out data as a T x P matrix 

% Get the hemodynamic response in micro-time resolution
SPM.xBF.UNITS    = 'secs'; % units of the hrf
SPM.xBF.T        = 16; % microtime resolution: number of time samples per scan
SPM.xBF.dt       = SPM.xY.RT/SPM.xBF.T;  % Delta-t per sample 
SPM.xBF.T0       = 1; % first time bin
SPM.xBF.name     = 'fitted_hrf';
SPM.xBF.order    = 1;
SPM.xBF.Volterra = 1;  % volterra expansion order?
SPM.xBF.bf = spm_hrf(SPM.xBF.dt, hrf_params);
SPM.xBF.length = size(SPM.xBF.bf,1)*SPM.xBF.dt; % support in seconds 

% Reconvolve the design matrix with new HRF
SPM = spmj_glm_convolve(SPM);

% Restimate the betas and get predicted and residual response
[beta, Yhat, Yres] = spmj_glm_fit(SPM,Yraw);

% Diagnostic plots
figure;
t=[SPM.xBF.dt*SPM.xBF.T0:SPM.xBF.dt:SPM.xBF.length]; 
plot(t,SPM.xBF.bf);
drawline([0],'dir','horz','linestyle','-');
title(sprintf('params = %s',num2str(hrf_params)));

% Run, convolved design matrix, sum of regressors of interest
figure;
X=SPM.xX.X(SPM.Sess(run).row,SPM.Sess(run).col);
plot(sum(X,2));
title(sprintf('Run %d - overall response', run));

figure;
Yhat_run1 = mean(Yhat(SPM.Sess(run).row,:),2);
Yres_run1 = mean(Yres(SPM.Sess(run).row,:),2);
Yadj_run1 = Yres_run1+Yres_run1; 
t= SPM.Sess(run).row;
hold on;
plot(t, Yhat_run1, 'r', 'LineWidth', 1.5)
plot(t, Yadj_run1, 'k', 'LineWidth', 1)
title(sprintf('Run %d - overall response', run));

% Get onset structure, cut-out the trials of choice, and plot evoked
% response
figure;
D = spmj_get_ons_struct(SPM);
Yadj = Yres+Yhat; 
for i=1:size(D.block,1)
    D.y_adj(i,:)=cut(mean(Yadj,2),pre,round(D.ons(i))-1,post,'padding','nan')';
    D.y_hat(i,:)=cut(mean(Yhat,2),pre,round(D.ons(i))-1,post,'padding','nan')';
    D.y_res(i,:)=cut(mean(Yres,2),pre,round(D.ons(i))-1,post,'padding','nan')';
end

T = getrow(D,mod(D.num,2)==1); % Get the first onset for each double 
traceplot([-pre:post],T.y_adj,'errorfcn','stderr'); % ,
hold on;
traceplot([-pre:post],T.y_hat,'linestyle',':',...
        'linewidth',3); %
% drawline([-7,0,7,14,21],'dir','vert','linestyle',':');

drawline([0],'dir','vert','linestyle','--');
drawline([0],'dir','horz','linestyle','-');
hold off;
xlabel('TR');
ylabel('activation');
















