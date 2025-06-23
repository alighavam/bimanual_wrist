%% Fit hrf params - redo GLM:
usr_path = userpath;
usr_path = usr_path(1:end-17);
baseDir = fullfile(usr_path,'Desktop','Projects','bimanual_wrist','data','fMRI');

% UCL:
% baseDir = fullfile(usr_path,'Desktop','Projects','bimanual_wrist','data','UCL');
sn = 102;
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
r = 3; % M1_L
% r = 12; % M1_R
hrf_params = [5 10 0.6 1 3 0 32];
% hrf_params = [4 10 0.6 0.1 2 0 32];
pre = 8;
post = 20;
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
regressor_names = {SPM.Vbeta.descrip};
lhand_idx = find(contains(regressor_names, 'lhand'));
rhand_idx = find(contains(regressor_names, 'rhand'));
bi_idx = find(contains(regressor_names, 'bi'));

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
    D.y_adj(i,:)=cut(mean(Yadj,2),pre,round(D.ons(i)),post,'padding','nan')';
    D.y_hat(i,:)=cut(mean(Yhat,2),pre,round(D.ons(i)),post,'padding','nan')';
    D.y_res(i,:)=cut(mean(Yres,2),pre,round(D.ons(i)),post,'padding','nan')';
end

T = D; % Get the first onset for each double
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


figure;
% FIND GAP TRIALS:
D = spmj_get_ons_struct(SPM);
% D.ons = D.ons - 1;
D.iti = zeros(size(D.ons));
% sort based on onsets:
blocks = unique(D.block)';
for b = blocks
    % sorting:
    rows = D.block==b;

    ons = D.ons(rows);
    event = D.event(rows);
    eventname = D.eventname(rows);
    num = D.num(rows);
    
    % sorting based on onset:
    [~, ix] = sort(ons);
    ons = ons(ix);
    event = event(ix);
    eventname = eventname(ix);
    num = num(ix);
    iti = diff(ons);

    % adding to dataframe:
    D.ons(rows) = ons;
    D.event(rows) = event;
    D.eventname(rows) = eventname;
    D.num(rows) = num;
    idx = find(rows);
    D.iti(idx(2:end)) = iti;
end

Yadj = Yres+Yhat;
for i=1:size(D.block,1)
    D.y_adj(i,:)=cut(mean(Yadj,2),pre,round(D.ons(i)),post,'padding','nan')';
    D.y_hat(i,:)=cut(mean(Yhat,2),pre,round(D.ons(i)),post,'padding','nan')';
    D.y_res(i,:)=cut(mean(Yres,2),pre,round(D.ons(i)),post,'padding','nan')';
end

% find gap trials:
idx = find(D.iti>10)-1;
Tgap = getrow(D, idx);

% select only specific condition:
idx = find(contains(Tgap.eventname, 'lhand'));
T = getrow(Tgap,idx);

subplot(3,1,1)
traceplot([-pre:post],T.y_adj,'errorfcn','stderr'); % ,
hold on;
traceplot([-pre:post],T.y_hat,'linestyle',':',...
        'linewidth',3); %
drawline([-7,0,12,19],'dir','vert','linestyle',':');

drawline([0],'dir','vert','linestyle','--');
drawline([0],'dir','horz','linestyle','-');
hold off;
xlabel('TR');
ylabel('activation');
title('lhand')

% select only specific condition:
idx = find(contains(Tgap.eventname, 'rhand'));
T = getrow(Tgap,idx);

subplot(3,1,2)
traceplot([-pre:post],T.y_adj,'errorfcn','stderr'); % ,
hold on;
traceplot([-pre:post],T.y_hat,'linestyle',':',...
        'linewidth',3); %
drawline([-7,0,12,19],'dir','vert','linestyle',':');

drawline([0],'dir','vert','linestyle','--');
drawline([0],'dir','horz','linestyle','-');
hold off;
xlabel('TR');
ylabel('activation');
title('rhand')

% select only specific condition:
idx = find(contains(Tgap.eventname, 'bi'));
T = getrow(Tgap,idx);

subplot(3,1,3)
traceplot([-pre:post],T.y_adj,'errorfcn','stderr'); % ,
hold on;
traceplot([-pre:post],T.y_hat,'linestyle',':',...
        'linewidth',3); %
drawline([-7,0,12,19],'dir','vert','linestyle',':');

drawline([0],'dir','vert','linestyle','--');
drawline([0],'dir','horz','linestyle','-');
hold off;
xlabel('TR');
ylabel('activation');
title('bi')

