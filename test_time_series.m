%% Fit hrf params - redo GLM:
usr_path = userpath;
usr_path = usr_path(1:end-17);
baseDir = fullfile(usr_path,'Desktop','Projects','bimanual_wrist','data','fMRI');

% UCL:
% baseDir = fullfile(usr_path,'Desktop','Projects','bimanual_wrist','data','UCL');
sn = 115;
glm = 99;

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

regions = {'', 'S1_l', 'M1_l', 'PMd_l', 'PMv_l', 'SMA_l', 'V1_l', 'SPLa_l', 'SPLp_l',...
           '', 'S1_r', 'M1_r', 'PMd_r', 'PMv_r', 'SMA_r', 'V1_r', 'SPLa_r', 'SPLp_r'};

r = 3;
pre = 8;
post = 20;
avg_raw = mean([region_data{2},region_data{3}],2);
% Yraw = [region_data{2},region_data{3}]; % Pick out data as a T x P matrix
Yraw = avg_raw;


%% ============== FIT HRF ==============
function [p_opt, hrf_fit, fit_stats] = fit_spm_hrf_to_target(SPM, Yraw, pre, post)
    % Initial guess from SPM defaults:
    % SPM default p ~ [6 16 1 1 6 0 32];
    p0 = [6,10,1,1,6,0,32];
    
    % Reasonable bounds (physiological and numeric):
    % [delay1, delay2, disp1, disp2, ratio, onset, length]
    lb = [ 3.5;  8; 0.1; 0.1; 1; 0; 32];
    ub = [12; 20;  3;   3; 10; 3; 32];
    
    % Optimize only first 6 parameters;
    p0_free = p0(1:6);
    lb_free = lb(1:6);
    ub_free = ub(1:6);
    
    % Residual function: normalize both target and model for shape matching
    function r = resid(p_free)
        p = [p_free(:); 32];
        % fprintf('params = %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %d\n\n',p(1),p(2),p(3),p(4),p(5),p(6),p(7))
        
        % Get the hemodynamic response in micro-time resolution
        SPM.xBF.UNITS    = 'secs'; % units of the hrf
        SPM.xBF.T        = 16; % microtime resolution: number of time samples per scan
        SPM.xBF.dt       = SPM.xY.RT/SPM.xBF.T;  % Delta-t per sample 
        SPM.xBF.T0       = 1; % first time bin
        SPM.xBF.name     = 'fitted_hrf';
        SPM.xBF.order    = 1;
        SPM.xBF.Volterra = 1;  % volterra expansion order?
        SPM.xBF.bf = spm_hrf(SPM.xBF.dt, p);
        SPM.xBF.length = size(SPM.xBF.bf,1)*SPM.xBF.dt; % support in seconds 
        
        % Reconvolve the design matrix with new HRF
        SPM = spmj_glm_convolve(SPM);
        
        % Restimate the betas and get predicted and residual response
        [beta, Yhat, Yres] = spmj_glm_fit(SPM,Yraw);

        D = spmj_get_ons_struct(SPM);
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
            D.iti(idx(2:end)) = iti*SPM.xY.RT;
        end
        
        Yadj = Yres+Yhat;
        for i=1:size(D.block,1)
            D.y_adj(i,:)=cut(mean(Yadj,2,"omitmissing"),pre,round(D.ons(i)),post,'padding','nan')';
            D.y_hat(i,:)=cut(mean(Yhat,2,"omitmissing"),pre,round(D.ons(i)),post,'padding','nan')';
            D.y_res(i,:)=cut(mean(Yres,2,"omitmissing"),pre,round(D.ons(i)),post,'padding','nan')';
        end
        
        % select trials that come just before a rest gap:
        idx = find(D.iti>10)-1;
        Tgap = getrow(D, idx);
        % ASSUMING M1 LEFT HEM:
        idx1 = find(contains(Tgap.eventname, 'bi'));
        idx2 = find(contains(Tgap.eventname, 'rhand'));
        T = getrow(Tgap,[idx1;idx2]);
        fit_signal = mean(T.y_hat,1,'omitmissing');
        target_signal = mean(T.y_adj,1,'omitmissing');
        
        % compare the fit with target:
        % r = fit_signal - target_signal;
        r = mean(Yres,2);
        % r = mean(Yres,2);
    end
    
    % Options
    opts = optimoptions('lsqnonlin', ...
        'Display','iter', ...
        'MaxFunctionEvaluations',5e3, ...
        'MaxIterations',1e3, ...
        'StepTolerance',1e-12, ...
        'FunctionTolerance',1e-10, ...
        'MaxIterations',30);

    % Run optimization
    [p_free_opt, ~, residual, exitflag, output] = lsqnonlin(@resid, p0_free, lb_free, ub_free, opts);

    % Assemble full parameter vector
    p_opt = [p_free_opt(:); 32];
    
    % Compute fitted HRF
    hrf_fit = spm_hrf(SPM.xY.RT, p_opt);
    hrf_fit = hrf_fit(:);

    fit_stats = 0;
end

[p_opt, hrf_fit, fit_stats] = fit_spm_hrf_to_target(SPM, Yraw, pre, post);

disp('Optimized parameters [delay1, delay2, disp1, disp2, ratio, onset, length:');
fprintf('s%d:\n',sn)
fprintf('%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n',p_opt(1),p_opt(2),p_opt(3),p_opt(4),p_opt(5),p_opt(6),p_opt(7))




%%
r = 3;
% hrf_params = [6,10,1,1,6,0,32];
hrf_params = [5.0435,10.9274,0.9961,0.3557,2.7909,0.8963,32.0000];
pre = 8;
post = 20;
run = 4;

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

% =========================================================================
% =========================================================================
% figure;
% % plot(mean(Yhat,2)); 
% % hold on;
% ytmp = Yhat+Yres;
% plot(mean(ytmp,2))

% =========================================================================
% =========================================================================
% Diagnostic plots
figure;
SPM.xBF.bf = SPM.xBF.bf/max(SPM.xBF.bf);
t=[SPM.xBF.dt*SPM.xBF.T0:SPM.xBF.dt:SPM.xBF.length];
plot(t,SPM.xBF.bf);
drawline([0],'dir','horz','linestyle','-');
title(sprintf('params = %s',num2str(hrf_params)));

% =========================================================================
% =========================================================================
% Run, convolved design matrix, sum of regressors of interest
% figure;
% X=SPM.xX.X(SPM.Sess(run).row,SPM.Sess(run).col);
% plot(sum(X,2));
% title(sprintf('Run %d - overall response', run));

% =========================================================================
% =========================================================================
% figure;
% Yhat_run1 = mean(Yhat(SPM.Sess(run).row,:),2);
% Yres_run1 = mean(Yres(SPM.Sess(run).row,:),2);
% Yadj_run1 = Yres_run1+Yres_run1; 
% t= SPM.Sess(run).row;
% hold on;
% plot(t, Yhat_run1, 'r', 'LineWidth', 1.5)
% plot(t, Yadj_run1, 'k', 'LineWidth', 1)
% title(sprintf('Run %d - overall response', run));

% =========================================================================
% =========================================================================
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
traceplot([-pre:post]*SPM.xY.RT,T.y_adj,'errorfcn','stderr'); % ,
hold on;
traceplot([-pre:post]*SPM.xY.RT,T.y_hat,'linestyle',':',...
        'linewidth',3); %
drawline([-7,0,7,14,21],'dir','vert','linestyle',':');

drawline([0],'dir','vert','linestyle','--');
drawline([0],'dir','horz','linestyle','-');
hold off;
xlabel('Time');
ylabel('activation');
title(sprintf('region %s',regions{r}))

% =========================================================================
% =========================================================================
% FIND GAP TRIALS:
D = spmj_get_ons_struct(SPM);
D.iti = zeros(size(D.ons));
% sort based on onsets:
blocks = unique(D.block)';
for b = blocks
    % sorting:
    rows = D.block==b;

    ons = D.ons(rows); % these are the onset frames (volumes) not times (seconds)
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
    D.iti(idx(2:end)) = iti*SPM.xY.RT; % turn iti to time (seconds)
end

Yadj = Yres+Yhat;
for i=1:size(D.block,1)
    D.y_adj(i,:)=cut(mean(Yadj,2),pre,round(D.ons(i)),post,'padding','nan')';
    D.y_hat(i,:)=cut(mean(Yhat,2),pre,round(D.ons(i)),post,'padding','nan')';
    D.y_res(i,:)=cut(mean(Yres,2),pre,round(D.ons(i)),post,'padding','nan')';
end

% =========================================================================
% ============================ trials that come before gap:
idx = find(D.iti>10)-1;
Tgap = getrow(D, idx);

% select only specific condition:
idx = find(contains(Tgap.eventname, 'lhand'));
T = getrow(Tgap,idx);
figure;
subplot(3,1,1)
traceplot([-pre:post]*SPM.xY.RT,T.y_adj,'errorfcn','stderr'); % ,
hold on;
traceplot([-pre:post]*SPM.xY.RT,T.y_hat,'linestyle',':',...
        'linewidth',3); %
drawline([-7,0,12,19],'dir','vert','linestyle',':');

drawline([0],'dir','vert','linestyle','--');
drawline([0],'dir','horz','linestyle','-');
hold off;
xlabel('Time (sec)');
ylabel('activation');
title(sprintf('lhand, region %s',regions{r}))

% select only specific condition:
idx = find(contains(Tgap.eventname, 'rhand'));
T = getrow(Tgap,idx);

subplot(3,1,2)
traceplot([-pre:post]*SPM.xY.RT,T.y_adj,'errorfcn','stderr'); % ,
hold on;
traceplot([-pre:post]*SPM.xY.RT,T.y_hat,'linestyle',':',...
        'linewidth',3); %
drawline([-7,0,12,19],'dir','vert','linestyle',':');

drawline([0],'dir','vert','linestyle','--');
drawline([0],'dir','horz','linestyle','-');
hold off;
xlabel('TR');
ylabel('activation');
title(sprintf('rhand, region %s',regions{r}))

% select only specific condition:
idx = find(contains(Tgap.eventname, 'bi'));
T = getrow(Tgap,idx);
subplot(3,1,3)
traceplot([-pre:post]*SPM.xY.RT,T.y_adj,'errorfcn','stderr'); % ,
hold on;
traceplot([-pre:post]*SPM.xY.RT,T.y_hat,'linestyle',':',...
        'linewidth',3); %
drawline([-7,0,12,19],'dir','vert','linestyle',':');
drawline([0],'dir','vert','linestyle','--');
drawline([0],'dir','horz','linestyle','-');
hold off;
xlabel('TR');
ylabel('activation');
title(sprintf('bi, region %s',regions{r}))
hold on;
% plot(linspace(0,32,length(SPM.xBF.bf))-7,SPM.xBF.bf);
% plot(linspace(0,32,length(SPM.xBF.bf)),SPM.xBF.bf);
% xlim([-pre,post])

%=========================================================================
%  ============================ the trials that come after the gap:
idx = find(D.iti>10);
T_aftergap = getrow(D, idx);

idx = find(contains(T_aftergap.eventname, 'lhand'));
T = getrow(T_aftergap,idx);
figure;
subplot(3,1,1)
traceplot([-pre:post]*SPM.xY.RT,T.y_adj,'errorfcn','stderr'); % ,
hold on;
traceplot([-pre:post]*SPM.xY.RT,T.y_hat,'linestyle',':',...
        'linewidth',3); %
drawline([-12,0,7,14],'dir','vert','linestyle',':');

drawline([0],'dir','vert','linestyle','--');
drawline([0],'dir','horz','linestyle','-');
hold off;
xlabel('TR');
ylabel('activation');
title(sprintf('lhand, region %s',regions{r}))


idx = find(contains(T_aftergap.eventname, 'rhand'));
T = getrow(T_aftergap,idx);
subplot(3,1,2)
traceplot([-pre:post]*SPM.xY.RT,T.y_adj,'errorfcn','stderr'); % ,
hold on;
traceplot([-pre:post]*SPM.xY.RT,T.y_hat,'linestyle',':',...
        'linewidth',3); %
drawline([-12,0,7,14],'dir','vert','linestyle',':');

drawline([0],'dir','vert','linestyle','--');
drawline([0],'dir','horz','linestyle','-');
hold off;
xlabel('TR');
ylabel('activation');
title(sprintf('rhand, region %s',regions{r}))

% select only specific condition:
idx = find(contains(T_aftergap.eventname, 'bi'));
T = getrow(T_aftergap,idx);
subplot(3,1,3)
traceplot([-pre:post]*SPM.xY.RT,T.y_adj,'errorfcn','stderr'); % ,
hold on;
traceplot([-pre:post]*SPM.xY.RT,T.y_hat,'linestyle',':',...
        'linewidth',3); %
drawline([-12,0,7,14],'dir','vert','linestyle',':');
drawline([0],'dir','vert','linestyle','--');
drawline([0],'dir','horz','linestyle','-');
hold off;
xlabel('TR');
ylabel('activation');
title(sprintf('bi, region %s',regions{r}))
% hold on;
% plot(linspace(0,32,length(SPM.xBF.bf))-7,SPM.xBF.bf);
% plot(linspace(0,32,length(SPM.xBF.bf)),SPM.xBF.bf);
% xlim([-pre,post])


