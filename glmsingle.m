% Start fresh
clear
clc
close all

% Add GLMsingle to the MATLAB path (in case the user has not already done so).
GLMsingle_dir = '/Users/aghavamp/Desktop/Projects/GLMsingle';

addpath(fullfile(GLMsingle_dir, 'matlab'));
addpath(fullfile(GLMsingle_dir, 'matlab', 'utilities'));

% if the submodules were installed we try to add their code to the path
addpath(fullfile(GLMsingle_dir, 'matlab', 'fracridge', 'matlab'));

% check that the dependencies are in the path
tmp = which('fracridge.m');
if isempty(tmp)
  error('fracridge is missing. Please install from: https://github.com/nrdg/fracridge.git')
end

clear GLMsingle_dir;

usr_path = userpath;
usr_path = usr_path(1:end-17);
baseDir = fullfile(usr_path,'Desktop','Projects','bimanual_wrist','data','fMRI');

glmEstDir = 'glm';
behavDir = 'behavioural';
anatomicalDir = 'anatomicals';
imagingDir = 'imaging_data';
wbDir = 'surfaceWB';
regDir = 'ROI';

pinfo = dload(fullfile(baseDir,'participants.tsv'));

sn_list = [101,102,103,104,106,107,108];
sn_list = [101]
for sn = sn_list
    % setup
    % get participant row from participant.tsv
    participant_row = getrow(pinfo, pinfo.sn==sn);
    
    % get participant_id
    participant_id = participant_row.participant_id{1};
    
    SPM_folder  = fullfile(baseDir,'glm1',participant_id);
    TR          = 1;
    stimdur     = 0.1;
    
    % Name of directory to which outputs will be saved
    outputdir   = fullfile(baseDir, 'glmsingle'); 
    
    %
    SPM = load(fullfile(SPM_folder,'SPM.mat'));
    
    mask = niftiread(fullfile(baseDir,'glm1',participant_id,'mask.nii'));
    
    % load design matrix
    tr = SPM.xY.RT;
    design = cell(1,length(SPM.Sess));
    for zz=1:length(SPM.Sess)  % for each run
    
      ncond = length(SPM.Sess(zz).U);    % number of conditions
      nvol = length(SPM.Sess(zz).row);   % number of volumes
      
      design{zz} = zeros(nvol,ncond);
    
      for yy=1:length(SPM.Sess(zz).U)    % for each condition
        design{zz}(round(SPM.Sess(zz).U(yy).ons/tr)+1,yy) = 1;  % set all of the onsets
      end
    end

    % make session indicator
    sessionindicator = ones(1,length(SPM.Sess))+1;
    run_ses1 = str2double(split(participant_row.run_ses1{1},'.'));
    sessionindicator(run_ses1) = 1;
    
    % load fMRI data
    data = cell(1,length(SPM.Sess));
    fname = unique(struct2table(SPM.xY.VY).fname);
    for zz=1:length(fname)
        tmp = niftiread(fname{zz});
        % mask4d = repmat(mask, [1,1,1,size(tmp,4)]);
        % tmp = tmp .* uint16(mask4d);
        data{zz} = single(tmp);
    end
    
    % 
    opt = struct('wantmemoryoutputs',[0 0 0 1],'wantlibrary',1,'wantglmdenoise',0,'wantfracridge',0,'sessionindicator',sessionindicator,'chunknum',100000);
    [results] = GLMestimatesingletrial(design,data,stimdur,tr,fullfile(outputdir, participant_id),opt);

    % =========== Ali's stuff after glm single fit
    % copy subject glm mask to glmsingle direcotry:
    copyfile(fullfile(baseDir,'glm1',participant_id,'mask.nii'),fullfile(outputdir,participant_id,'mask.nii'));

    % Save betas as nifti:
    % load glmsingle model
    modelname = 'TYPEB_FITHRF.mat';
    m = load(fullfile(outputdir, participant_id, modelname));
    
    % get event onsets and sort chronological:
    D = spmj_get_ons_struct(SPM);
    D.ons = D.ons-1;
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
        D.iti(idx(2:end),1) = iti;
    end
    
    info_base = niftiinfo(fullfile(baseDir,'glm1',participant_id,'beta_0001.nii'));
    sz = size(m.modelmd);
    
    % outside of mask voxels:
    outside_of_mask = find(mask==0); 
    for i = 1:sz(4)
        % make nifti:
        info = info_base;
        info.Filename = [];
        info.Filemoddate = [];
        info.Filesize = [];
        descrip = sprintf('glmsingle:beta (%.4d) - Sn(%d) %s', i, D.block(i), D.eventname{i});
        info.Description = descrip;
        info.raw.descrip = descrip;
        nii = m.modelmd(:,:,:,i);
        nii(outside_of_mask) = nan;
        % nii(nii==0) = NaN;
        % save nifti:
        niftidir = fullfile(outputdir,participant_id);
        niftiwrite(nii,fullfile(niftidir,sprintf('beta_%.4d.nii',i)), info);
    end
    
    % make reginfo.tsv:
    reginfo = {};
    reginfo.sn = sn * ones(size(D.block));
    reginfo.run = D.block;
    reginfo.name = D.eventname;
    reginfo.ons = D.ons;
    dsave(fullfile(outputdir, participant_id, 'reginfo.tsv'),reginfo);
    
    % save R2:
    R2 = m.R2;
    info = info_base;
    info.Filename = [];
    info.Filemoddate = [];
    info.Filesize = [];
    descrip = 'glmsingle:R2 percent';
    info.Description = descrip;
    info.raw.descrip = descrip;
    niftidir = fullfile(outputdir,participant_id);
    niftiwrite(R2,fullfile(niftidir,'R2.nii'), info);
    
    % Make t-maps:
    % load reginfo:
    reginfo = dload(fullfile(outputdir, participant_id, 'reginfo.tsv'));
    
    % load betas:
    betafiles = dir(fullfile(outputdir, participant_id, 'beta*.nii'));
    beta = {};
    info = {};
    for i = 1:length(betafiles)
        beta{i,1} = niftiread(fullfile(betafiles(i).folder, betafiles(i).name));
        % beta{i,1} = beta{i,1} .* single(mask);
        info{i,1} = niftiinfo(fullfile(betafiles(i).folder, betafiles(i).name));
    end
    
    % estimate t-maps:
    conditions = unique(reginfo.name);
    for i = 1:length(conditions)
        c = conditions{i};
        idx = find(strcmp(reginfo.name,c));
        select_betas = beta(idx);
        cat4d = cat(4, select_betas{:});
        tstats = nanmean(cat4d,4) ./ (std(cat4d,[],4)./sqrt(length(idx)));
        tstats(isnan(tstats)) = 0;
        % save tstats:
        infotmp = info{1};
        infotmp.Filename = [];
        infotmp.Filemoddate = [];
        infotmp.Filesize = [];
        descrip = sprintf('t-stats:%s',conditions{i});
        infotmp.Description = descrip;
        infotmp.raw.descrip = descrip;
        niftidir = fullfile(outputdir,participant_id);
        niftiwrite(tstats,fullfile(niftidir,sprintf('tmap_%s.nii',replace(conditions{i},":", "-"))), infotmp);
    end 
end



%% Extract and save time series as cifti
clc
clear
close

usr_path = userpath;
usr_path = usr_path(1:end-17);
baseDir = fullfile(usr_path,'Desktop','Projects','bimanual_wrist','data','fMRI');

% Add GLMsingle to the MATLAB path (in case the user has not already done so).
GLMsingle_dir = '/Users/aghavamp/Desktop/Projects/GLMsingle';

addpath(fullfile(GLMsingle_dir, 'matlab'));
addpath(fullfile(GLMsingle_dir, 'matlab', 'utilities'));

% if the submodules were installed we try to add their code to the path
addpath(fullfile(GLMsingle_dir, 'matlab', 'fracridge', 'matlab'));

clear GLMsingle_dir;

glmEstDir = 'glm';
behavDir = 'behavioural';
anatomicalDir = 'anatomicals';
imagingDir = 'imaging_data';
wbDir = 'surfaceWB';
regDir = 'ROI';
outputdir = fullfile(baseDir, 'glmsingle');
pinfo = dload(fullfile(baseDir,'participants.tsv'));

sn_list = [101, 102, 103, 104, 106, 107, 108];
T = [];
for sn = sn_list
    participant_row = getrow(pinfo, pinfo.sn==sn);
    % get participant_id
    participant_id = participant_row.participant_id{1};
    
    SPM_folder  = fullfile(baseDir,'glm1',participant_id);
    SPM = load(fullfile(SPM_folder,'SPM.mat'));
    reginfo = dload(fullfile(outputdir, participant_id, 'reginfo.tsv'));
    modelname = 'TYPEB_FITHRF.mat';
    m = load(fullfile(outputdir, participant_id, modelname));
    mask = niftiread(fullfile(baseDir,'glm1',participant_id,'mask.nii'));
    designinfo = load(fullfile(outputdir, participant_id, 'DESIGNINFO.mat'));
    R = load(fullfile(baseDir, regDir, participant_id, sprintf('%s_%s_glm%d_region.mat', participant_id, 'ROI', 1))); R=R.R;
    
    % get region voxels, "assuming all images have the SAME AFFINE":
    X=[];Y=[];Z=[];
    for r=1:length(R)
        if (~isempty(R{r}))
            [x,y,z]=spmj_affine_transform(R{r}.data(:,1),R{r}.data(:,2),R{r}.data(:,3),inv(SPM.xY.VY(1).mat));
            from(r)=size(X,1)+1;
            X=[X;x];Y=[Y;y];Z=[Z;z];
            to(r)=size(X,1);
        else 
            from(r)=size(X,1)+1;
            to(r)=size(X,1);
        end 
    end
    
    stimdur = 0.1;
    tr = 1;
    typeA = load(fullfile(outputdir, participant_id, 'TYPEA_ONOFF.mat'));
    typeD = load(fullfile(outputdir, participant_id, 'TYPED_FITHRF_GLMDENOISE_RR.mat'));
    temp = load(fullfile(outputdir, participant_id, 'DESIGNINFO.mat'));
    designSINGLE = temp.designSINGLE;
    
    % load fMRI data
    mask = niftiread(fullfile(outputdir,participant_id,'mask.nii'));
    SPM_folder  = fullfile(baseDir,'glm1',participant_id);
    SPM = load(fullfile(SPM_folder,'SPM.mat'));
    data = cell(1,length(SPM.Sess));
    fname = unique(struct2table(SPM.xY.VY).fname);
    for zz=1:length(fname)
        tmp = niftiread(fname{zz});
        mask4d = repmat(mask, [1,1,1,size(tmp,4)]);
        tmp = tmp .* uint16(mask4d);
        data{zz} = tmp;
    end
    
    % ============= MAKE THE Y_Adj =============
    y_adj = zeros(length(X),7400);
    for ii=1:length(designSINGLE)
        data_run = data{ii};
    
        % --- A. Get the GLMdenoise Principal Components for this run ---
        % 'pcregressors' is a cell array, one cell per run
        noise_pcs = typeD.pcregressors{ii}; 
        
        % The number of PCs used in the final model is stored in 'pcnum'
        num_pcs_to_use = typeD.pcnum;
        
        % Keep only the PCs that were actually used in the model
        if num_pcs_to_use > 0
            selected_pcs = noise_pcs(:, 1:num_pcs_to_use);
        else
            selected_pcs = []; % No PCs were selected
        end
        
        % --- B. Generate the Polynomial Regressors ---
        num_timepoints = size(data_run, 4);
        run_duration_min = (num_timepoints * tr) / 60;
        poly_degree = round(run_duration_min / 2);
        
        % Create the polynomial regressors (including the constant term)
        polynomial_regressors = zeros(num_timepoints, poly_degree + 1);
        for p = 0:poly_degree
            polynomial_regressors(:, p+1) = (linspace(-1, 1, num_timepoints)').^p;
        end
        
        % --- C. Combine into a single nuisance regressor matrix ---
        % Also add any 'extraregressors' you might have used
        extra_regs = []; % Populate this if you used opt.extraregressors
        nuisance_regressors = [selected_pcs, polynomial_regressors, extra_regs];
        
        volume_size = [size(data_run, 1), size(data_run, 2), size(data_run, 3)];
        linear_indices = sub2ind(volume_size, round(X), round(Y), round(Z));
        
        num_voxels_total = size(data_run, 1) * size(data_run, 2) * size(data_run, 3); % 90*90*32 = 259200
        num_timepoints = size(data_run, 4); % 740
        data_reshaped = reshape(data_run, num_voxels_total, num_timepoints);
        
        selected_data = double(data_reshaped(linear_indices, :)); % voxels x time data
        
        betas_nuisance = nuisance_regressors \ selected_data'; 
        y_noise = nuisance_regressors * betas_nuisance;
        % y_adj = [y_adj (selected_data' - y_noise)'];
        y_adj(:,(ii-1)*740+1:ii*740) = (selected_data' - y_noise)';
    end
    
    % ============= MAKE THE PREDICTIONS =============
    % --- 1. Vectorized Pre-computation ---
    fprintf('Starting vectorized processing...\n');
    
    fprintf('Generating HRF library once...\n');
    hrflibrary = getcanonicalhrflibrary(stimdur, tr)'; % timepoints x HRFs
    
    % Get volume dimensions.
    vol_size = size(typeD.modelmd);
    vol_size = vol_size(1:3);
    
    % Convert (X,Y,Z) coordinates into linear indices
    fprintf('Converting voxel coordinates to linear indices...\n');
    linear_indices = sub2ind(vol_size, round(X), round(Y), round(Z));
    
    % Extract all necessary data for ALL voxels at once using these indices:
    fprintf('Extracting model data for all voxels...\n');
    num_voxels = length(X);
    num_betas = size(typeD.modelmd, 4);
    
    % Reshape the 4D beta matrix into 2D (all_voxels_in_volume x betas)
    % and then select the voxels we want:
    all_betas_flat = reshape(typeD.modelmd, [], num_betas);
    all_betas = all_betas_flat(linear_indices, :); % [num_voxels, num_betas]
    
    % Extract the other per-voxel information.
    all_hrf_indices = typeD.HRFindex(linear_indices);
    all_meansignal = typeD.meanvol(linear_indices);
    
    % Pre-allocate memory for the final prediction matrix.
    num_sessions = length(designSINGLE);
    timepoints_per_session = 740; % Assuming this is fixed
    num_timepoints_total = num_sessions * timepoints_per_session;
    pred = zeros(num_voxels, num_timepoints_total);
    
    % --- 2. Grouped and Vectorized Processing ---
    % Find the unique HRF indices that are actually used by oir set of voxels.
    unique_hrfs = unique(all_hrf_indices);
    fprintf('Found %d unique HRFs to process. Starting main computation...\n', length(unique_hrfs));
    
    % Loop over the list of unique HRFs
    for h_idx = 1:length(unique_hrfs)
        current_hrf_index = unique_hrfs(h_idx);
        fprintf('Processing group for HRF index %d/%d...\n', h_idx, length(unique_hrfs));
        
        % Create a boolean mask to identify all voxels that use the current HRF.
        voxel_mask_for_hrf = (all_hrf_indices == current_hrf_index);
        
        % If no voxels use this HRF for some reason, skip.
        if ~any(voxel_mask_for_hrf)
            continue;
        end
        
        % Select the data for this specific group of voxels using the mask.
        betas_group = all_betas(voxel_mask_for_hrf, :);
        meansignal_group = all_meansignal(voxel_mask_for_hrf);
        
        % Loop through the sessions
        for ii = 1:num_sessions
            
            % Convolve the design with the HRF for this group
            design0 = conv2(designSINGLE{ii}, hrflibrary(:, current_hrf_index));
            design0 = design0(1:size(designSINGLE{ii}, 1), :); % Truncate
            
            % --- VECTORIZED PREDICTION CALCULATION FOR THE ENTIRE GROUP ---
            
            % Scale betas by mean signal. This is now a matrix operation.
            % The result is a (num_betas x num_voxels_in_group) matrix.
            betatemp_group = (betas_group' / 100) .* meansignal_group';
            
            % calculating the predicted time series for ALL voxels in the group at once.
            Xb_group = design0 * betatemp_group; % (timepoints x voxels_in_group)
            
            % Center the data for all voxels in the group
            Xb_group_centered = Xb_group - mean(Xb_group, 1);
            
            % store in final 'pred' matrix.
            time_indices = (ii-1)*timepoints_per_session + 1 : ii*timepoints_per_session;
            pred(voxel_mask_for_hrf, time_indices) = Xb_group_centered'; % Transpose to fit (voxels x time)
        end
    end
    
    % ROI summary:
    T_tmp = [];
    roi_name = ["...._L","..S1_L","..M1_L",".PMd_L",".PMv_L",".SMA_L","..V1_L","SPLa_L","SPLp_L","...._R","..S1_R","..M1_R",".PMd_R",".PMv_R",".SMA_R","..V1_R","SPLa_R","SPLp_R"];
    for i = 1:length(from)
        y_adj_roi = y_adj(from(i):to(i),:)';
        y_hat_roi = pred(from(i):to(i),:)';
        tmp = [];
        tmp.region = repmat(char(roi_name(i)),timepoints_per_session*length(designSINGLE),1);
        tmp.run = repelem(1:length(designSINGLE),timepoints_per_session)';
        tmp.y_adj = mean(y_adj_roi,2);
        tmp.y_hat = mean(y_hat_roi,2);
        tmp.y_res = mean(y_adj_roi - y_hat_roi,2);
        tmp.sn = repmat(sn,timepoints_per_session*length(designSINGLE),1);
        T_tmp = addstruct(T_tmp,tmp,'row','force');
    end
    
    T = addstruct(T,T_tmp,'row','force');
end
    
dsave('./analysis/timeseries_uwo.tsv',T);


%% find good hrf
clc
clear
close

usr_path = userpath;
usr_path = usr_path(1:end-17);
baseDir = fullfile(usr_path,'Desktop','Projects','bimanual_wrist','data','fMRI');

% Add GLMsingle to the MATLAB path (in case the user has not already done so).
GLMsingle_dir = '/Users/aghavamp/Desktop/Projects/GLMsingle';

addpath(fullfile(GLMsingle_dir, 'matlab'));
addpath(fullfile(GLMsingle_dir, 'matlab', 'utilities'));

% if the submodules were installed we try to add their code to the path
addpath(fullfile(GLMsingle_dir, 'matlab', 'fracridge', 'matlab'));

clear GLMsingle_dir;

glmEstDir = 'glm';
behavDir = 'behavioural';
anatomicalDir = 'anatomicals';
imagingDir = 'imaging_data';
wbDir = 'surfaceWB';
regDir = 'ROI';
outputdir = fullfile(baseDir, 'glmsingle');
pinfo = dload(fullfile(baseDir,'participants.tsv'));

tr = 1;
stim_duration = 0.1;
sn_list = [101, 102, 103, 104, 106, 107, 108];
sn_list = [115];
T = [];
for sn = sn_list
    participant_row = getrow(pinfo, pinfo.sn==sn);
    % get participant_id
    participant_id = participant_row.participant_id{1};
    
    SPM_folder  = fullfile(baseDir,'glm1',participant_id);
    SPM = load(fullfile(SPM_folder,'SPM.mat'));
    modelname = 'TYPEB_FITHRF.mat';
    m = load(fullfile(outputdir, participant_id, modelname));
    mask = niftiread(fullfile(baseDir,'glm1',participant_id,'mask.nii'));
    designinfo = load(fullfile(outputdir, participant_id, 'DESIGNINFO.mat'));
    R = load(fullfile(baseDir, regDir, participant_id, sprintf('%s_%s_glm%d_region.mat', participant_id, 'ROI', 1))); R=R.R;
    
    % get region voxels, "assuming all images have the SAME AFFINE":
    X=[];Y=[];Z=[];
    for r=1:length(R)
        if (~isempty(R{r}))
            [x,y,z]=spmj_affine_transform(R{r}.data(:,1),R{r}.data(:,2),R{r}.data(:,3),inv(SPM.xY.VY(1).mat));
            from(r)=size(X,1)+1;
            X=[X;x];Y=[Y;y];Z=[Z;z];
            to(r)=size(X,1);
        else 
            from(r)=size(X,1)+1;
            to(r)=size(X,1);
        end 
    end
end


%% visualize mask:
region_names = ['','S1_l','M1_l','PMd_l','PMv_l','SMA_l','V1_l','SPLa_l','SPLp_l',...
    '','S1_r','M1_r','PMd_r','PMv_r','SMA_r','V1_r','SPLa_r','SPLp_r'];
roi = 3;
x = round(X(from(roi):to(roi)));
y = round(Y(from(roi):to(roi)));
z = round(Z(from(roi):to(roi)));

mask_smooth = smooth3(double(mask), 'box', 3);

% ROI:
indices = sub2ind(size(mask), x, y, z);
ROI = zeros(size(mask));
ROI(indices) = 1;
ROI_smooth = smooth3(double(ROI), 'box', 3);

% 2. Create the figure and plot the first surface
figure;
% Generate isosurface data for the first volume at an isovalue of 2
s1 = isosurface(mask_smooth, 0.5); 
% Plot the first surface as a patch and set its color to red
p1 = patch(s1, 'FaceColor', 'red', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
isonormals(mask_smooth, p1); % Compute normals for lighting

% 3. Hold the plot to add the second surface
hold on;

% 4. Plot the second surface
% Generate isosurface data for the second volume at an isovalue of 1
s2 = isosurface(ROI_smooth, 0.5);
% Plot the second surface as a patch and set its color to blue
p2 = patch(s2, 'FaceColor', 'blue', 'EdgeColor', 'none');
isonormals(ROI_smooth, p2); % Compute normals for lighting

% 5. Finalize the plot appearance
title('Multiple Isosurfaces with Different Colors');
axis equal;    % Set aspect ratio
grid on;       % Add a grid
view(3);       % Set to 3D view

% Add lighting to make the surfaces look 3D
camlight;
lighting gouraud;

% Release the hold on the figure
hold off;


% Get HRF index for the ROI voxels:
hrfidx = m.HRFindex;
roi_hrfs = hrfidx(indices);
figure;
hist(roi_hrfs)

% visualize the hrf shape:
idx = mode(roi_hrfs);
hrflibrary = getcanonicalhrflibrary(stim_duration, tr);
target_hrf = hrflibrary(idx,:);

time_points_sec = (0:size(target_hrf, 2)-1) * tr;
% 5. Create the plot
figure;
plot(time_points_sec, target_hrf, 'k-', 'LineWidth', 2);
title(['HRF from Library (Index ', num2str(idx), ')']);
xlabel('Time (seconds)');
ylabel('Normalized Amplitude');
grid on;
axis tight; % Adjust axes to fit the data

% Add a vertical line at time = 0 for reference
hold on;
plot([0 0], ylim, 'k--'); 
hold off;

%% find good hrf


% Fit SPM HRF parameters to target_hrf (length 51) with TR=1, fixing p(7)=50.
% Requires: SPM on path for spm_hrf, and Optimization Toolbox for lsqnonlin.

function [p_opt, hrf_fit, fit_stats] = fit_spm_hrf_to_target(target_hrf)
    % Validate input
    if ~isvector(target_hrf) || length(target_hrf) ~= 51
        error('target_hrf must be a vector of length 51.');
    end
    target_hrf = target_hrf(:);  % column

    TR = 1;

    % Initial guess from SPM defaults:
    % SPM default p ~ [6 16 1 1 6 0 32]; here we fix p(7)=50 (length)
    p0 = [6; 16; 1; 1; 6; 0; 50];

    % Reasonable bounds (physiological and numeric):
    % [delay1, delay2, disp1, disp2, ratio, onset, length]
    lb = [ 1;  8; 0.1; 0.1; 1; 0; 50];
    ub = [12; 100;  4;   4; 10; 5; 50];

    % Optimize only first 6 parameters; p(7) fixed at 50.
    p0_free = p0(1:6);
    lb_free = lb(1:6);
    ub_free = ub(1:6);

    % Residual function: normalize both target and model for shape matching
    function r = resid(p_free)
        p = [p_free(:); 50];
        h = spm_hrf(TR, p);           % requires SPM
        h = h(:);
        % Ensure exactly 51 samples (SPM returns length p(7)+1 with TR=1)
        if length(h) ~= 51
            % Adjust by trimming/padding if needed (defensive)
            if length(h) > 51
                h = h(1:51);
            else
                h = [h; zeros(51 - length(h), 1)];
            end
        end
        % Normalize to unit L2 norm to compare shapes independent of scale
        h = h / max(norm(h), eps);
        t = target_hrf / max(norm(target_hrf), eps);
        r = h - t;
    end

    % Options
    opts = optimoptions('lsqnonlin', ...
        'Display','iter', ...
        'MaxFunctionEvaluations',5e3, ...
        'MaxIterations',1e3, ...
        'StepTolerance',1e-12, ...
        'FunctionTolerance',1e-10);

    % Run optimization
    [p_free_opt, ~, residual, exitflag, output] = lsqnonlin(@resid, p0_free, lb_free, ub_free, opts);

    % Assemble full parameter vector
    p_opt = [p_free_opt(:); 50];
    
    % Compute fitted HRF
    hrf_fit = spm_hrf(TR, p_opt);
    hrf_fit = hrf_fit(:);
    if length(hrf_fit) > 51, hrf_fit = hrf_fit(1:51); end

    % Return simple fit stats
    t = target_hrf / max(norm(target_hrf), eps);
    h = hrf_fit / max(norm(hrf_fit), eps);
    fit_stats.RMSE = sqrt(mean((h - t).^2));
    fit_stats.Corr = corr(h, t);
    fit_stats.exitflag = exitflag;
    fit_stats.output = output;
    fit_stats.initial_p = p0;
end

% 5.2744-32-1.2192-4-10-2-50

[p_opt, hrf_fit, fit_stats] = fit_spm_hrf_to_target(target_hrf);

disp('Optimized parameters [delay1, delay2, disp1, disp2, ratio, onset, length]:');
disp(p_opt.');

% figure;
% Plot the fitted HRF for visual comparison
figure;
plot(time_points_sec, hrf_fit/max(hrf_fit), 'k', 'LineWidth', 2, 'DisplayName', 'Fitted HRF'); hold on;
plot(time_points_sec, target_hrf, 'r--', 'LineWidth', 2, 'DisplayName', 'Target HRF');
legend('show', 'Location', 'best');

fprintf('%.4f-%.4f-%.4f-%.4f-%.4f-%.4f-%.4f\n',p_opt(1),p_opt(2),p_opt(3),p_opt(4),p_opt(5),p_opt(6),p_opt(7))


