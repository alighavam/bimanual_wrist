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
    tr = SPM.xY.RT;
    
    mask = niftiread(fullfile(baseDir,'glm1',participant_id,'mask.nii'));
    
    % load design matrix
    design = cell(1,length(SPM.Sess));
    for zz=1:length(SPM.Sess)  % for each run
    
      ncond = length(SPM.Sess(zz).U);    % number of conditions
      nvol = length(SPM.Sess(zz).row);   % number of volumes
      
      design{zz} = zeros(nvol,ncond);
    
      for yy=1:length(SPM.Sess(zz).U)    % for each condition
        design{zz}(round(SPM.Sess(zz).U(yy).ons/tr)+1,yy) = 1;  % set all of the onsets
      end
    end
    
    % load fMRI data
    data = cell(1,length(SPM.Sess));
    fname = unique(struct2table(SPM.xY.VY).fname);
    for zz=1:length(fname)
        tmp = niftiread(fname{zz});
        mask4d = repmat(mask, [1,1,1,size(tmp,4)]);
        tmp = tmp .* uint16(mask4d);
        data{zz} = tmp;
    end
    
    % 
    opt = struct('wantmemoryoutputs',[0 0 0 1]);
    [results] = GLMestimatesingletrial(design,data,stimdur,tr,fullfile(outputdir, participant_id),opt);

    % copy subject glm mask to glmsingle direcotry:
    copyfile(fullfile(baseDir,'glm1',participant_id,'mask.nii'),fullfile(outputdir,participant_id,'mask.nii'));

    
    % Save betas as nifti:
    % load glmsingle model
    modelname = 'TYPED_FITHRF_GLMDENOISE_RR.mat';
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
        beta{i,1} = beta{i,1} .* single(mask);
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


%% Define tmaps
