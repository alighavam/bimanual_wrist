betafiles = dir(fullfile('/Users/aghavamp/Desktop/Projects/bimanual_wrist/data/fMRI/glmsingle/s101', 'beta_*.nii'));
mask = niftiread(fullfile('/Users/aghavamp/Desktop/Projects/bimanual_wrist/data/fMRI/glmsingle/s101','mask.nii'));

beta_single = {};
info = {};
for i = 1:length(betafiles)
    beta_single{i,1} = niftiread(fullfile(betafiles(i).folder, betafiles(i).name));
    info{i,1} = niftiinfo(fullfile(betafiles(i).folder, betafiles(i).name));
end

beta = zeros(90,90,32,960);
for i = 1:length(beta_single)
    beta(:,:,:,i) = beta_single{i};
end
beta_flat = reshape(beta,[],960);


betafiles = dir(fullfile('/Users/aghavamp/Desktop/Projects/bimanual_wrist/data/fMRI/glm1/s101', 'beta_*.nii'));
beta_spm = {};

for i = 1:length(betafiles)
    beta_spm{i,1} = niftiread(fullfile(betafiles(i).folder, betafiles(i).name));
    % beta{i,1} = beta{i,1} .* single(mask);
end

spm = zeros(90,90,32,480);
for i = 1:480
    spm(:,:,:,i) = beta_spm{i};
end
spm_flat = reshape(spm,[],480);





%%

idx1 = find(~isnan(beta_flat(:,1)));
idx2 = find(~isnan(spm_flat(:,1)));

% voxel cov:
tmp1 = beta_flat(idx1); % P by 960
vox = randi(size(tmp1,1),1000,1);
tmp1 = tmp1(vox,:);
figure; imagesc(tmp1 * tmp1'); clim([-5,5])


tmp2 = spm_flat(idx2); % P by 960
tmp2 = tmp2(vox,:);
figure; imagesc(tmp2 * tmp2'); clim([-10000,10000])

%%

affine = info{1}.Transform.T;

% Define patch size (cubic side length for ~1000 voxels)
side_length = 20;  % 10x10x10 = 1000 voxels

% Example center voxel (adjust based on your brain region of interest)
center = [56, 47, 19];

% Calculate start/end indices (1-based in MATLAB), clamping to array bounds
sz = size(beta);
start_idx = max(1, center - floor(side_length / 2));
end_idx = min(sz(1:3), start_idx + side_length - 1);

% Extract the patch (shape e.g., [10 10 10 960] if center is interior)
glm_patch = beta(start_idx(1):end_idx(1), ...
                  start_idx(2):end_idx(2), ...
                  start_idx(3):end_idx(3), :);

sz = size(spm);
start_idx = max(1, center - floor(side_length / 2));
end_idx = min(sz(1:3), start_idx + side_length - 1);
spm_patch = spm(start_idx(1):end_idx(1), ...
                  start_idx(2):end_idx(2), ...
                  start_idx(3):end_idx(3), :);

glm_flat = reshape(glm_patch,[],960);
spm_flat = reshape(spm_patch,[],480);


% voxel cov:
% tmp1 = glm_flat; % P by 960
% figure; imagesc(tmp1 * tmp1'); clim([-500,500])
% 
% 
% tmp2 = spm_flat; % P by 480
% figure; imagesc(tmp2 * tmp2'); clim([-500000,500000])
% 
% reginfo = dload('/Users/aghavamp/Desktop/Projects/bimanual_wrist/data/fMRI/glmsingle/s101/reginfo.tsv');
% name = reginfo.name;
% idx_reg = contains(name,'lhand');
% 
% % beta cov
% figure; imagesc(tmp1' * tmp1); clim([-500,500]);
% hold on;
% xline([12,24,36])
% 
% tmp2 = spm_flat; % P by 480
% figure; imagesc(tmp2' * tmp2); clim([-5,5]*1e2);



figure; imagesc(glm_flat); clim([-4,4])
xline([96:96:960])

figure; imagesc(spm_flat); clim([-2,2]*1e2)
xline([48:48:480])

