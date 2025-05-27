old = [2,4,5,6,7,9,10];
new = [109,110,111,112,113,114,115];
glm = 1;
for i = 1:length(old)
    glmDir = fullfile(baseDir, [glmEstDir num2str(glm)]);
    spm_file_path = fullfile(glmDir, sprintf('s%d',new(i)), 'SPM.mat');
    old_path = ['/Users/aghavamp/Desktop/Projects/bimanual_wrist/data/UCL/glm1/',num2str(old(i),'%.2d')];
    new_path = fullfile(glmDir, sprintf('s%d',new(i)));
    spm_changepath(spm_file_path,old_path,new_path);
end

