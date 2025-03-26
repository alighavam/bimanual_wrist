sn = 7;
T = load(['/Users/alighavampour/Desktop/Projects/bimanual_wrist/data/fMRI/ROI/s' num2str(sn,'%.2d') '/time_series_glm1.mat']);
spm_file = load(['/Users/alighavampour/Desktop/Projects/bimanual_wrist/data/fMRI/glm1/s' num2str(sn,'%.2d') '/SPM.mat']);
SPM = spm_file.SPM;

%%
runs = 1:10;
roi = 2;

for run = runs
    x = 1:740;
    y_raw = T.y_raw(SPM.Sess(run).row,roi);
    y_adj = T.y_adj(SPM.Sess(run).row,roi);
    y_hat = T.y_hat(SPM.Sess(run).row,roi);
    y_res = T.y_res(SPM.Sess(run).row,roi);

    figure;
    % plot(x, y_raw, 'k', 'linewidth', 1.5); 
    % plot(x, y_res, 'k', 'linewidth', 1.5); 
    hold on;
    plot(x, y_adj, 'k', 'linewidth', 1.5); 
    plot(x, y_hat, 'r', 'linewidth', 1.5);
    title('run %d',run)
end