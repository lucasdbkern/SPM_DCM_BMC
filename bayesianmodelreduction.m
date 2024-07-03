% function bayesianmodelreduction(GCM)
% 
% % GCM - NxM structure array or filename of a .mat file containing it
% % field - parameter fields in DCM{i}.Ep to plot, or the fields to search
% % if only one DCM is provided per subject [default: {'A','B'}]
% 
% % Check if GCM is a filename, if so, load the data
% if ischar(GCM) || isstring(GCM)
%     loaded_data = load(GCM);
%     GCM = loaded_data.GCM;
% end
% 
% n = size(GCM, 1);
% 
% % Run spm_dcm_bmr for each subject
% for i = 1:n
%     fprintf('Running spm_dcm_bmr for subject %d...\n', i);
%     G = GCM(i, :);
% 
% 
%     % Run spm_dcm_bmr for this subject
%     [RCM, BMC, BMA] = spm_dcm_bmr(G);
% 
%     % Save the results
%     save_path = fullfile('subject_folder', sprintf('subject_%d', i), 'reduced');
%     if ~exist(save_path, 'dir')
%         mkdir(save_path);
%     end
%     save(fullfile(save_path, 'bmr_results.mat'), 'RCM', 'BMC', 'BMA');
% end
% 


function bayesianmodelreduction(GCM)
% GCM - NxM structure array or filename of a .mat file containing it
% field - parameter fields in DCM{i}.Ep to plot, or the fields to search
% if only one DCM is provided per subject [default: {'A','B'}]

% Check if GCM is a filename, if so, load the data
if ischar(GCM) || isstring(GCM)
    loaded_data = load(GCM);
    GCM = loaded_data.GCM;
    %GCM = cell2struct(GCM);
end

n = size(GCM, 1);

% Run spm_dcm_bmr for each subject
for i = 1:n
    fprintf('Running spm_dcm_bmr for subject %d...\n', i);
    G = GCM(i, :);
  
    

    % Run spm_dcm_bmr for this subject
    [RCM, BMC, BMA] = spm_dcm_bmr(G);

    % Save the results
    save_path = fullfile('subject_folder', sprintf('subject_%d', i), 'reduced');
    if ~exist(save_path, 'dir')
        mkdir(save_path);
    end
    save(fullfile(save_path, 'bmr_results.mat'), 'RCM', 'BMC', 'BMA');
end


fprintf('Bayesian Model Reduction completed for all subjects.\n');
end
