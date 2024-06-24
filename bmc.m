function bmc(param_key, subject_folder)
%==========================================================================
%This function performs Bayesian Model Comparison on the generated subjects.

% Ensure to wrap file paths in brackets {file, file2, file3} for proper functionality.

% Inputs:
%   param_key           - Parameter key for naming files.
%   subject_folder      - Folder containing the generated subjects.
 
%==========================================================================


spm('defaults', 'eeg')

% Define the subfolder path
subfolder = fullfile(subject_folder, sprintf('subject_%s', param_key));

% Check if the subfolder exists
if ~exist(subfolder, 'dir')
    error('The specified subfolder does not exist: %s', subfolder);
end

% List .mat files in the subfolder
file_list = dir(fullfile(subfolder, '*.mat'));
file_paths = fullfile({file_list.folder}, {file_list.name});


% Number of files/models
n = length(file_paths);
if n == 0
    error('No .mat files found in subject %s', param_key);
end

DCMs = {};
for i = 1:n
    loaded_data = load(file_paths{i});
    DCMs{i} = loaded_data.DCMi;
    DCMs{i}.M.Nmax = 64;
end

spm_plot_ci(DCMs{1}.Ep, DCMs{i}.Cp, DCMs{i}.x, '', 'log');
spm_plot_ci(DCMs{2}.Ep, DCMs{i}.Cp, DCMs{i}.x, '', 'log');
spm_plot_ci(DCMs{3}.Ep, DCMs{i}.Cp, DCMs{i}.x, '', 'log');





subplot(3, 1, 1);
plot(DCMs{1}.xY.y{1});
title('Full Model - Subject 1');

subplot(3, 1, 2);
plot(DCMs{2}.xY.y{1});
title('ORA Model - Subject 1');

subplot(3, 1, 3);
plot(DCMs{3}.xY.y{1});
title('TRA Model - Subject 1');



% Compare timeseries for i=1 and j=1, 2, 3
%compare_timeseries(DCMs, 1, [1, 2, 3]);


DCMij = cell(n, n);
for i = 1:n
    for j = 1:n
        DCMtemp = DCMs{i};
        DCMtemp.xY.y = DCMs{j}.xY.y;
        DCMtemp.name = spm_file(DCMtemp.name, 'prefix', 'bmc_');
        DCMtemp.assignment_tag = sprintf('Model %d, Timeseries %d', i, j); % Add tracking tag
        DCMij{i,j} = spm_dcm_erp(DCMtemp);
    end
end

% for i =1:n
%     for j =1:n
%         spm_plot_ci(DCMij{i,j}.Ep, DCMij{i,j}.Cp, DCMij{i,j}.x, '', '--r');
%     end
% end
% 
% for i = 1:n
%     for j=1:n
%         imagesc(DCMij{i,j}.M.pE.A{i}) 
%     end
% end


% Create the subfolder path
subfolder = fullfile(subject_folder, sprintf('subject_%s', param_key), 'derivatives');

% Check if the subfolder exists, and create it if it doesn't
if ~exist(subfolder, 'dir')
    mkdir(subfolder);
end

% Define the full file path for saving DCMij.mat
filename = fullfile(subfolder, 'DCMij.mat');

% Save the DCMij cell array to the specified file
save(filename, 'DCMij', '-v7.3');

fprintf('DCMij saved to %s\n', filename);


subplot(3, 1, 1);
plot(DCMij{1,1}.xY.y{1});
title('Full Model ');

subplot(3, 1, 2);
plot(DCMij{2,2}.xY.y{1});
title('ORA Model ');

subplot(3, 1, 3);
plot(DCMij{3,3}.xY.y{1});
title('TRA Model ');

% Check assignments
for i = 1:n
    for j = 1:n
        fprintf('DCMij{%d,%d} -> %s\n', i, j, DCMij{i,j}.assignment_tag);
    end
end

% Compare timeseries for i=1 and j=1, 2, 3
%    compare_timeseries(DCMij, 1, [1, 2, 3]);
end



function compare_timeseries(DCMij, i, j_values)
    threshold = 10^-1;
    for j = j_values
        timeseries_i = DCMij{i,i}.xY.y{1};
        timeseries_j = DCMij{i,j}.xY.y{1};

        % Initialize the matrix to store absolute distances
        abs_distance_matrix = abs(timeseries_i - timeseries_j);

        % Check if distances satisfy the threshold
        threshold_matrix = abs_distance_matrix > threshold;

        % Calculate the sum of distances per channel
        sum_distances_per_channel = sum(abs_distance_matrix, 1);
        % Calculate the average distance per channel
        avg_distance_per_channel = sum_distances_per_channel / size(abs_distance_matrix, 1);

        % Display the results
        fprintf('Comparison of timeseries for i=%d and j=%d:\n', i, j);
        disp('Absolute distance matrix:');
        disp(abs_distance_matrix);
        disp('Threshold matrix (>10^-1):');
        disp(threshold_matrix);
        disp('Sum of distances per channel:');
        disp(sum_distances_per_channel);
        disp('Average distance per channel:');
        disp(avg_distance_per_channel);
    end
end



