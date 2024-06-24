function multiple_subject_bmc(dcm_file_path, n)
% ===========================================================================
% This script generates multiple subjects and performs BMC for each set.

% Inputs:

% dcm_file_path - File path for the DCM files.

% n - Number of subjects to simulate.
 

%==========================================================================

Ep = 0.1  %0.1; %   noise_magnitude_Ep   - Magnitude of noise for Ep
u = 0.25 %0.45; %   noise_magnitude_u    - Magnitude of noise for u
erp = 0.015 %0.035; %   noise_magnitude_erp  - Magnitude of noise for erp

% Loop to create files for each subject triple
%--------------------------------------------------------------------------
for i = 1:n
    subjectcreator(dcm_file_path, 'ERP', Ep, u, erp, 1, 'FullA1.mat', 'FullA2.mat', 'FullA3.mat', i,'/Users/lucaskern/Desktop/Desktop/UCL/Research_Project/roving-oddball-dataset/subject_folder');
    subjectcreator(dcm_file_path, 'ERP', Ep, u, erp, 1, 'ORA1.mat', 'ORA2.mat', 'ORA3.mat', i, '/Users/lucaskern/Desktop/Desktop/UCL/Research_Project/roving-oddball-dataset/subject_folder');
    subjectcreator(dcm_file_path, 'ERP', Ep, u, erp, 1, 'TRA1.mat', 'TRA2.mat', 'TRA3.mat', i, '/Users/lucaskern/Desktop/Desktop/UCL/Research_Project/roving-oddball-dataset/subject_folder');
   

    % run the Bayesian Model Comparison
    bmc(num2str(i),'/Users/lucaskern/Desktop/Desktop/UCL/Research_Project/roving-oddball-dataset/subject_folder')
end
 