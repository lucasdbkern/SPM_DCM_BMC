function subjectcreator(dcm_file_path, model_type, noise_magnitude_Ep, noise_magnitude_u, noise_magnitude_erp, Nm, customA1, customA2, customA3, param_key, save_folder)
%==========================================================================
%This function generates simulated subjects based on a given DCM file 
% and custom connectivity matrices.

% Inputs:
%   dcm_file_path        - File path for the DCM file
%   model_type           - Model type for simulation
%   noise_magnitude_Ep   - Magnitude of noise for Ep
%   noise_magnitude_u    - Magnitude of noise for u
%   noise_magnitude_erp  - Magnitude of noise for erp
%   Nm                   - Number of subjects to simulate
%   customA1, customA2, customA3 - Custom matrices to replace A{1,1}, A{1,2}, A{1,3}
%   param_key            - Parameter key for naming files
%   save_folder          - Folder path to save the generated files
%==========================================================================

% Load the DCM file
load(dcm_file_path, 'DCM');
M = DCM.M;
spm('defaults', 'eeg');

% Configure DCM
DCM = spm_dcm_erp_dipfit(DCM);
DCM.M.Nmax =1;
DCM.options.Nmodes = 16;
M.dipfit.type = model_type;
M.dipfit.modality = model_type;
M.dipfit.Ns = 1;
M.dipfit.Nc = 1;
U.dt = DCM.xU.dt;
U.X = DCM.xU.X;


% Load custom matrices if provided

if ~isfile(customA1)
    error('file not found')
end

if ~isfile(customA2)
    error('file not found')
end

if ~isfile(customA3)
    error('file not found')
end

if ~isempty(customA1)
    loadedA1 = load(customA1);
    fields = fieldnames(loadedA1);
    A1 = loadedA1.(fields{1});
    DCM.A{1, 1} = A1;
end
if ~isempty(customA2)
    loadedA2 = load(customA2);
    fields = fieldnames(loadedA2);
    A2 = loadedA2.(fields{1});
    DCM.A{1, 2} = A2;
end
if ~isempty(customA3)
    loadedA3 = load(customA3);
    fields = fieldnames(loadedA3);
    A3 = loadedA3.(fields{1});
    DCM.A{1, 3} = A3;
end

% Update DCM.B if all custom matrices are provided
if ~isempty(customA1) && ~isempty(customA2) && ~isempty(customA3)
    DCM.B{1, 1} = (~~DCM.B{1, 1})&(DCM.A{1, 1}|DCM.A{1, 2}|DCM.A{1, 3});
end


for i=1:3
    DCM.Ep.A{1, i} = DCM.Ep.A{i} .*  DCM.A{1, i} - (1 - DCM.A{1, i}) .* 32;
    DCM.M.pE.A{1, i} = DCM.M.pE.A{i} .*  DCM.A{1, i} - (1 - DCM.A{1, i}) .* 32;
    DCM.M.pC.A{1, i} = DCM.M.pC.A{i} .*  DCM.A{1, i} - (1 - DCM.A{1, i}) .* 32; 
end 




i = 1;

% Creates the subfolder within the save_folder
subfolder = fullfile(save_folder, sprintf('subject_%d', param_key));
if ~exist(subfolder, 'dir')
    mkdir(subfolder);
end

while i <= Nm
    
    Pvec = spm_vec(DCM.Ep); 
    sqrtpC = sqrt(spm_vec(DCM.M.pC));
    noise = normrnd(0, 1, [158, 1]) .* sqrtpC; % ERP: 158
    Pvec_noisy = Pvec + (noise_magnitude_Ep * noise); 
    P_noisy = spm_unvec(Pvec_noisy, DCM.Ep);
    
    U.u = spm_erp_u((1:M.ns)*U.dt, P_noisy, M);
    noise_u = noise_magnitude_u * randn(size(U.u));
    U.u = U.u + noise_u; 

    DCMtemp = DCM;
    DCMtemp.M.hE = 8;
    DCMtemp.xU.u = U.u;
    DCMtemp.M.P = P_noisy;
    DCMtemp.M.pE = P_noisy;
    DCMtemp.M.pC = spm_unvec(0*spm_vec(DCMtemp.M.pC),DCMtemp.M.pC);
    DCMtemp.M.gC = spm_unvec(0*spm_vec(DCMtemp.M.gC),DCMtemp.M.gC);
    DCMtemp.M.gE = DCMtemp.Eg;
    DCMtemp.name = strcat(DCM.name, '_temp');
    DCMtemp = spm_dcm_erp(DCMtemp);
    erp = DCMtemp.H;
    pst = DCMtemp.xY.pst;

    all_erp = {};
    all_pst = [];
    
    for cond = 1:length(erp)
        y = erp{cond} + (noise_magnitude_erp * randn(size(erp{cond})));
        
        if all(all(y == 0)) || any(any(isnan(y)))
            continue; 
        end

        all_erp = [all_erp, y]; % Collecting all erp outputs
        all_pst = [all_pst, pst]; % Collecting time stamps
        % % % Plotting ERP for the current participant and condition
        % figure; 
        % plot(pst, y, 'b');
        % title(sprintf('ERP: Participant %d, Condition %d', i, cond), 'FontSize', 16);
        % 
        % % % Plotting U.u time series for the current participant
        % figure; 
        % plot((1:M.ns)*U.dt, U.u);
        % title(sprintf('U.u Time Series: Participant %d', i), 'FontSize', 16);
        % xlabel('Time (s)');
        % ylabel('Input Signal');
    end
    
    if ~isempty(all_erp)
        % Update DCM structure 
        DCMi = DCM;
        DCMi.xY.y = all_erp;
        DCMi.xY.pst = all_pst;
        DCMi.xU = U;
        DCMi.M = M;
        DCMi.Ep = P_noisy;
        %DCMi.Ep.A = DCM.Ep.A;
        DCMi.M.pE.A = DCM.M.pE.A;
        DCMi.M.pC.A = DCM.M.pC.A;
        DCMi.A = DCM.A;
        DCMi.B =DCM.B;
        DCMi.M.pE.C= DCM.M.pE.C * 4;
        

        % Initialize parts of the filename based on custom matrices
        part1 = '';
        
        % Check if customA1 is not empty and adjust filename part
        if ~isempty(customA1)
            [~, name, ~] = fileparts(customA1); % Extract the name without extension
            part1 = name; % Use the extracted name directly in the filename part
        end

        filename = sprintf('subject_%s_%d.mat', part1, param_key);
        
        % Change DCM name within the file
        DCMi.name = filename;

        % Create the full file path
        fullFilePath = fullfile(subfolder, filename);

        % Save the updated DCM structure
        save(fullFilePath, 'DCMi', '-v7.3');

        i = i + 1;
    end
end
