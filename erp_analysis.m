%% Import data
% -------------------------------------------------------------------------
datapath = fullfile(pwd, 'data', 'AEF');

spm('defaults', 'eeg');

dataset = 2; 
model   = 'ERP'; 
switch dataset
    case 1
        sub = 'sub-001'; 
        ses = 'ses-002'; 
        runs = {'run-002', 'run-003'}; 
    case 2
        sub  = 'sub-002'; 
        ses  = 'ses-001'; 
        runs = {'run-001', 'run-003'}; 
end
task = 'task-tone'; 
spath = fullfile(datapath, sub, ses, 'meg'); 

fname = ['be_hdfff' sub '_' ses '_' task '_' runs{1} '_meg.mat']; 
fpath = fullfile(spath, fname);

D = spm_eeg_load(fpath);


%% Define ROIs and coordinates
% -------------------------------------------------------------------------

Sname = {
    'lA1'
    'rA1'
    'lSTG'
    'rSTG'
    'rIFG'
}; 

Lpos = [
   -42 -21 7; % lA1
    46 -13 8; % rA1
   -61 -31 9; % lSTG
    58 -23 9; % rSTG
    46  20 6; % rIFGr
]'; 

lA1  = 1;
rA1  = 2;
lSTG = 3; 
rSTG = 4;
rIFG = 5;

nr = numel(Sname);

%% Define the connectivity structure
% -------------------------------------------------------------------------

A = []; 

% Forward connections
A{1} = zeros(nr, nr);
A{1}(lSTG,  lA1) = 1; 
A{1}(rSTG,  rA1) = 1; 
A{1}(rIFG, rSTG) = 1; 

% Backward connections
A{2} = A{1}'; 

% Lateral connections
A{3} = zeros(nr, nr);
A{3}(rA1,  lA1) = 1; 
A{3}(rSTG, lSTG) = 1; 
A{3} = A{3} + A{3}';  

% Modulated forward - backward connections and self inhibition
B = {A{1} + A{2} + eye(nr)};

% Inputs
C = zeros(nr, 1);
C(lA1) = 1; 
C(rA1) = 1; 


%%
DCM = []; 

% Specify the data file 
% ---------------------------------------------------------------------
DCM.xY.Dfile    = fullfile(D); 
DCM.xY.modality = 'MEG';

    
% Specify options
% ---------------------------------------------------------------------
DCM.options.Tdcm     = [-50 450]; % twin in ms
DCM.options.spatial  = 'ECD';      % Equivalent Current Dipole
DCM.options.model    = 'ERP';      % Canonical microcircuit model 
DCM.options.analysis = 'ERP';      % ERP analysis
DCM.options.D        =  2;         % 
DCM.options.Nmodes = 16;

DCM.Sname = Sname;                  % Sources names
DCM.Lpos  = Lpos;                   % Sources positions

DCM.options.trials = [1 2];        % Model conditions 1 and 2
DCM.xU.X           = [1; 0]; 
DCM.xU.name        = {'Deviant'}; 

DCM.A = A; 
DCM.B = B; 
DCM.C = C; 

DCM = spm_dcm_erp(DCM);

%% 
spm_dcm_erp_results(DCM)
