function process_subjects()

% Define the path to the main folder containing subject subfolders
subject_folder = '/Users/lucaskern/Desktop/Desktop/UCL/Research_Project/roving-oddball-dataset/subject_folder';

folderPattern = fullfile(subject_folder, 'subject_*');
subfolders = dir(folderPattern);
numSubjects = length(subfolders);

% Extract numbers from folder names and sort them
subjectNumbers = zeros(numSubjects, 1);
for i = 1:numSubjects
    num = regexp(subfolders(i).name, '\d+', 'match');
    subjectNumbers(i) = str2double(num{1});
end
% Sort numbers and use them to sort subfolders array
[~, sortIndex] = sort(subjectNumbers);
sortedSubfolders = subfolders(sortIndex);

disp(sortedSubfolders);

GCM = cell(numSubjects, 3);

% Loop through each sorted subject folder
for i = 1:numSubjects
    subjectFolder = fullfile(subject_folder, sortedSubfolders(i).name);
    subjectNum = regexp(sortedSubfolders(i).name, '\d+', 'match');
    subjectNum = subjectNum{1}; 

 
    entry1 = load(fullfile(subjectFolder, sprintf('subject_FullA1_%s.mat', subjectNum)));
    entry2 = load(fullfile(subjectFolder, sprintf('subject_ORA1_%s.mat', subjectNum)));
    entry3 = load(fullfile(subjectFolder, sprintf('subject_TRA1_%s.mat', subjectNum)));
   
    % Assign the entries to the corresponding row in GCM
    GCM{i, 1} = entry1;
    GCM{i, 2} = entry2;
    GCM{i, 3} = entry3;
end
disp(GCM);

save(fullfile(subject_folder, 'GCM.mat'), 'GCM', '-v7.3');
end
