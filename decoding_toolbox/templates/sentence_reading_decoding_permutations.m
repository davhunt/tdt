function [] = sentence_reading_decoding_permutations(subject_num, input_dir, spreadsheetFile, window_type, contrast, num_permutations)

% This function does a permutation analysis on the sentence reading data, generating decoding accuracy maps from 
% betas with randomly permuted labels, to provide a null distribution of maps with which to compare the real-labelled
% results.
% Four runs/blocks of data were collected for this experiment, with trial types balanced across all blocks, and so a
% leave-one-run-out (LORO) cross-validation is performed where a classifier is trained on three blocks, and tested on
% held out block.

% subject_num: Subject number of subject, will be converted to string
% input_dir: Directory under which we can find subject folders
% spreadsheetFile: Excel spreadsheet with subject trial data (onsets + durations)
% window_type: Directory under subject folders, what type of window was used to calculate betas
% contrast: The contrast to run permutation analysis on, can be "AP," "CAOR," or "GUG"
% num_permutations: Number of permutations to do, takes 5-20 minutes per permutation, ~5 for AP/CAOR and ~20 for GUG

% Notes:
% - For contrasts "AP" or "CAOR" only use grammatical trials, while "GUG" uses all trials
% - Each block has balanced numbers of G, UG, A, P, CA, OR trials, 6 of each except one of AG, PUG, PG, or AUG which has 7
%    (49 trials per block, 196 total)

tStart = tic;

if exist('/N/slate/brainevo/Implicit_Learning') == 7
    base_folder = '/N/slate/brainevo/Implicit_Learning'; % If on slate
    ants_registrations_dir = '/N/slate/brainevo/Implicit_Learning/ants_registrations';
elseif exist('/Volumes/data1') == 7
    base_folder = '/Volumes/data1'; % Or on habilis
    ants_registrations_dir = '/Users/lab/Downloads/ants_registrations';
else
    error('Can''t find base folder, are we on habilis or the HPC?');
end

if exist([base_folder '/tdt_3.999I/decoding_toolbox']) == 7 % slate
    addpath(genpath([base_folder '/tdt_3.999I/decoding_toolbox']));
else
    addpath(genpath('/Users/lab/Downloads/tdt_3.999I/decoding_toolbox')) % habilis
end
if exist([base_folder '/spm12']) == 7
    addpath(genpath([base_folder '/spm12']));
else
    addpath(genpath('/Users/lab/Downloads/spm12'))
end

cfg = decoding_defaults;

sheetName = 'SPM_events';
subjChar = num2str(subject_num, '%02d'); % later filepaths require two digit subject IDs
T = readtable(spreadsheetFile, "Sheet", sheetName);
T = T(T.subject == subject_num, :);
blocks = unique(T.scanner_sequence_number, 'stable');
if numel(blocks) ~= 4
    error('Number of blocks should be 4, seeing %d.', numel(blocks));
end


%labelvalue1 = 1; % value for active, conjoined active, or grammatical
%labelvalue2 = -1; % value for passive, object relative, or non-grammatical
if contrast == "AP"
    % keep grammatical, active/passive trials
    T = T(matches(T.grammar_type, ["A","P"]),:);
    T = T(T.grammaticality == "G", :);
    labels = double(T.grammar_type == "A");
    labels(labels == 0) = -1;
elseif contrast == "CAOR"
    % keep grammatical, CA/OR trials
    T = T(matches(T.grammar_type, ["CA","OR"]),:);
    T = T(T.grammaticality == "G", :);
    labels = double(T.grammar_type == "CA");
    labels(labels == 0) = -1;
elseif contrast == "GUG"
    % keep all trials
    labels = double(T.grammaticality == "G");
    labels(labels == 0) = -1;
else
    error('Contrast must be "AP", "CAOR", or "GUG".')
end

num_betas = size(T,1);
names = cell(num_betas,1); chunks = zeros(num_betas,1);
for i = 1:size(T,1)
    sentence_id = char(T.sentence_ID(i));
    names{i} = fullfile(input_dir, subjChar, window_type, ['beta_' sentence_id '.nii']);
    chunks(i) = find(blocks == T.scanner_sequence_number(i));
    if ~isfile(names{i})
        error('Can''t find file %s.', names{i});
    end
end

cfg.files.name = names;
cfg.files.label = labels;
cfg.files.chunk = chunks;
cfg.files.mask = [ants_registrations_dir '/sentence_reading/' subjChar '/mask_epi.nii.gz'];

%cfg.design = make_design_cv(cfg);
cfg.design = make_design_cv(cfg);
perm_designs = make_design_permutation(cfg,num_permutations);

cfg.plot_design = 0; % no plot in parfor
cfg.verbose = 0;

% Set the analysis that should be performed (default is 'searchlight')
cfg.analysis = 'searchlight'; % standard alternatives: 'wholebrain', 'ROI' (pass ROIs in cfg.files.mask, see below)
cfg.searchlight.radius = 3; % use searchlight of radius 3 (by default in voxels), see more details below
cfg.searchlight.spherical = 1;
cfg.scale.method = 'min0max1';
cfg.scale.estimation = 'all'; % scaling across all data is equivalent to no scaling (i.e. will yield the same results), it only changes the data range which allows libsvm to compute faster
cfg.results.output = {'accuracy_minus_chance','confusion_matrix'};
% Set the output directory where data will be saved, e.g. 'c:\exp\results\buttonpress'

myCluster = parcluster('local');
parpool(num_workers);
sc = parallel.pool.Constant(RandStream('Threefry', 'Seed', 'shuffle'));

%perm_results_dir = [base_folder '/Complex_seq_analysis/' subject_num '/MVPA/' subfolder '/perm/perm' num2str(perm_num, '%.3d')];
perm_results_dir = char(fullfile(input_dir,subjChar,window_type,strcat(contrast,'_perm')));
if ~isfolder(perm_results_dir)
    mkdir(perm_results_dir)
end

% Check on perms already done
rx = '^perm([0-9]+)$';
file_list = dir(perm_results_dir);
perms_done = {};
for i = 1:length(file_list)
    %if ~isempty(regexp(file_list(i).name, rx, 'once'))
    tokens = regexp(file_list(i).name, rx, 'tokens');
    if ~isempty(tokens)
        perms_done{end+1} = tokens{1}{1};
    else
        start_perm = 1;
    end
end
if size(perms_done, 2) > 0
    perms_done = str2double(perms_done);
    start_perm = max(perms_done) + 1;
end

parfor perm_num=1:num_permutations
    cfg.results.dir = char(fullfile(perm_results_dir,['perm' num2str(start_perm+perm_num-1, '%.3d')]));
    if ~isfolder(cfg.results.dir)
        mkdir(cfg.results.dir)
    end
    cfg.design = perm_designs{perm_num};
    cfg.design.unbalanced_data = 'ok'; % number of training/testing trials of a type (A, CA, P, OR) may be off by one to number of contrasting trials in a block (12 vs 13)
    cfg.files.label = cfg.design.label(:,1);
    results = decoding(cfg);
    gzip([cfg.results.dir '/res_accuracy_minus_chance.nii']);
    delete([cfg.results.dir '/res_accuracy_minus_chance.nii']);
end

delete(gcp('nocreate'));

tEnd = toc(tStart);
disp(['Running LORO CV on subject ' subjChar ' took ' num2str(tEnd, '%.2f') ' seconds. Done.'])

end
