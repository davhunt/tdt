function [] = gng_decoding_permutations_sixtoten(subject_num, num_permutations, num_cv_splits, subfolder)

tStart = tic;
% Script for filling in CVs 6 through 10 when the first 5 have been computed
% Load label permutations from the res_cfg.mat file from one of the first 5 CVs and use these labels for 6 - 10

% set global parameters
vox_radius = 3;
%subfolder = 'durations_unsmoothed_native';
%subfolder = 'breakpoint_rts_unsmoothed_native';
num_workers = 34; % could pass this in to function?

if exist('/N/slate/brainevo/Implicit_Learning') == 7
    base_folder = ['/N/slate/brainevo/Implicit_Learning']; % If on slate
elseif exist('/Volumes/data1') == 7
    base_folder = ['/Volumes/data1']; % Or on habilis
else
    error('Can''t find base folder, are we on habilis or the HPC?');
end

if exist([base_folder '/tdt_3.999I/decoding_toolbox']) == 7 % slate
    addpath(genpath([base_folder '/tdt_3.999I/decoding_toolbox']));
else
    addpath(genpath('/Users/lab/Downloads/tdt_3.999I/decoding_toolbox')); % habilis
end
if exist([base_folder '/spm12']) == 7
    addpath(genpath([base_folder '/spm12']));
else
    addpath(genpath('/Users/lab/Downloads/spm12'));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
base_folder = ['/N/project/brainevo']; % overwrite base_folder to write in /N/project/brainevo so we have enough space
%%%%%%%%%%%%%%%%%%%%%%%%%%

% create perm dir, write to log
if ~isfolder([base_folder '/Complex_seq_analysis/0' num2str(subject_num, '%02d') '/MVPA/' subfolder '/perm'])
    mkdir([base_folder '/Complex_seq_analysis/0' num2str(subject_num, '%02d') '/MVPA/' subfolder '/perm']);
end
info_txt = [base_folder '/Complex_seq_analysis/0' num2str(subject_num, '%02d') '/MVPA/' subfolder '/perm/info.txt'];
%parameters_info = sprintf(['Permutations for subject %d: %d CV splits, searchlight radius %d voxels, %s\n'], subject_num, num_cv_splits, vox_radius, subfolder);
parameters_info = sprintf(['Permutations for subject %d: 6 through 10 CV splits, searchlight radius %d voxels, %s\n'], subject_num, vox_radius, subfolder);
writelines(parameters_info, info_txt, WriteMode="append");

% Set up unpermuted labels
labelname1 = 'HGa'; labelname2 = 'HGb'; labelname3 = 'HGc'; labelname4 = 'HGd'; labelname5 = 'HGe'; labelname6 = 'HGf'; labelname7 = 'HGg'; labelname8 = 'HGh'; labelname9 = 'HGi'; labelname10 = 'HGj'; labelname11 = 'HGk'; labelname12 = 'HGl'; labelname13 = 'HGm'; labelname14 = 'HGn';
% 14 high chunk-strength grammatical
labelname15 = 'HNGa'; labelname16 = 'HNGb'; labelname17 = 'HNGc'; labelname18 = 'HNGd'; labelname19 = 'HNGe'; labelname20 = 'HNGf'; labelname21 = 'HNGg'; labelname22 = 'HNGh'; labelname23 = 'HNGi'; labelname24 = 'HNGj'; labelname25 = 'HNGk'; labelname26 = 'HNGl'; labelname27 = 'HNGm'; labelname28 = 'HNGn';
% 14 high chunk-strength nongrammatical
labelname29 = 'LGa'; labelname30 = 'LGb'; labelname31 = 'LGc'; labelname32 = 'LGd'; labelname33 = 'LGe'; labelname34 = 'LGf'; labelname35 = 'LGg'; labelname36 = 'LGh'; labelname37 = 'LGi'; labelname38 = 'LGj'; labelname39 = 'LGk'; labelname40 = 'LGl'; labelname41 = 'LGm'; labelname42 = 'LGn';
% 14 low chunk-strength grammatical
labelname43 = 'LNGa'; labelname44 = 'LNGb'; labelname45 = 'LNGc'; labelname46 = 'LNGd'; labelname47 = 'LNGe'; labelname48 = 'LNGf'; labelname49 = 'LNGg'; labelname50 = 'LNGh'; labelname51 = 'LNGi'; labelname52 = 'LNGj'; labelname53 = 'LNGk'; labelname54 = 'LNGl'; labelname55 = 'LNGm'; labelname56 = 'LNGn';
% 14 low chunk-strength nongrammatical

labelvalue1 = 1; % value for grammatical
labelvalue2 = -1; % value for non-grammatical

labels_arr = zeros(1,56); labels_arr(1:14) = labelvalue1; labels_arr(15:28) = labelvalue2; labels_arr(29:42) = labelvalue1; labels_arr(43:56) = labelvalue2;
labelnames_arr = {};
for i = 1:56
    labelnames_arr{i} = eval(['labelname' num2str(i)]);
end

% Initialize cfg struct and fill in decoding parameters
cfg = decoding_defaults;
cfg.plot_design = 0; % no plot with parfor
cfg.verbose = 0;
cfg.analysis = 'searchlight'; % standard alternatives: 'wholebrain', 'ROI' (pass ROIs in cfg.files.mask, see below)
cfg.searchlight.radius = vox_radius;
cfg.searchlight.spherical = 1;
cfg.scale.method = 'min0max1';
cfg.scale.estimation = 'all'; % scaling across all data is equivalent to no scaling (i.e. will yield the same results), it only changes the data range which allows libsvm to compute faster
cfg.results.output = {'accuracy_minus_chance','confusion_matrix'};

% Set beta_loc and mask
beta_loc = [base_folder '/Complex_seq_analysis/' num2str(subject_num, '%.3d') '/MVPA/' subfolder];
regressor_names = design_from_spm(beta_loc);
cfg.files.mask = [base_folder '/Complex_seq_analysis/' num2str(subject_num, '%.3d') '/MVPA/' subfolder '/mask.nii'];

% Extract all information for the cfg.files structure (labels will be [1 -1] if not changed above)
cfg = decoding_describe_data(cfg,labelnames_arr,labels_arr,regressor_names,beta_loc);

rng('shuffle') % make sure randsample in generating the design below is truly random
cfg.design = make_design_custom_GNG(cfg, 0); % need to seed RNG in function % not permuted labels yet
cfg.files.chunk = cfg.design.chunk; cfg.design = rmfield(cfg.design, 'chunk');
cfg.design.function.name = 'make_design_cv';
%perm_designs = make_design_permutation(cfg,num_permutations,0);
% Only want permuted labels from permuted designs, not train or test, will replace those later


% See how many permutations have already written been written, continue from there (i.e. if there are 100 already start from 101)
% rx = '^perm([0-9]+)$';
% perm_folder = ['/N/slate/brainevo/Implicit_Learning/Complex_seq_analysis/' num2str(subject_num, '%.3d') '/MVPA/' subfolder '/perm'];
% file_list = dir(perm_folder);
% perms_done = {};
% for i = 1:length(file_list)
%     %if ~isempty(regexp(file_list(i).name, rx, 'once'))
%     tokens = regexp(file_list(i).name, rx, 'tokens');
%     if ~isempty(tokens)
%         perms_done{end+1} = tokens{1}{1};
%     else
%         start_perm = 1;
%     end
% end
% if size(perms_done, 2) > 0
%     perms_done = str2double(perms_done);
%     start_perm = max(perms_done) + 1;
% end
% log_info = sprintf(['%s Proceeding with permutations, starting with %d and ending with %d\n'], datetime, start_perm, start_perm+num_permutations-1);
% writelines(log_info, info_txt, WriteMode="append");

%%%% Find minimum permutation number with less than 10 files (assuming all of these only 5 have been done, so do 6-10)
just_five_perms = [];
for i = 1:500
    contents = dir([base_folder '/Complex_seq_analysis/' num2str(subject_num, '%.3d') '/MVPA/' subfolder '/perm/perm' num2str(i, '%03d')]);
    if numel(contents) < 10
        just_five_perms = [just_five_perms, i];
    end
end
start_perm = min(just_five_perms);
%%%%









myCluster = parcluster('local');
parpool(num_workers);
sc = parallel.pool.Constant(RandStream('Threefry', 'Seed', 'shuffle'));

%parfor i_perm = start_perm:start_perm+num_permutations-1
parfor i_perm = 1:num_permutations

    stream = sc.Value;
    stream.Substream = i_perm;

    perm_num = i_perm + start_perm - 1;
    %perm_results_dir = [base_folder '/Complex_seq_analysis/' num2str(subject_num, '%.3d') '/MVPA/' subfolder '/perm/perm' num2str(i_perm, '%.3d')];
    perm_results_dir = [base_folder '/Complex_seq_analysis/' num2str(subject_num, '%.3d') '/MVPA/' subfolder '/perm/perm' num2str(perm_num, '%.3d')];

    if ~isfolder(perm_results_dir)
        mkdir(perm_results_dir)
    end
    dispv(1, 'Permutation %i/%i', perm_num, num_permutations)

    A = load([base_folder '/Complex_seq_analysis/' num2str(subject_num, '%.3d') '/MVPA/' subfolder '/perm/perm' num2str(perm_num, '%.3d') '/results_rep01/res_cfg.mat']);
    % use results_rep01 but could use any results_repXX since the labels should be the same for all of them
    cfg_copy = cfg;
    %cfg_copy.design.label = perm_designs{i_perm}.label; % get permuted labels from perm_designs
    cfg_copy.design.label = A.cfg.design.label; % get permuted labels from res_cfg.mat

    %for spl = 1:num_cv_splits
    for spl = 6:10

        spl_design = make_design_custom_GNG(cfg_copy, 1, stream); % need to seed RNG in function

        cfg_copy.results.dir = [perm_results_dir '/results_rep' num2str(spl, '%.2d')];
        if ~isfolder(cfg_copy.results.dir)
            mkdir(cfg_copy.results.dir)
        end
        cfg_copy.design = spl_design; % get train and test sets from make_design_custom_GNG, generate new ones each split

        results = decoding(cfg_copy);
        gzip([cfg_copy.results.dir '/res_accuracy_minus_chance.nii']);
        delete([cfg_copy.results.dir '/res_accuracy_minus_chance.nii']);
    end
end

delete(gcp('nocreate'));

tEnd = toc(tStart);
log_info = sprintf(['%s Done with %d permutations\n'], datetime, num_permutations);
writelines(log_info, info_txt, WriteMode="append");

end
