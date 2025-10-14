function [] = decoding_template_gng_permutations_parallel_all_splits(subject_num, num_permutations, cv_splits)

tStart = tic;
%% Demo how to create a permutation design
%
% This demo shows how to calculate a permutation analysis.
%
% For this, we need the cfg of the normal design.
% This can 
% - either be created new (not shown below, see e.g. the decoding_tutorial)
% - or loaded from the cfg.mat of an already calculated analysis (thats
%   what we do here)
%
% Data that fits to this script can be computed with
% demo8_demodata_decoding_tutorial_motion_direction.m
%
% Kai, 2016/07/25
 
%% Check that SPM and TDT are available on the path
addpath(genpath('/Users/lab/Downloads/tdt_3.999I/decoding_toolbox'));
addpath(genpath('/Users/lab/Downloads/spm12'));

%%%% cv_splits can be row vector i.e. [1, 2, 3, 4, 5, 6, 7, 8, 9, 10], or just a number (i.e. 4)

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
cfg.searchlight.radius = 3; % use searchlight of radius 3 (by default in voxels), see more details below
cfg.searchlight.spherical = 1;
cfg.scale.method = 'min0max1';
cfg.scale.estimation = 'all'; % scaling across all data is equivalent to no scaling (i.e. will yield the same results), it only changes the data range which allows libsvm to compute faster
cfg.results.output = {'accuracy_minus_chance','confusion_matrix'};

% Set beta_loc and mask
%beta_loc = ['/Volumes/data1/Complex_seq_analysis/' num2str(subject_num, '%.3d') '/MVPA/durations_unsmoothed'];
%regressor_names = design_from_spm(beta_loc);
%cfg.files.mask = ['/Volumes/data1/Complex_seq_analysis/' num2str(subject_num, '%.3d') '/MVPA/durations_unsmoothed/mask.nii'];
beta_loc = ['/Volumes/data1/Complex_seq_analysis/' num2str(subject_num, '%.3d') '/MVPA/onsets_only_unsmoothed'];
regressor_names = design_from_spm(beta_loc);
cfg.files.mask = ['/Volumes/data1/Complex_seq_analysis/' num2str(subject_num, '%.3d') '/MVPA/onsets_only_unsmoothed/mask.nii'];

% Extract all information for the cfg.files structure (labels will be [1 -1] if not changed above)
cfg = decoding_describe_data(cfg,labelnames_arr,labels_arr,regressor_names,beta_loc);

rng('shuffle') % make sure randsample in generating the design below is truly random
%cfg.design = make_design_custom_GNG(cfg, 1, stream); % need to seed RNG in function 
cfg.design = make_design_custom_GNG(cfg, 0); % need to seed RNG in function % not permuted labels yet
cfg.files.chunk = cfg.design.chunk; cfg.design = rmfield(cfg.design, 'chunk');
cfg.design.function.name = 'make_design_cv';
perm_designs = make_design_permutation(cfg,num_permutations,0);
% Only want permuted labels from permuted designs, not train or test, will replace those later

myCluster = parcluster('local');
parpool(myCluster.NumWorkers);
sc = parallel.pool.Constant(RandStream('Threefry', 'Seed', 'shuffle'));

parfor i_perm = 1:num_permutations
%for i_perm = 1:num_permutations
%%%%%%%%
    stream = sc.Value;
    stream.Substream = i_perm;

    %perm_results_dir = ['/Volumes/data1/Complex_seq_analysis/' num2str(subject_num, '%.3d') '/MVPA/durations_unsmoothed/perm/perm' num2str(i_perm, '%.3d')];
    perm_results_dir = ['/Volumes/data1/Complex_seq_analysis/' num2str(subject_num, '%.3d') '/MVPA/onsets_only_unsmoothed/perm/perm' num2str(i_perm, '%.3d')];
    if ~isfolder(perm_results_dir)
        mkdir(perm_results_dir)
    end
    dispv(1, 'Permutation %i/%i', i_perm, num_permutations)


    %// perm_cfg = cfg;
    %// %perm_cfg.design = make_design_custom_GNG(cfg, stream); % need to seed RNG in function 
    %// perm_cfg.design = make_design_custom_GNG(cfg); % need to seed RNG in function 
    %// perm_cfg.files.chunk = cfg.design.chunk; rmfield(cfg.design, 'chunk');
    %// perm_cfg.design.function.name = 'make_design_cv';
    %// permuted_designs = make_design_permutation(cfg,num_permutations,0);
    %// % Only want permuted labels from permuted_designs, not train or test

    cfg_copy = cfg;
    cfg_copy.design.label = perm_designs{i_perm}.label; % get permuted labels from perm_designs
    %cfg_copy.design.label = perm_designs{i_perm-100}.label; % get permuted labels from perm_designs

    for spl = 1:length(cv_splits)

        spl_design = make_design_custom_GNG(cfg_copy, 1, stream); % need to seed RNG in function
        %spl_design = cfg_copy.design;

        %cfg_copy = cfg;
        cfg_copy.results.dir = [perm_results_dir '/results_rep' num2str(spl, '%.2d')];
        if ~isfolder(cfg_copy.results.dir)
            mkdir(cfg_copy.results.dir)
        end
        cfg_copy.design = spl_design; % get train and test sets from make_design_custom_GNG, generate new ones each split
        %cfg_copy.design.label = perm_designs{i_perm}.label; % get permuted labels from perm_designs

        results = decoding(cfg_copy);
    end
end

delete(gcp('nocreate'));

tEnd = toc(tStart);
disp(['Running ' num2str(length(cv_splits)) ' splits on subject ' num2str(subject_num) ' took ' num2str(tEnd, '%.2f') ' seconds. Done.'])

end