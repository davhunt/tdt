function [] = gng_decoding_roi(subject_num, num_permutations, num_cv_splits, subfolder, roi_file, outdir)

%% Performs MVPA decoding analysis with TDT using an ROI mask instead of multiple "searchlights"

% subject_num: Subject number of subject, will be converted to string
% num_permutations: Number of permutations to do, each permutation will have the same random shuffling of labels (but maybe different assignments of data to folds)
% num_cv_splits: Number of different fold assignments/CVs to do per permutation, should be ~ 5 - 10. https://www.sciencedirect.com/science/article/pii/S1053811921004225
% subfolder: Folder under "MVPA" to set as working dir, save results in (e.g. "durations_unsmoothed")
% roi_file: ROI mask file (.nii)
% outdir: Directory under "subfolder" to save results in

tStart = tic;

% set global parameters
subject_num = num2str(subject_num, '%03d');
vox_radius = 3;
%subfolder = 'breakpoint_rts_unsmoothed_native';
%num_workers = 34; % could pass this in to function?

if exist('/N/slate/brainevo/Implicit_Learning') == 7
    base_folder = ['/N/slate/brainevo/Implicit_Learning']; % If on slate
elseif exist('/Volumes/data1') == 7
    base_folder = ['/Volumes/data1']; % Or on habilis
else
    error('Can''t find base folder, are we on habilis or the HPC?');
end

if exist([base_folder '/tdt_3.999I/decoding_toolbox']) == 7 % slate
    addpath(genpath([base_folder '/tdt_3.999I/decoding_toolbox']));
    addpath(genpath([base_folder '/spm12']));
    addpath(genpath([base_folder '/matlab_nifti_tools'));
else
    addpath(genpath('/Users/lab/Downloads/tdt_3.999I/decoding_toolbox')); % habilis
    addpath(genpath('/Users/lab/Downloads/spm12'));
    addpath(genpath('/Users/lab/Downloads/matlab_nifti_tools'));
end

% create perm dir, write to log
if ~isfolder([base_folder '/Complex_seq_analysis/' subject_num '/MVPA/' subfolder '/perm'])
    mkdir([base_folder '/Complex_seq_analysis/' subject_num '/MVPA/' subfolder '/perm']);
end
info_txt = [base_folder '/Complex_seq_analysis/' subject_num '/MVPA/' subfolder '/perm/info.txt'];
parameters_info = sprintf(['Permutations for subject %s: %d CV splits, searchlight radius %d voxels, %s\n'], subject_num, num_cv_splits, vox_radius, subfolder);
writelines(parameters_info, info_txt, WriteMode="append");

% Or on habilis...
%addpath(genpath('/Users/lab/Downloads/tdt_3.999I/decoding_toolbox'));
%addpath(genpath('/Users/lab/Downloads/spm12'));

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

results_dir = [base_folder '/Complex_seq_analysis/' subject_num '/MVPA/' subfolder '/' outdir '/TDT_results'];
if ~isfolder(results_dir)
    mkdir(results_dir);
end

% Intersect ROI w/ mask
wholebrain_mask = [base_folder '/Complex_seq_analysis/' subject_num '/MVPA/' subfolder '/mask.nii'];
wholebrain_mask = load_untouch_nii(wholebrain_mask)
roi_file_nii = load_untouch_nii(roi_file);
wholebrain_mask_roi_intersect = uint8(wholebrain_mask.img) .* uint8(roi_file_nii.img);
roi_file_nii.img = wholebrain_mask_roi_intersect;
wholebrain_mask_roi_intersect_out = [base_folder '/Complex_seq_analysis/' subject_num '/MVPA/' subfolder '/brocas_rois/intersect_mask.nii.gz']; % can overwrite this
save_untouch_nii(roi_file_nii, wholebrain_mask_roi_intersect_out);

myCluster = parcluster('local');
parpool(myCluster.NumWorkers); % On habilis should be "6" could be much more on HPC but certainly don't need more than 6 for real-labelled decoding
sc = parallel.pool.Constant(RandStream('Threefry', 'Seed', 'shuffle'));

parfor spl = 1:num_cv_splits
    % set up RNG seed
    stream = sc.Value;
    stream.Substream = spl;

    % Initialize cfg struct and fill in decoding parameters
    cfg = decoding_defaults;
    cfg.plot_design = 0; % no plot with parfor
    cfg.verbose = 0;
    cfg.analysis = 'ROI'; % standard alternatives: 'wholebrain', 'ROI' (pass ROIs in cfg.files.mask, see below)
    cfg.searchlight.radius = vox_radius;
    cfg.searchlight.spherical = 1;
    cfg.scale.method = 'min0max1';
    cfg.scale.estimation = 'all'; % scaling across all data is equivalent to no scaling (i.e. will yield the same results), it only changes the data range which allows libsvm to compute faster
    cfg.results.output = {'accuracy_minus_chance','confusion_matrix'};
    % Set the output directory where data will be saved, e.g. 'c:\exp\results\buttonpress'
    cfg.results.dir = [results_dir '/TDT_results_rep' num2str(spl, '%.2d')];
    if ~isfolder(cfg.results.dir)
        mkdir(cfg.results.dir);
    end

    % Set beta_loc and mask
    beta_loc = [base_folder '/Complex_seq_analysis/' subject_num '/MVPA/' subfolder];
    regressor_names = design_from_spm(beta_loc);
    cfg.files.mask = wholebrain_mask_roi_intersect_out;
    %cfg.files.mask = [base_folder '/Complex_seq_analysis/' subject_num '/MVPA/' subfolder '/mask.nii'];
    %cfg.files.mask = [base_folder '/Complex_seq_analysis/' subject_num '/MVPA/' subfolder '/mask_ants.nii']; % ANTs gray matter mask from antsBrainExtraction.sh

    % Extract all information for the cfg.files structure (labels will be [1 -1] if not changed above)
    cfg = decoding_describe_data(cfg,labelnames_arr,labels_arr,regressor_names,beta_loc);

    %rng('shuffle') % make sure randsample in generating the design below is truly random
    cfg.design = make_design_custom_GNG(cfg, 0, stream); % need to seed RNG in function % not permuted labels yet
    cfg.files.chunk = cfg.design.chunk;
    cfg.design = rmfield(cfg.design, 'chunk');
    %cfg.design.function.name = 'make_design_cv'; % ?
    results = decoding(cfg);
    gzip([cfg.results.dir '/res_accuracy_minus_chance.nii']); % really just one value in whole image
    delete([cfg.results.dir '/res_accuracy_minus_chance.nii']);
end

delete(gcp('nocreate'))

tEnd = toc(tStart);
disp(['Running ' num2str(num_cv_splits) ' CV splits on subject ' num2str(subject_num) ' took ' num2str(tEnd, '%.2f') ' seconds. Done.'])



perm_designs = make_design_permutation(cfg,num_permutations,0);
% Only want permuted labels from permuted designs, not train or test, will replace those later





end
