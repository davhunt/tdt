function [] = gng_decoding_allparams_real_labels(subject_num, num_cv_splits, subfolder)

% Formally called this "initial CV splits", decoding for real-labelled data

% subject_num: Subject number of subject, will be converted to string
% num_cv_splits: Number of different fold assignments/CVs to do per permutation, should be ~ 5 - 10. https://www.sciencedirect.com/science/article/pii/S1053811921004225
% subfolder: Folder under "MVPA" to set as working dir, save results in (e.g. "durations_unsmoothed")
tStart = tic;

% from decoding_template_grammatical_ungrammatical, from decoding_template_Haxby

% This script is a template that can be used for a decoding analysis on 
% brain image data. It is for people who have betas available from an 
% SPM.mat and want to automatically extract the relevant images used for
% classification, as well as corresponding labels and decoding chunk numbers
% (e.g. run numbers). If you don't have this available, then use
% decoding_template_nobetas.m

% Make sure the decoding toolbox and your favorite software (SPM or AFNI)
% are on the Matlab path (e.g. addpath('/home/decoding_toolbox') )
% Are we on Slate or Habilis?
if exist('/N/slate/brainevo/Implicit_Learning') == 7
    base_folder = ['/N/slate/brainevo/Implicit_Learning']; % If on slate
elseif exist('/Volumes/data1') == 7
    base_folder = ['/Volumes/data1']; % Or on habilis
else
    error('Can''t find base folder, are we on habilis or the HPC?');
end

% add TDT and SPM
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

% Set the label names to the regressor names which you want to use for 
% decoding, e.g. 'button left' and 'button right'
% don't remember the names? -> run display_regressor_names(beta_loc)
% infos on '*' (wildcard) or regexp -> help decoding_describe_data
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

myCluster = parcluster('local');
parpool(myCluster.NumWorkers);

sc = parallel.pool.Constant(RandStream('Threefry', 'Seed', 'shuffle'));

parfor spl = 1:num_cv_splits % 10 CV splits reasonable?
    % set up RNG seed
    stream = sc.Value;
    stream.Substream = spl;

    % Set defaults
    cfg = decoding_defaults;
    cfg.plot_design = 0; % no plot in parfor
    cfg.verbose = 0;

    % Set the analysis that should be performed (default is 'searchlight')
    cfg.analysis = 'searchlight'; % standard alternatives: 'wholebrain', 'ROI' (pass ROIs in cfg.files.mask, see below)
    cfg.searchlight.radius = 5; % use searchlight of radius 3 (by default in voxels), see more details below % use 5
    cfg.searchlight.spherical = 1;
    cfg.scale.method = 'min0max1';
    cfg.scale.estimation = 'all'; % scaling across all data is equivalent to no scaling (i.e. will yield the same results), it only changes the data range which allows libsvm to compute faster
    cfg.results.output = {'accuracy_minus_chance','confusion_matrix'};
    % Set the output directory where data will be saved, e.g. 'c:\exp\results\buttonpress'
    cfg.results.dir = [base_folder '/Complex_seq_analysis/' num2str(subject_num, '%.3d') '/MVPA/' subfolder '/TDT_results/TDT_results_rep' num2str(spl, '%.2d')];
    if ~isfolder(cfg.results.dir)
        mkdir(cfg.results.dir)
    end

    % Set the filepath where your SPM.mat and all related betas are, e.g. 'c:\exp\glm\model_button'
    beta_loc = [base_folder '/Complex_seq_analysis/' num2str(subject_num, '%.3d') '/MVPA/' subfolder];

    % Set the filename of your brain mask (or your ROI masks as cell matrix) 
    % for searchlight or wholebrain e.g. 'c:\exp\glm\model_button\mask.img' OR 
    % for ROI e.g. {'c:\exp\roi\roimaskleft.img', 'c:\exp\roi\roimaskright.img'}
    % You can also use a mask file with multiple masks inside that are
    % separated by different integer values (a "multi-mask")
    cfg.files.mask = [base_folder '/Complex_seq_analysis/' num2str(subject_num, '%.3d') '/MVPA/' subfolder '/mask.nii'];

    %% Set additional parameters
    % Set additional parameters manually if you want (see decoding.m or
    % decoding_defaults.m). Below some example parameters that you might want 
    % to use a searchlight with radius 12 mm that is spherical:

    % cfg.searchlight.unit = 'mm';
    % cfg.searchlight.radius = 12; % if you use this, delete the other searchlight radius row at the top!
    %cfg.searchlight.spherical = 1;
    % cfg.verbose = 2; % you want all information to be printed on screen
    % cfg.decoding.train.classification.model_parameters = '-s 0 -t 0 -c 1 -b 0 -q'; 

    % Enable scaling min0max1 (otherwise libsvm can get VERY slow)
    % if you dont need model parameters, and if you use libsvm, use:
    %cfg.scale.method = 'min0max1';
    %cfg.scale.estimation = 'all'; % scaling across all data is equivalent to no scaling (i.e. will yield the same results), it only changes the data range which allows libsvm to compute faster

    % if you like to change the decoding software (default: libsvm):
    % cfg.decoding.software = 'liblinear'; % for more, see decoding_toolbox\decoding_software\. 
    % Note: cfg.decoding.software and cfg.software are easy to confuse.
    % cfg.decoding.software contains the decoding software (standard: libsvm)
    % cfg.software contains the data reading software (standard: SPM/AFNI)

    % Some other cool stuff
    % Check out 
    %   combine_designs(cfg, cfg2)
    % if you like to combine multiple designs in one cfg.

    %% Decide whether you want to see the searchlight/ROI/... during decoding
    %cfg.plot_selected_voxels = 500; % 0: no plotting, 1: every step, 2: every second step, 100: every hundredth step...

    %% Add additional output measures if you like
    % See help decoding_transform_results for possible measures

    % cfg.results.output = {'accuracy_minus_chance', 'AUC'}; % 'accuracy_minus_chance' by default
    %cfg.results.output = {'accuracy_minus_chance','confusion_matrix'};

    % You can also use all methods that start with "transres_", e.g. use
    %   cfg.results.output = {'SVM_pattern'};
    % will use the function transres_SVM_pattern.m to get the pattern from 
    % linear svm weights (see Haufe et al, 2015, Neuroimage)

    %% Nothing needs to be changed below for a standard leave-one-run out cross
    %% validation analysis.

    % The following function extracts all beta names and corresponding run
    % numbers from the SPM.mat
    regressor_names = design_from_spm(beta_loc);

    % Extract all information for the cfg.files structure (labels will be [1 -1] if not changed above)
    cfg = decoding_describe_data(cfg,labelnames_arr,labels_arr,regressor_names,beta_loc);

    % This creates the leave-one-run-out cross validation design:
    %cfg.design = make_design_custom_GNG(cfg);
    cfg.design = make_design_custom_GNG(cfg, 0, stream); % need to seed RNG in function 

    cfg.files.chunk = cfg.design.chunk;
    cfg.design = rmfield(cfg.design,'chunk');

    %cfg.plot_design = 0; % no plot in parfor
    %cfg.verbose = 0;
    % Run decoding
    results = decoding(cfg);
end
delete(gcp('nocreate'))

tEnd = toc(tStart);
disp(['Running ' num2str(num_cv_splits) ' CV splits on subject ' num2str(subject_num) ' took ' num2str(tEnd, '%.2f') ' seconds. Done.'])

end