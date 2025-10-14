function design = make_design_custom_GNG(cfg, varargin)

%
% This function randomly generates a cfg.design for a 7-fold (k-fold) cross-validation for the grammatical non-grammatical non-linguistic Reber grammar task data (beta map data)
% 7 folds of 8 conditions each, chosen according to grammaticality and chunk strength, with 2 HG, 2 LG, 2 HNG, and 2 LNG
% 
% No longer paired up by grammatical and its nongrammatical counterpart, was originally thinking that was a good idea but maybe not
% 
% e.g. a subset might be ["HGa", "HNGb", "LGc", "LNGl", "HGg", "HNGm", "LGn", "LNGm"]
%
% OUT
%   design.label: matrix with one column for each CV step, containing a
%       label for each image used for decoding (a replication of the vector
%       cfg.files.label across CV steps)
%       56 x 7
%   design.train: binary matrix with one column for each CV step, containing
%       a 1 for each image used for training in this CV step and 0 for all
%       images not used
%       56 x 7
%   design.test: same as in design.train, but this time for all test images (should be inverses, i.e. for each fold all beta maps not used to train classifier should be in test set)
%       56 x 7
%   design.set: 1xn vector, describing the set number of each CV step
%       1 x 7
%   design.function: Information about function used to create design
%       make_design_custom_GNG.m

%addpath(genpath('/Applications/MATLAB_R2023b.app/toolbox/stats'));
if length(varargin) >= 1
    permuted = varargin{1};
else
    permuted = 0;
end
if length(varargin) >= 2
    stream = varargin{2};
end


design.function.name = mfilename;
design.function.ver = 'v20140107';

%%%% design.set
design.set = ones(1, 7);
%%%%

if ~permuted
    conds_H_arr = [1:28]; conds_L_arr = [29:56];
    conds_G_arr = [1:14, 29:42]; conds_NG_arr = [15:28, 43:56];
    conds_HG_arr = [1:14]; conds_HNG_arr = [15:28]; conds_LG_arr = [29:42]; conds_LNG_arr = [43:56];
else
    conds_G_arr = find(cfg.design.label(:,1) == 1)'; conds_NG_arr = find(cfg.design.label(:,1) == -1)'; % uses first column but should all be the same
    conds_HG_arr = randsample(conds_G_arr, 14, false); conds_LG_arr = conds_G_arr(~ismember(conds_G_arr, conds_HG_arr));
    conds_HNG_arr = randsample(conds_NG_arr, 14, false); conds_LNG_arr = conds_NG_arr(~ismember(conds_NG_arr, conds_HNG_arr));
    conds_H_arr = [conds_HG_arr conds_HNG_arr]; conds_L_arr = [conds_LG_arr conds_LNG_arr];
    conds_H_arr = sort(conds_H_arr); conds_L_arr = sort(conds_L_arr);
    %conds_HG_arr_tmp = conds_HG_arr; conds_LG_arr_tmp = conds_LG_arr; conds_HNG_arr_tmp = conds_HNG_arr; conds_LNG_arr_tmp = conds_LNG_arr;
    % ^ need modify-able arrays to remove elements from as they are allocated to subsets %nvm we do not
end

%%%%%%%%%%%%% Now, always have 2 HG, 2 HNG, 2 LG, 2 LNG in each partition, can change later if we want

subsets = {}; % will be 1 x 7 cell array with each subset in the partition
for i = 1:7
    if exist('stream', 'var')
        subsets{i} = randsample(stream, conds_HG_arr, 2, false);
        subsets{i} = [subsets{i} randsample(stream, conds_LG_arr, 2, false)];
        subsets{i} = [subsets{i} randsample(stream, conds_HNG_arr, 2, false)];
        subsets{i} = [subsets{i} randsample(stream, conds_LNG_arr, 2, false)];
    else
        subsets{i} = randsample(conds_HG_arr, 2, false);
        subsets{i} = [subsets{i} randsample(conds_LG_arr, 2, false)];
        subsets{i} = [subsets{i} randsample(conds_HNG_arr, 2, false)];
        subsets{i} = [subsets{i} randsample(conds_LNG_arr, 2, false)];
    end
    conds_HG_arr = conds_HG_arr(~ismember(conds_HG_arr,subsets{i}));
    conds_LG_arr = conds_LG_arr(~ismember(conds_LG_arr,subsets{i}));
    conds_HNG_arr = conds_HNG_arr(~ismember(conds_HNG_arr,subsets{i}));
    conds_LNG_arr = conds_LNG_arr(~ismember(conds_LNG_arr,subsets{i}));
end





%conds_arr = [1:14 29:42];
%conds_H_arr = [1:14]; conds_L_arr = [29:42];
%conds_H_arr = []; conds_L_arr = [];
%% create arrays of high chunk strength and low chunk strength trials
%% never mind
%subsets = {}; % will be 1 x 7 cell array with each subset in the partition
%for i = 1:7
%    %subsets{i} = randsample(conds_H_arr, 2, false);
%    %subsets{i} = [subsets{i} randsample(conds_L_arr, 2, false)];
%    if exist('stream', 'var')
%        subsets{i} = randsample(stream, conds_H_arr, 2, false);
%        subsets{i} = [subsets{i} randsample(stream, conds_L_arr, 2, false)];
%    else
%        subsets{i} = randsample(conds_H_arr, 2, false);
%        subsets{i} = [subsets{i} randsample(conds_L_arr, 2, false)];
%    end
%    conds_H_arr = conds_H_arr(~ismember(conds_H_arr,subsets{i}));
%    conds_L_arr = conds_L_arr(~ismember(conds_L_arr,subsets{i}));
%    subsets{i} = [subsets{i} subsets{i}+14];
%    subsets{i} = sort(subsets{i});
%end

%%%%%%%% design.label
design.label = zeros(56, 7);
if ~permuted
    design.label(1:14,:) = 1;
    design.label(15:28,:) = -1;
    design.label(29:42,:) = 1;
    design.label(43:56,:) = -1;
else
    design.label(conds_G_arr, :) = 1;
    design.label(conds_NG_arr, :) = -1;
end
%%%%%%%

%%%%%%% design.train, design.test
design.train = ones(56,7); design.test = zeros(56,7);
for i = 1:7
    design.train(subsets{i}, i) = 0;
    design.test(subsets{i}, i) = 1;
end
%%%%%%%%

%%%%%% design.files.chunk
chunks = [];
for cond = 1:56
    for i = 1:7
        if ismember(cond,subsets{i})
            chunks = [chunks; i];
        end
    end
end
design.chunk = chunks;
%%%%%%


end