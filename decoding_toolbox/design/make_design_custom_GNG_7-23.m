function design = make_design_custom_GNG(cfg, varargin)

%
% This function randomly generates a cfg.design for a 7-fold (k-fold) cross-validation for the grammatical non-grammatical non-linguistic Reber grammar task data (beta map data)
% 7 folds of 8 conditions each, paired up by grammaticality, and chosen to have 4 H and 4 L chunk strength
% e.g. a subset might be ["HGa", "HNGa", "LGc", "LNGc", "HGg", "HNGg", "LGn", "LNGn"]
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
    stream = varargin{1};
end


design.function.name = mfilename;
design.function.ver = 'v20140107';

%%%% design.set
design.set = ones(1, 7);
%%%%

%conds_arr = [1:14 29:42];
conds_H_arr = [1:14]; conds_L_arr = [29:42];
subsets = {}; % will be 1 x 7 cell array with each subset in the partition
for i = 1:7
    %subsets{i} = randsample(conds_H_arr, 2, false);
    %subsets{i} = [subsets{i} randsample(conds_L_arr, 2, false)];
    if exist('stream', 'var')
        subsets{i} = randsample(stream, conds_H_arr, 2, false);
        subsets{i} = [subsets{i} randsample(stream, conds_L_arr, 2, false)];
    else
        subsets{i} = randsample(conds_H_arr, 2, false);
        subsets{i} = [subsets{i} randsample(conds_L_arr, 2, false)];
    end
    conds_H_arr = conds_H_arr(~ismember(conds_H_arr,subsets{i}));
    conds_L_arr = conds_L_arr(~ismember(conds_L_arr,subsets{i}));
    subsets{i} = [subsets{i} subsets{i}+14];
    subsets{i} = sort(subsets{i});
end

%%%%%%%% design.label
design.label = zeros(56, 7);
design.label(1:14,:) = 1;
design.label(15:28,:) = -1;
design.label(29:42,:) = 1;
design.label(43:56,:) = -1;
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