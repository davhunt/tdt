function design = make_design_custom_crossval(cfg_combined, varargin)

% This function builds the design struct for the implicit grammar -> sentence reading task cross classification
% TODO: parameterize all the 56's 32's and 88's
% But all subjects have 56 G/UG implicit reading betas and 32 G/UG sentence reading betas I believe (?)

% Won't ever need these arguments I don't think?
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

% Just one "chunk" (iteration), all of the complex seq trials in training data and all of the sentence reading trials
% in the testing data
design.label = cfg_combined.files.label;
design.train = [ones(56,1); zeros(32,1)];
design.test = [zeros(56,1); ones(32,1)];
design.chunk = ones(88,1);
design.set = ones(88,1);


end