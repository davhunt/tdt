%%%%%%%%%%%% Modified by David for use with our Brainevo lab grammatical/ungrammatical data
%%%%%%%%%%%% save outputs to .nii files
% function p = stats_permutation(n_correct,reference,tail)
% 
% Function to calculate p-value for a permutation test.
% The permutation sample needs to have been calculated already, e.g. using
% make_design_permutation.m
% The minimal p-value provided by this function is 1/(n+1) where n is the
% number of existing permutations. For example, for 1000 permutations, the
% minimal p-value is 0.00999.
% The results of this function are valid only if dependencies that existed
% in the original results also exist in the permutation. Example: When
% permuting labels and having multiple chunks (e.g. runs), permutation is
% only valid within chunk, not between chunks (otherwise, you get inflated
% p-values). Such dependencies need to remain in the permutations.
%
%   INPUT:
%       n_correct: measured decoding value (e.g. how many were judged correctly, accuracy, etc.) (nx1 vector)
%       reference: same measure from permutation samples (nxm matrix)
%       tail: 'left','right','both' ('both' will pick the smaller p-value)
%
%   OUTPUT:
%       p: probability of finding such a value
%       z: p-values converted to z-values (for external reference)

% 14/07/30 Martin Hebart

% TODO: if numerically imprecise numbers are entered, we might end up with
% wrong results when n_correct = reference often, because it might falsely
% be perceived as too large or too small. Consider rounding results to a
% certain number of digits

%function [p,z] = stats_permutation(n_correct,reference,tail)

function p = stats_permutation(results_folder, perms_folder, nreps, nperms, tail)

% Ex:
% results_folder: /Volumes/data1/Complex_seq_analysis/014/MVPA/durations_unsmoothed/TDT_results
% perms_folder: /Volumes/data1/Complex_seq_analysis/014/MVPA/durations_unsmoothed/perm
% nreps: 10, number of repetitions for each permutation, same as number of repetitions as found in "TDT_results"
% nperms: 100, total permutations in perm folder, matches "m" in reference
% tail: 'right'
%%%%%%%%%% build n_correct, reference matricces from input dirs

addpath(genpath('/Users/lab/Downloads/matlab_nifti_tools'));

avgs_nii = load_untouch_nii([results_folder '/TDT_results_avg_allreps/res_accuracy_minus_chance.nii']); % from stats_average_repetitions
orig_sz = size(avgs_nii.img);
n_correct = avgs_nii.img(:);

n = size(n_correct, 1);
reference = zeros(n, nperms);
for i = 1:nperms
    p_avg = zeros(n, 1);
    for j = 1:nreps
        nif = load_untouch_nii([perms_folder '/perm' sprintf('%.3d',i) '/results_rep' sprintf('%.2d',j) '/res_accuracy_minus_chance.nii']);
        p_avg(:) = p_avg(:) + nif.img(:);
    end
    p_avg = p_avg / nreps;
    reference(:,i) = p_avg;
end

%%%%%%%%%%


sz = size(n_correct);
sz_ref = size(reference);

if ~any(sz==1) || ndims(sz)>2
    error('Wrong dimensionality of input variable n_correct. Please enter an nx1 vector.')
end

if sz(1)==1 && sz(2)~=1
    n_correct = n_correct';
    reference = reference';
    warningv('STATS_PERMUTATION:FLIP','Input variable n_correct was an 1xn vector. Flipping n_correct and reference...')
end

if sz_ref(1) ~= sz(1)
    error('Input variable reference needs to be nxm variable, and input variable n_correct an nx1 variable.')
end

% get percentage of cases more extreme than the measured value
switch lower(tail)
    
    case 'left'

        if exist('bsxfun','builtin') % New method for Matlab 7.4+ (fast)
            p = (1/(sz_ref(2)+1))*(sum(bsxfun(@ge,n_correct,reference),2)+1);
        else
            p = (1/(sz_ref(2)+1))*(sum(repmat(n_correct,1,sz_ref(2))>=reference,2)+1);
        end
        
    case 'right'
        
        if exist('bsxfun','builtin') % New method for Matlab 7.4+ (fast)
            p = (1/(sz_ref(2)+1))*(sum(bsxfun(@le,n_correct,reference),2)+1);
        else
            p = (1/(sz_ref(2)+1))*(sum(repmat(n_correct,1,sz_ref(2))<=reference,2)+1);
        end

    case 'both'

        if exist('bsxfun','builtin') % New method for Matlab 7.4+ (fast)
            p1 = (1/(sz_ref(2)+1))*(sum(bsxfun(@ge,n_correct,reference),2)+1);
            p2 = (1/(sz_ref(2)+1))*(sum(bsxfun(@le,n_correct,reference),2)+1);
            p = 2*min(p1,p2); % pick the smaller of both and multiply by 2 (because test is two-sided)
            p = min(p,1); % since it's a Fisher p-value it can lead to values larger than 1, this corrects for this fact.
        else
            n_correctmat = repmat(n_correct,1,sz_ref(2));
            p1 = (1/(sz_ref(2)+1))*(sum(n_correctmat>=reference,2)+1);
            p2 = (1/(sz_ref(2)+1))*(sum(n_correctmat<=reference,2)+1);
            p = 2*min(p1,p2);
            p = min(p,1); % since it's a Fisher p-value it can lead to values larger than 1, this corrects for this fact.
        end
        
    otherwise
        error('unknown method ''%s'' for input variable ''tail''',tail)
end

p_vals_img = reshape(p, orig_sz);
out_nii = avgs_nii;
out_nii.img = p_vals_img;
save_untouch_nii(out_nii, [results_folder '/pvals_map_post_permutation.nii.gz']);
