% function p = stats_binomial(n_correct,n_samples,chancelevel,tail)
% 
function [] = stats_average_repetitions(reps_folder)

addpath(genpath('/Users/lab/Downloads/matlab_nifti_tools'));

all_dirs = dir(reps_folder);

if exist([reps_folder '/TDT_results_avg_allreps']) ~= 7;
    mkdir([reps_folder '/TDT_results_avg_allreps']);
end

firstFlag = 1;
n = 0; % total number of niftis/CV reps
for i = 1:length(all_dirs)
    if regexp(all_dirs(i).name, '^TDT_results_rep[0-9]+$')
        nif = load_untouch_nii([reps_folder '/' all_dirs(i).name '/res_accuracy_minus_chance.nii']);
        if firstFlag
            sz = size(nif.img);
            firstFlag = 0;
            sum_all_niis = zeros(sz);
        end
        sum_all_niis(:,:,:) = sum_all_niis(:,:,:) + nif.img(:,:,:);
        n = n+1;
    end
end

sum_all_niis(:,:,:) = sum_all_niis(:,:,:) / n;

nif.img = sum_all_niis;

save_untouch_nii(nif, [reps_folder '/TDT_results_avg_allreps/res_accuracy_minus_chance.nii']);