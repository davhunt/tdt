function [] = stats_average_repetitions(reps_folder, type)

% Function to average all the CV splits seen in a folder ("TDT_results") or nested set of folders (each permutation folder within an overall "perm" folder)
% Modified for "ROI" analysis, where there is one accuracy value per repetition for the whole ROI, instead of "searchlight"

% reps_folder should be the folder that contains the "permXXX" folders, in the case of average repetitions for each permutation, or
% the folder "TDT_results" that contains the "TDT_results_repXX" folders
% type can be "real" (real data, just one directory) or "perm" (permuted data, many directories)

if exist('/Users/lab/Downloads/matlab_nifti_tools') == 7
    addpath(genpath('/Users/lab/Downloads/matlab_nifti_tools')); % On Habilis
elseif exist('/N/slate/brainevo/Implicit_Learning/matlab_nifti_tools') == 7
    addpath(genpath('/N/slate/brainevo/Implicit_Learning/matlab_nifti_tools')); % On Quartz/HPC
else
    error('Cannot find matlab_nifti_tools dir, exiting.');
end

if strcmp(type,'real')
    all_dirs = dir(reps_folder);

    if exist([reps_folder '/TDT_results_avg_allreps']) ~= 7
        mkdir([reps_folder '/TDT_results_avg_allreps']);
    end

    firstFlag = 1;
    n = 0; % total number of niftis/CV reps
    for i = 1:length(all_dirs)
        if regexp(all_dirs(i).name, '^TDT_results_rep[0-9]+$')
            %nif = load_untouch_nii([reps_folder '/' all_dirs(i).name '/res_accuracy_minus_chance.nii.gz']);
            load([reps_folder '/' all_dirs(i).name '/res_accuracy_minus_chance.mat']);
            val = results.accuracy_minus_chance.output;
            if firstFlag
                %sz = size(nif.img);
                firstFlag = 0;
                %sum_all_niis = zeros(sz);
                sum_all_vals = 0;
            end
            %sum_all_niis(:,:,:) = sum_all_niis(:,:,:) + nif.img(:,:,:);
            sum_all_vals = sum_all_vals + val;
            n = n+1;
        end
    end

    %sum_all_niis(:,:,:) = sum_all_niis(:,:,:) / n;
    avg_all_reps = sum_all_vals / n;

    %nif.img = sum_all_niis;

    %save_untouch_nii(nif, [reps_folder '/TDT_results_avg_allreps/res_accuracy_minus_chance.nii.gz']);
    save([reps_folder '/TDT_results_avg_allreps/res_accuracy_minus_chance.mat'], 'avg_all_reps');

elseif strcmp(type,'perm')
    all_dirs = dir(reps_folder);
    n_perms = 0;
    for i = 1:size(all_dirs,1)
        if startsWith(all_dirs(i).name, 'perm')
            n_perms = n_perms + 1;
        end
    end
    for i = 1:n_perms
        if exist([reps_folder '/perm' num2str(i, '%03d')]) ~= 7
            fprintf(['Don''t see a folder ' reps_folder '/perm%03d, check your perm folder input\n'], i);
        end
        %if exist([reps_folder '/perm' num2str(i, '%03d') '/perm' num2str(i, '%03d') '_avg_allreps']) ~= 7
        if exist([reps_folder '/perm' num2str(i, '%03d') '/results_avg_allreps']) ~= 7
            mkdir([reps_folder '/perm' num2str(i, '%03d') '/results_avg_allreps']);
        end
        firstFlag = 1;
        n = 0; % total number of niftis/CV reps
        all_reps = dir([reps_folder '/perm' num2str(i, '%03d')]);
        for j = 1:length(all_reps)
            if regexp(all_reps(j).name, '^results_rep[0-9]+.*$')
                %nif = load_untouch_nii([reps_folder '/perm' num2str(i, '%03d') '/' all_reps(j).name '/res_accuracy_minus_chance.nii.gz']);
                load([reps_folder '/perm' num2str(i, '%03d') '/' all_reps(j).name '/res_accuracy_minus_chance.mat']);
                val = results.accuracy_minus_chance.output;
                if firstFlag
                    %sz = size(nif.img);
                    firstFlag = 0;
                    %sum_all_niis = zeros(sz);
                    sum_all_vals = 0;
                end
                %sum_all_niis(:,:,:) = sum_all_niis(:,:,:) + nif.img(:,:,:);
                sum_all_vals = sum_all_vals + val;
                n = n+1;
            end
        end
        %sum_all_niis(:,:,:) = sum_all_niis(:,:,:) / n;
        avg_all_reps = sum_all_vals / n;

        %nif.img = sum_all_niis;

        %save_untouch_nii(nif, [reps_folder '/perm' num2str(i, '%03d') '/results_avg_allreps/res_accuracy_minus_chance.nii.gz']);
        save([reps_folder '/perm' num2str(i, '%03d') '/results_avg_allreps/res_accuracy_minus_chance.mat'], 'avg_all_reps');


    end


else
    fprintf('Must specify either "real" or "perm," exiting.\n')

end

end
