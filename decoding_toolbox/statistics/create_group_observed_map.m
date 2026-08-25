% for future reference
% also, should this be t-values?

>> subj_masks = {};
>> for sub = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24]
sub = sprintf('%02d', sub);
subj_masks{end+1} = load_untouch_nii(['/N/project/brainevo/Complex_seq_analysis/0' sub '/MVPA/breakpoint_rts_unsmoothed_native_antsmask_new_bps/mask_mni.nii']);
end
out_img = zeros(91,109,91);

>> for i = 1:91
for j = 1:109
for k = 1:91
val = 0;
if mask_nii.img(i,j,k) == 1
n = 0;
for s = 1:22
if subj_masks{s}.img(i,j,k) == 1
val = val + subjs{s}.img(i,j,k);
n = n+1;
end
end
out_img(i,j,k) = val / n;
end
end
end
end
>> tmp_nii = mask_nii;
>> tmp_nii.img = out_img;
>> save_untouch_nii(tmp_nii, '/N/slate/brainevo/Implicit_Learning/Complex_seq_analysis/breakpoint_rts_unsmoothed_native_antsmask_new_bps_group_results/group_observed_decoding_accuracies.nii.gz')
>> 
