function runPLS_fmri(ANALYSIS)
% runPLS_fmri
%
% ANALYSIS = 'ELO' or 'A3D'
%

addpath /Users/mlombardo/Dropbox/matlab/spm12;
addpath /Users/mlombardo/Dropbox/matlab/Pls/plsgui;
addpath /Users/mlombardo/Dropbox/matlab/Pls/plscmd;

rootpath = '/Users/mlombardo/Dropbox/Manuscripts/AUTISMS_3D_MSEL_VABS_Toddlers/Autisms3D_fMRI/pls';
codepath = fullfile(rootpath,'code');
datapath = fullfile(rootpath,'data');
resultpath = fullfile(rootpath,'results','fmri');

%% Prepare input arguments for pls_analysis.m
cd(codepath);
[datamat_lst, data_fname, sublist2use, brainmask_idx] = preparefMRIData(ANALYSIS);

wgcna_data = readtable(fullfile(datapath,'module_eigengenes.csv'));
wgcna_data.Properties.RowNames=wgcna_data.subjectId;
stackedMEdata = [];

for icell = 1:length(datamat_lst)
    num_subj_lst(icell) = size(datamat_lst{icell},1); % initialize num_subj_lst
    stackedMEdata = [stackedMEdata; table2array(wgcna_data(sublist2use{icell},2:end))];
end % for icell

num_cond = 1; %{ones(1,length(datamat_lst))};  % initialize num_cond

% specify option structure
option.method = 3; %[1] | 2 | 3 | 4 | 5 | 6
option.num_perm = 10000; %( single non-negative integer )
option.is_struct = 0;%[0] | 1
option.num_split = 0; %( single non-negative integer )
option.num_boot = 10000; %( single non-negative integer )
option.clim = 95; %( [95] single number between 0 and 100 )
option.stacked_behavdata = stackedMEdata;
option.cormode = 0; %[0] | 2 | 4 | 6
option.boot_type = 'strat'; %['strat'] | 'nonstrat'


%% run pls_analysis.m
result = pls_analysis(datamat_lst, num_subj_lst, num_cond, option);

%% Compute percentage of cross-block covariance
result.crossblockCovPercent = result.s.^2/sum(result.s.^2);

%% find the significant LVs
result.sigLVs = find(result.perm_result.sprob<=0.05);
result.sigLVs_pvals = result.perm_result.sprob(result.sigLVs);
disp(sprintf('Significant LVs: LV %d',result.sigLVs));


%% compute brain BSR
for i = 1:size(result.boot_result.compare_u,2)
    result.brain_bsr(:,i) = result.boot_result.compare_u(:,i)./result.boot_result.u_se(:,i);
end

%% save result
fname = fullfile(resultpath,sprintf('pls_%s.mat',ANALYSIS));
save(fname,'result');


%% save out brain BSR image
brain_mask = spm_read_vols(spm_vol(fullfile(datapath,'fmri','whole_brain_maps','finalmask4pls.img')));
brain_mask = brain_mask==1;
BSRimg = zeros(size(brain_mask));
LVnum = 1;
[x,y,z] = ind2sub(size(brain_mask),brainmask_idx);
for i = 1:length(x)
    BSRimg(x(i),y(i),z(i)) = result.brain_bsr(i,LVnum);
end
BSRimg_rev = BSRimg.*-1;
BSRimg_thr = BSRimg; BSRimg_thr(BSRimg_thr<1.96 & BSRimg_thr>-1.96) = 0;
BSRimg_thr_rev = BSRimg_rev; BSRimg_thr_rev(BSRimg_thr_rev<1.96 & BSRimg_thr_rev>-1.96) = 0;

V = spm_vol(fullfile(datapath,'fmri','whole_brain_maps','Good_A5V2H.img'));
V.fname = sprintf('%s_bsrimg_LV%d.img',ANALYSIS,LVnum);
V.private.dat.fname = V.fname;
cd(resultpath);
spm_write_vol(V,BSRimg);

V.fname = sprintf('%s_bsrimg_LV%d_REV.img',ANALYSIS,LVnum);
V.private.dat.fname = V.fname;
spm_write_vol(V,BSRimg_rev);

V.fname = sprintf('%s_bsrimg_LV%d_thresh1.96.img',ANALYSIS,LVnum);
V.private.dat.fname = V.fname;
spm_write_vol(V,BSRimg_thr);

V.fname = sprintf('%s_bsrimg_LV%d_REV_thresh1.96.img',ANALYSIS,LVnum);
V.private.dat.fname = V.fname;
spm_write_vol(V,BSRimg_thr_rev);


%% compute bootstrap CIs
cd(codepath);
LVnum = 1;
resultmat = fname;
[tab2write] = compute_bootstrap_ci_pls(resultmat, LVnum, ANALYSIS);

end % function runPLS_fmri



%%
function [datamat_lst, data_fname, sublist2use, brainmask_idx] = preparefMRIData(ANALYSIS)

addpath /Users/mlombardo/Dropbox/matlab/spm12;

rootpath = '/Users/mlombardo/Dropbox/Manuscripts/AUTISMS_3D_MSEL_VABS_Toddlers/Autisms3D_fMRI';
codepath = fullfile(rootpath,'pls','code');
fmripath = fullfile(rootpath,'pls','data','fmri','whole_brain_maps');
pheno_data = fullfile(rootpath,'FilteringData_UCSDE_MRI','results','psc_data_nn_neu.csv');
pheno_data = readtable(pheno_data);

mask = spm_read_vols(spm_vol(fullfile(fmripath,'mask.img')));
mask = mask==1;

if strcmp(ANALYSIS,'ELO')

    grps2use = {'Good','Poor','TD'};
    datamat_lst = cell(1,length(grps2use));
    sublist2use = cell(1,length(grps2use));

    mask2use = ismember(pheno_data.subtype_lang,grps2use);
    pheno_data_sub = pheno_data(mask2use,:);
    
    for igrp = 1:length(grps2use)
        
        tmp_sub = pheno_data_sub(ismember(pheno_data_sub.subtype_lang,grps2use(igrp)),:);
        sublist2use{igrp} = sort(tmp_sub.subjectId);

%         datamat_lst{igrp} = zeros(length(sublist2use{igrp}), sum(mask(:)));
        datamat_lst{igrp} = zeros(length(sublist2use{igrp}), length(mask(:)));
        data_fname{igrp} = cell(length(sublist2use{igrp}), 1);

        for isub = 1:length(sublist2use{igrp})
            
            tmp_fname = dir(fullfile(fmripath,sprintf('*%s.img',sublist2use{igrp}{isub})));
            data = spm_read_vols(spm_vol(fullfile(fmripath, tmp_fname.name)));
%             datamat_lst{igrp}(isub,:) = data(mask);
            datamat_lst{igrp}(isub,:) = data(:);
            data_fname{igrp}{isub,1} = tmp_fname.name;        
        
        end % for isub

        nan_mask(igrp,:) = sum(isnan(datamat_lst{igrp}),1)>0;
    end % for igrp

elseif strcmp(ANALYSIS,'A3D')

    grps2use = {'HIGH','LOW','TD'};
    datamat_lst = cell(1,length(grps2use));
    sublist2use = cell(1,length(grps2use));

    mask2use = ismember(pheno_data.subtype3d,grps2use);
    pheno_data_sub = pheno_data(mask2use,:);

    for igrp = 1:length(grps2use)

        tmp_sub = pheno_data_sub(ismember(pheno_data_sub.subtype3d,grps2use(igrp)),:);
        sublist2use{igrp} = sort(tmp_sub.subjectId);

%         datamat_lst{igrp} = zeros(length(sublist2use{igrp}), sum(mask(:)));
        datamat_lst{igrp} = zeros(length(sublist2use{igrp}), length(mask(:)));
        data_fname{igrp} = cell(length(sublist2use{igrp}), 1);

        for isub = 1:length(sublist2use{igrp})
        
            tmp_fname = dir(fullfile(fmripath,sprintf('*%s.img',sublist2use{igrp}{isub})));
            data = spm_read_vols(spm_vol(fullfile(fmripath, tmp_fname.name)));
%             datamat_lst{igrp}(isub,:) = data(mask);
            datamat_lst{igrp}(isub,:) = data(:);
            data_fname{igrp}{isub,1} = tmp_fname.name;        
        
        end % for isub
        
        nan_mask(igrp,:) = sum(isnan(datamat_lst{igrp}),1)>0;
    end % for igrp

end % if strcmp(ANALYSIS,'ELO')

% final constraint to grab brain voxels that are in the brain mask or are
% not NaN in any of the groups
nobrain_mask = mask(:)==0;
all_masks = [nan_mask; nobrain_mask'];
final_nobrain_mask = sum(all_masks,1)>0;

final_brain_mask = ~final_nobrain_mask;

datamat_lst{1} = datamat_lst{1}(:,final_brain_mask);
datamat_lst{2} = datamat_lst{2}(:,final_brain_mask);
datamat_lst{3} = datamat_lst{3}(:,final_brain_mask);

final_brain_mask2export = reshape(final_brain_mask, size(mask));
brainmask_idx = find(final_brain_mask2export);
V = spm_vol(fullfile(fmripath,'mask.img'));
V.fname = 'finalmask4pls.img';
V.private.dat.fname = V.fname;
cd(fmripath);
spm_write_vol(V, final_brain_mask2export);
cd(rootpath);

end % function preparefMRIData