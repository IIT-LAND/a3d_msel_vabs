function result = runPLS_smri(ANALYSIS,data_type)
% runPLS_smri
%
% ANALYSIS = 'ELO' or 'A3D'
% data_type = 'ct' or 'sa'
%

addpath /Users/mlombardo/Dropbox/matlab/spm12;
addpath /Users/mlombardo/Dropbox/matlab/Pls/plsgui;
addpath /Users/mlombardo/Dropbox/matlab/Pls/plscmd;

rootpath = '/Users/mlombardo/Dropbox/Manuscripts/AUTISMS_3D_MSEL_VABS_Toddlers/Autisms3D_fMRI/pls';
datapath = fullfile(rootpath,'data','smri');
resultpath = fullfile(rootpath,'results','smri');

if strcmp(ANALYSIS,'ELO')
    MEfiles = {fullfile(datapath,'ASDGood_all_MEdata.txt'), ...
        fullfile(datapath,'ASDPoor_all_MEdata.txt'), ...
        fullfile(datapath,'TD_all_MEdata.txt')};

    Brain_files = {fullfile(datapath,sprintf('ASDGood_all_cortical_%s_sexAdj_HaglerGeneticTemplate.txt',data_type)), ...
        fullfile(datapath,sprintf('ASDPoor_all_cortical_%s_sexAdj_HaglerGeneticTemplate.txt',data_type)), ...
        fullfile(datapath,sprintf('TD_all_cortical_%s_sexAdj_HaglerGeneticTemplate.txt',data_type))};
elseif strcmp(ANALYSIS,'A3D')
    MEfiles = {fullfile(datapath,'HIGH_all_MEdata.txt'), ...
        fullfile(datapath,'LOW_all_MEdata.txt'), ...
        fullfile(datapath,'TD_all_MEdata.txt')};

    Brain_files = {fullfile(datapath,sprintf('HIGH_all_cortical_%s_sexAdj_HaglerGeneticTemplate.txt',data_type)), ...
        fullfile(datapath,sprintf('LOW_all_cortical_%s_sexAdj_HaglerGeneticTemplate.txt',data_type)), ...
        fullfile(datapath,sprintf('TD_all_cortical_%s_sexAdj_HaglerGeneticTemplate.txt',data_type))};
end % if strcmp(ANALYSIS,'ELO')

nregions_per_hemi = 12; nregions_total = nregions_per_hemi*2; 

fname2save = fullfile(resultpath,sprintf('%s_%s.mat',ANALYSIS,data_type));


%% read in text files
% read in ME files
for ifile = 1:length(MEfiles)
    fname = MEfiles{ifile};
    tab2use = readtable(fname);
    ME{ifile} = table2array(tab2use(:,2:end));
    me_subids{ifile} = tab2use(:,1);
end % for ifile

me_names = tab2use.Properties.VariableNames(2:end);

% read in brain vol files
for ifile = 1:length(Brain_files)
    fname = Brain_files{ifile};
    tab2use = readtable(fname);
    BrainVols{ifile} = table2array(tab2use(:,2:end));
    brainvols_subids{ifile} = tab2use(:,1);
end % for ifile

brainreg_names = tab2use.Properties.VariableNames(2:end);

tmp_array = [];
for ime = 1:length(ME)
    tmp_array = [tmp_array;ME{ime}];
end % for ime
MEstacked = tmp_array;

%% Prepare input arguments for pls_analysis.m
datamat_lst = BrainVols; % initialize datamat_lst
for icell = 1:length(datamat_lst)
    num_subj_lst(icell) = size(datamat_lst{icell},1); % initialize num_subj_lst
end % for icell

num_cond = 1; %{ones(1,length(datamat_lst))};  % initialize num_cond

% specify option structure
option.method = 3; %[1] | 2 | 3 | 4 | 5 | 6
option.num_perm = 10000; %( single non-negative integer )
option.is_struct = 0;%[0] | 1
option.num_split = 0; %( single non-negative integer )
option.num_boot = 10000; %( single non-negative integer )
option.clim = 95; %( [95] single number between 0 and 100 )
option.stacked_behavdata = MEstacked;
option.cormode = 0; %[0] | 2 | 4 | 6
option.boot_type = 'strat'; %['strat'] | 'nonstrat'


%% run pls_analysis.m
result = pls_analysis(datamat_lst, num_subj_lst, num_cond, option);
result.me_names = me_names;

%% Compute percentage of cross-block covariance
result.crossblockCovPercent = result.s.^2/sum(result.s.^2);

%% fix p-values
result.perm_result.sprob = (result.perm_result.sp+1)./(result.perm_result.num_perm+1);

%% find the significant LVs
result.sigLVs = find(result.perm_result.sprob<=0.05);
result.sigLVs_pvals = result.perm_result.sprob(result.sigLVs);
disp(sprintf('Significant LVs: LV %d',result.sigLVs));


%% compute brain BSR
for i = 1:size(result.boot_result.compare_u,2)
    result.brain_bsr(:,i) = result.boot_result.compare_u(:,i)./result.boot_result.u_se(:,i);
end
result.brainreg_names = brainreg_names';


%% compute bootstrap CIs
cis2use = {[2.5,97.5],[0.5,99.5]};
for i = 1:length(result.sigLVs)
    LVnum = result.sigLVs(i);
    for j = 1:length(cis2use)
        ci_bounds = cis2use{j};
        clear asd_good_bootres asd_poor_bootres td_bootres asd_good_ci asd_poor_ci td_ci asd_good_corr asd_poor_corr td_corr asd_good asd_poor td all_data;
        ORDERBY = 'TD';
        mod_names_orig = me_names;
        
        mod_names = repmat(mod_names_orig',length(MEfiles),1);
        
        asd_good_idx = 1:length(mod_names_orig);
        asd_poor_idx = (length(mod_names_orig)+1):(length(mod_names_orig)*2);
        td_idx = ((length(mod_names_orig)*2)+1):(length(mod_names_orig)*3);
        
        asd_good_bootres = squeeze(result.boot_result.distrib(asd_good_idx,LVnum,:));
        asd_poor_bootres = squeeze(result.boot_result.distrib(asd_poor_idx,LVnum,:));
        td_bootres = squeeze(result.boot_result.distrib(td_idx,LVnum,:));
        
        asd_good_ci = prctile(asd_good_bootres',ci_bounds)';
        asd_poor_ci = prctile(asd_poor_bootres',ci_bounds)';
        td_ci = prctile(td_bootres',ci_bounds)';
        
        asd_good_corr = result.boot_result.orig_corr(asd_good_idx,LVnum);
        asd_poor_corr = result.boot_result.orig_corr(asd_poor_idx,LVnum);
        td_corr = result.boot_result.orig_corr(td_idx,LVnum);
        
        if strcmp(ORDERBY,'TD')
            [idx, plot_order] = sort(td_corr,'ascend');
            plot_order(plot_order) = 1:length(mod_names_orig);
        elseif strcmp(ORDERBY,'Good')
            [idx, plot_order] = sort(asd_good_corr,'ascend');
            plot_order(plot_order) = 1:length(mod_names_orig);
        elseif strcmp(ORDERBY,'Poor')
            [idx, plot_order] = sort(asd_poor_corr,'ascend');
            plot_order(plot_order) = 1:length(mod_names_orig);
        end
        
        asd_good = [asd_good_corr asd_good_ci];
        asd_poor = [asd_poor_corr asd_poor_ci];
        td = [td_corr td_ci];
        
        all_data = [asd_good; asd_poor; td];
        all_data(:,4) = (sign(all_data(:,1)) == sign(all_data(:,2))) & (sign(all_data(:,1)) == sign(all_data(:,3)));
        all_data(:,5) = repmat(plot_order,3,1);
        
        if strcmp(ANALYSIS,'ELO')
            group_labels = [repmat({'Good'},length(mod_names_orig),1);repmat({'Poor'},length(mod_names_orig),1);repmat({'TD'},length(mod_names_orig),1)];
        elseif strcmp(ANALYSIS,'A3D')
            group_labels = [repmat({'Type II'},length(mod_names_orig),1);repmat({'Type I'},length(mod_names_orig),1);repmat({'TD'},length(mod_names_orig),1)];
        end % if strcmp(ANALYSIS,'ELO')
        tab2write = cell2table([group_labels, mod_names num2cell(all_data)], ...
            'VariableNames',{'Grp','ModName','corr','lo_lim','up_lim','nonzero','plot_order'});
        
        [f_path, f_name, f_ext] = fileparts(fname2save);
        file2save = fullfile(f_path,sprintf('%s_%s_MEcorr_bootCI4plotting_LV%d_ci%d.csv',ANALYSIS,data_type,LVnum,ci_bounds(2)-ci_bounds(1)));
        writetable(tab2write,file2save,'FileType','text','delimiter',',');
    end
end

%% output brain BSR
for i = 1:length(result.sigLVs)
    bsr_names{i} = sprintf('BSR_LV%d',result.sigLVs(i));
end

if strcmp(data_type, 'ct')
    feature_names = [repmat({'Cortical Thickness'},nregions_per_hemi*2,1)];
elseif strcmp(data_type, 'sa')
    feature_names = [repmat({'Surface Area'},nregions_per_hemi*2,1)];
end % if strcmp(data_type, 'ct')

hemi_names = [repmat({'LH'},nregions_per_hemi,1);repmat({'RH'},nregions_per_hemi,1)];
tab2write = cell2table([brainreg_names', feature_names, hemi_names, num2cell(result.brain_bsr(:,result.sigLVs))],'VariableNames',[{'brainreg','feature','hemisphere'},bsr_names]);
writetable(tab2write, fullfile(f_path,sprintf('%s_%s_brainBSR4plotting.csv',ANALYSIS,data_type)),'FileType','text','Delimiter',',');

%% save result
save(fname2save,'result');

end % function runPLS_smri




