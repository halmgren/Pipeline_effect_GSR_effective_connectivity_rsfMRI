function Group_mean_PEB_GSR_regr_hemodyna_paper_GSR(name_ROI_def,ntwrk_name,procedure,N_prec_comp,SPM_dir,Work_dir)

for network_number=1:length(ntwrk_name)
    
    tmp=0;
    
    disp(network_number);
    for number_dataset=1:4
        
        [dataset,number_subject,single_band,slice_time_seconds]=Dataset_info_paper_GSR(number_dataset);
        for subject=1:number_subject
            clear PEB PEB_A_Basic PEB_hemodyna_Basic PEB_fluct_Basic PEB_A_GSR_regr PEB_hemodyna_GSR_regr PEB_fluct_GSR_regr;
            
            %%%%%%%%%%%%%
            %If any of the parts of the DCM parameters shows NaNs (ill-conditioned matrix), then this subject will be discarded from these particular analyses
            %for both analyses types
            %%%%%%%%%%%%%
            try
                load([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/Basic/' name_ROI_def '/Full_model/PEB_A_mean_comp_GSR_regr_single_' ntwrk_name{network_number} '.mat']);
                PEB_A_Basic=PEB;
            catch
                PEB_A_Basic.Ep=[];
            end
            
            try
                load([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/Basic/' name_ROI_def '/Full_model/PEB_hemodyna_mean_comp_GSR_regr_single_' ntwrk_name{network_number} '.mat']);
                PEB_hemodyna_Basic=PEB;
            catch
                PEB_hemodyna_Basic.Ep=[];
            end
            
            try
                load([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/Basic/' name_ROI_def '/Full_model/PEB_fluct_mean_comp_GSR_regr_single_' ntwrk_name{network_number} '.mat']);
                PEB_fluct_Basic=PEB;
            catch
                PEB_fluct_Basic.Ep=[];
            end
            
            try
                load([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/GSR_regr/' name_ROI_def '/Full_model/PEB_A_mean_comp_GSR_regr_single_' ntwrk_name{network_number} '.mat']);
                PEB_A_GSR_regr=PEB;
            catch
                PEB_A_GSR_regr.Ep=[];
            end
            
            try
                load([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/GSR_regr/' name_ROI_def '/Full_model/PEB_hemodyna_mean_comp_GSR_regr_single_' ntwrk_name{network_number} '.mat']);
                PEB_hemodyna_GSR_regr=PEB;
            catch
                PEB_hemodyna_GSR_regr.Ep=[];
            end
            
            try
                load([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/GSR_regr/' name_ROI_def '/Full_model/PEB_fluct_mean_comp_GSR_regr_single_' ntwrk_name{network_number} '.mat']);
                PEB_fluct_GSR_regr=PEB;
            catch
                PEB_fluct_GSR_regr.Ep=[];
            end
            
            if any(isnan(PEB_A_Basic.Ep))||any(isnan(PEB_hemodyna_Basic.Ep))||any(isnan(PEB_fluct_Basic.Ep))||any(isnan(PEB_A_GSR_regr.Ep))||any(isnan(PEB_hemodyna_GSR_regr.Ep))||any(isnan(PEB_fluct_GSR_regr.Ep))
                disp(['Dataset: ' dataset '; Subject: ' num2str(subject) '; Procedure: ' procedure]);
                disp(['Network: ' ntwrk_name{network_number}]);
                disp(['Parameter type: Hemodyna']);
                continue;
            end
            
            
            tmp=tmp+1;
            
            try
                load([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_hemodyna_mean_comp_GSR_regr_single_' ntwrk_name{network_number} '.mat']);
            catch ME
                tmp=tmp-1;
                disp('PEB group')
                disp(['Dataset: ' dataset '; Subject: ' num2str(subject) '; Procedure: ' procedure]);
                disp(['Network: ' ntwrk_name{network_number}]);
                disp(['Parameter type: Hemodyna']);
                if strcmp(ME.identifier,'MATLAB:load:couldNotReadFile')
                    continue;
                else
                    save([Work_dir '/Results_paper_GSR/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_group/ERROR_hemodyna_mean_' dataset '_subject_' num2str(subject) '_single_' ntwrk_name{network_number} '_comp_GSR_regr.mat'],'ME');
                    continue;
                end
            end
            
            PEB_group{tmp}=PEB;
            
            
            clear PEB;
        end
    end
    
    M = struct();
    M.alpha = 1;
    M.beta  = 16;
    M.hE    = 0;
    M.hC    = 1/16;
    M.Q     = cell2mat(N_prec_comp);
    M.X=ones(length(PEB_group),1);
    
    [PEB DCM]=spm_dcm_peb_of_peb(PEB_group',M,[]);
    save([Work_dir '/Results_paper_GSR/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_group/PEB_hemodyna_mean_comp_GSR_regr_group_' cell2mat(N_prec_comp) '_' ntwrk_name{network_number} '.mat'],'PEB','DCM','PEB_group')
    clear PEB_group;
    
end

end