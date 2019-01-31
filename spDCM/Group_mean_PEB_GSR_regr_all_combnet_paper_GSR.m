function Group_mean_PEB_GSR_regr_all_combnet_paper_GSR(name_ROI_def,ntwrk_abbrev2,procedure,N_prec_comp,SPM_dir,Work_dir)

for network_number=1:length(ntwrk_abbrev2)
    
    tmp=0;
    
    disp(network_number);
    for number_dataset=1:4
        
        [dataset,number_subject,single_band,slice_time_seconds]=Dataset_info_paper_GSR(number_dataset);
        for subject=1:number_subject
            clear PEB PEB_all_Basic PEB_all_GSR_regr;
            
            try
                load([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/Basic/' name_ROI_def '/Full_model/PEB_all_mean_comp_GSR_regr_single_' ntwrk_abbrev2{network_number} '.mat']);
                PEB_all_Basic=PEB;
            catch
                PEB_all_Basic.Ep=[];
            end
            
            try
                load([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/GSR_regr/' name_ROI_def '/Full_model/PEB_all_mean_comp_GSR_regr_single_' ntwrk_abbrev2{network_number} '.mat']);
                PEB_all_GSR_regr=PEB;
            catch
                PEB_all_GSR_regr.Ep=[];
            end
            
            if any(isnan(PEB_all_Basic.Ep))||any(isnan(PEB_all_GSR_regr.Ep))
                disp(['Dataset: ' dataset '; Subject: ' num2str(subject) '; Procedure: ' procedure]);
                disp(['Network: ' ntwrk_abbrev2{network_number}]);
                disp(['Parameter type: All']);
                continue;
            end
            
            tmp=tmp+1;
            
            try
                load([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_all_mean_comp_GSR_regr_single_' ntwrk_abbrev2{network_number} '.mat']);
            catch ME
                
                tmp=tmp-1;
                disp('PEB group');
                disp(['Dataset: ' dataset '; Subject: ' num2str(subject) '; Procedure: ' procedure]);
                disp(['Network: ' ntwrk_abbrev2{network_number}]);
                disp(['Parameter type: All']);
                if strcmp(ME.identifier,'MATLAB:load:couldNotReadFile')
                    continue;
                else
                    save([Work_dir '/Results_paper_GSR/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_group/ERROR_all_mean_' dataset '_subject_' num2str(subject) '_single_' ntwrk_abbrev2{network_number} '_comp_GSR_regr.mat'],'ME');
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
    save([Work_dir '/Results_paper_GSR/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_group/PEB_all_mean_comp_GSR_regr_group_' cell2mat(N_prec_comp) '_' ntwrk_abbrev2{network_number} '.mat'],'PEB','DCM','PEB_group')
    clear PEB_group;
    
end

end