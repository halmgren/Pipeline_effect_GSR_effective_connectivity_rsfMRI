function PEB_group_mean_paper_GSR_DatasetHCP(SPM_dir,Work_dir)

all_ROI_defs={'Smith'};
all_procedure_names={'Basic','GSR_regr'};                        %all_procedure_names={'Basic','GSR','GSR_regr'};

for number_ROI_def=1:length(all_ROI_defs)
    
    name_ROI_def=all_ROI_defs{number_ROI_def};
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Give regions name and coördinates
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [ROI_list]=Define_ROIs_paper_GSR(name_ROI_def);
    
    [ROI_list2]=Define_comb_ROIs_paper_GSR(name_ROI_def);
    
    tmp1=0;
    
    for VOI_number=1:size(ROI_list,1)
        ntwrk=ROI_list{VOI_number,1}(1:3);
        
        if VOI_number>1 && strcmp(ROI_list{VOI_number,1}(1:3),ROI_list{VOI_number-1,1}(1:3))
            ntwrk_size(tmp1)=ntwrk_size(tmp1)+1;
            continue
            
        else
            tmp1=tmp1+1;
            ntwrk_size(tmp1)=1;
            ntwrk_name{tmp1}=ROI_list{VOI_number,1}(1:3);
        end
    end
    
    for number_procedure=1:length(all_procedure_names)
        procedure=all_procedure_names{number_procedure};
        
        for N_prec_comp={'single'}    %N_prec_comp={'single','fields','all'}
            
            %%%%%%%%%%%%%%%%%%%
            %Compare GSR-Basic
            %%%%%%%%%%%%%%%%%%%
            if strcmp(procedure,'Basic')||strcmp(procedure,'GSR_regr')
                
                mkdir([Work_dir '/Results_paper_GSR/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_group/']);
                
                %%%%%%%%%%
                %A-matrix
                %%%%%%%%%%
                for network_number=1:length(ntwrk_name)
                    
                    tmp=0;
                    tmp2=0;
                    
                    disp(network_number);
                    for number_dataset=5
                        
                        [dataset,number_subject,single_band,slice_time_seconds]=Dataset_info_paper_GSR(number_dataset);
                        for subject=1:number_subject
                            clear GCM flag flag_comp;
                            clear Posterior_estimates_var Posterior_estimates_max Posterior_estimates_par isnan(Posterior_estimates_mot Posterior_estimates_thr
                            
                            load([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_summary/DCM/'  procedure '/' name_ROI_def '/Full_model/GCM_' ntwrk_name{network_number} '_full_estim.mat']);
                            load([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/QC/Above_treshold_marks_' ntwrk_name{network_number} '.mat']);
                            
                            flag=zeros(1,length(GCM));
                            for diagn=1:length(GCM)
                                if ~isnan(Posterior_estimates_var(1,1,diagn))||~isnan(Posterior_estimates_max(1,1,diagn))||~isnan(Posterior_estimates_par(1,1,diagn))||~isnan(Posterior_estimates_mot(1,1,diagn))||~isnan(Posterior_estimates_thr(1,1,diagn))
                                    flag(diagn)=1;
                                end
                            end
                            
                            if strcmp(procedure,'Basic')
                                procedure_comp='GSR_regr';
                                clear Posterior_estimates_var Posterior_estimates_max Posterior_estimates_par isnan(Posterior_estimates_mot Posterior_estimates_thr;
                                load([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/'  procedure_comp '/' name_ROI_def '/Full_model/QC/Above_treshold_marks_' ntwrk_name{network_number} '.mat']);
                            end
                            
                            if strcmp(procedure,'GSR_regr')
                                procedure_comp='Basic';
                                clear Posterior_estimates_var Posterior_estimates_max Posterior_estimates_par isnan(Posterior_estimates_mot Posterior_estimates_thr;
                                load([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/'  procedure_comp '/' name_ROI_def '/Full_model/QC/Above_treshold_marks_' ntwrk_name{network_number} '.mat']);
                            end
                            
                            flag_comp=zeros(1,length(GCM));
                            for diagn=1:length(GCM)
                                if ~isnan(Posterior_estimates_var(1,1,diagn))||~isnan(Posterior_estimates_max(1,1,diagn))||~isnan(Posterior_estimates_par(1,1,diagn))||~isnan(Posterior_estimates_mot(1,1,diagn))||~isnan(Posterior_estimates_thr(1,1,diagn))
                                    flag_comp(diagn)=1;
                                end
                            end
                            
                            if ~(flag==1||flag_comp==1)
                                tmp=tmp+1;
                                GCM_full{tmp}=GCM{1};
                                subjects_included(tmp)=subject;
                            else
                                tmp2=tmp2+1;
                                Omitted_subjects(tmp2)=subject;
                            end
                        end
                    end
                    
                    clear GCM
                    
                    save([Work_dir '/Results_paper_GSR/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_group/Omitted_subjects_HCP_' ntwrk_name{network_number} '_comp_GSR_regr.mat'],'Omitted_subjects');
                    save([Work_dir '/Results_paper_GSR/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_group/Subjects_included_datasetHCP_GSR_regr_' cell2mat(N_prec_comp) '_' ntwrk_name{network_number} '.mat'],'subjects_included');
                    
                    GCM_full=GCM_full';
                    if length(GCM_full)<8
                        continue
                    end
                    
                    %PEB settings
                    M = struct();
                    M.alpha = 1;
                    M.beta  = 16;
                    M.hE    = 0;
                    M.hC    = 1/16;
                    M.Q     = cell2mat(N_prec_comp);
                    M.X=ones(length(GCM_full),1);
                    
%                     lastwarn('');
%                     try
%                         %All parameters together
%                         disp('all param in same model');
%                         [PEB DCM]=spm_dcm_peb(GCM_full,M,'all');
%                         
%                         save([Work_dir '/Results_paper_GSR/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_group/PEB_all_mean_comp_GSR_regr_HCP_' cell2mat(N_prec_comp) '_' ntwrk_name{network_number} '.mat'],'PEB','DCM','GCM_full');
%                     catch ME
%                         disp(['Dataset: ' dataset '; Subject: ' num2str(subject) '; Procedure: ' procedure]);
%                         disp('all');
%                         save([Work_dir '/Results_paper_GSR/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_group/ERROR_PEB_all_mean_comp_GSR_regr_HCP_' cell2mat(N_prec_comp) '_' ntwrk_name{network_number} '.mat'],'ME');
%                     end
%                     
%                     if ~isempty(lastwarn)
%                         [msgstr,msgid]=lastwarn;
%                         save([Work_dir '/Results_paper_GSR/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_group/WARNING_PEB_all_mean_comp_GSR_regr_HCP_' cell2mat(N_prec_comp) '_' ntwrk_name{network_number} '.mat'],'msgstr','msgid');
%                     end
                    
                    clear PEB DCM;
                    
                    lastwarn('');
                    try
                        %Connectivity
                        disp('Connectivity matrix only');
                        [PEB DCM]=spm_dcm_peb(GCM_full,M,'A');
                        
                        save([Work_dir '/Results_paper_GSR/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_group/PEB_A_mean_comp_GSR_regr_HCP_' cell2mat(N_prec_comp) '_' ntwrk_name{network_number} '.mat'],'PEB','DCM','GCM_full');
                    catch ME
                        disp(['Dataset: ' dataset '; Subject: ' num2str(subject) '; Procedure: ' procedure]);
                        disp('A');
                        save([Work_dir '/Results_paper_GSR/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_group/ERROR_PEB_A_mean_comp_GSR_regr_HCP_' cell2mat(N_prec_comp) '_' ntwrk_name{network_number} '.mat'],'ME');
                    end
                    
                    if ~isempty(lastwarn)
                        [msgstr,msgid]=lastwarn;
                        save([Work_dir '/Results_paper_GSR/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_group/WARNING_PEB_A_mean_comp_GSR_regr_HCP_' cell2mat(N_prec_comp) '_' ntwrk_name{network_number} '.mat'],'msgstr','msgid');
                    end
                    
                    clear PEB DCM;
%                     
%                     lastwarn('');
%                     try
%                         %Hemodynamics
%                         disp('Hemodynamics only');
%                         [PEB DCM]=spm_dcm_peb(GCM_full,M,{'transit','decay','epsilon'});
%                         
%                         save([Work_dir '/Results_paper_GSR/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_group/PEB_hemodyna_mean_comp_GSR_regr_HCP_' cell2mat(N_prec_comp) '_' ntwrk_name{network_number} '.mat'],'PEB','DCM','GCM_full');
%                     catch ME
%                         disp(['Dataset: ' dataset '; Subject: ' num2str(subject) '; Procedure: ' procedure]);
%                         disp('Hemodyna');
%                         save([Work_dir '/Results_paper_GSR/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_group/ERROR_PEB_hemodyna_mean_comp_GSR_regr_HCP_' cell2mat(N_prec_comp) '_' ntwrk_name{network_number} '.mat'],'ME');
%                     end
%                     
%                     if ~isempty(lastwarn)
%                         [msgstr,msgid]=lastwarn;
%                         save([Work_dir '/Results_paper_GSR/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_group/WARNING_PEB_hemodyna_mean_comp_GSR_regr_HCP_' cell2mat(N_prec_comp) '_' ntwrk_name{network_number} '.mat'],'msgstr','msgid');
%                     end
%                     
%                     clear PEB DCM;
                    
                    lastwarn('');
                    try
                        %Spectral noise
                        disp('Spectral noise only');
                        [PEB DCM]=spm_dcm_peb(GCM_full,M,{'a','b','c'});
                        
                        save([Work_dir '/Results_paper_GSR/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_group/PEB_fluct_mean_comp_GSR_regr_HCP_' cell2mat(N_prec_comp) '_' ntwrk_name{network_number} '.mat'],'PEB','DCM','GCM_full');
                    catch ME
                        disp(['Dataset: ' dataset '; Subject: ' num2str(subject) '; Procedure: ' procedure]);
                        disp('fluct');
                        save([Work_dir '/Results_paper_GSR/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_group/ERROR_PEB_fluct_mean_comp_GSR_regr_HCP_' cell2mat(N_prec_comp) '_' ntwrk_name{network_number} '.mat'],'ME');
                    end
                    
                    if ~isempty(lastwarn)
                        [msgstr,msgid]=lastwarn;
                        save([Work_dir '/Results_paper_GSR/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_group/WARNING_PEB_fluct_mean_comp_GSR_regr_HCP_' cell2mat(N_prec_comp) '_' ntwrk_name{network_number} '.mat'],'msgstr','msgid');
                    end
                    
                    clear PEB DCM Omitted_subjects subjects_included;
                    
                    clear GCM_full Posterior_estimates_var Lower_bound_ci_var Higher_bound_ci_var Length_error_bar_var Posterior_estimates_max Lower_bound_ci_max Higher_bound_ci_max Length_error_bar_max Posterior_estimates_par Lower_bound_ci_par Higher_bound_ci_par Length_error_bar_par Posterior_estimates_mot Lower_bound_ci_mot Higher_bound_ci_mot Length_error_bar_mot Posterior_estimates_thr Lower_bound_ci_thr Higher_bound_ci_thr Length_error_bar_thr Posterior_estimates;
                    
                end
                
            end
        end
    end
    
    if exist('ROI_list2')
        
        tmp=1;
        tmp2=1;
        clear ntwrk_name2 ntwrk_name3;
        for VOI_number=1:size(ROI_list2,1)
            ntwrk=ROI_list2{VOI_number,1}(1:3);
            region=ROI_list2{VOI_number,1}(5:end);
            
            if VOI_number==1
                tmp2=1;
                ntwrk_abbrev2{tmp}=ROI_list2{VOI_number,1}(1:3);
            end
            
            if VOI_number>1 && ~strcmp(ROI_list2{VOI_number,1}(1:3),ROI_list2{VOI_number-1,1}(1:3))
                tmp2=1;
                tmp=tmp+1;
                ntwrk_abbrev2{tmp}=ROI_list2{VOI_number,1}(1:3);
            end
            
        end
        
        for number_procedure=1:length(all_procedure_names)
            procedure=all_procedure_names{number_procedure};
            
            
            for N_prec_comp={'single'}    %N_prec_comp={'single','fields'}
                
                %%%%%%%%%%%%%%%%%%%
                %Compare GSR-Basic
                %%%%%%%%%%%%%%%%%%%
                if strcmp(procedure,'Basic')||strcmp(procedure,'GSR_regr')
                    
                    
                    mkdir([Work_dir '/Results_paper_GSR/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_group/']);
                    
                    
                    %%%%%%%%%%
                    %A-matrix
                    %%%%%%%%%%
                    for network_number=1:length(ntwrk_abbrev2)
                        
                        tmp=0;
                        tmp2=0;
                        
                        disp(network_number);
                        for number_dataset=5
                            
                            [dataset,number_subject,single_band,slice_time_seconds]=Dataset_info_paper_GSR(number_dataset);
                            for subject=1:number_subject
                                clear GCM flag flag_comp;
                                clear Posterior_estimates_var Posterior_estimates_max Posterior_estimates_par isnan(Posterior_estimates_mot Posterior_estimates_thr
                                
                                load([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_summary/DCM/'  procedure '/' name_ROI_def '/Full_model/GCM_' ntwrk_abbrev2{network_number} '_full_estim.mat']);
                                load([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/QC/Above_treshold_marks_' ntwrk_abbrev2{network_number} '.mat']);
                                
                                flag=zeros(1,length(GCM));
                                for diagn=1:length(GCM)
                                    if ~isnan(Posterior_estimates_var(1,1,diagn))||~isnan(Posterior_estimates_max(1,1,diagn))||~isnan(Posterior_estimates_par(1,1,diagn))||~isnan(Posterior_estimates_mot(1,1,diagn))||~isnan(Posterior_estimates_thr(1,1,diagn))
                                        flag(diagn)=1;
                                    end
                                end
                                
                                if strcmp(procedure,'Basic')
                                    procedure_comp='GSR_regr';
                                    clear Posterior_estimates_var Posterior_estimates_max Posterior_estimates_par isnan(Posterior_estimates_mot Posterior_estimates_thr;
                                    load([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/'  procedure_comp '/' name_ROI_def '/Full_model/QC/Above_treshold_marks_' ntwrk_abbrev2{network_number} '.mat']);
                                end
                                
                                if strcmp(procedure,'GSR_regr')
                                    procedure_comp='Basic';
                                    clear Posterior_estimates_var Posterior_estimates_max Posterior_estimates_par isnan(Posterior_estimates_mot Posterior_estimates_thr;
                                    load([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/'  procedure_comp '/' name_ROI_def '/Full_model/QC/Above_treshold_marks_' ntwrk_abbrev2{network_number} '.mat']);
                                end
                                
                                flag_comp=zeros(1,length(GCM));
                                for diagn=1:length(GCM)
                                    if ~isnan(Posterior_estimates_var(1,1,diagn))||~isnan(Posterior_estimates_max(1,1,diagn))||~isnan(Posterior_estimates_par(1,1,diagn))||~isnan(Posterior_estimates_mot(1,1,diagn))||~isnan(Posterior_estimates_thr(1,1,diagn))
                                        flag_comp(diagn)=1;
                                    end
                                end
                                
                                if ~(flag==1||flag_comp==1)
                                    tmp=tmp+1;
                                    GCM_full{tmp}=GCM{1};
                                    subjects_included(tmp)=subject;
                                else
                                    tmp2=tmp2+1;
                                    Omitted_subjects(tmp2)=subject;
                                end
                            end
                        end
                        
                        save([Work_dir '/Results_paper_GSR/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_group/Omitted_subjects_HCP_' ntwrk_abbrev2{network_number} '_comp_GSR_regr.mat'],'Omitted_subjects');
                        save([Work_dir '/Results_paper_GSR/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_group/Subjects_included_datasetHCP_GSR_regr_' cell2mat(N_prec_comp) '_' ntwrk_abbrev2{network_number} '.mat'],'subjects_included');
                        
                        GCM_full=GCM_full';
                        
                        %PEB settings
                        M = struct();
                        M.alpha = 1;
                        M.beta  = 16;
                        M.hE    = 0;
                        M.hC    = 1/16;
                        M.Q     = cell2mat(N_prec_comp);
                        M.X=ones(length(GCM_full),1);
                        
                        %                 lastwarn('');
                        %                 try
                        %                     %All parameters together
                        %                     disp('all param in same model');
                        %                     [PEB DCM]=spm_dcm_peb(GCM_full,M,'all');
                        %
                        %                     save([Work_dir '/Results_paper_GSR/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_group/PEB_all_mean_comp_GSR_regr_HCP_' cell2mat(N_prec_comp) '_' ntwrk_abbrev2{network_number} '.mat'],'PEB','DCM','GCM_full','-v7.3');
                        %                 catch ME
                        %                     disp(['Dataset: ' dataset '; Subject: ' num2str(subject) '; Procedure: ' procedure]);
                        %                     disp('all');
                        %                     save([Work_dir '/Results_paper_GSR/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_group/ERROR_PEB_all_mean_comp_GSR_regr_HCP_' cell2mat(N_prec_comp) '_' ntwrk_abbrev2{network_number} '.mat'],'ME');
                        %                 end
                        %
                        %                 if ~isempty(lastwarn)
                        %                     [msgstr,msgid]=lastwarn;
                        %                     save([Work_dir '/Results_paper_GSR/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_group/WARNING_PEB_all_mean_comp_GSR_regr_HCP_' cell2mat(N_prec_comp) '_' ntwrk_abbrev2{network_number} '.mat'],'msgstr','msgid');
                        %                 end
                        %
                        %                 clear PEB DCM;
                        
                        lastwarn('');
                        try
                            %Connectivity
                            disp('Connectivity matrix only');
                            [PEB DCM]=spm_dcm_peb(GCM_full,M,'A');
                            
                            save([Work_dir '/Results_paper_GSR/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_group/PEB_A_mean_comp_GSR_regr_HCP_' cell2mat(N_prec_comp) '_' ntwrk_abbrev2{network_number} '.mat'],'PEB','DCM','GCM_full','-v7.3');
                        catch ME
                            disp(['Dataset: ' dataset '; Subject: ' num2str(subject) '; Procedure: ' procedure]);
                            disp('A');
                            save([Work_dir '/Results_paper_GSR/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_group/ERROR_PEB_A_mean_comp_GSR_regr_HCP_' cell2mat(N_prec_comp) '_' ntwrk_abbrev2{network_number} '.mat'],'ME');
                        end
                        
                        if ~isempty(lastwarn)
                            [msgstr,msgid]=lastwarn;
                            save([Work_dir '/Results_paper_GSR/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_group/WARNING_PEB_A_mean_comp_GSR_regr_HCP_' cell2mat(N_prec_comp) '_' ntwrk_abbrev2{network_number} '.mat'],'msgstr','msgid');
                        end
                        
                        clear PEB DCM;
%                         
%                         lastwarn('');
%                         try
%                             %Hemodynamics
%                             disp('Hemodynamics only');
%                             [PEB DCM]=spm_dcm_peb(GCM_full,M,{'transit','decay','epsilon'});
%                             
%                             save([Work_dir '/Results_paper_GSR/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_group/PEB_hemodyna_mean_comp_GSR_regr_HCP_' cell2mat(N_prec_comp) '_' ntwrk_abbrev2{network_number} '.mat'],'PEB','DCM','GCM_full','-v7.3');
%                         catch ME
%                             disp(['Dataset: ' dataset '; Subject: ' num2str(subject) '; Procedure: ' procedure]);
%                             disp('Hemodyna');
%                             save([Work_dir '/Results_paper_GSR/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_group/ERROR_PEB_hemodyna_mean_comp_GSR_regr_HCP_' cell2mat(N_prec_comp) '_' ntwrk_abbrev2{network_number} '.mat'],'ME');
%                         end
%                         
%                         if ~isempty(lastwarn)
%                             [msgstr,msgid]=lastwarn;
%                             save([Work_dir '/Results_paper_GSR/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_group/WARNING_PEB_hemodyna_mean_comp_GSR_regr_HCP_' cell2mat(N_prec_comp) '_' ntwrk_abbrev2{network_number} '.mat'],'msgstr','msgid');
%                         end
%                         
%                         clear PEB DCM;
                        
                        lastwarn('');
                        try
                            %Spectral noise
                            disp('Spectral noise only');
                            [PEB DCM]=spm_dcm_peb(GCM_full,M,{'a','b','c'});
                            
                            save([Work_dir '/Results_paper_GSR/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_group/PEB_fluct_mean_comp_GSR_regr_HCP_' cell2mat(N_prec_comp) '_' ntwrk_abbrev2{network_number} '.mat'],'PEB','DCM','GCM_full','-v7.3');
                        catch ME
                            disp(['Dataset: ' dataset '; Subject: ' num2str(subject) '; Procedure: ' procedure]);
                            disp('fluct');
                            save([Work_dir '/Results_paper_GSR/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_group/ERROR_PEB_fluct_mean_comp_GSR_regr_HCP_' cell2mat(N_prec_comp) '_' ntwrk_abbrev2{network_number} '.mat'],'ME');
                        end
                        
                        if ~isempty(lastwarn)
                            [msgstr,msgid]=lastwarn;
                            save([Work_dir '/Results_paper_GSR/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_group/WARNING_PEB_fluct_mean_comp_GSR_regr_HCP_' cell2mat(N_prec_comp) '_' ntwrk_abbrev2{network_number} '.mat'],'msgstr','msgid');
                        end
                        
                        clear PEB DCM Omitted_subjects subjects_included;
                        
                        clear GCM_full Posterior_estimates_var Lower_bound_ci_var Higher_bound_ci_var Length_error_bar_var Posterior_estimates_max Lower_bound_ci_max Higher_bound_ci_max Length_error_bar_max Posterior_estimates_par Lower_bound_ci_par Higher_bound_ci_par Length_error_bar_par Posterior_estimates_mot Lower_bound_ci_mot Higher_bound_ci_mot Length_error_bar_mot Posterior_estimates_thr Lower_bound_ci_thr Higher_bound_ci_thr Length_error_bar_thr Posterior_estimates;
                        
                        
                    end
                end
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Compute hemispheric assymetry: Connectivity-matrix only (GSR_regr)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for number_ROI_def=1:length(all_ROI_defs)
    
    name_ROI_def=all_ROI_defs{number_ROI_def};
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Give regions name and coördinates
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [ROI_list]=Define_ROIs_paper_GSR(name_ROI_def);
    
    [ROI_list2]=Define_comb_ROIs_paper_GSR(name_ROI_def);
    
    tmp1=0;
    
    for VOI_number=1:size(ROI_list,1)
        ntwrk=ROI_list{VOI_number,1}(1:3);
        
        if VOI_number>1 && strcmp(ROI_list{VOI_number,1}(1:3),ROI_list{VOI_number-1,1}(1:3))
            ntwrk_size(tmp1)=ntwrk_size(tmp1)+1;
            continue
            
        else
            tmp1=tmp1+1;
            ntwrk_size(tmp1)=1;
            ntwrk_name{tmp1}=ROI_list{VOI_number,1}(1:3);
        end
    end
    
    for number_procedure=1:length(all_procedure_names)
        procedure=all_procedure_names{number_procedure};
        if strcmp(procedure,'Basic')||strcmp(procedure,'GSR_regr')
            
            for N_prec_comp={'single'}    %N_prec_comp={'single','fields','all'}
                tmp=0;
                
                for VOI_number=1:size(ROI_list,1)
                    ntwrk=ROI_list{VOI_number,1}(1:3);
                    
                    if VOI_number>1 && strcmp(ROI_list{VOI_number,1}(1:3),ROI_list{VOI_number-1,1}(1:3))
                        ntwrk_size(tmp)=ntwrk_size(tmp)+1;
                        continue
                        
                    else
                        tmp=tmp+1;
                        ntwrk_size(tmp)=1;
                        ntwrk_name{tmp}=ROI_list{VOI_number,1}(1:3);
                    end
                end
                
                for network_number=1:length(ntwrk_name)
                    disp(ntwrk_name{network_number});
                    try
                        load([Work_dir '/Results_paper_GSR/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_group/PEB_A_mean_comp_GSR_regr_HCP_' cell2mat(N_prec_comp) '_' ntwrk_name{network_number} '.mat'],'PEB','DCM','GCM_full');
                        
                        %Salience network
                        if strcmp(ntwrk_name{network_number},'SAL')
                            T = 0;
                            C = [0 0 0 0 0 1 0 0 1 1 1 0 0 1 1 -1 -1 -1 0 0 -1 -1 -1 0 0 zeros(1,length(PEB.Ep)-25)]'/6;
                            c = C'*spm_vec(PEB.Ep);
                            v = C'*PEB.Cp*C;
                            PP   = 1-spm_Ncdf(T,c,v);
                        end
                        
                        %Somatomotor network
                        if strcmp(ntwrk_name{network_number},'SMR')
                            T = 0;
                            C = [0 0 0 1 0 1 -1 -1 0 zeros(1,length(PEB.Ep)-9)]'/2;
                            c = C'*spm_vec(PEB.Ep);
                            v = C'*PEB.Cp*C;
                            PP   = 1-spm_Ncdf(T,c,v);
                        end
                        
                        %DMN
                        if strcmp(ntwrk_name{network_number},'DMN')
                            T = 0;
                            C = [0 0 0 0 1 0 1 1 0 0 0 0 -1 -1 -1 0 zeros(1,length(PEB.Ep)-16)]'/3;
                            c = C'*spm_vec(PEB.Ep);
                            v = C'*PEB.Cp*C;
                            PP   = 1-spm_Ncdf(T,c,v);
                        end
                        
                        mean_diff=c;
                        var_of_sum=v;
                        posterior_probability=PP;
                        
                        save([Work_dir '/Results_paper_GSR/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_group/Lateralization_index_A_comp_GSR_regr_HCP_' cell2mat(N_prec_comp) '_' ntwrk_name{network_number} '.mat'],'mean_diff','var_of_sum','posterior_probability');
                        clear mean_diff var_of_sum c v C PP T PEB DCM posterior_probability;
                    catch ME
                        disp('Hemispheric Assymetry')
                        disp(['Dataset: ' dataset '; Subject: ' num2str(subject) '; Procedure: ' procedure]);
                        if strcmp(ME.identifier,'MATLAB:load:couldNotReadFile')
                            continue
                        else
                            save([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/ERROR_Lateralization_index_A_comp_GSR_regr_' cell2mat(N_prec_comp) '_' ntwrk_name{network_number} '.mat'],'ME');
                        end
                    end
                end
            end
        end
        
        if exist('ROI_list2')
       
        
            tmp=1;
            tmp2=1;
            clear ntwrk_name2 ntwrk_name3;
            for VOI_number=1:size(ROI_list2,1)
                ntwrk=ROI_list2{VOI_number,1}(1:3);
                region=ROI_list2{VOI_number,1}(5:end);
                
                if VOI_number==1
                    tmp2=1;
                    ntwrk_abbrev2{tmp}=ROI_list2{VOI_number,1}(1:3);
                end
                
                if VOI_number>1 && ~strcmp(ROI_list2{VOI_number,1}(1:3),ROI_list2{VOI_number-1,1}(1:3))
                    tmp2=1;
                    tmp=tmp+1;
                    ntwrk_abbrev2{tmp}=ROI_list2{VOI_number,1}(1:3);
                end
                
            end
            
            %%%%%%%%%%%%%%%%%%%%
            %GSR_regr - basic
            %%%%%%%%%%%%%%%%%%%
            if strcmp(procedure,'Basic')||strcmp(procedure,'GSR_regr')
                
                %%%%%%%%%%
                %A-matrix
                %%%%%%%%%%
                for N_prec_comp={'single'}                % Also possible with different number of precision components using N_prec_comp={'single','fields','all'}
                    for network_number=1:length(ntwrk_abbrev2)
                        load([Work_dir '/Results_paper_GSR/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_group/PEB_A_mean_comp_GSR_regr_HCP_' cell2mat(N_prec_comp) '_' ntwrk_abbrev2{network_number} '.mat'],'PEB')
                        
                        %Salience network
                        if strcmp(ntwrk_abbrev2{network_number},'TG2')
                            %from SAL to SMR
                            T1 = 0;
                            C1 = [zeros(1,5) ones(1,3) zeros(1,3) zeros(1,5) ones(1,3) zeros(1,3) zeros(1,5) ones(1,3) zeros(1,3) zeros(1,5) ones(1,3) zeros(1,3) zeros(1,5) ones(1,3) zeros(1,3) zeros(1,66)]';
                            C1 = C1./sum(C1);
                            c1 = C1'*spm_vec(PEB.Ep);
                            v1 = C1'*PEB.Cp*C1;
                            PP1   = 1-spm_Ncdf(T1,c1,v1);
                            
                            %from SAL to DMN
                            T2 = 0;
                            C2 = [zeros(1,5) zeros(1,3) ones(1,3) zeros(1,5) zeros(1,3) ones(1,3) zeros(1,5) zeros(1,3) ones(1,3) zeros(1,5) zeros(1,3) ones(1,3) zeros(1,5) zeros(1,3) ones(1,3) zeros(1,66)]';
                            C2 = C2./sum(C2);
                            c2 = C2'*spm_vec(PEB.Ep);
                            v2 = C2'*PEB.Cp*C2;
                            PP2   = 1-spm_Ncdf(T2,c2,v2);
                            
                            %from SMR to SAL
                            T3 = 0;
                            C3 = [zeros(1,55) ones(1,5) zeros(1,3) zeros(1,3) ones(1,5) zeros(1,3) zeros(1,3) ones(1,5) zeros(1,3) zeros(1,3) zeros(1,33)]';
                            C3 = C3./sum(C3);
                            c3 = C3'*spm_vec(PEB.Ep);
                            v3 = C3'*PEB.Cp*C3;
                            PP3   = 1-spm_Ncdf(T3,c3,v3);
                            
                            %from SMR to DMN
                            T4 = 0;
                            C4 = [zeros(1,55) zeros(1,5) zeros(1,3) ones(1,3) zeros(1,5) zeros(1,3) ones(1,3) zeros(1,5) zeros(1,3) ones(1,3) zeros(1,33)]';
                            C4 = C4./sum(C4);
                            c4 = C4'*spm_vec(PEB.Ep);
                            v4 = C4'*PEB.Cp*C4;
                            PP4   = 1-spm_Ncdf(T4,c4,v4);
                            
                            %from DMN to SAL
                            T5 = 0;
                            C5= [zeros(1,88) ones(1,5) zeros(1,3) zeros(1,3) ones(1,5) zeros(1,3) zeros(1,3) ones(1,5) zeros(1,3) zeros(1,3)]';
                            C5 = C5./sum(C5);
                            c5 = C5'*spm_vec(PEB.Ep);
                            v5 = C5'*PEB.Cp*C5;
                            PP5   = 1-spm_Ncdf(T5,c5,v5);
                            
                            %from DMN to SMR
                            T6 = 0;
                            C6 = [zeros(1,88) zeros(1,5) ones(1,3) zeros(1,3) zeros(1,5) ones(1,3) zeros(1,3) zeros(1,5) ones(1,3) zeros(1,3)]';
                            C6 = C6./sum(C6);
                            c6 = C6'*spm_vec(PEB.Ep);
                            v6 = C6'*PEB.Cp*C6;
                            PP6   = 1-spm_Ncdf(T6,c6,v6);
                            
                            save([Work_dir '/Results_paper_GSR/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_group/BN_effect_A_comp_GSR_regr_HCP_' cell2mat(N_prec_comp) '_' ntwrk_abbrev2{network_number} '.mat'],'c1','v1','PP1','c2','v2','PP2','c3','v3','PP3','c4','v4','PP4','c5','v5','PP5','c6','v6','PP6');
                            clear mean_diff var_of_sum c v C PP T PEB posterior_probability;
                        end
                        
                        
                    end
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %All parameters in one PEB model: NOT CORRECT
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %
                    %             for network_number=1:length(ntwrk_abbrev2)
                    %                 load([Work_dir '/Results_paper_GSR/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_group/PEB_all_mean_comp_GSR_group_' ntwrk_abbrev2e{network_number} '.mat'],'PEB')
                    %
                    %                 %Salience network
                    %                 if strcmp(ntwrk_name{network_number},'SAL')
                    %                     T = 0;
                    %                     C = [0 0 0 0 0 1 0 0 1 1 1 0 0 1 1 -1 -1 -1 0 0 -1 -1 -1 0 0 zeros(1,length(PEB.Ep)-25)]'/6;
                    %                     c = C'*spm_vec(PEB.Ep);
                    %                     v = C'*PEB.Cp*C;
                    %                     PP   = 1-spm_Ncdf(T,c,v);
                    %                 end
                    %
                    %                 mean_diff=c;
                    %                 var_of_sum=v;
                    %                 posterior_probability=PP;
                    %
                    %                 save([Work_dir '/Results_paper_GSR/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_group/BN_effect_all_comp_GSR_group_' ntwrk_abbrev2{network_number} '.mat'],'mean_diff','var_of_sum','posterior_probability');
                    %                 clear mean_diff var_of_sum c v C PP T PEB posterior_probability;
                    %             end
                end
            end
        end
            
    end
end

end
