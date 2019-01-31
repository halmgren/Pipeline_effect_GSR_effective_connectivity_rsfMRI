function PEB_subject_mean_paper_GSR(dataset,subject,name_ROI_def,ROI_list,ROI_list2,procedure,all_procedure_names,SPM_dir,Work_dir)

for N_prec_comp={'single'}    %for different number of different precision components using N_prec_comp={'single','fields','all'}

    disp(dataset);
    disp(num2str(subject));
    disp(procedure);
    
    %%%%%%%%%%%%%%%%%%%%%%
    %GSR-Basic comparison
    %%%%%%%%%%%%%%%%%%%%%%
    if (strcmp(procedure,'Basic')||strcmp(procedure,'GSR'))&&strcmp(all_procedure_names{2},'GSR') %only do this if GSR is in the loop
        
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
            cd([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_summary/DCM/'  procedure '/' name_ROI_def '/Full_model']);
            load([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_summary/DCM/'  procedure '/' name_ROI_def '/Full_model/GCM_' ntwrk_name{network_number} '_full_estim.mat']);
            load([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/QC/Above_treshold_marks_' ntwrk_name{network_number} '.mat']);
        
            flag=zeros(1,length(GCM));
            for diagn=1:length(GCM)
                if ~isnan(Posterior_estimates_var(1,1,diagn))||~isnan(Posterior_estimates_max(1,1,diagn))||~isnan(Posterior_estimates_par(1,1,diagn))||~isnan(Posterior_estimates_mot(1,1,diagn))||~isnan(Posterior_estimates_thr(1,1,diagn))
                    flag(diagn)=1;
                end
            end
            
            if strcmp(procedure,'Basic')
                procedure_comp='GSR';
                load([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/'  procedure_comp '/' name_ROI_def '/Full_model/QC/Above_treshold_marks_' ntwrk_name{network_number} '.mat']);
            end
            
            if strcmp(procedure,'GSR')
                procedure_comp='Basic';
                load([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/'  procedure_comp '/' name_ROI_def '/Full_model/QC/Above_treshold_marks_' ntwrk_name{network_number} '.mat']);
            end
            
            flag_comp=zeros(1,length(GCM));
            for diagn=1:length(GCM)
                if ~isnan(Posterior_estimates_var(1,1,diagn))||~isnan(Posterior_estimates_max(1,1,diagn))||~isnan(Posterior_estimates_par(1,1,diagn))||~isnan(Posterior_estimates_mot(1,1,diagn))||~isnan(Posterior_estimates_thr(1,1,diagn))
                    flag_comp(diagn)=1;
                end
            end
            
            Omitted_sessions=find(flag==1|flag_comp==1);
            save([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/QC/Omitted_sessions_' ntwrk_name{network_number} '_comp_GSR.mat'],'Omitted_sessions');
            
            GCM(find(flag==1|flag_comp==1))=[];
            
            if length(GCM)<8
                continue
            end
            
            %PEB settings
            M = struct();
            M.alpha = 1;
            M.beta  = 16;
            M.hE    = 0;
            M.hC    = 1/16;
            M.Q     = cell2mat(N_prec_comp);
            M.X=ones(length(GCM),1);
            
%             lastwarn('');
%             try
%                 
%                 %All parameters together
%                 disp('all param in same model');
%                 [PEB DCM]=spm_dcm_peb(GCM,M,'all');
%                 %test warning
%                 %             inv([0 0;0 0]);
%                 save([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/PEB_all_mean_comp_GSR_' cell2mat(N_prec_comp) '_' ntwrk_name{network_number} '.mat'],'PEB','DCM','GCM');
%             catch ME
%                 disp(['Dataset: ' dataset '; Subject: ' num2str(subject) '; Procedure: ' procedure]);
%                 disp('all');
%                 save([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/ERROR_PEB_all_mean_comp_GSR_' cell2mat(N_prec_comp) '_' ntwrk_name{network_number} '.mat'],'ME');
%             end
%             
%             if ~isempty(lastwarn)
%                 [msgstr,msgid]=lastwarn;
%                 save([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/WARNING_PEB_all_mean_comp_GSR_' cell2mat(N_prec_comp) '_' ntwrk_name{network_number} '.mat'],'msgstr','msgid');
%             end
            
%             clear PEB DCM;
            
            lastwarn('');
            try
                %Connectivity
                disp('Connectivity matrix only');
                [PEB DCM]=spm_dcm_peb(GCM,M,'A');
                
                save([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/PEB_A_mean_comp_GSR_' cell2mat(N_prec_comp) '_' ntwrk_name{network_number} '.mat'],'PEB','DCM','GCM');
            catch ME
                disp(['Dataset: ' dataset '; Subject: ' num2str(subject) '; Procedure: ' procedure]);
                disp('A');
                save([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/ERROR_PEB_A_mean_comp_GSR_' cell2mat(N_prec_comp) '_' ntwrk_name{network_number} '.mat'],'ME');
            end
            
            if ~isempty(lastwarn)
                [msgstr,msgid]=lastwarn;
                save([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/WARNING_PEB_A_mean_comp_GSR_' cell2mat(N_prec_comp) '_' ntwrk_name{network_number} '.mat'],'msgstr','msgid');
            end
            
            clear PEB DCM;
            
%             lastwarn('');
%             try
%                 %Hemodynamics
%                 disp('Hemodynamics only');
%                 [PEB DCM]=spm_dcm_peb(GCM,M,{'transit','decay','epsilon'});
%                 
%                 save([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/PEB_hemodyna_mean_comp_GSR_' cell2mat(N_prec_comp) '_' ntwrk_name{network_number} '.mat'],'PEB','DCM','GCM');
%             catch ME
%                 disp(['Dataset: ' dataset '; Subject: ' num2str(subject) '; Procedure: ' procedure]);
%                 disp('Hemodyna');
%                 save([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/ERROR_PEB_hemodyna_mean_comp_GSR_' cell2mat(N_prec_comp) '_' ntwrk_name{network_number} '.mat'],'ME');
%             end
%             
%             if ~isempty(lastwarn)
%                 [msgstr,msgid]=lastwarn;
%                 save([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/WARNING_PEB_hemodyna_mean_comp_GSR_' cell2mat(N_prec_comp) '_' ntwrk_name{network_number} '.mat'],'msgstr','msgid');
%             end
%             
%             clear PEB DCM;
            
            lastwarn('');
            try
                %Spectral noise
                disp('Spectral noise only');
                [PEB DCM]=spm_dcm_peb(GCM,M,{'a','b','c'});
                
                save([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/PEB_fluct_mean_comp_GSR_' cell2mat(N_prec_comp) '_' ntwrk_name{network_number} '.mat'],'PEB','DCM','GCM');
            catch ME
                disp(['Dataset: ' dataset '; Subject: ' num2str(subject) '; Procedure: ' procedure]);
                disp('fluct');
                save([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/ERROR_PEB_fluct_mean_comp_GSR_' cell2mat(N_prec_comp) '_' ntwrk_name{network_number} '.mat'],'ME');
            end
            
            if ~isempty(lastwarn)
                [msgstr,msgid]=lastwarn;
                save([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/WARNING_PEB_fluct_mean_comp_GSR_' cell2mat(N_prec_comp) '_' ntwrk_name{network_number} '.mat'],'msgstr','msgid');
            end
            
            clear PEB DCM;
            
            clear GCM Posterior_estimates_var Lower_bound_ci_var Higher_bound_ci_var Length_error_bar_var Posterior_estimates_max Lower_bound_ci_max Higher_bound_ci_max Length_error_bar_max Posterior_estimates_par Lower_bound_ci_par Higher_bound_ci_par Length_error_bar_par Posterior_estimates_mot Lower_bound_ci_mot Higher_bound_ci_mot Length_error_bar_mot Posterior_estimates_thr Lower_bound_ci_thr Higher_bound_ci_thr Length_error_bar_thr Posterior_estimates;
            
        end
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%
    %GSR-Basic comparison: GSR_regr
    %%%%%%%%%%%%%%%%%%%%%%
    if strcmp(procedure,'Basic')||strcmp(procedure,'GSR_regr')
        
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
            cd([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_summary/DCM/'  procedure '/' name_ROI_def '/Full_model']);
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
                load([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/'  procedure_comp '/' name_ROI_def '/Full_model/QC/Above_treshold_marks_' ntwrk_name{network_number} '.mat']);
            end
            
            if strcmp(procedure,'GSR_regr')
                procedure_comp='Basic';
                load([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/'  procedure_comp '/' name_ROI_def '/Full_model/QC/Above_treshold_marks_' ntwrk_name{network_number} '.mat']);
            end
            
            flag_comp=zeros(1,length(GCM));
            for diagn=1:length(GCM)
                if ~isnan(Posterior_estimates_var(1,1,diagn))||~isnan(Posterior_estimates_max(1,1,diagn))||~isnan(Posterior_estimates_par(1,1,diagn))||~isnan(Posterior_estimates_mot(1,1,diagn))||~isnan(Posterior_estimates_thr(1,1,diagn))
                    flag_comp(diagn)=1;
                end
            end
            
            Omitted_sessions=find(flag==1|flag_comp==1);
            save([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/QC/Omitted_sessions_' ntwrk_name{network_number} '_comp_GSR_regr.mat'],'Omitted_sessions');
            
            GCM(find(flag==1|flag_comp==1))=[];
            
            if length(GCM)<8
                continue
            end
            
            %PEB settings
            M = struct();
            M.alpha = 1;
            M.beta  = 16;
            M.hE    = 0;
            M.hC    = 1/16;
            M.Q     = cell2mat(N_prec_comp);
            M.X=ones(length(GCM),1);
            
%             lastwarn('');
%             try
%                 %All parameters together
%                 disp('all param in same model');
%                 [PEB DCM]=spm_dcm_peb(GCM,M,'all');
%                 
%                 save([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/PEB_all_mean_comp_GSR_regr_' cell2mat(N_prec_comp) '_' ntwrk_name{network_number} '.mat'],'PEB','DCM','GCM');
%             catch ME
%                 disp(['Dataset: ' dataset '; Subject: ' num2str(subject) '; Procedure: ' procedure]);
%                 disp('all');
%                 save([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/ERROR_PEB_all_mean_comp_GSR_regr_' cell2mat(N_prec_comp) '_' ntwrk_name{network_number} '.mat'],'ME');
%             end
%             
%             if ~isempty(lastwarn)
%                 [msgstr,msgid]=lastwarn;
%                 save([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/WARNING_PEB_all_mean_comp_GSR_regr_' cell2mat(N_prec_comp) '_' ntwrk_name{network_number} '.mat'],'msgstr','msgid');
%             end
%             
%             clear PEB DCM;
            
            lastwarn('');
            try
                %Connectivity
                disp('Connectivity matrix only');
                [PEB DCM]=spm_dcm_peb(GCM,M,'A');
                
                save([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/PEB_A_mean_comp_GSR_regr_' cell2mat(N_prec_comp) '_' ntwrk_name{network_number} '.mat'],'PEB','DCM','GCM');
            catch ME
                disp(['Dataset: ' dataset '; Subject: ' num2str(subject) '; Procedure: ' procedure]);
                disp('A');
                save([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/ERROR_PEB_A_mean_comp_GSR_regr_' cell2mat(N_prec_comp) '_' ntwrk_name{network_number} '.mat'],'ME');
            end
            
            if ~isempty(lastwarn)
                [msgstr,msgid]=lastwarn;
                save([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/WARNING_PEB_A_mean_comp_GSR_regr_' cell2mat(N_prec_comp) '_' ntwrk_name{network_number} '.mat'],'msgstr','msgid');
            end
            
            clear PEB DCM;
            
%             lastwarn('');
%             try
%                 %Hemodynamics
%                 disp('Hemodynamics only');
%                 [PEB DCM]=spm_dcm_peb(GCM,M,{'transit','decay','epsilon'});
%                 
%                 save([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/PEB_hemodyna_mean_comp_GSR_regr_' cell2mat(N_prec_comp) '_' ntwrk_name{network_number} '.mat'],'PEB','DCM','GCM');
%             catch ME
%                 disp(['Dataset: ' dataset '; Subject: ' num2str(subject) '; Procedure: ' procedure]);
%                 disp('Hemodyna');
%                 save([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/ERROR_PEB_hemodyna_mean_comp_GSR_regr_' cell2mat(N_prec_comp) '_' ntwrk_name{network_number} '.mat'],'ME');
%             end
%             
%             if ~isempty(lastwarn)
%                 [msgstr,msgid]=lastwarn;
%                 save([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/WARNING_PEB_hemodyna_mean_comp_GSR_regr_' cell2mat(N_prec_comp) '_' ntwrk_name{network_number} '.mat'],'msgstr','msgid');
%             end
%             
%             clear PEB DCM;
            
            lastwarn('');
            try
                %Spectral noise
                disp('Spectral noise only');
                [PEB DCM]=spm_dcm_peb(GCM,M,{'a','b','c'});
                
                save([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/PEB_fluct_mean_comp_GSR_regr_' cell2mat(N_prec_comp) '_' ntwrk_name{network_number} '.mat'],'PEB','DCM','GCM');
            catch ME
                disp(['Dataset: ' dataset '; Subject: ' num2str(subject) '; Procedure: ' procedure]);
                disp('fluct');
                save([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/ERROR_PEB_fluct_mean_comp_GSR_regr_' cell2mat(N_prec_comp) '_' ntwrk_name{network_number} '.mat'],'ME');
            end
            
            if ~isempty(lastwarn)
                [msgstr,msgid]=lastwarn;
                save([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/WARNING_PEB_fluct_mean_comp_GSR_regr_' cell2mat(N_prec_comp) '_' ntwrk_name{network_number} '.mat'],'msgstr','msgid');
            end
            
            clear PEB DCM;
            
            clear GCM Posterior_estimates_var Lower_bound_ci_var Higher_bound_ci_var Length_error_bar_var Posterior_estimates_max Lower_bound_ci_max Higher_bound_ci_max Length_error_bar_max Posterior_estimates_par Lower_bound_ci_par Higher_bound_ci_par Length_error_bar_par Posterior_estimates_mot Lower_bound_ci_mot Higher_bound_ci_mot Length_error_bar_mot Posterior_estimates_thr Lower_bound_ci_thr Higher_bound_ci_thr Length_error_bar_thr Posterior_estimates;
            
        end
        
    end
end

for N_prec_comp={'single'}  %for different number of different precision components using N_prec_comp={'single','fields'}
    if exist('ROI_list2')
        
        if (strcmp(procedure,'Basic')||strcmp(procedure,'GSR'))&&strcmp(all_procedure_names{2},'GSR') %only do this is GSR is in the loop
            
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
            
            for network_number=1:length(ntwrk_abbrev2)
                disp(ntwrk_abbrev2{network_number});
                cd([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_summary/DCM/'  procedure '/' name_ROI_def '/Full_model']);
                load([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_summary/DCM/'  procedure '/' name_ROI_def '/Full_model/GCM_' ntwrk_abbrev2{network_number} '_full_estim.mat']);
                load([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/QC/Above_treshold_marks_' ntwrk_abbrev2{network_number} '.mat']);
        
                flag=zeros(1,length(GCM));
                for diagn=1:length(GCM)
                    if ~isnan(Posterior_estimates_var(1,1,diagn))||~isnan(Posterior_estimates_max(1,1,diagn))||~isnan(Posterior_estimates_par(1,1,diagn))||~isnan(Posterior_estimates_mot(1,1,diagn))||~isnan(Posterior_estimates_thr(1,1,diagn))
                        flag(diagn)=1;
                    end
                end
                
                if strcmp(procedure,'Basic')
                    procedure_comp='GSR';
                    load([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/'  procedure_comp '/' name_ROI_def '/Full_model/QC/Above_treshold_marks_' ntwrk_abbrev2{network_number} '.mat']);
                end
                
                if strcmp(procedure,'GSR')
                    procedure_comp='Basic';
                    load([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/'  procedure_comp '/' name_ROI_def '/Full_model/QC/Above_treshold_marks_' ntwrk_abbrev2{network_number} '.mat']);
                end
                
                flag_comp=zeros(1,length(GCM));
                for diagn=1:length(GCM)
                    if ~isnan(Posterior_estimates_var(1,1,diagn))||~isnan(Posterior_estimates_max(1,1,diagn))||~isnan(Posterior_estimates_par(1,1,diagn))||~isnan(Posterior_estimates_mot(1,1,diagn))||~isnan(Posterior_estimates_thr(1,1,diagn))
                        flag_comp(diagn)=1;
                    end
                end
                
                Omitted_sessions=find(flag==1|flag_comp==1);
                save([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/QC/Omitted_sessions_' ntwrk_abbrev2{network_number} '_comp_GSR.mat'],'Omitted_sessions');
                
                GCM(find(flag==1|flag_comp==1))=[];
                
                if length(GCM)<8
                    continue
                end
                
                %PEB settings
                M = struct();
                M.alpha = 1;
                M.beta  = 16;
                M.hE    = 0;
                M.hC    = 1/16;
                M.Q     = cell2mat(N_prec_comp);
                M.X=ones(length(GCM),1);
                
%                 lastwarn('');
%                 try
%                     
%                     %All parameters together
%                     disp('all param in same model');
%                     [PEB DCM]=spm_dcm_peb(GCM,M,'all');
%                     %test warning
%                     %             inv([0 0;0 0]);
%                     save([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/PEB_all_mean_comp_GSR_' cell2mat(N_prec_comp) '_' ntwrk_abbrev2{network_number} '.mat'],'PEB','DCM','GCM');
%                 catch ME
%                     disp(['Dataset: ' dataset '; Subject: ' num2str(subject) '; Procedure: ' procedure]);
%                     disp('all');
%                     save([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/ERROR_PEB_all_mean_comp_GSR_' cell2mat(N_prec_comp) '_' ntwrk_abbrev2{network_number} '.mat'],'ME');
%                 end
%                 
%                 if ~isempty(lastwarn)
%                     [msgstr,msgid]=lastwarn;
%                     save([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/WARNING_PEB_all_mean_comp_GSR_' cell2mat(N_prec_comp) '_' ntwrk_abbrev2{network_number} '.mat'],'msgstr','msgid');
%                 end
%                 
%                 clear PEB DCM;
                
                lastwarn('');
                try
                    %Connectivity
                    disp('Connectivity matrix only');
                    [PEB DCM]=spm_dcm_peb(GCM,M,'A');
                    
                    save([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/PEB_A_mean_comp_GSR_' cell2mat(N_prec_comp) '_' ntwrk_abbrev2{network_number} '.mat'],'PEB','DCM','GCM');
                catch ME
                    disp(['Dataset: ' dataset '; Subject: ' num2str(subject) '; Procedure: ' procedure]);
                    disp('A');
                    save([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/ERROR_PEB_A_mean_comp_GSR_' cell2mat(N_prec_comp) '_' ntwrk_abbrev2{network_number} '.mat'],'ME');
                end
                
                if ~isempty(lastwarn)
                    [msgstr,msgid]=lastwarn;
                    save([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/WARNING_PEB_A_mean_comp_GSR_' cell2mat(N_prec_comp) '_' ntwrk_abbrev2{network_number} '.mat'],'msgstr','msgid');
                end
                
                clear PEB DCM;
                
%                 lastwarn('');
%                 try
%                     %Hemodynamics
%                     disp('Hemodynamics only');
%                     [PEB DCM]=spm_dcm_peb(GCM,M,{'transit','decay','epsilon'});
%                     
%                     save([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/PEB_hemodyna_mean_comp_GSR_' cell2mat(N_prec_comp) '_' ntwrk_abbrev2{network_number} '.mat'],'PEB','DCM','GCM');
%                 catch ME
%                     disp(['Dataset: ' dataset '; Subject: ' num2str(subject) '; Procedure: ' procedure]);
%                     disp('Hemodyna');
%                     save([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/ERROR_PEB_hemodyna_mean_comp_GSR_' cell2mat(N_prec_comp) '_' ntwrk_abbrev2{network_number} '.mat'],'ME');
%                 end
%                 
%                 if ~isempty(lastwarn)
%                     [msgstr,msgid]=lastwarn;
%                     save([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/WARNING_PEB_hemodyna_mean_comp_GSR_' cell2mat(N_prec_comp) '_' ntwrk_abbrev2{network_number} '.mat'],'msgstr','msgid');
%                 end
%                 
%                 clear PEB DCM;
                
                lastwarn('');
                try
                    %Spectral noise
                    disp('Spectral noise only');
                    [PEB DCM]=spm_dcm_peb(GCM,M,{'a','b','c'});
                    
                    save([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/PEB_fluct_mean_comp_GSR_' cell2mat(N_prec_comp) '_' ntwrk_abbrev2{network_number} '.mat'],'PEB','DCM','GCM');
                catch ME
                    disp(['Dataset: ' dataset '; Subject: ' num2str(subject) '; Procedure: ' procedure]);
                    disp('fluct');
                    save([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/ERROR_PEB_fluct_mean_comp_GSR_' cell2mat(N_prec_comp) '_' ntwrk_abbrev2{network_number} '.mat'],'ME');
                end
                
                if ~isempty(lastwarn)
                    [msgstr,msgid]=lastwarn;
                    save([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/WARNING_PEB_fluct_mean_comp_GSR_' cell2mat(N_prec_comp) '_' ntwrk_abbrev2{network_number} '.mat'],'msgstr','msgid');
                end
                
                clear PEB DCM;
                
                clear GCM Posterior_estimates_var Lower_bound_ci_var Higher_bound_ci_var Length_error_bar_var Posterior_estimates_max Lower_bound_ci_max Higher_bound_ci_max Length_error_bar_max Posterior_estimates_par Lower_bound_ci_par Higher_bound_ci_par Length_error_bar_par Posterior_estimates_mot Lower_bound_ci_mot Higher_bound_ci_mot Length_error_bar_mot Posterior_estimates_thr Lower_bound_ci_thr Higher_bound_ci_thr Length_error_bar_thr Posterior_estimates;
                
            end
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%
        %GSR-Basic comparison: GSR_regr
        %%%%%%%%%%%%%%%%%%%%%%
        if strcmp(procedure,'Basic')||strcmp(procedure,'GSR_regr')
            
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
            
            for network_number=1:length(ntwrk_abbrev2)
                disp(ntwrk_abbrev2{network_number});
                cd([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_summary/DCM/'  procedure '/' name_ROI_def '/Full_model']);
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
                    load([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/'  procedure_comp '/' name_ROI_def '/Full_model/QC/Above_treshold_marks_' ntwrk_abbrev2{network_number} '.mat']);
                end
                
                if strcmp(procedure,'GSR_regr')
                    procedure_comp='Basic';
                    load([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/'  procedure_comp '/' name_ROI_def '/Full_model/QC/Above_treshold_marks_' ntwrk_abbrev2{network_number} '.mat']);
                end
                
                flag_comp=zeros(1,length(GCM));
                for diagn=1:length(GCM)
                    if ~isnan(Posterior_estimates_var(1,1,diagn))||~isnan(Posterior_estimates_max(1,1,diagn))||~isnan(Posterior_estimates_par(1,1,diagn))||~isnan(Posterior_estimates_mot(1,1,diagn))||~isnan(Posterior_estimates_thr(1,1,diagn))
                        flag_comp(diagn)=1;
                    end
                end
                
                Omitted_sessions=find(flag==1|flag_comp==1);
                save([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/QC/Omitted_sessions_' ntwrk_abbrev2{network_number} '_comp_GSR_regr.mat'],'Omitted_sessions');
                
                GCM(find(flag==1|flag_comp==1))=[];
                
                if length(GCM)<8
                    continue
                end
                
                %PEB settings
                M = struct();
                M.alpha = 1;
                M.beta  = 16;
                M.hE    = 0;
                M.hC    = 1/16;
                M.Q     = cell2mat(N_prec_comp);
                M.X=ones(length(GCM),1);
                
%                 lastwarn('');
%                 try
%                     %All parameters together
%                     disp('all param in same model');
%                     [PEB DCM]=spm_dcm_peb(GCM,M,'all');
%                     
%                     save([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/PEB_all_mean_comp_GSR_regr_' cell2mat(N_prec_comp) '_' ntwrk_abbrev2{network_number} '.mat'],'PEB','DCM','GCM');
%                 catch ME
%                     disp(['Dataset: ' dataset '; Subject: ' num2str(subject) '; Procedure: ' procedure]);
%                     disp('all');
%                     save([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/ERROR_PEB_all_mean_comp_GSR_regr_' cell2mat(N_prec_comp) '_' ntwrk_abbrev2{network_number} '.mat'],'ME');
%                 end
%                 
%                 if ~isempty(lastwarn)
%                     [msgstr,msgid]=lastwarn;
%                     save([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/WARNING_PEB_all_mean_comp_GSR_regr_' cell2mat(N_prec_comp) '_' ntwrk_abbrev2{network_number} '.mat'],'msgstr','msgid');
%                 end
%                 
%                 clear PEB DCM;
                
                lastwarn('');
                try
                    %Connectivity
                    disp('Connectivity matrix only');
                    [PEB DCM]=spm_dcm_peb(GCM,M,'A');
                    
                    save([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/PEB_A_mean_comp_GSR_regr_' cell2mat(N_prec_comp) '_' ntwrk_abbrev2{network_number} '.mat'],'PEB','DCM','GCM');
                catch ME
                    disp(['Dataset: ' dataset '; Subject: ' num2str(subject) '; Procedure: ' procedure]);
                    disp('A');
                    save([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/ERROR_PEB_A_mean_comp_GSR_regr_' cell2mat(N_prec_comp) '_' ntwrk_abbrev2{network_number} '.mat'],'ME');
                end
                
                if ~isempty(lastwarn)
                    [msgstr,msgid]=lastwarn;
                    save([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/WARNING_PEB_A_mean_comp_GSR_regr_' cell2mat(N_prec_comp) '_' ntwrk_abbrev2{network_number} '.mat'],'msgstr','msgid');
                end
                
                clear PEB DCM;
                
%                 lastwarn('');
%                 try
%                     %Hemodynamics
%                     disp('Hemodynamics only');
%                     [PEB DCM]=spm_dcm_peb(GCM,M,{'transit','decay','epsilon'});
%                     
%                     save([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/PEB_hemodyna_mean_comp_GSR_regr_' cell2mat(N_prec_comp) '_' ntwrk_abbrev2{network_number} '.mat'],'PEB','DCM','GCM');
%                 catch ME
%                     disp(['Dataset: ' dataset '; Subject: ' num2str(subject) '; Procedure: ' procedure]);
%                     disp('Hemodyna');
%                     save([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/ERROR_PEB_hemodyna_mean_comp_GSR_regr_' cell2mat(N_prec_comp) '_' ntwrk_abbrev2{network_number} '.mat'],'ME');
%                 end
%                 
%                 if ~isempty(lastwarn)
%                     [msgstr,msgid]=lastwarn;
%                     save([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/WARNING_PEB_hemodyna_mean_comp_GSR_regr_' cell2mat(N_prec_comp) '_' ntwrk_abbrev2{network_number} '.mat'],'msgstr','msgid');
%                 end
%                 
%                 clear PEB DCM;
                
                lastwarn('');
                try
                    %Spectral noise
                    disp('Spectral noise only');
                    [PEB DCM]=spm_dcm_peb(GCM,M,{'a','b','c'});
                    
                    save([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/PEB_fluct_mean_comp_GSR_regr_' cell2mat(N_prec_comp) '_' ntwrk_abbrev2{network_number} '.mat'],'PEB','DCM','GCM');
                catch ME
                    disp(['Dataset: ' dataset '; Subject: ' num2str(subject) '; Procedure: ' procedure]);
                    disp('fluct');
                    save([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/ERROR_PEB_fluct_mean_comp_GSR_regr_' cell2mat(N_prec_comp) '_' ntwrk_abbrev2{network_number} '.mat'],'ME');
                end
                
                if ~isempty(lastwarn)
                    [msgstr,msgid]=lastwarn;
                    save([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/WARNING_PEB_fluct_mean_comp_GSR_regr_' cell2mat(N_prec_comp) '_' ntwrk_abbrev2{network_number} '.mat'],'msgstr','msgid');
                end
                
                clear PEB DCM;
                
                clear GCM Posterior_estimates_var Lower_bound_ci_var Higher_bound_ci_var Length_error_bar_var Posterior_estimates_max Lower_bound_ci_max Higher_bound_ci_max Length_error_bar_max Posterior_estimates_par Lower_bound_ci_par Higher_bound_ci_par Length_error_bar_par Posterior_estimates_mot Lower_bound_ci_mot Higher_bound_ci_mot Length_error_bar_mot Posterior_estimates_thr Lower_bound_ci_thr Higher_bound_ci_thr Length_error_bar_thr Posterior_estimates;
                
            end
            
        end
        
    end
end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Compute hemispheric assymetry: Connectivity-matrix only
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     if (strcmp(procedure,'Basic')||strcmp(procedure,'GSR'))&&strcmp(all_procedure_names{2},'GSR') %only do this is GSR is in the loop
%         
%         tmp=0;
%         
%         for VOI_number=1:size(ROI_list,1)
%             ntwrk=ROI_list{VOI_number,1}(1:3);
%             
%             if VOI_number>1 && strcmp(ROI_list{VOI_number,1}(1:3),ROI_list{VOI_number-1,1}(1:3))
%                 ntwrk_size(tmp)=ntwrk_size(tmp)+1;
%                 continue
%                 
%             else
%                 tmp=tmp+1;
%                 ntwrk_size(tmp)=1;
%                 ntwrk_name{tmp}=ROI_list{VOI_number,1}(1:3);
%             end
%         end
%         
%         for network_number=1:length(ntwrk_name)
%             disp(ntwrk_name{network_number});
%             try
%                 load([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/PEB_A_mean_comp_GSR_' cell2mat(N_prec_comp) '_' ntwrk_name{network_number} '.mat'],'PEB','DCM');
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
%                 %Somatomotor network
%                 if strcmp(ntwrk_name{network_number},'SMR')
%                     T = 0;
%                     C = [0 0 0 1 0 1 -1 -1 0 zeros(1,length(PEB.Ep)-9)]'/2;
%                     c = C'*spm_vec(PEB.Ep);
%                     v = C'*PEB.Cp*C;
%                     PP   = 1-spm_Ncdf(T,c,v);
%                 end
%                 
%                 %DMN
%                 if strcmp(ntwrk_name{network_number},'DMN')
%                     T = 0;
%                     C = [0 0 0 0 1 0 1 1 0 0 0 0 -1 -1 -1 0 zeros(1,length(PEB.Ep)-16)]'/3;
%                     c = C'*spm_vec(PEB.Ep);
%                     v = C'*PEB.Cp*C;
%                     PP   = 1-spm_Ncdf(T,c,v);
%                 end
%                 
%                 mean_diff=c;
%                 var_of_sum=v;
%                 posterior_probability=PP;
%                 
%                 save([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/Lateralization_index_A_comp_GSR_' cell2mat(N_prec_comp) '_' ntwrk_name{network_number} '.mat'],'mean_diff','var_of_sum','posterior_probability');
%                 clear mean_diff var_of_sum c v C PP T PEB DCM posterior_probability;
%             catch ME
%                 disp('Hemispheric Assymetry')
%                 disp(['Dataset: ' dataset '; Subject: ' num2str(subject) '; Procedure: ' procedure]);
%                 if strcmp(ME.identifier,'MATLAB:load:couldNotReadFile')
%                     continue
%                 else
%                     save([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/ERROR_Lateralization_index_A_comp_GSR_' cell2mat(N_prec_comp) '_' ntwrk_name{network_number} '.mat'],'ME');
%                 end
%             end
%         end
%     end
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %Compute hemispheric assymetry: All parameters in same model
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     if (strcmp(procedure,'Basic')||strcmp(procedure,'GSR'))&&strcmp(all_procedure_names{2},'GSR') %only do this is GSR is in the loop
%         
%         tmp=0;
%         
%         for VOI_number=1:size(ROI_list,1)
%             ntwrk=ROI_list{VOI_number,1}(1:3);
%             
%             if VOI_number>1 && strcmp(ROI_list{VOI_number,1}(1:3),ROI_list{VOI_number-1,1}(1:3))
%                 ntwrk_size(tmp)=ntwrk_size(tmp)+1;
%                 continue
%                 
%             else
%                 tmp=tmp+1;
%                 ntwrk_size(tmp)=1;
%                 ntwrk_name{tmp}=ROI_list{VOI_number,1}(1:3);
%             end
%         end
%         
%         for network_number=1:length(ntwrk_name)
%             disp(ntwrk_name{network_number});
%             try
%                 load([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/PEB_all_mean_comp_GSR_' cell2mat(N_prec_comp) '_' ntwrk_name{network_number} '.mat'],'PEB','DCM');
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
%                 %Somatomotor network
%                 if strcmp(ntwrk_name{network_number},'SMR')
%                     T = 0;
%                     C = [0 0 0 1 0 1 -1 -1 0 zeros(1,length(PEB.Ep)-9)]'/2;
%                     c = C'*spm_vec(PEB.Ep);
%                     v = C'*PEB.Cp*C;
%                     PP   = 1-spm_Ncdf(T,c,v);
%                 end
%                 
%                 %DMN
%                 if strcmp(ntwrk_name{network_number},'DMN')
%                     T = 0;
%                     C = [0 0 0 0 1 0 1 1 0 0 0 0 -1 -1 -1 0 zeros(1,length(PEB.Ep)-16)]'/3;
%                     c = C'*spm_vec(PEB.Ep);
%                     v = C'*PEB.Cp*C;
%                     PP   = 1-spm_Ncdf(T,c,v);
%                 end
%                 
%                 mean_diff=c;
%                 var_of_sum=v;
%                 posterior_probability=PP;
%                 
%                 save([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/Lateralization_index_all_comp_GSR_' cell2mat(N_prec_comp) '_' ntwrk_name{network_number} '.mat'],'mean_diff','var_of_sum','posterior_probability');
%                 clear mean_diff var_of_sum c v C PP T PEB DCM posterior_probability;
%             catch ME
%                 disp('Hemispheric Assymetry')
%                 disp(['Dataset: ' dataset '; Subject: ' num2str(subject) '; Procedure: ' procedure]);
%                 if strcmp(ME.identifier,'MATLAB:load:couldNotReadFile')
%                     continue
%                 else
%                     save([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/ERROR_Lateralization_index_all_comp_GSR_' cell2mat(N_prec_comp) '_' ntwrk_name{network_number} '.mat'],'ME');
%                 end
%             end
%         end
%     end
%     
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute hemispheric assymetry: Connectivity-matrix only (GSR_regr)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if strcmp(procedure,'Basic')||strcmp(procedure,'GSR_regr')
%     for N_prec_comp={'single'}              %N_prec_comp={'single','fields','all'}
%         tmp=0;
%         
%         for VOI_number=1:size(ROI_list,1)
%             ntwrk=ROI_list{VOI_number,1}(1:3);
%             
%             if VOI_number>1 && strcmp(ROI_list{VOI_number,1}(1:3),ROI_list{VOI_number-1,1}(1:3))
%                 ntwrk_size(tmp)=ntwrk_size(tmp)+1;
%                 continue
%                 
%             else
%                 tmp=tmp+1;
%                 ntwrk_size(tmp)=1;
%                 ntwrk_name{tmp}=ROI_list{VOI_number,1}(1:3);
%             end
%         end
%         
%         for network_number=1:length(ntwrk_name)
%             disp(ntwrk_name{network_number});
%             try
%                 load([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/PEB_A_mean_comp_GSR_regr_' cell2mat(N_prec_comp) '_' ntwrk_name{network_number} '.mat'],'PEB','DCM');
%                 
%                 Salience network
%                 if strcmp(ntwrk_name{network_number},'SAL')
%                     T = 0;
%                     C = [0 0 0 0 0 1 0 0 1 1 1 0 0 1 1 -1 -1 -1 0 0 -1 -1 -1 0 0 zeros(1,length(PEB.Ep)-25)]'/6;
%                     c = C'*spm_vec(PEB.Ep);
%                     v = C'*PEB.Cp*C;
%                     PP   = 1-spm_Ncdf(T,c,v);
%                 end
%                 
%                 Somatomotor network
%                 if strcmp(ntwrk_name{network_number},'SMR')
%                     T = 0;
%                     C = [0 0 0 1 0 1 -1 -1 0 zeros(1,length(PEB.Ep)-9)]'/2;
%                     c = C'*spm_vec(PEB.Ep);
%                     v = C'*PEB.Cp*C;
%                     PP   = 1-spm_Ncdf(T,c,v);
%                 end
%                 
%                 DMN
%                 if strcmp(ntwrk_name{network_number},'DMN')
%                     T = 0;
%                     C = [0 0 0 0 1 0 1 1 0 0 0 0 -1 -1 -1 0 zeros(1,length(PEB.Ep)-16)]'/3;
%                     c = C'*spm_vec(PEB.Ep);
%                     v = C'*PEB.Cp*C;
%                     PP   = 1-spm_Ncdf(T,c,v);
%                 end
%                 
%                 mean_diff=c;
%                 var_of_sum=v;
%                 posterior_probability=PP;
%                 
%                 save([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/Lateralization_index_A_comp_GSR_regr_' cell2mat(N_prec_comp) '_' ntwrk_name{network_number} '.mat'],'mean_diff','var_of_sum','posterior_probability');
%                 clear mean_diff var_of_sum c v C PP T PEB DCM posterior_probability;
%             catch ME
%                 disp('Hemispheric Assymetry')
%                 disp(['Dataset: ' dataset '; Subject: ' num2str(subject) '; Procedure: ' procedure]);
%                 if strcmp(ME.identifier,'MATLAB:load:couldNotReadFile')
%                     continue
%                 else
%                     save([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/ERROR_Lateralization_index_A_comp_GSR_regr_' cell2mat(N_prec_comp) '_' ntwrk_name{network_number} '.mat'],'ME');
%                 end
%             end
%         end
%     end
% end
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %Compute hemispheric assymetry: All parameters in same model
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     if strcmp(procedure,'Basic')||strcmp(procedure,'GSR_regr')
%         
%         tmp=0;
%         
%         for VOI_number=1:size(ROI_list,1)
%             ntwrk=ROI_list{VOI_number,1}(1:3);
%             
%             if VOI_number>1 && strcmp(ROI_list{VOI_number,1}(1:3),ROI_list{VOI_number-1,1}(1:3))
%                 ntwrk_size(tmp)=ntwrk_size(tmp)+1;
%                 continue
%                 
%             else
%                 tmp=tmp+1;
%                 ntwrk_size(tmp)=1;
%                 ntwrk_name{tmp}=ROI_list{VOI_number,1}(1:3);
%             end
%         end
%         
%         for network_number=1:length(ntwrk_name)
%             disp(ntwrk_name{network_number});
%             
%             try
%                 load([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/PEB_all_mean_comp_GSR_regr_' cell2mat(N_prec_comp) '_' ntwrk_name{network_number} '.mat'],'PEB','DCM');
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
%                 %Somatomotor network
%                 if strcmp(ntwrk_name{network_number},'SMR')
%                     T = 0;
%                     C = [0 0 0 1 0 1 -1 -1 0 zeros(1,length(PEB.Ep)-9)]'/2;
%                     c = C'*spm_vec(PEB.Ep);
%                     v = C'*PEB.Cp*C;
%                     PP   = 1-spm_Ncdf(T,c,v);
%                 end
%                 
%                 %DMN
%                 if strcmp(ntwrk_name{network_number},'DMN')
%                     T = 0;
%                     C = [0 0 0 0 1 0 1 1 0 0 0 0 -1 -1 -1 0 zeros(1,length(PEB.Ep)-16)]'/3;
%                     c = C'*spm_vec(PEB.Ep);
%                     v = C'*PEB.Cp*C;
%                     PP   = 1-spm_Ncdf(T,c,v);
%                 end
%                 
%                 mean_diff=c;
%                 var_of_sum=v;
%                 posterior_probability=PP;
%                 
%                 save([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/Lateralization_index_all_comp_GSR_regr_' cell2mat(N_prec_comp) '_' ntwrk_name{network_number} '.mat'],'mean_diff','var_of_sum','posterior_probability');
%                 clear mean_diff var_of_sum c v C PP T PEB DCM posterior_probability;
%             catch ME
%                 disp('Hemispheric Assymetry')
%                 disp(['Dataset: ' dataset '; Subject: ' num2str(subject) '; Procedure: ' procedure]);
%                 if strcmp(ME.identifier,'MATLAB:load:couldNotReadFile')
%                     continue
%                 else
%                     save([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/ERROR_Lateralization_index_all_comp_GSR_regr_' cell2mat(N_prec_comp) '_' ntwrk_name{network_number} '.mat'],'ME');
%                 end
%             end
%         end
%     end
end