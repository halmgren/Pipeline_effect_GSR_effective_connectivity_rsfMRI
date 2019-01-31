function Number_included_sessions_paper_GSR(SPM_dir,Work_dir)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Only written for GSR_regr, not for GSR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

all_ROI_defs={'Smith'}; %ROI definition based on Smith et al., 2009
all_procedure_names={'Basic','GSR_regr'};                    %all_procedure_names={'Basic','GSR','GSR_regr'};

%%%%%%%%%%%%%%%%%%%%%%
%GSR-Basic comparison: GSR_regr
%%%%%%%%%%%%%%%%%%%%%%

for number_ROI_def=1:length(all_ROI_defs)
    
    name_ROI_def=all_ROI_defs{number_ROI_def};
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Give regions name and coördinates
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [ROI_list]=Define_ROIs_paper_GSR(name_ROI_def);
    
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
            for network_number=1:length(ntwrk_name)
                
                N_incl_ses=0;
                N_excl_ses=0;
                
                subjects_included=[];
                Omitted_subjects=[];
                
                tmp=0;
                
                for number_dataset=1:4
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %Extract information about datasets
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    [dataset,number_subject,single_band,slice_time_seconds]=Dataset_info_paper_GSR(number_dataset);
                    %%%%%%%%%%%%%%%%%%%%%%
                    %GSR-Basic comparison
                    %%%%%%%%%%%%%%%%%%%%%%
                    for subject=1:number_subject
                        tmp=tmp+1;
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
                        
                        GCM(find(flag==1|flag_comp==1))=[];
                        
                        if length(GCM)<8
                            disp(subject);
                            Omitted_subjects = [Omitted_subjects tmp];
                            continue
                        end
                        
                        subjects_included = [subjects_included tmp];
                        
                        Used_sessions=length(GCM);
                        Excl_sessions=length(find(flag==1|flag_comp==1));
                        
                        N_incl_ses=N_incl_ses+Used_sessions;
                        N_excl_ses=N_excl_ses+Excl_sessions;
                        
                    end
                end
                direct=[Work_dir '/Results_paper_GSR/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_group/N_incl_sessions/'];
                mkdir(direct);
                cd(direct);
                save(['N_incl_ses_' ntwrk_name{network_number} '_GSR_regr.mat'],'N_incl_ses','N_excl_ses');
                save(['Omitted_subjects_A_GSR_regr_single_' ntwrk_name{network_number} '.mat'],'Omitted_subjects');
                save(['Subjects_included_A_GSR_regr_single_' ntwrk_name{network_number} '.mat'],'subjects_included');
                
            end
            
            
        end
    end
end

for number_ROI_def=1:length(all_ROI_defs)
    
    name_ROI_def=all_ROI_defs{number_ROI_def};
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Give regions name and coördinates
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %DMN, salience, and somatomotor network were studied for this paper
    [ROI_list]=Define_ROIs_paper_GSR(name_ROI_def);
    [ROI_list2]=Define_comb_ROIs_paper_GSR(name_ROI_def);
    
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
            if strcmp(procedure,'Basic')||strcmp(procedure,'GSR_regr')
                for network_number=1:length(ntwrk_abbrev2)
                    
                    N_incl_ses=0;
                    N_excl_ses=0;
                    
                    subjects_included=[];
                    Omitted_subjects=[];
                    
                    tmp=0;
                    for number_dataset=1:4
                        
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %Extract information about datasets
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        
                        [dataset,number_subject,single_band,slice_time_seconds]=Dataset_info_paper_GSR(number_dataset);
                        %%%%%%%%%%%%%%%%%%%%%%
                        %GSR-Basic comparison
                        %%%%%%%%%%%%%%%%%%%%%%
                        for subject=1:number_subject
                            tmp=tmp+1;
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
                            
                            GCM(find(flag==1|flag_comp==1))=[];
                            
                            if length(GCM)<8
                                disp(subject);
                                Omitted_subjects = [Omitted_subjects tmp];
                                continue
                            end
                            
                            subjects_included = [subjects_included tmp];
                            
                            Used_sessions=length(GCM);
                            Excl_sessions=length(find(flag==1|flag_comp==1));
                            
                            N_incl_ses=N_incl_ses+Used_sessions;
                            N_excl_ses=N_excl_ses+Excl_sessions;
                            
                        end
                    end
                    direct=[Work_dir '/Results_paper_GSR/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_group/N_incl_sessions/'];
                    mkdir(direct);
                    cd(direct);
                    save(['N_incl_ses_' ntwrk_abbrev2{network_number} '_GSR_regr.mat'],'N_incl_ses','N_excl_ses');
                    save(['Omitted_subjects_A_GSR_regr_single_' ntwrk_abbrev2{network_number} '.mat'],'Omitted_subjects');
                    save(['Subjects_included_A_GSR_regr_single_' ntwrk_abbrev2{network_number} '.mat'],'subjects_included');
                    
                end
                
                
            end
        end
    end
end

%%%%%%%%%%%%%%
%Dataset HCP
%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%
%GSR-Basic comparison: GSR_regr
%%%%%%%%%%%%%%%%%%%%%%

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
                direct=[Work_dir '/Results_paper_GSR/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_group/N_incl_sessions/'];
                mkdir(direct);
                cd(direct);
                N_excl_ses=length(Omitted_subjects);
                N_incl_ses=length(subjects_included);
                
                save(['Omitted_subjects_HCP_A_GSR_regr_single_' ntwrk_name{network_number} '.mat'],'Omitted_subjects');
                save(['Subjects_included_HCP_A_GSR_regr_single_' ntwrk_name{network_number} '.mat'],'subjects_included');
                save(['N_incl_ses_HCP_' ntwrk_name{network_number} '_GSR_regr.mat'],'N_incl_ses','N_excl_ses');
                clear Omitted_subjects subjects_included N_incl_ses N_excl_ses
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
                    direct=[Work_dir '/Results_paper_GSR/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_group/N_incl_sessions/'];
                    mkdir(direct);
                    cd(direct);
                    N_excl_ses=length(Omitted_subjects);
                    N_incl_ses=length(subjects_included);
                    
                    save(['Omitted_subjects_HCP_A_GSR_regr_single_' ntwrk_abbrev2{network_number} '.mat'],'Omitted_subjects');
                    save(['Subjects_included_HCP_A_GSR_regr_single_' ntwrk_abbrev2{network_number} '.mat'],'subjects_included');
                    
                    save(['N_incl_ses_HCP_' ntwrk_abbrev2{network_number} '_GSR_regr.mat'],'N_incl_ses','N_excl_ses');
                    clear Omitted_subjects subjects_included N_incl_ses N_excl_ses
                end
                
            end
        end
        
        
    end
    
end