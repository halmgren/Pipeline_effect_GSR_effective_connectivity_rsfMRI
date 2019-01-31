function Above_threshold_sessions_paper_GSR(dataset,number_subject,ROI_list,ROI_list2,name_ROI_def,threshold_expl_var,threshold_max_conn,threshold_n_par_est,threshold_FD,threshold_FD_HCP,threshold_threshold_VOIs,procedure,SPM_dir,Work_dir)

if ~strcmp(dataset,'DatasetHCP')
    for subject=1:number_subject
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Extract folder (session) names
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        cd([Work_dir '/' dataset '/']);
        session_number=dir(['sub-' sprintf('%03d',subject)]);
        session_number(1:2)=[];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Extract network name and size
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
        
        if strcmp(procedure,'Basic')||strcmp(procedure,'GSR')||strcmp(procedure,'GSR_regr')
            
            mkdir([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/' procedure '/' name_ROI_def '/Full_model/']);
            cd([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/' procedure '/' name_ROI_def '/Full_model/']);
        
            
            for network=1:length(ntwrk_name)
                
                [Positive_sign_prop,Negative_sign_prop,Ns_sign_prop,Number_HQ_sessions]=deal(NaN(ntwrk_size(network),ntwrk_size(network)));
                
                %%%%%%%%%%%%%%%%%%
                %Load diagnostics
                %%%%%%%%%%%%%%%%%%
                
                cd([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_summary/DCM/' procedure '/' name_ROI_def '/Full_model/QC']);
            
                load(['Diagnostics_' ntwrk_name{network} '.mat']);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %Check whether diagnostics reach criteria
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                clear 'Posterior_estimates_var' 'Lower_bound_ci_var' 'Higher_bound_ci_var' 'Length_error_bar_var' 'Posterior_estimates_max' 'Lower_bound_ci_max' 'Higher_bound_ci_max' 'Length_error_bar_max' 'Posterior_estimates_par' 'Lower_bound_ci_par' 'Higher_bound_ci_par' 'Length_error_bar_par' 'Posterior_estimates_mot' 'Lower_bound_ci_mot' 'Higher_bound_ci_mot' 'Length_error_bar_mot' 'Posterior_estimates_thr' 'Lower_bound_ci_thr' 'Higher_bound_ci_thr' 'Length_error_bar_thr' 'Posterior_estimates_zero_conf' 'Lower_bound_ci_zero_conf' 'Higher_bound_ci_zero_conf' 'Length_error_bar_zero_conf';
                
                [Posterior_estimates_var(:,:,:),Lower_bound_ci_var(:,:,:),Higher_bound_ci_var(:,:,:),Length_error_bar_var(:,:,:)]=deal(NaN(size(Posterior_estimates,1),size(Posterior_estimates,2),size(Posterior_estimates,3)));
                [Posterior_estimates_max(:,:,:),Lower_bound_ci_max(:,:,:),Higher_bound_ci_max(:,:,:),Length_error_bar_max(:,:,:)]=deal(NaN(size(Posterior_estimates,1),size(Posterior_estimates,2),size(Posterior_estimates,3)));
                [Posterior_estimates_par(:,:,:),Lower_bound_ci_par(:,:,:),Higher_bound_ci_par(:,:,:),Length_error_bar_par(:,:,:)]=deal(NaN(size(Posterior_estimates,1),size(Posterior_estimates,2),size(Posterior_estimates,3)));
                [Posterior_estimates_mot(:,:,:),Lower_bound_ci_mot(:,:,:),Higher_bound_ci_mot(:,:,:),Length_error_bar_mot(:,:,:)]=deal(NaN(size(Posterior_estimates,1),size(Posterior_estimates,2),size(Posterior_estimates,3)));
                [Posterior_estimates_thr(:,:,:),Lower_bound_ci_thr(:,:,:),Higher_bound_ci_thr(:,:,:),Length_error_bar_thr(:,:,:)]=deal(NaN(size(Posterior_estimates,1),size(Posterior_estimates,2),size(Posterior_estimates,3)));
                [Posterior_estimates_zero_conf(:,:,:),Lower_bound_ci_zero_conf(:,:,:),Higher_bound_ci_zero_conf(:,:,:),Length_error_bar_zero_conf(:,:,:)]=deal(NaN(size(Posterior_estimates,1),size(Posterior_estimates,2),size(Posterior_estimates,3)));
                
                %%%%%%%%%%%%%%%%%%
                %If session doesn't reach a criterium, its estimates are saved
                %in the corresponding file; initially done to display
                %longitudinal data
                %%%%%%%%%%%%%%%%%%
                for session=1:length(Expl_var)
                    if Expl_var(session)<threshold_expl_var
                        Posterior_estimates_var(:,:,session)=Posterior_estimates(:,:,session);
                        Lower_bound_ci_var(:,:,session)=Lower_bound_ci(:,:,session);
                        Higher_bound_ci_var(:,:,session)=Higher_bound_ci(:,:,session);
                        Length_error_bar_var(:,:,session)=Length_error_bar(:,:,session);
                    end
                    
                    if Max_conn(session)<threshold_max_conn
                        Posterior_estimates_max(:,:,session)=Posterior_estimates(:,:,session);
                        Lower_bound_ci_max(:,:,session)=Lower_bound_ci(:,:,session);
                        Higher_bound_ci_max(:,:,session)=Higher_bound_ci(:,:,session);
                        Length_error_bar_max(:,:,session)=Length_error_bar(:,:,session);
                    end
                    
                    if N_est_par(session)<threshold_n_par_est
                        Posterior_estimates_par(:,:,session)=Posterior_estimates(:,:,session);
                        Lower_bound_ci_par(:,:,session)=Lower_bound_ci(:,:,session);
                        Higher_bound_ci_par(:,:,session)=Higher_bound_ci(:,:,session);
                        Length_error_bar_par(:,:,session)=Length_error_bar(:,:,session);
                    end
                    
                    if max(FD(:,session))>threshold_FD
                        Posterior_estimates_mot(:,:,session)=Posterior_estimates(:,:,session);
                        Lower_bound_ci_mot(:,:,session)=Lower_bound_ci(:,:,session);
                        Higher_bound_ci_mot(:,:,session)=Higher_bound_ci(:,:,session);
                        Length_error_bar_mot(:,:,session)=Length_error_bar(:,:,session);
                    end
                    
                    if max(thresholds_VOIs(:,session))>threshold_threshold_VOIs
                        Posterior_estimates_thr(:,:,session)=Posterior_estimates(:,:,session);
                        Lower_bound_ci_thr(:,:,session)=Lower_bound_ci(:,:,session);
                        Higher_bound_ci_thr(:,:,session)=Higher_bound_ci(:,:,session);
                        Length_error_bar_thr(:,:,session)=Length_error_bar(:,:,session);
                    end
                    
                    %confidence interval includes zero (not used for
                    %any diagnostics)
                    
                    for connection1=1:ntwrk_size(network)
                        for connection2=1:ntwrk_size(network)
                            if connection1==connection2
                                continue
                            elseif sign(Lower_bound_ci(connection1,connection2,session))~=sign(Higher_bound_ci(connection1,connection2,session))
                                Posterior_estimates_zero_conf(connection1,connection2,session)=Posterior_estimates(connection1,connection2,session);
                                Lower_bound_ci_zero_conf(connection1,connection2,session)=Lower_bound_ci(connection1,connection2,session);
                                Higher_bound_ci_zero_conf(connection1,connection2,session)=Higher_bound_ci(connection1,connection2,session);
                                Length_error_bar_zero_conf(connection1,connection2,session)=Length_error_bar(connection1,connection2,session);
                            end
                        end
                    end
                    
                    if ~isnan(Posterior_estimates_var(1,1,session))||~isnan(Posterior_estimates_max(1,1,session))||~isnan(Posterior_estimates_par(1,1,session))||~isnan(Posterior_estimates_mot(1,1,session))||~isnan(Posterior_estimates_thr(1,1,session))
                        Posterior_estimates(:,:,session)=NaN;
                    end
                    
                    for connection1=1:ntwrk_size(network)
                        for connection2=1:ntwrk_size(network)
                            if ~isnan(Posterior_estimates_zero_conf(connection1,connection2,session))
                                Posterior_estimates(connection1,connection2,session)=NaN;
                            end
                        end
                    end
                    
                end
                
                mkdir([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/' procedure '/' name_ROI_def '/Full_model/QC']);
                cd([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/' procedure '/' name_ROI_def '/Full_model/QC']);
             
                save(['Above_treshold_marks_' ntwrk_name{network} '.mat'],'Posterior_estimates_var','Lower_bound_ci_var','Higher_bound_ci_var','Length_error_bar_var','Posterior_estimates_max','Lower_bound_ci_max','Higher_bound_ci_max','Length_error_bar_max','Posterior_estimates_par','Lower_bound_ci_par','Higher_bound_ci_par','Length_error_bar_par','Posterior_estimates_mot','Lower_bound_ci_mot','Higher_bound_ci_mot','Length_error_bar_mot','Posterior_estimates_thr','Lower_bound_ci_thr','Higher_bound_ci_thr','Length_error_bar_thr','Posterior_estimates');
                
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
            
            tmp=0;
            
            for VOI_number=1:size(ROI_list2,1)
                ntwrk=ROI_list2{VOI_number,1}(1:3);
                
                if VOI_number>1 && strcmp(ROI_list2{VOI_number,1}(1:3),ROI_list2{VOI_number-1,1}(1:3))
                    ntwrk_size2(tmp)=ntwrk_size2(tmp)+1;
                    continue
                    
                else
                    tmp=tmp+1;
                    ntwrk_size2(tmp)=1;
                    %                 ntwrk_name{tmp}=ROI_list2{VOI_number,1}(1:3);
                end
            end
            
            if strcmp(procedure,'Basic')||strcmp(procedure,'GSR')||strcmp(procedure,'GSR_regr')
                
                mkdir([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/' procedure '/' name_ROI_def '/Full_model/']);
                cd([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/' procedure '/' name_ROI_def '/Full_model/']);
        
                for network=1:length(ntwrk_abbrev2)
                    
                    [Positive_sign_prop,Negative_sign_prop,Ns_sign_prop,Number_HQ_sessions]=deal(NaN(ntwrk_size2(network),ntwrk_size2(network)));
                    
                    %%%%%%%%%%%%%%%%%%
                    %Load diagnostics
                    %%%%%%%%%%%%%%%%%%
                    
                    cd([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_summary/DCM/' procedure '/' name_ROI_def '/Full_model/QC']);
            
                    load(['Diagnostics_' ntwrk_abbrev2{network} '.mat']);
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %Check whether diagnostics reach criteria
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    clear 'Posterior_estimates_var' 'Lower_bound_ci_var' 'Higher_bound_ci_var' 'Length_error_bar_var' 'Posterior_estimates_max' 'Lower_bound_ci_max' 'Higher_bound_ci_max' 'Length_error_bar_max' 'Posterior_estimates_par' 'Lower_bound_ci_par' 'Higher_bound_ci_par' 'Length_error_bar_par' 'Posterior_estimates_mot' 'Lower_bound_ci_mot' 'Higher_bound_ci_mot' 'Length_error_bar_mot' 'Posterior_estimates_thr' 'Lower_bound_ci_thr' 'Higher_bound_ci_thr' 'Length_error_bar_thr' 'Posterior_estimates_zero_conf' 'Lower_bound_ci_zero_conf' 'Higher_bound_ci_zero_conf' 'Length_error_bar_zero_conf';
                    
                    [Posterior_estimates_var(:,:,:),Lower_bound_ci_var(:,:,:),Higher_bound_ci_var(:,:,:),Length_error_bar_var(:,:,:)]=deal(NaN(size(Posterior_estimates,1),size(Posterior_estimates,2),size(Posterior_estimates,3)));
                    [Posterior_estimates_max(:,:,:),Lower_bound_ci_max(:,:,:),Higher_bound_ci_max(:,:,:),Length_error_bar_max(:,:,:)]=deal(NaN(size(Posterior_estimates,1),size(Posterior_estimates,2),size(Posterior_estimates,3)));
                    [Posterior_estimates_par(:,:,:),Lower_bound_ci_par(:,:,:),Higher_bound_ci_par(:,:,:),Length_error_bar_par(:,:,:)]=deal(NaN(size(Posterior_estimates,1),size(Posterior_estimates,2),size(Posterior_estimates,3)));
                    [Posterior_estimates_mot(:,:,:),Lower_bound_ci_mot(:,:,:),Higher_bound_ci_mot(:,:,:),Length_error_bar_mot(:,:,:)]=deal(NaN(size(Posterior_estimates,1),size(Posterior_estimates,2),size(Posterior_estimates,3)));
                    [Posterior_estimates_thr(:,:,:),Lower_bound_ci_thr(:,:,:),Higher_bound_ci_thr(:,:,:),Length_error_bar_thr(:,:,:)]=deal(NaN(size(Posterior_estimates,1),size(Posterior_estimates,2),size(Posterior_estimates,3)));
                    [Posterior_estimates_zero_conf(:,:,:),Lower_bound_ci_zero_conf(:,:,:),Higher_bound_ci_zero_conf(:,:,:),Length_error_bar_zero_conf(:,:,:)]=deal(NaN(size(Posterior_estimates,1),size(Posterior_estimates,2),size(Posterior_estimates,3)));
                    
                    %%%%%%%%%%%%%%%%%%
                    %If session doesn't reach a criterium, its estimates are saved
                    %in the corresponding file; initially done to display
                    %longitudinal data
                    %%%%%%%%%%%%%%%%%%
                    for session=1:length(Expl_var)
                        if Expl_var(session)<threshold_expl_var
                            Posterior_estimates_var(:,:,session)=Posterior_estimates(:,:,session);
                            Lower_bound_ci_var(:,:,session)=Lower_bound_ci(:,:,session);
                            Higher_bound_ci_var(:,:,session)=Higher_bound_ci(:,:,session);
                            Length_error_bar_var(:,:,session)=Length_error_bar(:,:,session);
                        end
                        
                        if Max_conn(session)<threshold_max_conn
                            Posterior_estimates_max(:,:,session)=Posterior_estimates(:,:,session);
                            Lower_bound_ci_max(:,:,session)=Lower_bound_ci(:,:,session);
                            Higher_bound_ci_max(:,:,session)=Higher_bound_ci(:,:,session);
                            Length_error_bar_max(:,:,session)=Length_error_bar(:,:,session);
                        end
                        
                        if N_est_par(session)<threshold_n_par_est
                            Posterior_estimates_par(:,:,session)=Posterior_estimates(:,:,session);
                            Lower_bound_ci_par(:,:,session)=Lower_bound_ci(:,:,session);
                            Higher_bound_ci_par(:,:,session)=Higher_bound_ci(:,:,session);
                            Length_error_bar_par(:,:,session)=Length_error_bar(:,:,session);
                        end
                        
                        if max(FD(:,session))>threshold_FD
                            Posterior_estimates_mot(:,:,session)=Posterior_estimates(:,:,session);
                            Lower_bound_ci_mot(:,:,session)=Lower_bound_ci(:,:,session);
                            Higher_bound_ci_mot(:,:,session)=Higher_bound_ci(:,:,session);
                            Length_error_bar_mot(:,:,session)=Length_error_bar(:,:,session);
                        end
                        
                        if max(thresholds_VOIs(:,session))>threshold_threshold_VOIs
                            Posterior_estimates_thr(:,:,session)=Posterior_estimates(:,:,session);
                            Lower_bound_ci_thr(:,:,session)=Lower_bound_ci(:,:,session);
                            Higher_bound_ci_thr(:,:,session)=Higher_bound_ci(:,:,session);
                            Length_error_bar_thr(:,:,session)=Length_error_bar(:,:,session);
                        end
                        
                        %confidence interval includes zero (not used for
                        %any diagnostics)
                        
                        for connection1=1:ntwrk_size2(network)
                            for connection2=1:ntwrk_size2(network)
                                if connection1==connection2
                                    continue
                                elseif sign(Lower_bound_ci(connection1,connection2,session))~=sign(Higher_bound_ci(connection1,connection2,session))
                                    Posterior_estimates_zero_conf(connection1,connection2,session)=Posterior_estimates(connection1,connection2,session);
                                    Lower_bound_ci_zero_conf(connection1,connection2,session)=Lower_bound_ci(connection1,connection2,session);
                                    Higher_bound_ci_zero_conf(connection1,connection2,session)=Higher_bound_ci(connection1,connection2,session);
                                    Length_error_bar_zero_conf(connection1,connection2,session)=Length_error_bar(connection1,connection2,session);
                                end
                            end
                        end
                        
                        if ~isnan(Posterior_estimates_var(1,1,session))||~isnan(Posterior_estimates_max(1,1,session))||~isnan(Posterior_estimates_par(1,1,session))||~isnan(Posterior_estimates_mot(1,1,session))||~isnan(Posterior_estimates_thr(1,1,session))
                            Posterior_estimates(:,:,session)=NaN;
                        end
                        
                        for connection1=1:ntwrk_size2(network)
                            for connection2=1:ntwrk_size2(network)
                                if ~isnan(Posterior_estimates_zero_conf(connection1,connection2,session))
                                    Posterior_estimates(connection1,connection2,session)=NaN;
                                end
                            end
                        end
                        
                    end
                    
                    mkdir([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/' procedure '/' name_ROI_def '/Full_model/QC']);
                    cd([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/' procedure '/' name_ROI_def '/Full_model/QC']);
             
                    save(['Above_treshold_marks_' ntwrk_abbrev2{network} '.mat'],'Posterior_estimates_var','Lower_bound_ci_var','Higher_bound_ci_var','Length_error_bar_var','Posterior_estimates_max','Lower_bound_ci_max','Higher_bound_ci_max','Length_error_bar_max','Posterior_estimates_par','Lower_bound_ci_par','Higher_bound_ci_par','Length_error_bar_par','Posterior_estimates_mot','Lower_bound_ci_mot','Higher_bound_ci_mot','Length_error_bar_mot','Posterior_estimates_thr','Lower_bound_ci_thr','Higher_bound_ci_thr','Length_error_bar_thr','Posterior_estimates');
                    
                end
            end
        end
    end
end

if strcmp(dataset,'DatasetHCP')
    for subject=1:number_subject
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Extract folder (session) names
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        cd([Work_dir '/' dataset '/']);
        session_number=dir(['sub-' sprintf('%03d',subject)]);
        session_number(1:2)=[];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Extract network name and size
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
        
        if strcmp(procedure,'Basic')||strcmp(procedure,'GSR')||strcmp(procedure,'GSR_regr')
            
            mkdir([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/' procedure '/' name_ROI_def '/Full_model/']);
        	cd([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/' procedure '/' name_ROI_def '/Full_model/']);
        
            for network=1:length(ntwrk_name)
                
                [Positive_sign_prop,Negative_sign_prop,Ns_sign_prop,Number_HQ_sessions]=deal(NaN(ntwrk_size(network),ntwrk_size(network)));
                
                %%%%%%%%%%%%%%%%%%
                %Load diagnostics
                %%%%%%%%%%%%%%%%%%
                
                cd([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_summary/DCM/' procedure '/' name_ROI_def '/Full_model/QC']);
            
                load(['Diagnostics_' ntwrk_name{network} '.mat']);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %Check whether diagnostics reach criteria
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                clear 'Posterior_estimates_var' 'Lower_bound_ci_var' 'Higher_bound_ci_var' 'Length_error_bar_var' 'Posterior_estimates_max' 'Lower_bound_ci_max' 'Higher_bound_ci_max' 'Length_error_bar_max' 'Posterior_estimates_par' 'Lower_bound_ci_par' 'Higher_bound_ci_par' 'Length_error_bar_par' 'Posterior_estimates_mot' 'Lower_bound_ci_mot' 'Higher_bound_ci_mot' 'Length_error_bar_mot' 'Posterior_estimates_thr' 'Lower_bound_ci_thr' 'Higher_bound_ci_thr' 'Length_error_bar_thr' 'Posterior_estimates_zero_conf' 'Lower_bound_ci_zero_conf' 'Higher_bound_ci_zero_conf' 'Length_error_bar_zero_conf';
                
                [Posterior_estimates_var(:,:,:),Lower_bound_ci_var(:,:,:),Higher_bound_ci_var(:,:,:),Length_error_bar_var(:,:,:)]=deal(NaN(size(Posterior_estimates,1),size(Posterior_estimates,2),size(Posterior_estimates,3)));
                [Posterior_estimates_max(:,:,:),Lower_bound_ci_max(:,:,:),Higher_bound_ci_max(:,:,:),Length_error_bar_max(:,:,:)]=deal(NaN(size(Posterior_estimates,1),size(Posterior_estimates,2),size(Posterior_estimates,3)));
                [Posterior_estimates_par(:,:,:),Lower_bound_ci_par(:,:,:),Higher_bound_ci_par(:,:,:),Length_error_bar_par(:,:,:)]=deal(NaN(size(Posterior_estimates,1),size(Posterior_estimates,2),size(Posterior_estimates,3)));
                [Posterior_estimates_mot(:,:,:),Lower_bound_ci_mot(:,:,:),Higher_bound_ci_mot(:,:,:),Length_error_bar_mot(:,:,:)]=deal(NaN(size(Posterior_estimates,1),size(Posterior_estimates,2),size(Posterior_estimates,3)));
                [Posterior_estimates_thr(:,:,:),Lower_bound_ci_thr(:,:,:),Higher_bound_ci_thr(:,:,:),Length_error_bar_thr(:,:,:)]=deal(NaN(size(Posterior_estimates,1),size(Posterior_estimates,2),size(Posterior_estimates,3)));
                [Posterior_estimates_zero_conf(:,:,:),Lower_bound_ci_zero_conf(:,:,:),Higher_bound_ci_zero_conf(:,:,:),Length_error_bar_zero_conf(:,:,:)]=deal(NaN(size(Posterior_estimates,1),size(Posterior_estimates,2),size(Posterior_estimates,3)));
                
                %%%%%%%%%%%%%%%%%%
                %If session doesn't reach a criterium, its estimates are saved
                %in the corresponding file; initially done to display
                %longitudinal data
                %%%%%%%%%%%%%%%%%%
                for session=1:length(Expl_var)
                    if Expl_var(session)<threshold_expl_var
                        Posterior_estimates_var(:,:,session)=Posterior_estimates(:,:,session);
                        Lower_bound_ci_var(:,:,session)=Lower_bound_ci(:,:,session);
                        Higher_bound_ci_var(:,:,session)=Higher_bound_ci(:,:,session);
                        Length_error_bar_var(:,:,session)=Length_error_bar(:,:,session);
                    end
                    
                    if Max_conn(session)<threshold_max_conn
                        Posterior_estimates_max(:,:,session)=Posterior_estimates(:,:,session);
                        Lower_bound_ci_max(:,:,session)=Lower_bound_ci(:,:,session);
                        Higher_bound_ci_max(:,:,session)=Higher_bound_ci(:,:,session);
                        Length_error_bar_max(:,:,session)=Length_error_bar(:,:,session);
                    end
                    
                    if N_est_par(session)<threshold_n_par_est
                        Posterior_estimates_par(:,:,session)=Posterior_estimates(:,:,session);
                        Lower_bound_ci_par(:,:,session)=Lower_bound_ci(:,:,session);
                        Higher_bound_ci_par(:,:,session)=Higher_bound_ci(:,:,session);
                        Length_error_bar_par(:,:,session)=Length_error_bar(:,:,session);
                    end
                    
                    if max(FD(:,session))>threshold_FD_HCP
                        Posterior_estimates_mot(:,:,session)=Posterior_estimates(:,:,session);
                        Lower_bound_ci_mot(:,:,session)=Lower_bound_ci(:,:,session);
                        Higher_bound_ci_mot(:,:,session)=Higher_bound_ci(:,:,session);
                        Length_error_bar_mot(:,:,session)=Length_error_bar(:,:,session);
                    end
                    
                    if max(thresholds_VOIs(:,session))>threshold_threshold_VOIs
                        Posterior_estimates_thr(:,:,session)=Posterior_estimates(:,:,session);
                        Lower_bound_ci_thr(:,:,session)=Lower_bound_ci(:,:,session);
                        Higher_bound_ci_thr(:,:,session)=Higher_bound_ci(:,:,session);
                        Length_error_bar_thr(:,:,session)=Length_error_bar(:,:,session);
                    end
                    
                    %confidence interval includes zero (not used for
                    %any diagnostics)
                    
                    for connection1=1:ntwrk_size(network)
                        for connection2=1:ntwrk_size(network)
                            if connection1==connection2
                                continue
                            elseif sign(Lower_bound_ci(connection1,connection2,session))~=sign(Higher_bound_ci(connection1,connection2,session))
                                Posterior_estimates_zero_conf(connection1,connection2,session)=Posterior_estimates(connection1,connection2,session);
                                Lower_bound_ci_zero_conf(connection1,connection2,session)=Lower_bound_ci(connection1,connection2,session);
                                Higher_bound_ci_zero_conf(connection1,connection2,session)=Higher_bound_ci(connection1,connection2,session);
                                Length_error_bar_zero_conf(connection1,connection2,session)=Length_error_bar(connection1,connection2,session);
                            end
                        end
                    end
                    
                    if ~isnan(Posterior_estimates_var(1,1,session))||~isnan(Posterior_estimates_max(1,1,session))||~isnan(Posterior_estimates_par(1,1,session))||~isnan(Posterior_estimates_mot(1,1,session))||~isnan(Posterior_estimates_thr(1,1,session))
                        Posterior_estimates(:,:,session)=NaN;
                    end
                    
                    for connection1=1:ntwrk_size(network)
                        for connection2=1:ntwrk_size(network)
                            if ~isnan(Posterior_estimates_zero_conf(connection1,connection2,session))
                                Posterior_estimates(connection1,connection2,session)=NaN;
                            end
                        end
                    end
                    
                end
                
                mkdir([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/' procedure '/' name_ROI_def '/Full_model/QC']);
                cd([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/' procedure '/' name_ROI_def '/Full_model/QC']);
            
                save(['Above_treshold_marks_' ntwrk_name{network} '.mat'],'Posterior_estimates_var','Lower_bound_ci_var','Higher_bound_ci_var','Length_error_bar_var','Posterior_estimates_max','Lower_bound_ci_max','Higher_bound_ci_max','Length_error_bar_max','Posterior_estimates_par','Lower_bound_ci_par','Higher_bound_ci_par','Length_error_bar_par','Posterior_estimates_mot','Lower_bound_ci_mot','Higher_bound_ci_mot','Length_error_bar_mot','Posterior_estimates_thr','Lower_bound_ci_thr','Higher_bound_ci_thr','Length_error_bar_thr','Posterior_estimates');
                
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
            
            tmp=0;
            
            for VOI_number=1:size(ROI_list2,1)
                ntwrk=ROI_list2{VOI_number,1}(1:3);
                
                if VOI_number>1 && strcmp(ROI_list2{VOI_number,1}(1:3),ROI_list2{VOI_number-1,1}(1:3))
                    ntwrk_size2(tmp)=ntwrk_size2(tmp)+1;
                    continue
                    
                else
                    tmp=tmp+1;
                    ntwrk_size2(tmp)=1;
                    %                 ntwrk_name{tmp}=ROI_list2{VOI_number,1}(1:3);
                end
            end
            
            if strcmp(procedure,'Basic')||strcmp(procedure,'GSR')||strcmp(procedure,'GSR_regr')
                
                mkdir([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/' procedure '/' name_ROI_def '/Full_model/']);
                cd([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/' procedure '/' name_ROI_def '/Full_model/']);
        
                for network=1:length(ntwrk_abbrev2)
                    
                    [Positive_sign_prop,Negative_sign_prop,Ns_sign_prop,Number_HQ_sessions]=deal(NaN(ntwrk_size2(network),ntwrk_size2(network)));
                    
                    %%%%%%%%%%%%%%%%%%
                    %Load diagnostics
                    %%%%%%%%%%%%%%%%%%
                    
                    cd([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_summary/DCM/' procedure '/' name_ROI_def '/Full_model/QC']);
            
                    load(['Diagnostics_' ntwrk_abbrev2{network} '.mat']);
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %Check whether diagnostics reach criteria
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    clear 'Posterior_estimates_var' 'Lower_bound_ci_var' 'Higher_bound_ci_var' 'Length_error_bar_var' 'Posterior_estimates_max' 'Lower_bound_ci_max' 'Higher_bound_ci_max' 'Length_error_bar_max' 'Posterior_estimates_par' 'Lower_bound_ci_par' 'Higher_bound_ci_par' 'Length_error_bar_par' 'Posterior_estimates_mot' 'Lower_bound_ci_mot' 'Higher_bound_ci_mot' 'Length_error_bar_mot' 'Posterior_estimates_thr' 'Lower_bound_ci_thr' 'Higher_bound_ci_thr' 'Length_error_bar_thr' 'Posterior_estimates_zero_conf' 'Lower_bound_ci_zero_conf' 'Higher_bound_ci_zero_conf' 'Length_error_bar_zero_conf';
                    
                    [Posterior_estimates_var(:,:,:),Lower_bound_ci_var(:,:,:),Higher_bound_ci_var(:,:,:),Length_error_bar_var(:,:,:)]=deal(NaN(size(Posterior_estimates,1),size(Posterior_estimates,2),size(Posterior_estimates,3)));
                    [Posterior_estimates_max(:,:,:),Lower_bound_ci_max(:,:,:),Higher_bound_ci_max(:,:,:),Length_error_bar_max(:,:,:)]=deal(NaN(size(Posterior_estimates,1),size(Posterior_estimates,2),size(Posterior_estimates,3)));
                    [Posterior_estimates_par(:,:,:),Lower_bound_ci_par(:,:,:),Higher_bound_ci_par(:,:,:),Length_error_bar_par(:,:,:)]=deal(NaN(size(Posterior_estimates,1),size(Posterior_estimates,2),size(Posterior_estimates,3)));
                    [Posterior_estimates_mot(:,:,:),Lower_bound_ci_mot(:,:,:),Higher_bound_ci_mot(:,:,:),Length_error_bar_mot(:,:,:)]=deal(NaN(size(Posterior_estimates,1),size(Posterior_estimates,2),size(Posterior_estimates,3)));
                    [Posterior_estimates_thr(:,:,:),Lower_bound_ci_thr(:,:,:),Higher_bound_ci_thr(:,:,:),Length_error_bar_thr(:,:,:)]=deal(NaN(size(Posterior_estimates,1),size(Posterior_estimates,2),size(Posterior_estimates,3)));
                    [Posterior_estimates_zero_conf(:,:,:),Lower_bound_ci_zero_conf(:,:,:),Higher_bound_ci_zero_conf(:,:,:),Length_error_bar_zero_conf(:,:,:)]=deal(NaN(size(Posterior_estimates,1),size(Posterior_estimates,2),size(Posterior_estimates,3)));
                    
                    %%%%%%%%%%%%%%%%%%
                    %If session doesn't reach a criterium, its estimates are saved
                    %in the corresponding file; initially done to display
                    %longitudinal data
                    %%%%%%%%%%%%%%%%%%
                    for session=1:length(Expl_var)
                        if Expl_var(session)<threshold_expl_var
                            Posterior_estimates_var(:,:,session)=Posterior_estimates(:,:,session);
                            Lower_bound_ci_var(:,:,session)=Lower_bound_ci(:,:,session);
                            Higher_bound_ci_var(:,:,session)=Higher_bound_ci(:,:,session);
                            Length_error_bar_var(:,:,session)=Length_error_bar(:,:,session);
                        end
                        
                        if Max_conn(session)<threshold_max_conn
                            Posterior_estimates_max(:,:,session)=Posterior_estimates(:,:,session);
                            Lower_bound_ci_max(:,:,session)=Lower_bound_ci(:,:,session);
                            Higher_bound_ci_max(:,:,session)=Higher_bound_ci(:,:,session);
                            Length_error_bar_max(:,:,session)=Length_error_bar(:,:,session);
                        end
                        
                        if N_est_par(session)<threshold_n_par_est
                            Posterior_estimates_par(:,:,session)=Posterior_estimates(:,:,session);
                            Lower_bound_ci_par(:,:,session)=Lower_bound_ci(:,:,session);
                            Higher_bound_ci_par(:,:,session)=Higher_bound_ci(:,:,session);
                            Length_error_bar_par(:,:,session)=Length_error_bar(:,:,session);
                        end
                        
                        if max(FD(:,session))>threshold_FD_HCP
                            Posterior_estimates_mot(:,:,session)=Posterior_estimates(:,:,session);
                            Lower_bound_ci_mot(:,:,session)=Lower_bound_ci(:,:,session);
                            Higher_bound_ci_mot(:,:,session)=Higher_bound_ci(:,:,session);
                            Length_error_bar_mot(:,:,session)=Length_error_bar(:,:,session);
                        end
                        
                        if max(thresholds_VOIs(:,session))>threshold_threshold_VOIs
                            Posterior_estimates_thr(:,:,session)=Posterior_estimates(:,:,session);
                            Lower_bound_ci_thr(:,:,session)=Lower_bound_ci(:,:,session);
                            Higher_bound_ci_thr(:,:,session)=Higher_bound_ci(:,:,session);
                            Length_error_bar_thr(:,:,session)=Length_error_bar(:,:,session);
                        end
                        
                        %confidence interval includes zero (not used for
                        %any diagnostics)
                        
                        for connection1=1:ntwrk_size2(network)
                            for connection2=1:ntwrk_size2(network)
                                if connection1==connection2
                                    continue
                                elseif sign(Lower_bound_ci(connection1,connection2,session))~=sign(Higher_bound_ci(connection1,connection2,session))
                                    Posterior_estimates_zero_conf(connection1,connection2,session)=Posterior_estimates(connection1,connection2,session);
                                    Lower_bound_ci_zero_conf(connection1,connection2,session)=Lower_bound_ci(connection1,connection2,session);
                                    Higher_bound_ci_zero_conf(connection1,connection2,session)=Higher_bound_ci(connection1,connection2,session);
                                    Length_error_bar_zero_conf(connection1,connection2,session)=Length_error_bar(connection1,connection2,session);
                                end
                            end
                        end
                        
                        if ~isnan(Posterior_estimates_var(1,1,session))||~isnan(Posterior_estimates_max(1,1,session))||~isnan(Posterior_estimates_par(1,1,session))||~isnan(Posterior_estimates_mot(1,1,session))||~isnan(Posterior_estimates_thr(1,1,session))
                            Posterior_estimates(:,:,session)=NaN;
                        end
                        
                        for connection1=1:ntwrk_size2(network)
                            for connection2=1:ntwrk_size2(network)
                                if ~isnan(Posterior_estimates_zero_conf(connection1,connection2,session))
                                    Posterior_estimates(connection1,connection2,session)=NaN;
                                end
                            end
                        end
                        
                    end
                    
                    mkdir([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/' procedure '/' name_ROI_def '/Full_model/QC']);
                    cd([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/' procedure '/' name_ROI_def '/Full_model/QC']);
             
                    save(['Above_treshold_marks_' ntwrk_abbrev2{network} '.mat'],'Posterior_estimates_var','Lower_bound_ci_var','Higher_bound_ci_var','Length_error_bar_var','Posterior_estimates_max','Lower_bound_ci_max','Higher_bound_ci_max','Length_error_bar_max','Posterior_estimates_par','Lower_bound_ci_par','Higher_bound_ci_par','Length_error_bar_par','Posterior_estimates_mot','Lower_bound_ci_mot','Higher_bound_ci_mot','Length_error_bar_mot','Posterior_estimates_thr','Lower_bound_ci_thr','Higher_bound_ci_thr','Length_error_bar_thr','Posterior_estimates');
                    
                end
            end
        end
    end
end

end
