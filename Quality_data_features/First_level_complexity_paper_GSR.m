function First_level_complexity_paper_GSR(dataset,number_subject,name_ROI_def,ROI_list,ROI_list2,procedure,all_procedure_names,SPM_dir,Work_dir)

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%First-level complexity
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%First level complexity is independent of number of precision 
%components (since these are second- and third level parameters) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~strcmp(dataset,'DatasetHCP')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Within-network connectivity
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for subject=1:number_subject
        cd([Work_dir '/' dataset '/']);
        session_number=dir(['sub-' sprintf('%03d',subject)]);
        session_number(1:2)=[];  

        disp(dataset);
        disp(num2str(subject));
        disp(procedure);
        
        %%%%%%%%%%%%%%%%%%%%%%
        %GSR-Basic comparison
        %%%%%%%%%%%%%%%%%%%%%%
        
        if (strcmp(procedure,'Basic')||strcmp(procedure,'GSR'))&&strcmp(all_procedure_names{2},'GSR')
            
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
                
                lastwarn('');
                try
                    
                    load([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/PEB_A_mean_comp_GSR_single_' ntwrk_name{network_number} '.mat'],'PEB','DCM','GCM');
                    
                    %initiate variables
                    pC = NaN(ntwrk_size(network_number)^2,ntwrk_size(network_number)^2,length(GCM));
                    ipC = pC;
                    qC = pC;
                    pE = NaN(ntwrk_size(network_number)^2,length(GCM));
                    qE = pE;
                    
                    for session=1:length(GCM)
                        %negative entropy
                        FiLe_negent(session) =  -0.5 * spm_logdet(2 * pi * exp(1) * GCM{session}.Cp(1:ntwrk_size(network_number)^2,1:ntwrk_size(network_number)^2)); %adapted from https://github.com/spm/spm12/blob/master/spm_dcm_bdc.m
                        
                        pC(:,:,session) = GCM{session}.M.pC(1:ntwrk_size(network_number)^2,1:ntwrk_size(network_number)^2);
                        ipC(:,:,session) = pinv(pC(:,:,session));
                        k = rank(full(pC(:,:,session)));
                        pE(:,session) = spm_vec(GCM{session}.M.pE.A);
                        
                        qC(:,:,session) = GCM{session}.Cp(1:ntwrk_size(network_number)^2,1:ntwrk_size(network_number)^2);
                        qE(:,session) = spm_vec(GCM{session}.Ep.A);
                        
                        %divergence between prior and posterior as normal distributions
                        FiLe_complex_Normal_Distr(session) = 0.5 * (trace(ipC(:,:,session)*qC(:,:,session)) + (pE(:,session) - qE(:,session))'*ipC(:,:,session)*(pE(:,session) - qE(:,session)) - k + log(det(pC(:,:,session)) / det(qC(:,:,session))));
                        
                        %complexity term of Free energy
                        FiLe_complex_ind(session)=-sum(GCM{session}.L(2:3));
                    end
                    
                    FiLe_complex_avg=sum(FiLe_complex_ind);
                    save([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/First_level_complexity_GSR_' ntwrk_name{network_number} '.mat'],'FiLe_complex_ind','FiLe_complex_avg','FiLe_complex_Normal_Distr','FiLe_negent','pC','qC','pE','qE','ipC');
                    
                catch
                    continue
                end
                
                clear PEB DCM;
                
                clear GCM Posterior_estimates_var Lower_bound_ci_var Higher_bound_ci_var Length_error_bar_var Posterior_estimates_max Lower_bound_ci_max Higher_bound_ci_max Length_error_bar_max Posterior_estimates_par Lower_bound_ci_par Higher_bound_ci_par Length_error_bar_par Posterior_estimates_mot Lower_bound_ci_mot Higher_bound_ci_mot Length_error_bar_mot Posterior_estimates_thr Lower_bound_ci_thr Higher_bound_ci_thr Length_error_bar_thr Posterior_estimates FiLe_complex_avg FiLe_complex_ind;
                
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %GSR_regr - Basic comparison
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
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
                
                try
                    
                    load([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/PEB_A_mean_comp_GSR_regr_single_'  ntwrk_name{network_number} '.mat'],'PEB','DCM','GCM');
                    %initiate variables
                    pC = NaN(ntwrk_size(network_number)^2,ntwrk_size(network_number)^2,length(GCM));
                    ipC = pC;
                    qC = pC;
                    pE = NaN(ntwrk_size(network_number)^2,length(GCM));
                    qE = pE;
                    
                    for session=1:length(GCM)
                        %negative entropy
                        FiLe_negent(session) =  -0.5 * spm_logdet(2 * pi * exp(1) * GCM{session}.Cp(1:ntwrk_size(network_number)^2,1:ntwrk_size(network_number)^2)); %adapted from https://github.com/spm/spm12/blob/master/spm_dcm_bdc.m
                        
                        pC(:,:,session) = GCM{session}.M.pC(1:ntwrk_size(network_number)^2,1:ntwrk_size(network_number)^2);
                        ipC(:,:,session) = pinv(pC(:,:,session));
                        k = rank(full(pC(:,:,session)));
                        pE(:,session) = spm_vec(GCM{session}.M.pE.A);
                        
                        qC(:,:,session) = GCM{session}.Cp(1:ntwrk_size(network_number)^2,1:ntwrk_size(network_number)^2);
                        qE(:,session) = spm_vec(GCM{session}.Ep.A);
                        
                        %divergence between prior and posterior as normal distributions
                        FiLe_complex_Normal_Distr(session) = 0.5 * (trace(ipC(:,:,session)*qC(:,:,session)) + (pE(:,session) - qE(:,session))'*ipC(:,:,session)*(pE(:,session) - qE(:,session)) - k + log(det(pC(:,:,session)) / det(qC(:,:,session))));
                        
                        %complexity term of Free energy
                        FiLe_complex_ind(session)=-sum(GCM{session}.L(2:3));
                    end
                    FiLe_complex_avg=sum(FiLe_complex_ind);
                    save([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/First_level_complexity_GSR_regr_' ntwrk_name{network_number} '.mat'],'FiLe_complex_ind','FiLe_complex_avg','FiLe_complex_Normal_Distr','FiLe_negent','pC','qC','pE','qE','ipC');
                    
                catch
                    continue;
                end
                
                clear PEB DCM;
                
                clear GCM Posterior_estimates_var Lower_bound_ci_var Higher_bound_ci_var Length_error_bar_var Posterior_estimates_max Lower_bound_ci_max Higher_bound_ci_max Length_error_bar_max Posterior_estimates_par Lower_bound_ci_par Higher_bound_ci_par Length_error_bar_par Posterior_estimates_mot Lower_bound_ci_mot Higher_bound_ci_mot Length_error_bar_mot Posterior_estimates_thr Lower_bound_ci_thr Higher_bound_ci_thr Length_error_bar_thr Posterior_estimates FiLe_complex_avg FiLe_complex_ind 'FiLe_complex_Normal_Distr' 'FiLe_negent' 'pC' 'qC' 'pE' 'qE' 'ipC' k;
                
            end
            
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Between-network connectivity
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for subject=1:number_subject
        cd([Work_dir '/' dataset '/']);
        session_number=dir(['sub-' sprintf('%03d',subject)]);
        session_number(1:2)=[];

        disp(dataset);
        disp(num2str(subject));
        disp(procedure);
        
        if  exist('ROI_list2') %Special function if second ROI list exists
            
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
            
            %%%%%%%%%%%%%%%%%%%%%%
            %GSR-Basic comparison
            %%%%%%%%%%%%%%%%%%%%%%
            
            if (strcmp(procedure,'Basic')||strcmp(procedure,'GSR'))&&strcmp(all_procedure_names{2},'GSR')
                
                for network_number=1:length(ntwrk_abbrev2)
                    disp(ntwrk_abbrev2{network_number});
                    
                    lastwarn('');
                    try
                        
                        load([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/PEB_A_mean_comp_GSR_single_' ntwrk_abbrev2{network_number} '.mat'],'PEB','DCM','GCM');
                        NW_size = size(GCM{1}.Ep.A,1);
                        
                        %initiate variables
                        pC = NaN(NW_size^2,NW_size^2,length(GCM));
                        ipC = pC;
                        qC = pC;
                        pE = NaN(NW_size^2,length(GCM));
                        qE = pE;
                        
                        for session=1:length(GCM)
                            %negative entropy
                            FiLe_negent(session) =  -0.5 * spm_logdet(2 * pi * exp(1) * GCM{session}.Cp(1:NW_size^2,1:NW_size^2)); %adapted from https://github.com/spm/spm12/blob/master/spm_dcm_bdc.m
                            
                            pC(:,:,session) = GCM{session}.M.pC(1:NW_size^2,1:NW_size^2);
                            ipC(:,:,session) = pinv(pC(:,:,session));
                            k = rank(full(pC(:,:,session)));
                            pE(:,session) = spm_vec(GCM{session}.M.pE.A);
                            
                            qC(:,:,session) = GCM{session}.Cp(1:NW_size^2,1:NW_size^2);
                            qE(:,session) = spm_vec(GCM{session}.Ep.A);
                            
                            %divergence between prior and posterior as normal distributions
                            FiLe_complex_Normal_Distr(session) = 0.5 * (trace(ipC(:,:,session)*qC(:,:,session)) + (pE(:,session) - qE(:,session))'*ipC(:,:,session)*(pE(:,session) - qE(:,session)) - k + log(det(pC(:,:,session)/qC(:,:,session)/100))+log(100^(NW_size^2)));
                            
                            %complexity term of Free energy
                            FiLe_complex_ind(session)=-sum(GCM{session}.L(2:3));
                        end
                        FiLe_complex_avg=sum(FiLe_complex_ind);
                        save([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/First_level_complexity_GSR_' ntwrk_abbrev2{network_number} '.mat'],'FiLe_complex_ind','FiLe_complex_avg','FiLe_complex_Normal_Distr','FiLe_negent','pC','qC','pE','qE','ipC');
                        
                    catch
                        continue
                    end
                    
                    clear PEB DCM;
                    
                    clear GCM Posterior_estimates_var Lower_bound_ci_var Higher_bound_ci_var Length_error_bar_var Posterior_estimates_max Lower_bound_ci_max Higher_bound_ci_max Length_error_bar_max Posterior_estimates_par Lower_bound_ci_par Higher_bound_ci_par Length_error_bar_par Posterior_estimates_mot Lower_bound_ci_mot Higher_bound_ci_mot Length_error_bar_mot Posterior_estimates_thr Lower_bound_ci_thr Higher_bound_ci_thr Length_error_bar_thr Posterior_estimates FiLe_complex_avg FiLe_complex_ind 'FiLe_complex_Normal_Distr' 'FiLe_negent' 'pC' 'qC' 'pE' 'qE' 'ipC' k;
                    
                end
                
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %GSR_regr - Basic comparison
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if strcmp(procedure,'Basic')||strcmp(procedure,'GSR_regr')
                
                for network_number=1:length(ntwrk_abbrev2)
                    disp(ntwrk_abbrev2{network_number});
                    
                    try
                        
                        load([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/PEB_A_mean_comp_GSR_regr_single_' ntwrk_abbrev2{network_number} '.mat'],'PEB','DCM','GCM');
                        NW_size = size(GCM{1}.Ep.A,1);
                        
                        %initiate variables
                        pC = NaN(NW_size^2,NW_size^2,length(GCM));
                        ipC = pC;
                        qC = pC;
                        pE = NaN(NW_size^2,length(GCM));
                        qE = pE;
                        
                        for session=1:length(GCM)
                            %negative entropy
                            FiLe_negent(session) =  -0.5 * spm_logdet(2 * pi * exp(1) * GCM{session}.Cp(1:NW_size^2,1:NW_size^2)); %adapted from https://github.com/spm/spm12/blob/master/spm_dcm_bdc.m
                            
                            pC(:,:,session) = GCM{session}.M.pC(1:NW_size^2,1:NW_size^2);
                            ipC(:,:,session) = pinv(pC(:,:,session));
                            k = rank(full(pC(:,:,session)));
                            pE(:,session) = spm_vec(GCM{session}.M.pE.A);
                            
                            qC(:,:,session) = GCM{session}.Cp(1:NW_size^2,1:NW_size^2);
                            qE(:,session) = spm_vec(GCM{session}.Ep.A);
                            
                            %divergence between prior and posterior as normal distributions
                            FiLe_complex_Normal_Distr(session) = 0.5 * (trace(ipC(:,:,session)*qC(:,:,session)) + (pE(:,session) - qE(:,session))'*ipC(:,:,session)*(pE(:,session) - qE(:,session)) - k + log(det(pC(:,:,session)/qC(:,:,session)/100))+log(100^(NW_size^2)));
                            
                            %complexity term of Free energy
                            FiLe_complex_ind(session)=-sum(GCM{session}.L(2:3));
                        end
                        FiLe_complex_avg=sum(FiLe_complex_ind);
                        save([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/First_level_complexity_GSR_regr_' ntwrk_abbrev2{network_number} '.mat'],'FiLe_complex_ind','FiLe_complex_avg','FiLe_complex_Normal_Distr','FiLe_negent','pC','qC','pE','qE','ipC');
                        
                    catch
                        continue;
                    end
                    
                    clear PEB DCM;
                    
                    clear GCM Posterior_estimates_var Lower_bound_ci_var Higher_bound_ci_var Length_error_bar_var Posterior_estimates_max Lower_bound_ci_max Higher_bound_ci_max Length_error_bar_max Posterior_estimates_par Lower_bound_ci_par Higher_bound_ci_par Length_error_bar_par Posterior_estimates_mot Lower_bound_ci_mot Higher_bound_ci_mot Length_error_bar_mot Posterior_estimates_thr Lower_bound_ci_thr Higher_bound_ci_thr Length_error_bar_thr Posterior_estimates FiLe_complex_avg FiLe_complex_ind 'FiLe_complex_Normal_Distr' 'FiLe_negent' 'pC' 'qC' 'pE' 'qE' 'ipC' k;
                    
                end
            end
        end
    end
end


%%%%%%%%%%%%%
%Dataset HCP
%%%%%%%%%%%%%

if strcmp(dataset,'DatasetHCP')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Within network connectivity
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for subject=1:number_subject
        cd([Work_dir '/' dataset '/']);
        session_number=dir(['sub-' sprintf('%03d',subject)]);
        session_number(1:2)=[];

        disp(dataset);
        disp(num2str(subject));
        disp(procedure);
        
        %%%%%%%%%%%%%%%%%%%%%%
        %GSR-Basic comparison
        %%%%%%%%%%%%%%%%%%%%%%
        
        if (strcmp(procedure,'Basic')||strcmp(procedure,'GSR'))&&strcmp(all_procedure_names{2},'GSR')
            
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
                
                lastwarn('');
                try
                    
                    cd([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_summary/DCM/'  procedure '/' name_ROI_def '/Full_model/']);
                    load(['GCM_' ntwrk_name{network_number} '_full_estim.mat'],'GCM');
                    %initiate variables
                    pC = NaN(ntwrk_size(network_number)^2,ntwrk_size(network_number)^2,length(GCM));
                    ipC = pC;
                    qC = pC;
                    pE = NaN(ntwrk_size(network_number)^2,length(GCM));
                    qE = pE;
                    
                    for session=1:length(GCM)
                        %negative entropy
                        FiLe_negent(session) =  -0.5 * spm_logdet(2 * pi * exp(1) * GCM{session}.Cp(1:ntwrk_size(network_number)^2,1:ntwrk_size(network_number)^2)); %adapted from https://github.com/spm/spm12/blob/master/spm_dcm_bdc.m
                        
                        pC(:,:,session) = GCM{session}.M.pC(1:ntwrk_size(network_number)^2,1:ntwrk_size(network_number)^2);
                        ipC(:,:,session) = pinv(pC(:,:,session));
                        k = rank(full(pC(:,:,session)));
                        pE(:,session) = spm_vec(GCM{session}.M.pE.A);
                        
                        qC(:,:,session) = GCM{session}.Cp(1:ntwrk_size(network_number)^2,1:ntwrk_size(network_number)^2);
                        qE(:,session) = spm_vec(GCM{session}.Ep.A);
                        
                        %divergence between prior and posterior as normal distributions
                        FiLe_complex_Normal_Distr(session) = 0.5 * (trace(ipC(:,:,session)*qC(:,:,session)) + (pE(:,session) - qE(:,session))'*ipC(:,:,session)*(pE(:,session) - qE(:,session)) - k + log(det(pC(:,:,session)) / det(qC(:,:,session))));
                        
                        %complexity term of Free energy
                        FiLe_complex_ind(session)=-sum(GCM{session}.L(2:3));
                    end
                    
                    FiLe_complex_avg=sum(FiLe_complex_ind);
                    save([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/First_level_complexity_GSR_' ntwrk_name{network_number} '.mat'],'FiLe_complex_ind','FiLe_complex_avg','FiLe_complex_Normal_Distr','FiLe_negent','pC','qC','pE','qE','ipC');
                    
                catch
                    continue
                end
                
                clear PEB DCM;
                
                clear GCM Posterior_estimates_var Lower_bound_ci_var Higher_bound_ci_var Length_error_bar_var Posterior_estimates_max Lower_bound_ci_max Higher_bound_ci_max Length_error_bar_max Posterior_estimates_par Lower_bound_ci_par Higher_bound_ci_par Length_error_bar_par Posterior_estimates_mot Lower_bound_ci_mot Higher_bound_ci_mot Length_error_bar_mot Posterior_estimates_thr Lower_bound_ci_thr Higher_bound_ci_thr Length_error_bar_thr Posterior_estimates FiLe_complex_avg FiLe_complex_ind 'FiLe_complex_Normal_Distr' 'FiLe_negent' 'pC' 'qC' 'pE' 'qE' 'ipC' k;
                
            end
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %GSR_regr - Basic comparison
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
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
                
                try
                    
                    cd([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_summary/DCM/'  procedure '/' name_ROI_def '/Full_model/']);
                    load(['GCM_' ntwrk_name{network_number} '_full_estim.mat'],'GCM');
                    %initiate variables
                    pC = NaN(ntwrk_size(network_number)^2,ntwrk_size(network_number)^2,length(GCM));
                    ipC = pC;
                    qC = pC;
                    pE = NaN(ntwrk_size(network_number)^2,length(GCM));
                    qE = pE;
                    
                    for session=1:length(GCM)
                        %negative entropy
                        FiLe_negent(session) =  -0.5 * spm_logdet(2 * pi * exp(1) * GCM{session}.Cp(1:ntwrk_size(network_number)^2,1:ntwrk_size(network_number)^2)); %adapted from https://github.com/spm/spm12/blob/master/spm_dcm_bdc.m
                        
                        pC(:,:,session) = GCM{session}.M.pC(1:ntwrk_size(network_number)^2,1:ntwrk_size(network_number)^2);
                        ipC(:,:,session) = pinv(pC(:,:,session));
                        k = rank(full(pC(:,:,session)));
                        pE(:,session) = spm_vec(GCM{session}.M.pE.A);
                        
                        qC(:,:,session) = GCM{session}.Cp(1:ntwrk_size(network_number)^2,1:ntwrk_size(network_number)^2);
                        qE(:,session) = spm_vec(GCM{session}.Ep.A);
                        
                        %divergence between prior and posterior as normal distributions
                        FiLe_complex_Normal_Distr(session) = 0.5 * (trace(ipC(:,:,session)*qC(:,:,session)) + (pE(:,session) - qE(:,session))'*ipC(:,:,session)*(pE(:,session) - qE(:,session)) - k + log(det(pC(:,:,session)) / det(qC(:,:,session))));
                        
                        %complexity term of Free energy
                        FiLe_complex_ind(session)=-sum(GCM{session}.L(2:3));
                    end
                    
                    FiLe_complex_avg=sum(FiLe_complex_ind);
                    save([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/First_level_complexity_GSR_regr_' ntwrk_name{network_number} '.mat'],'FiLe_complex_ind','FiLe_complex_avg','FiLe_complex_Normal_Distr','FiLe_negent','pC','qC','pE','qE','ipC');
                    
                catch
                    continue;
                end
                
                clear PEB DCM;
                
                clear GCM Posterior_estimates_var Lower_bound_ci_var Higher_bound_ci_var Length_error_bar_var Posterior_estimates_max Lower_bound_ci_max Higher_bound_ci_max Length_error_bar_max Posterior_estimates_par Lower_bound_ci_par Higher_bound_ci_par Length_error_bar_par Posterior_estimates_mot Lower_bound_ci_mot Higher_bound_ci_mot Length_error_bar_mot Posterior_estimates_thr Lower_bound_ci_thr Higher_bound_ci_thr Length_error_bar_thr Posterior_estimates FiLe_complex_avg FiLe_complex_ind 'FiLe_complex_Normal_Distr' 'FiLe_negent' 'pC' 'qC' 'pE' 'qE' 'ipC' k;
                
            end
            
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Between network connectivity
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for subject=1:number_subject
        cd([Work_dir '/' dataset '/']);
        session_number=dir(['sub-' sprintf('%03d',subject)]);
        session_number(1:2)=[];

        disp(dataset);
        disp(num2str(subject));
        disp(procedure);
        
        if  exist('ROI_list2') %Special function if second ROI list exists
            
            %collect regions of interest
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
            
            %%%%%%%%%%%%%%%%%%%%%%
            %GSR-Basic comparison
            %%%%%%%%%%%%%%%%%%%%%%
            
            if (strcmp(procedure,'Basic')||strcmp(procedure,'GSR'))&&strcmp(all_procedure_names{2},'GSR')
                
                for network_number=1:length(ntwrk_abbrev2)
                    disp(ntwrk_abbrev2{network_number});
                    
                    lastwarn('');
                    try
                        
                        cd([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_summary/DCM/'  procedure '/' name_ROI_def '/Full_model/']);
                        load(['GCM_' ntwrk_abbrev2{network_number} '_full_estim.mat'],'GCM');
                        
                        NW_size = size(GCM{1}.Ep.A,1);
                        
                        %initiate variables
                        pC = NaN(NW_size^2,NW_size^2,length(GCM));
                        ipC = pC;
                        qC = pC;
                        pE = NaN(NW_size^2,length(GCM));
                        qE = pE;
                        
                        for session=1:length(GCM)
                            %negative entropy
                            FiLe_negent(session) =  -0.5 * spm_logdet(2 * pi * exp(1) * GCM{session}.Cp(1:NW_size^2,1:NW_size^2)); %adapted from https://github.com/spm/spm12/blob/master/spm_dcm_bdc.m
                            
                            pC(:,:,session) = GCM{session}.M.pC(1:NW_size^2,1:NW_size^2);
                            ipC(:,:,session) = pinv(pC(:,:,session));
                            k = rank(full(pC(:,:,session)));
                            pE(:,session) = spm_vec(GCM{session}.M.pE.A);
                            
                            qC(:,:,session) = GCM{session}.Cp(1:NW_size^2,1:NW_size^2);
                            qE(:,session) = spm_vec(GCM{session}.Ep.A);
                            
                            %divergence between prior and posterior as normal distributions
                            FiLe_complex_Normal_Distr(session) = 0.5 * (trace(ipC(:,:,session)*qC(:,:,session)) + (pE(:,session) - qE(:,session))'*ipC(:,:,session)*(pE(:,session) - qE(:,session)) - k + log(det(pC(:,:,session)/qC(:,:,session)/100))+log(100^(NW_size^2)));
                            
                            %complexity term of Free energy
                            FiLe_complex_ind(session)=-sum(GCM{session}.L(2:3));
                        end
                        
                        FiLe_complex_avg=sum(FiLe_complex_ind);
                        save([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/First_level_complexity_GSR_' ntwrk_abbrev2{network_number} '.mat'],'FiLe_complex_ind','FiLe_complex_avg','FiLe_complex_Normal_Distr','FiLe_negent','pC','qC','pE','qE','ipC');
                        
                    catch
                        continue
                    end
                    
                    clear PEB DCM;
                    
                    clear GCM Posterior_estimates_var Lower_bound_ci_var Higher_bound_ci_var Length_error_bar_var Posterior_estimates_max Lower_bound_ci_max Higher_bound_ci_max Length_error_bar_max Posterior_estimates_par Lower_bound_ci_par Higher_bound_ci_par Length_error_bar_par Posterior_estimates_mot Lower_bound_ci_mot Higher_bound_ci_mot Length_error_bar_mot Posterior_estimates_thr Lower_bound_ci_thr Higher_bound_ci_thr Length_error_bar_thr Posterior_estimates FiLe_complex_avg FiLe_complex_ind 'FiLe_complex_Normal_Distr' 'FiLe_negent' 'pC' 'qC' 'pE' 'qE' 'ipC' k;
                    
                end
                
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %GSR_regr - Basic Comparison
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if strcmp(procedure,'Basic')||strcmp(procedure,'GSR_regr')
                
                for network_number=1:length(ntwrk_abbrev2)
                    disp(ntwrk_abbrev2{network_number});
                    
                    try
                        
                        cd([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_summary/DCM/'  procedure '/' name_ROI_def '/Full_model/']);
                        load(['GCM_' ntwrk_abbrev2{network_number} '_full_estim.mat'],'GCM');
                        
                        NW_size = size(GCM{1}.Ep.A,1);
                        
                        %initiate variables
                        pC = NaN(NW_size^2,NW_size^2,length(GCM));
                        ipC = pC;
                        qC = pC;
                        pE = NaN(NW_size^2,length(GCM));
                        qE = pE;
                        
                        for session=1:length(GCM)
                            %negative entropy
                            FiLe_negent(session) =  -0.5 * spm_logdet(2 * pi * exp(1) * GCM{session}.Cp(1:NW_size^2,1:NW_size^2)); %adapted from https://github.com/spm/spm12/blob/master/spm_dcm_bdc.m
                            
                            pC(:,:,session) = GCM{session}.M.pC(1:NW_size^2,1:NW_size^2);
                            ipC(:,:,session) = pinv(pC(:,:,session));
                            k = rank(full(pC(:,:,session)));
                            pE(:,session) = spm_vec(GCM{session}.M.pE.A);
                            
                            qC(:,:,session) = GCM{session}.Cp(1:NW_size^2,1:NW_size^2);
                            qE(:,session) = spm_vec(GCM{session}.Ep.A);
                            
                            %divergence between prior and posterior as normal distributions
                            FiLe_complex_Normal_Distr(session) = 0.5 * (trace(ipC(:,:,session)*qC(:,:,session)) + (pE(:,session) - qE(:,session))'*ipC(:,:,session)*(pE(:,session) - qE(:,session)) - k + log(det(pC(:,:,session)/qC(:,:,session)/100))+log(100^(NW_size^2)));
                            
                            %complexity term of Free energy
                            FiLe_complex_ind(session)=-sum(GCM{session}.L(2:3));
                        end
                        
                        FiLe_complex_avg=sum(FiLe_complex_ind);
                        save([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_results/DCM/'  procedure '/' name_ROI_def '/Full_model/First_level_complexity_GSR_regr_' ntwrk_abbrev2{network_number} '.mat'],'FiLe_complex_ind','FiLe_complex_avg','FiLe_complex_Normal_Distr','FiLe_negent','pC','qC','pE','qE','ipC');
                        
                    catch
                        continue;
                    end
                    
                    clear PEB DCM;
                    
                    clear GCM Posterior_estimates_var Lower_bound_ci_var Higher_bound_ci_var Length_error_bar_var Posterior_estimates_max Lower_bound_ci_max Higher_bound_ci_max Length_error_bar_max Posterior_estimates_par Lower_bound_ci_par Higher_bound_ci_par Length_error_bar_par Posterior_estimates_mot Lower_bound_ci_mot Higher_bound_ci_mot Length_error_bar_mot Posterior_estimates_thr Lower_bound_ci_thr Higher_bound_ci_thr Length_error_bar_thr Posterior_estimates FiLe_complex_avg FiLe_complex_ind 'FiLe_complex_Normal_Distr' 'FiLe_negent' 'pC' 'qC' 'pE' 'qE' 'ipC' k;
                    
                end
            end
        end
    end
end

end