function PEB_group_mean_paper_GSR(SPM_dir,Work_dir)

%%%%%%%%%%%%%%%%%%%
%NOT HCP DATASETS!
%%%%%%%%%%%%%%%%%%%
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
        
        for N_prec_comp={'single'}                % Also possible with different number of precision components using N_prec_comp={'single','fields','all'}
            
            %%%%%%%%%%%%%%%%%%%
            %Compare GSR-Basic
            %%%%%%%%%%%%%%%%%%%
            if (strcmp(procedure,'Basic')||strcmp(procedure,'GSR'))&&strcmp(all_procedure_names{2},'GSR')
                
                %%%%%%%%%%
                %A-matrix
                %%%%%%%%%%
                %Group_mean_PEB_GSR_A_matrix_paper_GSR(name_ROI_def,ntwrk_name,procedure,N_prec_comp,SPM_dir,Work_dir);
                
                %%%%%%%%%%%%%%
                %hemodynamics
                %%%%%%%%%%%%%%
                
                %Group_mean_PEB_GSR_hemodyna_paper_GSR(name_ROI_def,ntwrk_name,procedure,N_prec_comp,SPM_dir,Work_dir)
                
                %%%%%%%%%%%%%%%%
                %Spectral noise
                %%%%%%%%%%%%%%%%
                
                %Group_mean_PEB_GSR_fluct_paper_GSR(name_ROI_def,ntwrk_name,procedure,N_prec_comp,SPM_dir,Work_dir)
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %All parameters in one PEB model
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                %Group_mean_PEB_GSR_all_paper_GSR(name_ROI_def,ntwrk_name,procedure,N_prec_comp,SPM_dir,Work_dir)
                
                
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%
            %Compare GSR_regr-Basic
            %%%%%%%%%%%%%%%%%%%%%%%%
            if strcmp(procedure,'Basic')||strcmp(procedure,'GSR_regr')
                
                %%%%%%%%%%
                %A-matrix
                %%%%%%%%%%
                Group_mean_PEB_GSR_regr_A_matrix_paper_GSR(name_ROI_def,ntwrk_name,procedure,N_prec_comp,SPM_dir,Work_dir)
                
                %%%%%%%%%%%%%%
                %hemodynamics
                %%%%%%%%%%%%%%
                %Group_mean_PEB_GSR_regr_hemodyna_paper_GSR(name_ROI_def,ntwrk_name,procedure,N_prec_comp,SPM_dir,Work_dir)
                
                %%%%%%%%%%%%%%%%
                %Spectral noise
                %%%%%%%%%%%%%%%%
                %Group_mean_PEB_GSR_regr_fluct_paper_GSR(name_ROI_def,ntwrk_name,procedure,N_prec_comp,SPM_dir,Work_dir)
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %All parameters in one PEB model
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %Group_mean_PEB_GSR_regr_all_paper_GSR(name_ROI_def,ntwrk_name,procedure,N_prec_comp,SPM_dir,Work_dir)
                
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
            
            for N_prec_comp={'single'}                % Also possible with different number of precision components using N_prec_comp={'single','fields'}
                if (strcmp(procedure,'Basic')||strcmp(procedure,'GSR'))&&strcmp(all_procedure_names{2},'GSR')
                    
                    %%%%%%%%%%
                    %A-matrix
                    %%%%%%%%%%
                    %Group_mean_PEB_GSR_A_matrix_combnet_paper_GSR(name_ROI_def,ntwrk_abbrev2,procedure,N_prec_comp,SPM_dir,Work_dir)
                    
                    %%%%%%%%%%%%%%
                    %hemodynamics
                    %%%%%%%%%%%%%%
                    %Group_mean_PEB_GSR_hemodyna_combnet_paper_GSR(name_ROI_def,ntwrk_abbrev2,procedure,N_prec_comp,SPM_dir,Work_dir)
                    
                    %%%%%%%%%%%%%%%%
                    %Spectral noise
                    %%%%%%%%%%%%%%%%
                    %Group_mean_PEB_GSR_fluct_combnet_paper_GSR(name_ROI_def,ntwrk_abbrev2,procedure,N_prec_comp,SPM_dir,Work_dir)
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %All parameters in one PEB model
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %Group_mean_PEB_GSR_all_combnet_paper_GSR(name_ROI_def,ntwrk_abbrev2,procedure,N_prec_comp,SPM_dir,Work_dir)
                    
                end
                
                if strcmp(procedure,'Basic')||strcmp(procedure,'GSR_regr')
                   
                    %%%%%%%%%%
                    %A-matrix
                    %%%%%%%%%%
                    Group_mean_PEB_GSR_regr_A_matrix_combnet_paper_GSR(name_ROI_def,ntwrk_abbrev2,procedure,N_prec_comp,SPM_dir,Work_dir)
                    
                    %%%%%%%%%%%%%%
                    %hemodynamics
                    %%%%%%%%%%%%%%
                    %Group_mean_PEB_GSR_regr_hemodyna_combnet_paper_GSR(name_ROI_def,ntwrk_abbrev2,procedure,N_prec_comp,SPM_dir,Work_dir)
                    
                    %%%%%%%%%%%%%%%%
                    %Spectral noise
                    %%%%%%%%%%%%%%%%
                    Group_mean_PEB_GSR_regr_fluct_combnet_paper_GSR(name_ROI_def,ntwrk_abbrev2,procedure,N_prec_comp,SPM_dir,Work_dir)
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %All parameters in one PEB model
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %Group_mean_PEB_GSR_regr_all_combnet_paper_GSR(name_ROI_def,ntwrk_abbrev2,procedure,N_prec_comp,SPM_dir,Work_dir)
                    
                end
            end
        end
        
        %%%%%%%%%%%%%%%%%
        %Lateralization
        %%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%
        %GSR-Basic
        %%%%%%%%%%%%%%%%%
        if (strcmp(procedure,'Basic')||strcmp(procedure,'GSR'))&&strcmp(all_procedure_names{2},'GSR')
            
            %%%%%%%%%%
            %A-matrix
            %%%%%%%%%%
            for N_prec_comp={'single'}                % Also possible with different number of precision components using N_prec_comp={'single','fields','all'}
                for network_number=1:length(ntwrk_name)
                    load([Work_dir '/Results_paper_GSR/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_group/PEB_A_mean_comp_GSR_group_' cell2mat(N_prec_comp) '_' ntwrk_name{network_number} '.mat'],'PEB')
                    
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
                    
                    save([Work_dir '/Results_paper_GSR/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_group/Lateralization_index_A_comp_GSR_group_' cell2mat(N_prec_comp) '_' ntwrk_name{network_number} '.mat'],'mean_diff','var_of_sum','posterior_probability');
                    clear mean_diff var_of_sum c v C PP T PEB posterior_probability;
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %All parameters in one PEB model
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %
                %             for network_number=1:length(ntwrk_name)
                %                 load([Work_dir '/Results_paper_GSR/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_group/PEB_all_mean_comp_GSR_group_' ntwrk_name{network_number} '.mat'],'PEB')
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
                %                 save([Work_dir '/Results_paper_GSR/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_group/Lateralization_index_all_comp_GSR_group_' ntwrk_name{network_number} '.mat'],'mean_diff','var_of_sum','posterior_probability');
                %                 clear mean_diff var_of_sum c v C PP T PEB posterior_probability;
                %             end
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
                for network_number=1:length(ntwrk_name)
                    load([Work_dir '/Results_paper_GSR/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_group/PEB_A_mean_comp_GSR_regr_group_' cell2mat(N_prec_comp) '_' ntwrk_name{network_number} '.mat'],'PEB')
                    
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
                    
                    save([Work_dir '/Results_paper_GSR/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_group/Lateralization_index_A_comp_GSR_regr_group_' cell2mat(N_prec_comp) '_' ntwrk_name{network_number} '.mat'],'mean_diff','var_of_sum','posterior_probability');
                    clear mean_diff var_of_sum c v C PP T PEB posterior_probability;
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %All parameters in one PEB model
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %
                %             for network_number=1:length(ntwrk_name)
                %                 load([Work_dir '/Results_paper_GSR/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_group/PEB_all_mean_comp_GSR_regr_group_' ntwrk_name{network_number} '.mat'],'PEB')
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
                %                 save([Work_dir '/Results_paper_GSR/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_group/Lateralization_index_all_comp_GSR_regr_group_' ntwrk_name{network_number} '.mat'],'mean_diff','var_of_sum','posterior_probability');
                %                 clear mean_diff var_of_sum c v C PP T PEB posterior_probability;
                %             end
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%
        %BN effects
        %%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%
        %GSR-Basic
        %%%%%%%%%%%%%%%%%
        if exist('ROI_list2')
            if (strcmp(procedure,'Basic')||strcmp(procedure,'GSR'))&&strcmp(all_procedure_names{2},'GSR')
                
                %%%%%%%%%%
                %A-matrix
                %%%%%%%%%%
                for N_prec_comp={'single'}                % Also possible with different number of precision components using N_prec_comp={'single','fields','all'}
                    for network_number=1:length(ntwrk_abbrev2)
                        load([Work_dir '/Results_paper_GSR/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_group/PEB_A_mean_comp_GSR_group_' cell2mat(N_prec_comp) '_' ntwrk_abbrev2{network_number} '.mat'],'PEB')
                        
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
                            
                            save([Work_dir '/Results_paper_GSR/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_group/BN_effect_A_comp_GSR_group_' cell2mat(N_prec_comp) '_' ntwrk_abbrev2{network_number} '.mat'],'c1','v1','PP1','c2','v2','PP2','c3','v3','PP3','c4','v4','PP4','c5','v5','PP5','c6','v6','PP6');
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
            
            
            %%%%%%%%%%%%%%%%%%%%
            %GSR_regr - basic
            %%%%%%%%%%%%%%%%%%%
            if strcmp(procedure,'Basic')||strcmp(procedure,'GSR_regr')
                
                %%%%%%%%%%
                %A-matrix
                %%%%%%%%%%
                for N_prec_comp={'single'}                % Also possible with different number of precision components using N_prec_comp={'single','fields','all'}
                    for network_number=1:length(ntwrk_abbrev2)
                        load([Work_dir '/Results_paper_GSR/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_group/PEB_A_mean_comp_GSR_regr_group_' cell2mat(N_prec_comp) '_' ntwrk_abbrev2{network_number} '.mat'],'PEB')
                        
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
                            
                            save([Work_dir '/Results_paper_GSR/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_group/BN_effect_A_comp_GSR_regr_group_' cell2mat(N_prec_comp) '_' ntwrk_abbrev2{network_number} '.mat'],'c1','v1','PP1','c2','v2','PP2','c3','v3','PP3','c4','v4','PP4','c5','v5','PP5','c6','v6','PP6');
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