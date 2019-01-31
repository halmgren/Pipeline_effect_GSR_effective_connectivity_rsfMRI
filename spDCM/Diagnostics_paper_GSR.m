function Diagnostics_paper_GSR(dataset,number_subject,ROI_list,ROI_list2,name_ROI_def,procedure,SPM_dir,Work_dir)

for subject=1:number_subject
    cd([Work_dir '/' dataset '/']);
    session_number=dir(['sub-' sprintf('%03d',subject)]);
    session_number(1:2)=[];
    
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
        mkdir([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_summary/DCM/' procedure '/' name_ROI_def '/Full_model/QC']);
        
        for network=1:length(ntwrk_name)
            
            cd([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_summary/DCM/' procedure '/' name_ROI_def '/Full_model']);
            load(['GCM_' ntwrk_name{network} '_full_estim.mat']);
            
            for session=1:length(session_number)
                session_name=session_number(session).name;
                
                disp(['Dataset: ' dataset '; subject: ' num2str(subject) '; session: ' session_name]);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %Number of estimable parameters
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                qE    = spm_vec(GCM{session}.Ep);
                pE    = spm_vec(GCM{session}.M.pE);
                qC    = GCM{session}.Cp;
                pC    = GCM{session}.M.pC;
                k     = rank(full(pC));
                pC    = pinv(pC);
                
                D  = trace(pC*qC) + (pE - qE)'*pC*(pE - qE) - spm_logdet(qC*pC) - k;
                N_est_par(session)  = D/log(GCM{session}.v);
                
                clear qE pE qC pC k D;
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %maximum extrinsic connection strength
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Max_conn(session)=max(max(abs(GCM{session}.Ep.A - diag(diag(GCM{session}.Ep.A)))));
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %Proportion explained variance
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                PSS   = sum(sum(sum(abs(GCM{session}.Hc).^2)));
                RSS   = sum(sum(sum(abs(GCM{session}.Rc).^2)));
                Expl_var(session) = 100*PSS/(PSS + RSS);
                
                clear RSS PSS;
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %Free energy
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Free_energy(session)=GCM{session}.F;
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %Error bars and MAP estimates
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                qC  = GCM{session}.Cp;
                qE  = spm_vec(GCM{session}.Ep);
                i=spm_fieldindices(GCM{session}.Ep,'A');
                E=qE(i);
                C=qC(i,i);
                
                j = 1:size(E,1);
                
                ci    = spm_invNcdf(1 - 0.05);
                C = diag(C);
                c = ci*sqrt(C(j,:));
                
                E=spm_unvec(E,zeros(size(GCM{session}.Ep.A,1),size(GCM{session}.Ep.A,1)));
                c=spm_unvec(c,zeros(size(GCM{session}.Ep.A,1),size(GCM{session}.Ep.A,1)));
                
                Posterior_estimates(:,:,session)=E;
                Length_error_bar(:,:,session)=c;   %length of half of error bar
                Lower_bound_ci(:,:,session)=E-c;
                Higher_bound_ci(:,:,session)=E+c;
                
                clear qC qE i E C j ci C c;
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %First level thresholds VOIs
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                load([Work_dir '/' dataset '/sub-' sprintf('%03d',subject) '/' session_name '/func/VOI/' procedure '/' name_ROI_def '/FL_thresholds_' ntwrk_name{network} '.mat']);
                
                thresholds_VOIs(:,session)=thresholds;
                
                clear thresholds;
                
                %%%%%%%%
                %Motion
                %%%%%%%%
                if ~strcmp(dataset,'DatasetHCP')
                    if strcmp(procedure,'Basic')||strcmp(procedure,'GSR')||strcmp(procedure,'GSR_regr')
                        cd([Work_dir '/' dataset '/sub-' sprintf('%03d',subject) '/' session_name '/func/fourD_files/Basic']);
                    end
                    rp_name='rp_afunctional_disc_4D.txt';
                    rp = load(deblank(rp_name));
                    radius=50;
                    
                    %Change rotation in radials to mm (arch=radial*radius)
                    rp(:,4:6)=rp(:,4:6)*radius;
                    
                    %Calculate change in position between current and previous volume
                    rp_diff=diff(rp);
                    
                    %add zero as first element (see Power et al., 2014)
                    rp_zero=[zeros(1,size(rp_diff,2)); rp_diff];
                    
                    %compute framewise displacement (FD)
                    %%This particular session showed one scan less, so we added
                    %%a zero to make it same size (only maximum FD was used, so
                    %%doesn't alter any results)
                    if strcmp(dataset,'DatasetGordon')&&strcmp(session_name,'ses-008')&&(subject==6)
                        FD(:,session)=[sum(abs(rp_zero),2);0];
                    else
                        FD(:,session)=sum(abs(rp_zero),2);
                    end
                    
                end
                
                if strcmp(dataset,'DatasetHCP')
                    if strcmp(procedure,'Basic')||strcmp(procedure,'GSR')||strcmp(procedure,'GSR_regr')
                        cd([Work_dir '/' dataset '/sub-' sprintf('%03d',subject) '/' session_name '/func/fourD_files/Basic']);
                    end
                    rp_name='rp_afunctional_disc_4D.txt';
                    rp = load(deblank(rp_name));
                    
                    rp=rp(:,1:6);
                    radius=50;
                    
                    %Change rotation in radials to mm (arch=radial*radius)
                    rp(:,4:6)=rp(:,4:6)*radius/180*pi; %from degrees -> radial
                    
                    %Calculate change in position between current and previous volume
                    rp_diff=diff(rp);
                    
                    %add zero as first element (see Power et al., 2014)
                    rp_zero=[zeros(1,size(rp_diff,2)); rp_diff];
                    
                    %compute framewise displacement (FD)
                    %%This particular session showed one scan less, so we added
                    %%a zero to make it same size (only maximum FD was used, so
                    %%doesn't alter any results)
                    
                    FD(:,session)=sum(abs(rp_zero),2);
                    
                end
                
                clear rp rp_diff rp_zero
            end
            
            save([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_summary/DCM/' procedure '/' name_ROI_def '/Full_model/QC/Diagnostics_' ntwrk_name{network} '.mat'],'N_est_par','Max_conn','Expl_var','Posterior_estimates','Length_error_bar','Lower_bound_ci','Higher_bound_ci','thresholds_VOIs','FD','Free_energy');
            clear N_est_par Max_conn Expl_var Posterior_estimates Length_error_bar Lower_bound_ci Higher_bound_ci thresholds_VOIs FD Free_energy;
        end
    end
%     
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
        
        if strcmp(procedure,'Basic')||strcmp(procedure,'GSR')||strcmp(procedure,'GSR_regr')
            mkdir([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_summary/DCM/' procedure '/' name_ROI_def '/Full_model/QC']);
            
            for network=1:length(ntwrk_abbrev2)
                
                cd([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_summary/DCM/' procedure '/' name_ROI_def '/Full_model']);
                load(['GCM_' ntwrk_abbrev2{network} '_full_estim.mat']);
                
                for session=1:length(session_number)
                    

                    session_name=session_number(session).name;
                    disp(['Dataset: ' dataset '; subject: ' num2str(subject) '; session: ' session_name]);
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %Number of estimable parameters
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    qE    = spm_vec(GCM{session}.Ep);
                    pE    = spm_vec(GCM{session}.M.pE);
                    qC    = GCM{session}.Cp;
                    pC    = GCM{session}.M.pC;
                    k     = rank(full(pC));
                    pC    = pinv(pC);
                    
                    D  = trace(pC*qC) + (pE - qE)'*pC*(pE - qE) - spm_logdet(qC*pC) - k;
                    N_est_par(session)  = D/log(GCM{session}.v);
                    
                    clear qE pE qC pC k D;
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %maximum extrinsic connection strength
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    Max_conn(session)=max(max(abs(GCM{session}.Ep.A - diag(diag(GCM{session}.Ep.A)))));
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %Proportion explained variance
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    PSS   = sum(sum(sum(abs(GCM{session}.Hc).^2)));
                    RSS   = sum(sum(sum(abs(GCM{session}.Rc).^2)));
                    Expl_var(session) = 100*PSS/(PSS + RSS);
                    
                    clear RSS PSS;
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %Free energy
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    Free_energy(session)=GCM{session}.F;
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %Error bars and MAP estimates
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    qC  = GCM{session}.Cp;
                    qE  = spm_vec(GCM{session}.Ep);
                    i=spm_fieldindices(GCM{session}.Ep,'A');
                    E=qE(i);
                    C=qC(i,i);
                    
                    j = 1:size(E,1);
                    
                    ci    = spm_invNcdf(1 - 0.05);
                    C = diag(C);
                    c = ci*sqrt(C(j,:));
                    
                    E=spm_unvec(E,zeros(size(GCM{session}.Ep.A,1),size(GCM{session}.Ep.A,1)));
                    c=spm_unvec(c,zeros(size(GCM{session}.Ep.A,1),size(GCM{session}.Ep.A,1)));
                    
                    Posterior_estimates(:,:,session)=E;
                    Length_error_bar(:,:,session)=c;        %length of half of error bar
                    Lower_bound_ci(:,:,session)=E-c;
                    Higher_bound_ci(:,:,session)=E+c;
                    
                    clear qC qE i E C j ci C c;
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %First level thresholds VOIs
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %Begin new
                    cd([Work_dir '/' dataset '/sub-' sprintf('%03d',subject) '/' session_name '/func/VOI/' procedure '/' name_ROI_def]);
                    [pos_network,orig_network_abbrev2]=Rearrange_ROI_list2_thresholds(ROI_list,ROI_list2);
                    
                    tmp=1;
                    tmp2=0;
                    for tmpX=1:size(orig_network_abbrev2{network},1)
                        
                        load([Work_dir '/' dataset '/sub-' sprintf('%03d',subject) '/' session_name '/func/VOI/' procedure '/' name_ROI_def '/FL_thresholds_' orig_network_abbrev2{network}(tmpX,:) '.mat']);
                        
                        
                        if tmpX~=1
                            if ~all(orig_network_abbrev2{network}(tmpX,:)==orig_network_abbrev2{network}(tmpX-1,:))
                                tmp=tmp+1;
                                tmp2=0;
                            end
                        end
                        tmp2=tmp2+1;
                        thresholds_VOIs(tmpX,session)=thresholds(pos_network{network}{tmp}(tmp2));
                        
                    end
                    
                    thresholds=thresholds_VOIs(:,session)';
                    
                    save([Work_dir '/' dataset '/sub-' sprintf('%03d',subject) '/' session_name '/func/VOI/' procedure '/' name_ROI_def '/FL_thresholds_' ntwrk_abbrev2{network} '.mat'],'thresholds');
                    
                    %end new
                    
                    clear thresholds;
                    
                    %%%%%%%%
                    %Motion
                    %%%%%%%%
                    if ~strcmp(dataset,'DatasetHCP')
                        if strcmp(procedure,'Basic')||strcmp(procedure,'GSR')||strcmp(procedure,'GSR_regr')
                            cd([Work_dir '/' dataset '/sub-' sprintf('%03d',subject) '/' session_name '/func/fourD_files/Basic']);
                        end
                        rp_name='rp_afunctional_disc_4D.txt';
                        rp = load(deblank(rp_name));
                        radius=50;
                        
                        %Change rotation in radials to mm (arch=radial*radius)
                        rp(:,4:6)=rp(:,4:6)*radius;
                        
                        %Calculate change in position between current and previous volume
                        rp_diff=diff(rp);
                        
                        %add zero as first element (see Power et al., 2014)
                        rp_zero=[zeros(1,size(rp_diff,2)); rp_diff];
                        
                        %compute framewise displacement (FD)
                        %%This particular session showed one scan less, so we added
                        %%a zero to make it same size (only maximum FD was used, so
                        %%doesn't alter any results)
                        if strcmp(dataset,'DatasetGordon')&&strcmp(session_name,'ses-008')&&(subject==6)
                            FD(:,session)=[sum(abs(rp_zero),2);0];
                        else
                            FD(:,session)=sum(abs(rp_zero),2);
                        end
                        
                    end
                    
                    if strcmp(dataset,'DatasetHCP')
                        if strcmp(procedure,'Basic')||strcmp(procedure,'GSR')||strcmp(procedure,'GSR_regr')
                            cd([Work_dir '/' dataset '/sub-' sprintf('%03d',subject) '/' session_name '/func/fourD_files/Basic']);
                        end
                        rp_name='rp_afunctional_disc_4D.txt';
                        rp = load(deblank(rp_name));
                        
                        rp=rp(:,1:6);
                        radius=50;
                        
                        %Change rotation in radials to mm (arch=radial*radius)
                        rp(:,4:6)=rp(:,4:6)*radius/180*pi; %from degrees -> radial
                        
                        %Calculate change in position between current and previous volume
                        rp_diff=diff(rp);
                        
                        %add zero as first element (see Power et al., 2014)
                        rp_zero=[zeros(1,size(rp_diff,2)); rp_diff];
                        
                        %compute framewise displacement (FD)
                        %%This particular session showed one scan less, so we added
                        %%a zero to make it same size (only maximum FD was used, so
                        %%doesn't alter any results)
                        
                        FD(:,session)=sum(abs(rp_zero),2);
                        
                    end
                    
                    clear rp rp_diff rp_zero
                end
                
                save([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_summary/DCM/' procedure '/' name_ROI_def '/Full_model/QC/Diagnostics_' ntwrk_abbrev2{network} '.mat'],'N_est_par','Max_conn','Expl_var','Posterior_estimates','Length_error_bar','Lower_bound_ci','Higher_bound_ci','thresholds_VOIs','FD','Free_energy');
                clear N_est_par Max_conn Expl_var Posterior_estimates Length_error_bar Lower_bound_ci Higher_bound_ci thresholds_VOIs FD Free_energy;
                
                
            end
            
        end
    end
end

end