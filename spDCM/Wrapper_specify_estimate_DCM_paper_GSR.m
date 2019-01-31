function Wrapper_specify_estimate_DCM_paper_GSR(dataset,number_subject,name_ROI_def,ROI_list,ROI_list2,procedure,SPM_dir,Work_dir)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%this wrapper function makes it possible to parallelly specify and estimate several sessions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%(1) longitudinal datasets
%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~strcmp(dataset,'DatasetHCP')
    for subject=1:number_subject
        
        %%%%%%%%%%%%%%%%%%%%%%%%%
        %Collect session names
        %%%%%%%%%%%%%%%%%%%%%%%%%
        cd([Work_dir '/' dataset '/']);
        session_number=dir(['sub-' sprintf('%03d',subject)]);
        session_number(1:2)=[];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %invert all sessions seperately
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        parfor session=1:length(session_number)
            session_name=session_number(session).name;
            Specify_estimate_DCM_paper_GSR(dataset,session_name,subject,name_ROI_def,ROI_list,ROI_list2,procedure,SPM_dir,Work_dir);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%put all sessions for small networks in a separate file
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tmp=0;
        
        for VOI_number=1:size(ROI_list,1)
            ntwrk=ROI_list{VOI_number,1}(1:3);
            
            if VOI_number>1 && strcmp(ROI_list{VOI_number,1}(1:3),ROI_list{VOI_number-1,1}(1:3))
                continue
                
            else
                tmp=tmp+1;
                ntwrk_name{tmp}=ROI_list{VOI_number,1}(1:3);
            end
        end
        
        if strcmp(procedure,'Basic')||strcmp(procedure,'GSR')||strcmp(procedure,'GSR_regr')
            %put all session of each network in a separate file
            mkdir([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_summary/DCM/' procedure '/' name_ROI_def '/Full_model']);
            
            for network=1:length(ntwrk_name)
                
                for session=1:length(session_number)
                    session_name=session_number(session).name;
                    cd([Work_dir '/' dataset '/sub-' sprintf('%03d',subject) '/' session_name '/func/DCM/'  procedure '/' name_ROI_def '/Full_model']);
                    
                    load(['GCM_full_estim_' ntwrk_name{network} '.mat']);
                    
                    GCM_full{session}=GCM{1};
                end
                
                %Make network-specific files
                GCM=GCM_full';
                cd([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_summary/DCM/'  procedure '/' name_ROI_def '/Full_model']);
                save(['GCM_' ntwrk_name{network} '_full_estim.mat'],'GCM');
                
                clear GCM GCM_full;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Combined network
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if  exist('ROI_list2') 
            
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
                %put all session of each network in a separate file & Bayesian parameter averaging (BPA)
                mkdir([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_summary/DCM/' procedure '/' name_ROI_def '/Full_model']);
                
                for network=1:length(ntwrk_abbrev2)
                    
                    for session=1:length(session_number)
                        session_name=session_number(session).name;
                        cd([Work_dir '/' dataset '/sub-' sprintf('%03d',subject) '/' session_name '/func/DCM/'  procedure '/' name_ROI_def '/Full_model']);
                        
                        load(['GCM_full_estim_' ntwrk_abbrev2{network} '.mat']);
                        
                        GCM_full{session}=GCM{1};
                    end
                    
                    %Make network-specific files
                    GCM=GCM_full';
                    cd([Work_dir '/' dataset '/sub-' sprintf('%03d', subject) '_summary/DCM/'  procedure '/' name_ROI_def '/Full_model']);
                    save(['GCM_' ntwrk_abbrev2{network} '_full_estim.mat'],'GCM');
                    
                    clear GCM GCM_full;
                end
            end
        end
        
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%(2) cross-sectional datasets
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(dataset,'DatasetHCP')
    parfor subject=1:number_subject
        
        Specify_estimate_DCM_paper_GSR_HCP(dataset,subject,name_ROI_def,ROI_list,ROI_list2,procedure,SPM_dir,Work_dir);
        
    end
    
end

end