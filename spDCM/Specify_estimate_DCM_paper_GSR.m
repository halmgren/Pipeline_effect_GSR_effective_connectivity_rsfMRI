function Specify_estimate_DCM_paper_GSR(dataset,session_name,subject,name_ROI_def,ROI_list,ROI_list2,procedure,SPM_dir,Work_dir)

%%%%%%%%%%%%%%%%%%%%%
%specify spDCM file
%%%%%%%%%%%%%%%%%%%%%
disp('Specify and estimate full model DCM')
disp(['Dataset: ' dataset '; subject: ' num2str(subject) '; session: ' session_name]);

if strcmp(procedure,'Basic')||strcmp(procedure,'GSR')||strcmp(procedure,'GSR_regr')
    cd([Work_dir '/' dataset '/sub-' sprintf('%03d',subject) '/' session_name '/func/VOI/' procedure '/' name_ROI_def]);
    
    tmp=0;
    
    for VOI_number=1:size(ROI_list,1)
        ntwrk=ROI_list{VOI_number,1}(1:3);
        
        if VOI_number>1 && strcmp(ROI_list{VOI_number,1}(1:3),ROI_list{VOI_number-1,1}(1:3))
            continue
        else
            tmp=tmp+1;
            
            ntwrk_name{tmp}=spm_select('FPList',pwd,['^VOI_' ntwrk '.*\.mat$']);
            ntwrk_abbrev{tmp}=ROI_list{VOI_number,1}(1:3);
        end
    end
    
    cd([Work_dir '/' dataset '/sub-' sprintf('%03d',subject) '/' session_name '/func/scan_info/']);
    load('General_information.mat');
    
    cd('../fourD_files/Basic');
    [nifti_images,~]=spm_select('ExtFPList',pwd,'^swrafunctional_disc.*\.nii$',inf);
    length_scan=size(nifti_images,1);
    
    for network_number=1:length(ntwrk_name)
        clear DCM;
        number_regions=size(ntwrk_name{network_number},1);
        
        DCM.a=ones(number_regions,number_regions);
        DCM.b=zeros(number_regions,number_regions);
        DCM.c=zeros(number_regions,1);
        DCM.d=zeros(number_regions,number_regions,0);
        
        DCM.U.u=zeros(length_scan,1);
        DCM.U.name={'null'};
        
        xY=[];
        for region_number=1:number_regions
            if isspace(ntwrk_name{network_number}(region_number,end))==1
                K=load(ntwrk_name{network_number}(region_number,1:end-1),'xY');
            else
                K=load(ntwrk_name{network_number}(region_number,1:end),'xY');
            end
            xY = spm_cat_struct(xY,K.xY);
        end
        
        DCM.xY=xY;
        
        DCM.Y.dt=TR;
        DCM.Y.X0=xY(1).X0;
        DCM.Y.Q=spm_Ce(ones(1,number_regions)*length_scan);
        
        for  region_number = 1:number_regions
            DCM.Y.y(:,region_number)  = xY(region_number).u;
            DCM.Y.name{region_number} = xY(region_number).name;
        end
        
        DCM.v=length_scan;
        DCM.n=number_regions;
        
        DCM.TE=TE;  %unused in current implementation of DCM
        DCM.delays=ones(number_regions,1)*TR/2;
        
        DCM.options.nonlinear=0;
        DCM.options.two_state=0;
        DCM.options.stochastic=0;
        DCM.options.centre=0;
        DCM.options.induced=1;
        disp(num2str(number_regions))
        GCM_full{1}=DCM;
        
        %%%%%%%%%%%%%%%%%%%%
        %Estimate DCM files
        %%%%%%%%%%%%%%%%%%%%
        mkdir([Work_dir '/' dataset '/sub-' sprintf('%03d',subject) '/' session_name '/func/DCM/' procedure '/' name_ROI_def '/Full_model']);
        cd([Work_dir '/' dataset '/sub-' sprintf('%03d',subject) '/' session_name '/func/DCM/' procedure '/' name_ROI_def '/Full_model']);
        
        GCM=spm_dcm_fit_extra_output(GCM_full);
        save(['GCM_full_estim_' ntwrk_abbrev{network_number} '.mat'],'GCM');
        GCM=GCM_full;
        save(['GCM_full_raw_' ntwrk_abbrev{network_number} '.mat'],'GCM');
        
        delete([sprintf('DCM_%s',date) '.mat']);
        clear GCM GCM_full;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Extra networks consisting of a combination of previously defined regions
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cd([Work_dir '/' dataset '/sub-' sprintf('%03d',subject) '/' session_name '/func/VOI/' procedure '/' name_ROI_def]);
    
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
                ntwrk_name2{tmp}(tmp2,:)=spm_select('FPList',pwd,['^VOI_.*\_' region '_1.mat$']);
                ntwrk_abbrev2{tmp}=ROI_list2{VOI_number,1}(1:3);
            end
            
            if VOI_number>1 && ~strcmp(ROI_list2{VOI_number,1}(1:3),ROI_list2{VOI_number-1,1}(1:3))
                tmp2=1;
                tmp=tmp+1;
                ntwrk_name2{tmp}(1,:)=spm_select('FPList',pwd,['^VOI_.*\_' region '_1.mat$']);
                ntwrk_abbrev2{tmp}=ROI_list2{VOI_number,1}(1:3);
            end
            
            if VOI_number>1 && strcmp(ROI_list2{VOI_number,1}(1:3),ROI_list2{VOI_number-1,1}(1:3))
                
                tmp2=tmp2+1;
                
                %             if tmp2~=1
                if length(spm_select('FPList',pwd,['^VOI_.*\_' region '_1.mat$']))<length(ntwrk_name2{tmp}(tmp2-1,:))
                    ntwrk_name2{tmp}(tmp2,:)=[spm_select('FPList',pwd,['^VOI_.*\_' region '_1.mat$']) ' '];
                elseif length(spm_select('FPList',pwd,['^VOI_.*\_' region '_1.mat$']))>length(ntwrk_name2{tmp}(tmp2-1,:))
                    for i=1:tmp2-1
                        ntwrk_name3{tmp}(i,:)=[ntwrk_name2{tmp}(i,:) ' '];
                    end
                    ntwrk_name2{tmp}=ntwrk_name3{tmp};
                    ntwrk_name2{tmp}(tmp2,:)=spm_select('FPList',pwd,['^VOI_.*\_' region '_1.mat$']);
                else
                    ntwrk_name2{tmp}(tmp2,:)=spm_select('FPList',pwd,['^VOI_.*\_' region '_1.mat$']);
                end
                %             end
                
            end
        end
        
        for network_number=1:length(ntwrk_name2)
            clear DCM;
            number_regions=size(ntwrk_name2{network_number},1);
            
            DCM.a=ones(number_regions,number_regions);
            DCM.b=zeros(number_regions,number_regions);
            DCM.c=zeros(number_regions,1);
            DCM.d=zeros(number_regions,number_regions,0);
            
            DCM.U.u=zeros(length_scan,1);
            DCM.U.name={'null'};
            
            xY=[];
            for region_number=1:number_regions
                if isspace(ntwrk_name2{network_number}(region_number,end))==1
                    K=load(ntwrk_name2{network_number}(region_number,1:end-1),'xY');
                else
                    K=load(ntwrk_name2{network_number}(region_number,1:end),'xY');
                end
                xY = spm_cat_struct(xY,K.xY);
            end
            
            DCM.xY=xY;
            
            DCM.Y.dt=TR;
            DCM.Y.X0=xY(1).X0;
            DCM.Y.Q=spm_Ce(ones(1,number_regions)*length_scan);
            
            for  region_number = 1:number_regions
                DCM.Y.y(:,region_number)  = xY(region_number).u;
                DCM.Y.name{region_number} = xY(region_number).name;
            end
            
            DCM.v=length_scan;
            DCM.n=number_regions;
            
            DCM.TE=TE; %unused in current implementation of DCM
            DCM.delays=ones(number_regions,1)*TR/2;
            
            DCM.options.nonlinear=0;
            DCM.options.two_state=0;
            DCM.options.stochastic=0;
            DCM.options.centre=0;
            DCM.options.induced=1;
            DCM.options.nmax=12;
            disp(num2str(number_regions));
            GCM_full{1}=DCM;
            
            %%%%%%%%%%%%%%%%%%%%
            %Estimate DCM files
            %%%%%%%%%%%%%%%%%%%%
            mkdir([Work_dir '/' dataset '/sub-' sprintf('%03d',subject) '/' session_name '/func/DCM/' procedure '/' name_ROI_def '/Full_model']);
            cd([Work_dir '/' dataset '/sub-' sprintf('%03d',subject) '/' session_name '/func/DCM/' procedure '/' name_ROI_def '/Full_model']);
        
            GCM=spm_dcm_fit_extra_output(GCM_full);
            save(['GCM_full_estim_' ntwrk_abbrev2{network_number} '.mat'],'GCM');
            GCM=GCM_full;
            save(['GCM_full_raw_' ntwrk_abbrev2{network_number} '.mat'],'GCM');
            
            delete([sprintf('DCM_%s',date) '.mat']);
            clear GCM GCM_full;
        end
        
        
    end
end

end