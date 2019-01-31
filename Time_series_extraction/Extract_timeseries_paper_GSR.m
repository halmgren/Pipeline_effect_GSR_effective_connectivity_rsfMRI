function Extract_timeseries_paper_GSR(dataset,session_name,subject,name_ROI_def,ROI_list,procedure,SPM_dir,Work_dir)

disp(['dataset: ' dataset ' subject: ' num2str(subject) ' session: ' session_name]);
disp('Extract_timeseries');
clear matlabbatch;

%Extract information regarding TR
cd([Work_dir '/' dataset '/']);
cd(['sub-' sprintf('%03d',subject) '/' session_name '/func']);
cd('scan_info');
load('General_information.mat');

%Extract information regarding number of scans
if strcmp(procedure,'Basic')||strcmp(procedure,'GSR')||strcmp(procedure,'GSR_regr')
    cd('../fourD_files/Basic');
    [nifti_images,~]=spm_select('ExtFPList',pwd,'^swrafunctional_disc.*\.nii$',inf);
    length_scan=size(nifti_images,1);
end

cd ../..;
if 7~=exist('SPM','dir')
    mkdir('SPM');
end

cd SPM;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Estimation of SPM with DCT regressor
%%If SPM already exists, then we shouldn't estimate it again, since it is the same for all ROIs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(procedure,'Basic')||strcmp(procedure,'GSR')||strcmp(procedure,'GSR_regr')
    if 7~=exist(procedure,'dir')
        clear matlabbatch;
        
        load(['../regressors/Basic/dct_set.mat']);
        
        matlabbatch{1}.spm.stats.fmri_spec.dir = {[Work_dir '/' dataset '/sub-' sprintf('%03d',subject) '/' session_name '/func/SPM']};
        matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'scans';
        matlabbatch{1}.spm.stats.fmri_spec.timing.RT = TR;
        matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = image_dim(3);
        
        if strcmp(dataset,'DatasetKuehn')
            matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 18;
        elseif strcmp(dataset,'DatasetPoldrack')
            if image_dim(3)==68
                matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 9;
            elseif image_dim(3)==64
                matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
            end
        elseif strcmp(dataset,'DatasetKirby')
            matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 19;
        elseif strcmp(dataset,'DatasetGordon')
            matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 18;
        end
        
        %%fMRI model specification
        matlabbatch{1}.spm.stats.fmri_spec.sess.scans = cellstr(nifti_images);
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
        matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {''};
        matlabbatch{1}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});
        if strcmp(procedure,'Basic')||strcmp(procedure,'GSR')
            matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {
                [Work_dir '/' dataset '/sub-' sprintf('%03d',subject) '/' session_name '/func/regressors/Basic/dct_set.mat']
                [Work_dir '/' dataset '/sub-' sprintf('%03d',subject) '/' session_name '/func/regressors/Basic/VOI_CSF_signal_1.mat']
                [Work_dir '/' dataset '/sub-' sprintf('%03d',subject) '/' session_name '/func/regressors/Basic/VOI_WM_signal_1.mat']
                [Work_dir '/' dataset '/sub-' sprintf('%03d',subject) '/' session_name '/func/fourD_files/Basic/rp_afunctional_disc_4D.txt']
                };
        end
        
        
        if strcmp(procedure,'GSR_regr')
            matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {
                [Work_dir '/' dataset '/sub-' sprintf('%03d',subject) '/' session_name '/func/regressors/Basic/dct_set.mat']
                [Work_dir '/' dataset '/sub-' sprintf('%03d',subject) '/' session_name '/func/regressors/Basic/VOI_CSF_signal_1.mat']
                [Work_dir '/' dataset '/sub-' sprintf('%03d',subject) '/' session_name '/func/regressors/Basic/VOI_WM_signal_1.mat']
                [Work_dir '/' dataset '/sub-' sprintf('%03d',subject) '/' session_name '/func/fourD_files/Basic/rp_afunctional_disc_4D.txt']
                [Work_dir '/' dataset '/sub-' sprintf('%03d',subject) '/' session_name '/func/regressors/GSR_regr/global_signal.mat']
                };
        end
        
        matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = 128;
        matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
        matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
        matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
        
        if strcmp(procedure,'Basic')||strcmp(procedure,'GSR_regr')
            matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
        elseif strcmp(procedure,'GSR')
            matlabbatch{1}.spm.stats.fmri_spec.global = 'Scaling';
        end
        
        matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;
        matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
        matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
        
        %%fMRI model estimation
        matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
        matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
        matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
        
        %%Define contrasts
        matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
        matlabbatch{3}.spm.stats.con.consess{1}.fcon.name = 'Effects_of_interest';
        matlabbatch{3}.spm.stats.con.consess{1}.fcon.weights = eye(size(R,2));
        matlabbatch{3}.spm.stats.con.consess{1}.fcon.sessrep = 'none';
        matlabbatch{3}.spm.stats.con.delete = 0;
        
        spm_jobman('run',matlabbatch);
        
        mkdir(procedure);
        
        folder_content=dir;
        
        for beta_image=1:(length(folder_content(not([folder_content.isdir])))-6)
            beta_file=sprintf('beta_%04d.nii',beta_image);
            movefile(beta_file,procedure);
        end
        
        movefile('ess_0001.nii',procedure);
        movefile('mask.nii',procedure);
        movefile('ResMS.nii',procedure);
        movefile('RPV.nii',procedure);
        movefile('SPM.mat',procedure);
        movefile('spmF_0001.nii',procedure);
        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Timeseries Extraction for each ROI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mkdir([Work_dir '/' dataset '/sub-' sprintf('%03d',subject) '/' session_name '/func/VOI/' procedure '/' name_ROI_def]);

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

for ntwrk_number=1:length(ntwrk_name)
    if ntwrk_number==1
        cumul_ntwrk_size(ntwrk_number)=0;
    else
        cumul_ntwrk_size(ntwrk_number)=sum(ntwrk_size(1:ntwrk_number-1));
    end
end

for ntwrk_number=1:length(ntwrk_name)
    clear thresholds;
    thresholds=zeros(1,ntwrk_size(ntwrk_number));

    if strcmp(procedure,'Basic')||strcmp(procedure,'GSR')||strcmp(procedure,'GSR_regr')
        
        if strcmp(dataset,'DatasetPoldrack')||strcmp(dataset,'DatasetGordon')
            clear matlabbatch;
            for VOI_number=1:ntwrk_size(ntwrk_number)
            
                matlabbatch{1}.spm.util.voi.spmmat(1) = {[Work_dir '/' dataset '/sub-' sprintf('%03d',subject) '/' session_name '/func/SPM/' procedure '/SPM.mat']};
                matlabbatch{1}.spm.util.voi.adjust = 1;
                matlabbatch{1}.spm.util.voi.session = 1;
                matlabbatch{1}.spm.util.voi.name = ROI_list{cumul_ntwrk_size(ntwrk_number)+VOI_number,1};
                matlabbatch{1}.spm.util.voi.roi{1}.spm.spmmat(1) = {[Work_dir '/' dataset '/sub-' sprintf('%03d',subject) '/' session_name '/func/SPM/' procedure '/SPM.mat']};
                matlabbatch{1}.spm.util.voi.roi{1}.spm.contrast = 1;
                matlabbatch{1}.spm.util.voi.roi{1}.spm.conjunction = 1;
                matlabbatch{1}.spm.util.voi.roi{1}.spm.threshdesc = 'none';
                matlabbatch{1}.spm.util.voi.roi{1}.spm.thresh = 0.001;
                matlabbatch{1}.spm.util.voi.roi{1}.spm.extent = 0;
                matlabbatch{1}.spm.util.voi.roi{1}.spm.mask = struct('contrast', {}, 'thresh', {}, 'mtype', {});
                matlabbatch{1}.spm.util.voi.roi{2}.sphere.centre = ROI_list{cumul_ntwrk_size(ntwrk_number)+VOI_number,2};
                matlabbatch{1}.spm.util.voi.roi{2}.sphere.radius = 10;
                matlabbatch{1}.spm.util.voi.roi{2}.sphere.move.fixed = 1;
                matlabbatch{1}.spm.util.voi.roi{3}.sphere.centre = ROI_list{cumul_ntwrk_size(ntwrk_number)+VOI_number,2};
                matlabbatch{1}.spm.util.voi.roi{3}.sphere.radius = 8;
                matlabbatch{1}.spm.util.voi.roi{3}.sphere.move.global.spm = 1;
                matlabbatch{1}.spm.util.voi.roi{3}.sphere.move.global.mask = 'i2';
                matlabbatch{1}.spm.util.voi.expression = 'i1&i2&i3';
                
                count=0;
                err_count=0;
                while count==err_count
                    try
                        if count==0
                            matlabbatch{1}.spm.util.voi.roi{1}.spm.thresh = matlabbatch{1}.spm.util.voi.roi{1}.spm.thresh;
                        elseif count==1
                            matlabbatch{1}.spm.util.voi.roi{1}.spm.thresh = 0.01;
                        elseif count==2
                            matlabbatch{1}.spm.util.voi.roi{1}.spm.thresh = 0.05;
                        elseif count>=3
                            matlabbatch{1}.spm.util.voi.roi{1}.spm.thresh = 0.05+0.05*(count-2); %These sessions are excluded
                        end
                        spm_jobman('run',matlabbatch);
                    catch
                        err_count=err_count+1;
                    end
                    if count==22
                        break
                    end
                    count=count+1;
                    
                end
                
                disp(matlabbatch{1}.spm.util.voi.roi{1}.spm.thresh)
                
                thresholds(VOI_number)=matlabbatch{1}.spm.util.voi.roi{1}.spm.thresh;
                
                if err_count~=0
                    disp(['subject: ' num2str(subject)]);
                    disp(['session: ' session_name]);
                    
                    
                end
            end
        end
        
        %Kirby and Kuhn dataset contain less scans (150-200), therefore threshold
        %was set lower (0.05), and was increased with 0.05 if no significant voxels
        %appeared
        if strcmp(dataset,'DatasetKuehn')||strcmp(dataset,'DatasetKirby')
            clear matlabbatch;
            
            for VOI_number=1:ntwrk_size(ntwrk_number)
            
                matlabbatch{1}.spm.util.voi.spmmat(1) = {[Work_dir '/' dataset '/sub-' sprintf('%03d',subject) '/' session_name '/func/SPM/' procedure '/SPM.mat']};
                matlabbatch{1}.spm.util.voi.adjust = 1;
                matlabbatch{1}.spm.util.voi.session = 1;
                matlabbatch{1}.spm.util.voi.name = ROI_list{cumul_ntwrk_size(ntwrk_number)+VOI_number,1};
                matlabbatch{1}.spm.util.voi.roi{1}.spm.spmmat(1) = {[Work_dir '/' dataset '/sub-' sprintf('%03d',subject) '/' session_name '/func/SPM/' procedure '/SPM.mat']};
                matlabbatch{1}.spm.util.voi.roi{1}.spm.contrast = 1;
                matlabbatch{1}.spm.util.voi.roi{1}.spm.conjunction = 1;
                matlabbatch{1}.spm.util.voi.roi{1}.spm.threshdesc = 'none';
                matlabbatch{1}.spm.util.voi.roi{1}.spm.thresh = 0;
                matlabbatch{1}.spm.util.voi.roi{1}.spm.extent = 0;
                matlabbatch{1}.spm.util.voi.roi{1}.spm.mask = struct('contrast', {}, 'thresh', {}, 'mtype', {});
                matlabbatch{1}.spm.util.voi.roi{2}.sphere.centre = ROI_list{cumul_ntwrk_size(ntwrk_number)+VOI_number,2};
                matlabbatch{1}.spm.util.voi.roi{2}.sphere.radius = 10;
                matlabbatch{1}.spm.util.voi.roi{2}.sphere.move.fixed = 1;
                matlabbatch{1}.spm.util.voi.roi{3}.sphere.centre = ROI_list{cumul_ntwrk_size(ntwrk_number)+VOI_number,2};
                matlabbatch{1}.spm.util.voi.roi{3}.sphere.radius = 8;
                matlabbatch{1}.spm.util.voi.roi{3}.sphere.move.global.spm = 1;
                matlabbatch{1}.spm.util.voi.roi{3}.sphere.move.global.mask = 'i2';
                matlabbatch{1}.spm.util.voi.expression = 'i1&i2&i3';
                
                count=0;
                err_count=0;
                while count==err_count
                    try
                        matlabbatch{1}.spm.util.voi.roi{1}.spm.thresh = matlabbatch{1}.spm.util.voi.roi{1}.spm.thresh+0.05; %if treshold > 0.05, session will be excluded
                        spm_jobman('run',matlabbatch);
                    catch
                        err_count=err_count+1;
                    end
                    count=count+1; %if treshold > 0.05, session will be excluded
                end
                
                disp(err_count*0.05+0.05)
                
                thresholds(VOI_number)=matlabbatch{1}.spm.util.voi.roi{1}.spm.thresh;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Move time-series to new folder
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for VOI_number=1:ntwrk_size(ntwrk_number)
            
            movefile([Work_dir '/' dataset '/sub-' sprintf('%03d',subject) '/' session_name '/func/SPM/' procedure '/VOI_' ROI_list{cumul_ntwrk_size(ntwrk_number)+VOI_number,1} '_1.mat'],[Work_dir '/' dataset '/sub-' sprintf('%03d',subject) '/' session_name '/func/VOI/' procedure '/' name_ROI_def])
            movefile([Work_dir '/' dataset '/sub-' sprintf('%03d',subject) '/' session_name '/func/SPM/' procedure '/VOI_' ROI_list{cumul_ntwrk_size(ntwrk_number)+VOI_number,1} '_eigen.nii'],[Work_dir '/' dataset '/sub-' sprintf('%03d',subject) '/' session_name '/func/VOI/' procedure '/' name_ROI_def])
            movefile([Work_dir '/' dataset '/sub-' sprintf('%03d',subject) '/' session_name '/func/SPM/' procedure '/VOI_' ROI_list{cumul_ntwrk_size(ntwrk_number)+VOI_number,1} '_mask.nii'],[Work_dir '/' dataset '/sub-' sprintf('%03d',subject) '/' session_name '/func/VOI/' procedure '/' name_ROI_def])
            
        end

        cd([Work_dir '/' dataset '/sub-' sprintf('%03d',subject) '/' session_name '/func/VOI/' procedure '/' name_ROI_def]);
        
        save(['FL_thresholds_' ntwrk_name{ntwrk_number} '.mat'],'thresholds');
        
    end
end

end
