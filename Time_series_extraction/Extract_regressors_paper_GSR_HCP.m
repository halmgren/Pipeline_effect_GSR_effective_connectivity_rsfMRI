function Extract_regressors_paper_GSR_HCP(dataset,subject,procedure,SPM_dir,Work_dir)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Build regressors for DCT, WM, CSF, and movement
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Construct discrete cosine set (code adapted from Peter Zeidman (shared via spm archives))

cd([Work_dir '/' dataset '/']);
session_number=dir(['sub-' sprintf('%03d',subject)]);
session_number(1:2)=[];

for session=1:length(session_number)
    session_name=session_number(session).name;
    
    disp(['subject: ' num2str(subject) ' session: ' session_name]);
    disp('Extract_regressors');
    
    if strcmp(procedure,'Basic') %Same for all procedures, so only done once
        
        clear matlabbatch;
        cd(['sub-' sprintf('%03d',subject) '/' session_name '/func']);
        
        %Extract information regarding TR
        cd('scan_info');
        load('General_information.mat');
        
        %Extract information regarding number of scans
        cd(['../fourD_files/' procedure]);
        [nifti_images,~]=spm_select('ExtFPList',pwd,'^swrafunctional_disc.*\.nii$',inf);
        length_scan=size(nifti_images,1);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Construction of DCT regressors (part adapted from code Peter Zeidman on spm archives)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        dct_set = spm_dctmtx(length_scan,length_scan);
    
        upper_limit       = 0.1;                
        lower_limit       = 1/128;
        
        lower_limit_comp   = fix(2*(length_scan*TR)*lower_limit+1);
        
        upper_limit_comp   = fix(2*(length_scan*TR)*upper_limit+1);
        dct_set(:,[1:lower_limit_comp upper_limit_comp:end]) = [];
        
        R=dct_set;
        
        mkdir('../../regressors');
        cd('../../regressors');
        
        save('dct_set.mat','R');
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %specify and estimate SPM to get compute CSF and WM signal!
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        clear matlabbatch;
        
        matlabbatch{1}.spm.stats.fmri_spec.dir = {[Work_dir '/' dataset '/sub-' sprintf('%03d',subject) '/' session_name '/func/regressors']};
        matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'scans';
        matlabbatch{1}.spm.stats.fmri_spec.timing.RT = TR;
        
        if strcmp(dataset,'DatasetHCP')    %no slice timing correction, so no need to change this
            matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
            matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
        end
        
        matlabbatch{1}.spm.stats.fmri_spec.sess.scans = cellstr(nifti_images);
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
        matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {''};
        matlabbatch{1}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});
        matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {''}; 
        matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = 128;
        matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
        matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
        matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
        matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
        matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0;
        matlabbatch{1}.spm.stats.fmri_spec.mask = {[Work_dir '/' dataset '/sub-' sprintf('%03d',subject) '/' session_name '/func/fourD_files/Basic/brainmask_fs.nii,1']};
        matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
        
        matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
        matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
        matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
        
        matlabbatch{3}.spm.util.voi.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
        matlabbatch{3}.spm.util.voi.adjust = NaN;
        matlabbatch{3}.spm.util.voi.session = 1;
        matlabbatch{3}.spm.util.voi.name = 'WM_signal';
        matlabbatch{3}.spm.util.voi.roi{1}.mask.image(1) = cfg_dep('Model estimation: Analysis Mask', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','mask'));
        matlabbatch{3}.spm.util.voi.roi{1}.mask.threshold = 0.5;
        matlabbatch{3}.spm.util.voi.roi{2}.sphere.centre = [0 -24 -33];
        matlabbatch{3}.spm.util.voi.roi{2}.sphere.radius = 7;
        matlabbatch{3}.spm.util.voi.roi{2}.sphere.move.fixed = 1;
        
        matlabbatch{3}.spm.util.voi.expression = 'i1&i2';
        
        matlabbatch{4}.spm.util.voi.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
        matlabbatch{4}.spm.util.voi.adjust = NaN;
        matlabbatch{4}.spm.util.voi.session = 1;
        matlabbatch{4}.spm.util.voi.name = 'CSF_signal';
        matlabbatch{4}.spm.util.voi.roi{1}.mask.image(1) = cfg_dep('Model estimation: Analysis Mask', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','mask'));
        matlabbatch{4}.spm.util.voi.roi{1}.mask.threshold = 0.5;
        matlabbatch{4}.spm.util.voi.roi{2}.sphere.centre = [0 -40 -5];
        matlabbatch{4}.spm.util.voi.roi{2}.sphere.radius = 5;
        matlabbatch{4}.spm.util.voi.roi{2}.sphere.move.fixed = 1;
        
        
        matlabbatch{4}.spm.util.voi.expression = 'i1&i2';
        spm_jobman('run',matlabbatch);
        
        mkdir(procedure);
        
        movefile('beta_0001.nii',procedure);
        movefile('dct_set.mat',procedure);
        movefile('mask.nii',procedure);
        movefile('ResMS.nii',procedure);
        movefile('RPV.nii',procedure);
        movefile('SPM.mat',procedure);
        movefile('VOI_CSF_signal_1.mat',procedure);
        movefile('VOI_CSF_signal_eigen.nii',procedure);
        movefile('VOI_CSF_signal_mask.nii',procedure);
        movefile('VOI_WM_signal_1.mat',procedure);
        movefile('VOI_WM_signal_eigen.nii',procedure);
        movefile('VOI_WM_signal_mask.nii',procedure);
    end
    
    if strcmp(procedure,'GSR_regr')
        
        %load scans
        cd([Work_dir '/' dataset '/sub-' sprintf('%03d',subject) '/' session_name '/func/fourD_files/Basic']);
        
        [nifti_images,~]=spm_select('ExtFPList',pwd,'^swrafunctional_disc.*\.nii$',inf);
        
        %%%%%%%%%%%%%%%%%%%%%%%%
        %Compute global signal
        %%%%%%%%%%%%%%%%%%%%%%%%
        VY=spm_data_hdr_read(nifti_images);
        q     = length(VY);
        R     = zeros(q,1);
        
        for i = 1:q
            R(i) = spm_global(VY(i)); %calculate global signal
        end
        
        mkdir(['../../regressors/' procedure]);
        cd(['../../regressors/' procedure]);
        
        save(['global_signal.mat'],'R');
        
    end
end
end
