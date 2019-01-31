function Preprocess_functional_paper_GSR_HCP(dataset,slice_time_seconds,subject,procedure,SPM_dir,Work_dir)

if strcmp(procedure,'Basic') %Same for all procedure, so only done once
    
    if strcmp(dataset,'DatasetHCP')
        
        cd([Work_dir '/' dataset '/']);
        session_number=dir(['sub-' sprintf('%03d',subject)]);
        session_number(1:2)=[];
        session_one=session_number(1).name;
        
        for session=1:length(session_number)
            
            session_name=session_number(session).name;
            
            disp(['Dataset: ' dataset ' subject: ' num2str(subject) ' session: ' session_name]);
            disp('Preprocess_functional');
            
            clear matlabbatch;
            cd([Work_dir '/' dataset '/']);
            
            cd(['sub-' sprintf('%03d',subject) '/' session_name '/func']);
            cd(['fourD_files/' procedure]);
            [nifti_images,~]=spm_select('ExtFPList',pwd,'^wrafunctional_disc.*\.nii$',inf);
            
            matlabbatch{1}.spm.spatial.smooth.data = cellstr(nifti_images);
            matlabbatch{1}.spm.spatial.smooth.fwhm = [6 6 6]; %see, also Taghia et al., 2018, Nature Communications
            matlabbatch{1}.spm.spatial.smooth.dtype = 0;
            matlabbatch{1}.spm.spatial.smooth.im = 0;
            matlabbatch{1}.spm.spatial.smooth.prefix = 's';
            
            spm_jobman('run',matlabbatch);
            
        end
    end
end
end


