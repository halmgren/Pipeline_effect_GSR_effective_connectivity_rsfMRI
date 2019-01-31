function Move_MRI_files_paper_GSR(Work_dir)

    %%%%%%%%%%%%%%%%%
    %Day2day-dataset
    %%%%%%%%%%%%%%%%%
    dataset='DatasetKuehn';
    for subject=1:8
        cd([Work_dir '/' dataset '/sub-' sprintf('%03d',subject)]);
        session_number=dir;
        session_number(1:2)=[];
        
        for session=1:length(session_number)
            session_name=session_number(session).name;
            cd([session_name '/func']);
            mkdir('fourD_files');
            movefile('functional.nii','fourD_files/functional_4D.nii');
            cd ../..;
        end
        
    end

    %%%%%%%%%%%%%%%%%
    %myconnectome-dataset
    %%%%%%%%%%%%%%%%%
    dataset='DatasetPoldrack';
    for subject=1
        cd([Work_dir '/' dataset '/sub-' sprintf('%03d',subject)]);
        session_number=dir;
        session_number(1:2)=[];
        
        for session=1:length(session_number)
            session_name=session_number(session).name;
            cd([session_name '/func']);
            mkdir('fourD_files');
            movefile('functional.nii','fourD_files/functional_4D.nii');
            cd ../..;
        end
        
    end

    %%%%%%%%%%%%%%%%%
    %Kirby-dataset
    %%%%%%%%%%%%%%%%%
    dataset='DatasetKirby';
    for subject=1
        cd([Work_dir '/' dataset '/sub-' sprintf('%03d',subject)]);
        session_number=dir;
        session_number(1:2)=[];
        
        for session=1:length(session_number)
            session_name=session_number(session).name;
            cd([session_name '/func']);
            mkdir('fourD_files');
            movefile('functional.nii','fourD_files/functional_4D.nii');
            cd ../..;
        end
        
    end
 
    %%%%%%%%%%%%%%%%%
    %MSC-dataset
    %%%%%%%%%%%%%%%%%
    dataset='DatasetGordon';
    for subject=1:10
        cd([Work_dir '/' dataset '/sub-' sprintf('%03d',subject)]);
        session_number=dir;
        session_number(1:2)=[];
        
        for session=1:length(session_number)
            session_name=session_number(session).name;
            cd([session_name '/func']);
            mkdir('fourD_files');
            movefile('functional.nii','fourD_files/functional_4D.nii');
            cd ../..;
        end
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %HCP-dataset: already minimally preprocessed
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dataset='DatasetHCP';
    for subject=1:361
        cd([Work_dir '/' dataset '/sub-' sprintf('%03d',subject)]);
        session_number=dir;
        session_number(1:2)=[];
        
        %functional
        for session=1:length(session_number)
            session_name=session_number(session).name;
            cd([session_name '/func']);
            mkdir('fourD_files/Basic');
            movefile('functional.nii','fourD_files/Basic/wrafunctional_disc_4D.nii');
            movefile('brainmask_fs.nii','fourD_files/Basic/brainmask_fs.nii');
            movefile('Movement_Regressors.txt','fourD_files/Basic/rp_afunctional_disc_4D.txt');
            cd ../..;
        end
        
        %structural
        for session=1
            session_name=session_number(session).name;
            cd([session_name '/anat']);
            mkdir('Basic');
            movefile('structural.nii','Basic/wSkullstr_biascor_structural.nii');
            cd ../..;
        end
    end
end