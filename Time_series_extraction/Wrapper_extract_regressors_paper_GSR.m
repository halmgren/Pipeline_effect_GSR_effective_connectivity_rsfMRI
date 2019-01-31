function Wrapper_extract_regressors_paper_GSR(dataset,number_subject,procedure,SPM_dir,Work_dir)

%%this wrapper function makes it possible to  extract regressors for several sessions parallelly
if ~strcmp(dataset,'DatasetHCP')
    for subject=1:number_subject
        cd([Work_dir '/' dataset '/']);
        session_number=dir(['sub-' sprintf('%03d',subject)]);
        session_number(1:2)=[];
        
        parfor session=1:length(session_number)
            session_name=session_number(session).name;
            Extract_regressors_paper_GSR(dataset,session_name,subject,procedure,SPM_dir,Work_dir);
        end
        
    end
else
    parfor subject=1:number_subject
        
        Extract_regressors_paper_GSR_HCP(dataset,subject,procedure,SPM_dir,Work_dir);
        
    end
end
end