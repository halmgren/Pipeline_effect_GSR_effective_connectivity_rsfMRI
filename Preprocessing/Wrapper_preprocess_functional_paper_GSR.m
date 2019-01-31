function Wrapper_preprocess_functional_paper_GSR(dataset,number_subject,slice_time_seconds,procedure,SPM_dir,Work_dir)

%%this wrapper function makes it possible to parallelly preprocess several sessions
if ~strcmp(dataset,'DatasetHCP')
    for subject=1:number_subject
        cd([Work_dir '/' dataset '/']);
        session_number=dir(['sub-' sprintf('%03d',subject)]);
        session_number(1:2)=[];
        session_one=session_number(1).name;
        
        parfor session=1:length(session_number)
            session_name=session_number(session).name;
            Preprocess_functional_paper_GSR(dataset,slice_time_seconds,session_name,subject,session_one,procedure,SPM_dir,Work_dir);
        end
    end
end

if strcmp(dataset,'DatasetHCP')
    parfor subject=1:number_subject
        
        Preprocess_functional_paper_GSR_HCP(dataset,slice_time_seconds,subject,procedure,SPM_dir,Work_dir);
    end
end

end