function Wrapper_extract_timeseries_paper_GSR(dataset,number_subject,name_ROI_def,ROI_list,procedure,SPM_dir,Work_dir)

%%this wrapper function makes it possible to parallelly extract VOIs for several sessions
if ~strcmp(dataset,'DatasetHCP')
    for subject=1:number_subject
        cd([Work_dir '/' dataset '/']);
        session_number=dir(['sub-' sprintf('%03d',subject)]);
        session_number(1:2)=[];
        
        parfor session=1:length(session_number)
            session_name=session_number(session).name;
            Extract_timeseries_paper_GSR(dataset,session_name,subject,name_ROI_def,ROI_list,procedure,SPM_dir,Work_dir);
        end
    end
else
    
    parfor subject=1:number_subject
        
        Extract_timeseries_paper_GSR_HCP(dataset,subject,name_ROI_def,ROI_list,procedure,SPM_dir,Work_dir);
        
    end
    
end

end