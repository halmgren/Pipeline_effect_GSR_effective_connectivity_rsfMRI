function Wrapper_PEB_subject_mean_paper_GSR(dataset,number_subject,name_ROI_def,ROI_list,ROI_list2,procedure,all_procedure_names,SPM_dir,Work_dir)

parfor subject=1:number_subject
    cd([Work_dir '/' dataset '/']);
    session_number=dir(['sub-' sprintf('%03d',subject)]);
    session_number(1:2)=[];
    
    PEB_subject_mean_paper_GSR(dataset,subject,name_ROI_def,ROI_list,ROI_list2,procedure,all_procedure_names,SPM_dir,Work_dir);

end

end