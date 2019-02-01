function [dataset,number_subject,single_band,slice_time_seconds]=Dataset_info_paper_GSR(number_dataset)

if number_dataset==1
    dataset='DatasetKuehn';
    number_subject=8;
    
    %options slice_timing
    single_band=1;
    slice_time_seconds=1;
    
elseif number_dataset==2
    dataset='DatasetPoldrack';
    number_subject=1;
    
    %options slice_timing
    single_band=0;
    slice_time_seconds=1;
    
elseif number_dataset==3
    dataset='DatasetKirby';
    number_subject=1;
    
    %options slice_timing
    single_band=1;
    slice_time_seconds=0;
elseif number_dataset==4
    dataset='DatasetGordon';
    number_subject=10;
    
    %options slice_timing
    single_band=1;
    slice_time_seconds=1;
elseif number_dataset==5
    dataset='DatasetHCP';
    number_subject=361;
    
    %options slice_timing
    single_band=0;
    slice_time_seconds='XX';
end

end