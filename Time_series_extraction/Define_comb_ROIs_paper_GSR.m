function [ROI_list2]=Define_comb_ROIs_paper_GSR(name_ROI_def)

%naming is very important: first three letters should correspond to
%network, further letters correspond to region
%regions should contain three or four letter, no more, no less
%!!!Important: region names should be in the same order as in original 'ROI_list'!!!

if strcmp(name_ROI_def,'Smith')
    
    ROI_list2={
        'TG2_ACC',[0 36 22];...    
        'TG2_lINS',[-36 12 -2];...
        'TG2_lMFG',[-28 54 14];...
        'TG2_rINS',[32 20 4];...
        'TG2_rMFG',[30 52 14];...
        'TG2_SMA',[0 -12 50];...       %%Paracentral lobe
        'TG2_lSMC',[-38 -26 56];...
        'TG2_rSMC',[44 -16 48];...
        'TG2_PRC',[2,-58,30];...      %%Precuneus/posteriorPCC
        'TG2_lIPC',[-44,-60,24];...    %inferior parietal
        'TG2_rIPC',[54,-62,28];...
        }; 
end

end