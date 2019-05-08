%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Code for the paper 'The effect of global signal regression on DCM estimates of noise and effective connectivity from resting state fMRI'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SOFTWARE VERSIONS
%
%OS: ubuntu 16.04 LTS; Intel® Core™ i7-6800K CPU @ 3.40GHz × 12
%MATLAB: 8.6.0 (R2015b)
%SPM: SPM12 (revision 6906)
%DCM: DCM12 (revision 6801)
%PEB (revision 6778)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DATASETS
%
%'MyConnectome',            Downloaded from https://openfmri.org/dataset/ds000031/
%'The Midnight Scan Club'   Downloaded from https://openfmri.org/dataset/ds000224/
%'Kirby'                    Downloaded from https://www.nitrc.org/projects/kirbyweekly
%'Day2day'                  For availability, see Filevich et al. (2017); section 'Availability of data and materials'
%'HCP'                      Downloadable from https://www.humanconnectome.org (900 subjects release); for more information see https://www.humanconnectome.org/study/hcp-young-adult/document/900-subjects-data-release
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DATASET PREPARATION
%
%Because different datasets were structured in different formats (which might change in future 
%releases), they should all be (re)structured according to the BIDS format (http://bids.neuroimaging.io/).
%In addition, 4D resting state files should be renamed as 'functional.nii' and T1-weighted images should be 
%renamed as 'structural.nii'. In future releases we will try to encorporate some of these (re)structuring steps in
%the code.
%
%The Longitudinal datasets should thus be structured as follows:
%
%Work_dir
%   my_dataset1
%       sub-001          (for subject two digits are required, e.g., sub-05 or sub-10)
%           ses-001     (for session three digits are required, e.g., ses-005 or ses-015)
%               anat
%                   structural.nii
%               func
%                   functional.nii
%           ses-002
%               func
%                   functional.nii
%           ...
%
%       sub-002
%           ses-001
%               anat
%                   structural.nii
%               func
%                   functional.nii
%           ...
%       ...
%
%   my_dataset2
%       sub-001
%           ses-001
%               anat
%                   structural.nii
%               func
%                   functional.nii
%           ...
%
%!!!!!!!!!!!!!!!!!!!!!!!!!!
%THE NAMES OF THE DATASET-FOLDERS ARE FIXED:
%'Myconnectome':            change 'my_dataset' to 'DatasetPoldrack' 
%'The Midnight Scan Club':  change 'my_dataset' to 'DatasetGordon' 
%'Kirby':                   change 'my_dataset' to 'DatasetKirby' 
%'Day2day':                 change 'my_dataset' to 'DatasetKuehn'
%'HCP':                     change 'my_dataset' to 'DatasetHCP'
%!!!!!!!!!!!!!!!!!!!!!!!!!
%
%Anatomical (T1) images for longitudinal data analyses that were used, were stored in the 'anat' folder of the FIRST session that was analyzed
%
%               - 'MyConnectome':         T1 image of session 12 was used and should be stored in 'anat' folder of session 13 (first fMRI session analyzed)
%               - 'Midnight Scan Club':   First T1 image acquired (for each subject) was used and should be stored in 'anat' folder of session 1
%               - 'Day2day' and 'Kirby':  T1 image was acquired on same day as first fMRI scan
%
%For human connectome project (900 subject release) data analyses, the following images were used: 
%                           Preprocessed anatomical (T1) images (originally called and stored at '(subject_code)/MNINonLinear/T1w_restore_brain.nii.gz')
%                               !!name should be changed to 'structural.nii', see tree structure above
%                           Skull-stripped brain mask (originally called and stored at '(subject_code)/MNINonLinear/brainmask_fs.nii.gz'
%                               !!each subjects' brainmask should be gunzipped AND added to the respective func-directory (using the original name, i.e., brainmask_fs.nii)
%                           Minimally preprocessed LR functional images of the first session were used (originally calledand stored at: ('(subject_code)/MNINonLinear/Results/rfMRI_REST1_LR/rfMRI_REST1_LR.nii.gz'))
%                               !!name should be changed to 'functional.nii', see tree structure above
%                           Movement regressors (originally called '(subject_code)/MNINonLinear/Results/rfMRI_REST1_LR/Movement_Regressors.txt')
%                               !!each subjects' movement regressor file should be added to the respective func-directory with original name (i.e., Movement_Regressors.txt'))
%
%The directory of HCP data should look like this:
%Work dir
%   DatasetHCP
%       sub-001
%           ses-001
%               anat
%                   structural.nii
%               func
%                   functional.nii
%                   brainmask_fs.nii
%                   Movement_Regressors.txt
%                           
%
%All neuroimaging data should be in 4D nifti format AND should be unpacked/extracted
%
%The number of iterations in spm_dcm_peb.m was increased to 128 (changed 'for n = 1:64' to 'for n = 1:128')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%We encountered two duplicates in the 'Kirby' dataset. These are deleted from the present release, so if 
%you have downloaded the dataset before 2018-01-02 you should redownload it before using this code.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%For comments and feedback, please contact Hannes_Almgren@ugent.be
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Specify SPM directory and working directory for datasets (do NOT include '/' at the end)
SPM_dir='/home/hannes/Desktop/spm12';
Work_dir='/media/hannes/Almgren_Disk4';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%From here on everything is automatic
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%2.1 Datasets and subjects
%%%%%%%%%%%%%%%%%%%%%%%%%%%
Produce_list_HCP_subjects_codes_paper_GSR(Work_dir); %Save list of (unrelated) subjects' codes used in the present study

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%2.2.1 Preprocessing (part 1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Move_artefact_and_pilot_scans_paper_GSR(Work_dir);

Scanning_parameters_paper_GSR(Work_dir); %Create files with information regarding scan parameters (TR etc); most of them are in fact not used but calculated from data or from header

Move_MRI_files_paper_GSR(Work_dir); %Move functional images to separate folder

Remove_first_scans_paper_GSR(Work_dir); %Discard first five scans of each session; !NOT DONE FOR HCP DATA BECAUSE ENOUGH DATAPOINTS!

tic

for number_dataset=[1:5]
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Extract information about datasets
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [dataset,number_subject,single_band,slice_time_seconds]=Dataset_info_paper_GSR(number_dataset);
    
    all_ROI_defs={'Smith'}; %ROI definition based on Smith et al., 2009
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Important, the script allows for two types of correction for the global signal.
    %if 'GSR_regr' is specified, than regular GS regression is performed (see also paper)
    %if 'GSR' is specified, than GS normalization is performed (default in SPM)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    all_procedure_names={'Basic','GSR_regr'};                %all_procedure_names={'Basic','GSR','GSR_regr'};

    for number_ROI_def=1:length(all_ROI_defs)
        
        name_ROI_def=all_ROI_defs{number_ROI_def};
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Give regions name and coördinates
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %DMN, salience, and somatomotor network were studied for this paper
        [ROI_list]=Define_ROIs_paper_GSR(name_ROI_def);
        
        [ROI_list2]=Define_comb_ROIs_paper_GSR(name_ROI_def);
        
        for number_procedure=1:length(all_procedure_names)
            
            procedure=all_procedure_names{number_procedure};
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %2.2.1 Preprocessing (part 2)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %structural
            Preprocess_structural_paper_GSR(dataset,number_subject,procedure,SPM_dir,Work_dir);
            
            %Functional
            Wrapper_preprocess_functional_paper_GSR(dataset,number_subject,slice_time_seconds,procedure,SPM_dir,Work_dir);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %2.2.2 Time-series extraction & 2.2.5 GSR
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Extract regressors for CSV, WM, and discrete cosine set
            Wrapper_extract_regressors_paper_GSR(dataset,number_subject,procedure,SPM_dir,Work_dir);
            
            %Extract timeseries for all regions
            Wrapper_extract_timeseries_paper_GSR(dataset,number_subject,name_ROI_def,ROI_list,procedure,SPM_dir,Work_dir);
                        
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %2.2.3 (Spectral) DCM for resting state fMRI (part 1)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%Estimation of DCMs
            Wrapper_specify_estimate_DCM_paper_GSR(dataset,number_subject,name_ROI_def,ROI_list,ROI_list2,procedure,SPM_dir,Work_dir)
        
            %%Save diagnostics such as number of estimable parameters, explained variance, max connection strength, excessive motion, alpha level ROI extraction, etc...
            Diagnostics_paper_GSR(dataset,number_subject,ROI_list,ROI_list2,name_ROI_def,procedure,SPM_dir,Work_dir)
            
            %Thresholds to exclude DCMs
            [threshold_expl_var,threshold_max_conn,threshold_n_par_est,threshold_FD,threshold_FD_HCP,threshold_threshold_VOIs]=Define_QC_tresholds_paper_GSR(procedure);
            
            %Compute and save whether session reach specific threshold
            Above_threshold_sessions_paper_GSR(dataset,number_subject,ROI_list,ROI_list2,name_ROI_def,threshold_expl_var,threshold_max_conn,threshold_n_par_est,threshold_FD,threshold_FD_HCP,threshold_threshold_VOIs,procedure,SPM_dir,Work_dir)
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %2.2.4 Parametric empirical Bayes (part 1)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if ~strcmp(dataset,'DatasetHCP')
               Wrapper_PEB_subject_mean_paper_GSR(dataset,number_subject,name_ROI_def,ROI_list,ROI_list2,procedure,all_procedure_names,SPM_dir,Work_dir) %PEB estimation over sessions
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %2.2.6 Assessing the quality of data features
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            First_level_complexity_paper_GSR(dataset,number_subject,name_ROI_def,ROI_list,ROI_list2,procedure,all_procedure_names,SPM_dir,Work_dir) %Sum complexity over sessions
            
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%2.2.3 (Spectral) DCM for resting state fMRI (part 2)
%Number of excluded sessions for longi datasets (and subjects for longi datasets)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Number_included_sessions_paper_GSR(SPM_dir,Work_dir)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%2.2.4 Parametric empirical Bayes (part 2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PEB_group_mean_paper_GSR(SPM_dir,Work_dir) %longitudinal datasets only

PEB_group_mean_paper_GSR_DatasetHCP(SPM_dir,Work_dir)

%%%%%%%%%%%%%%%%%%%
%Calculate Results
%%%%%%%%%%%%%%%%%%%
Compute_effect_GSR_paper_GSR(SPM_dir,Work_dir)

Compute_effect_GSR_HCP_paper_GSR(SPM_dir,Work_dir)

Compute_effect_GSR_FL_complex_paper_GSR(SPM_dir,Work_dir)

Compute_effect_GSR_FL_complex_HCP_paper_GSR(SPM_dir,Work_dir)

%%%%%%%%%%%%%%%%%
%Produce Figures
%%%%%%%%%%%%%%%%%
Figures_paper_GSR(SPM_dir,Work_dir);
toc
