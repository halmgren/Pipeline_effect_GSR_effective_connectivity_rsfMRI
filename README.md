This repository contains the code that was used for the analyses of the study concerning the effect of global signal regression on spectral DCM estimates

IMPORTANT: This is still the code for the preprint, the code for the published paper will be uploaded soon.

The main pipeline is located in the 'Main_code'-folder. The instructions for reproduction (and data availability) are included in the 'Pipeline_paper_GSR.m'-file, and the pipeline can be applied by executing the same file.

The input files are included in the five datasets ('Myconnectome', 'Midnight Scan Club', 'Day2day', 'Kirby', 'HCP') and should be converted to BIDS format (if not downloaded in this format). The code does all analyses in the paper, and creates figures (results section is not yet added). Differences in software versions can have an effect on results and inference.

Correspondence between code and paper goes as follows:

Folders on github:

- 'Subjects': Shows and produces comprehensive list of HCP subjects that were included (see paragraph 2.1 in paper)
- 'Preprocessing': All preprocessing steps for anatomical and functional images (see paragraph 2.2.1 in paper)
- 'Time_series_extraction': All steps for time-series extraction (see paragraph 2.2.2 in paper)
- 'spDCM': Steps and details concerning spectral DCM specification and estimation (see paragraph 2.2.3 in paper)
- 'PEB': Steps and details concerning PEB specification and estimation (see paragraph 2.2.4 in paper)
- 'Quality_data_features': Steps and details concerning quality of data features and precision of posterior estimates (see paragraph 2.2.6 in paper)
- 'Figures': Steps and details to reconstruct the figures (see figures in paragraph 3.1 and 3.2)
- 'Datasets_information': Information concerning the datasets (e.g., multi- vs single band)

This code is provided without any warranty.

For feedback and comments, please contact Hannes Almgren (Hannes.Almgren@ugent.be)
