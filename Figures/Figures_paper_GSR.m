function Figures_paper_GSR(SPM_dir,Work_dir)

%%%%%%%%%
%Figures
%%%%%%%%%
%Figure 1
Mesh_DMN_paper_GSR(SPM_dir,Work_dir);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%WITHIN-NETWORK CONNECTIVITY: Longitudinal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Figure 2
Imagesc_group_PEB_Basic_DMN(SPM_dir,Work_dir);
Imagesc_group_PEB_GSR_regr_DMN(SPM_dir,Work_dir);

Imagesc_group_PEB_Basic_SMR(SPM_dir,Work_dir);
Imagesc_group_PEB_GSR_regr_SMR(SPM_dir,Work_dir);

Imagesc_group_PEB_Basic_SAL(SPM_dir,Work_dir);
Imagesc_group_PEB_GSR_regr_SAL(SPM_dir,Work_dir);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%BETWEEN-NETWORK CONNECTIVITY: longitudinal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Figure 3
Imagesc_group_PEB_Basic(SPM_dir,Work_dir);
Imagesc_group_PEB_GSR_regr(SPM_dir,Work_dir)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NOISE PARAMETERS: longitudinal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Figure 4
Figure_influence_spectral_btwn_ntwrk_paper_GSR_option4(SPM_dir,Work_dir);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%WITHIN-NETWORK CONNECTIVITY: HCP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Figure 5
Imagesc_group_PEB_Basic_DMN_HCP(SPM_dir,Work_dir);
Imagesc_group_PEB_GSR_regr_DMN_HCP(SPM_dir,Work_dir);

Imagesc_group_PEB_Basic_SMR_HCP(SPM_dir,Work_dir);
Imagesc_group_PEB_GSR_regr_SMR_HCP(SPM_dir,Work_dir);

Imagesc_group_PEB_Basic_SAL_HCP(SPM_dir,Work_dir);
Imagesc_group_PEB_GSR_regr_SAL_HCP(SPM_dir,Work_dir);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%BETWEEN-NETWORK CONNECTIVITY: longitudinal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Figure 6
Imagesc_group_PEB_Basic_HCP(SPM_dir,Work_dir);
Imagesc_group_PEB_GSR_regr_HCP(SPM_dir,Work_dir);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NOISE PARAMETERS: longitudinal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Figure 7
Figure_influence_spectral_btwn_ntwrk_paper_GSR_option4_HCP(SPM_dir,Work_dir);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%OPTIONAL FIGURES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Figure with group and subject-level estimates: small RSNs
Figure_influence_connectivity_paper_GSR_option4(SPM_dir,Work_dir);
Figure_influence_connectivity_paper_GSR_option4_HCP(SPM_dir,Work_dir);

%Figure with group and subject-level estimates: between-network connectivity RSNs
Figure_influence_btwn_ntwrk_connectivity_paper_GSR_option4(SPM_dir,Work_dir);
Figure_influence_btwn_ntwrk_connectivity_paper_GSR_option4_HCP(SPM_dir,Work_dir);

end