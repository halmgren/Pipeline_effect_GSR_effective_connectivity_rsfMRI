function Figure_influence_btwn_ntwrk_connectivity_paper_GSR_option4_HCP(SPM_dir,Work_dir)

all_ROI_defs={'Smith'}; %ROI definition based on Smith et al., 2009

all_procedure_names={'GSR_regr'};

for number_ROI_def=1:length(all_ROI_defs)
    
    name_ROI_def=all_ROI_defs{number_ROI_def};
    
    [ROI_list]=Define_ROIs_paper_GSR(name_ROI_def);
    [ROI_list2]=Define_comb_ROIs_paper_GSR(name_ROI_def);
    
    cd([Work_dir '/DatasetKuehn/sub-001/ses-001/func/VOI/Basic/Smith']);
    
    [pos_network,orig_network_abbrev2]=Rearrange_ROI_list2_thresholds(ROI_list,ROI_list2);
    
    tmp=1;
    tmp2=1;
    clear ntwrk_name2 ntwrk_name3;
    for VOI_number=1:size(ROI_list2,1)
        ntwrk=ROI_list2{VOI_number,1}(1:3);
        region=ROI_list2{VOI_number,1}(5:end);
        
        if VOI_number==1
            tmp2=1;
            ntwrk_abbrev2{tmp}=ROI_list2{VOI_number,1}(1:3);
            
        end
        
        if VOI_number>1 && ~strcmp(ROI_list2{VOI_number,1}(1:3),ROI_list2{VOI_number-1,1}(1:3))
            tmp2=1;
            tmp=tmp+1;
            ntwrk_abbrev2{tmp}=ROI_list2{VOI_number,1}(1:3);
        end
    end
    
    for number_procedure=1:length(all_procedure_names)
        procedure=all_procedure_names{number_procedure};
        
        direct=[Work_dir '/Figures_paper_GSR/PEB_Figures/' procedure '/' name_ROI_def '/Full_model/Figure_influence_connect'];
        
        mkdir(direct);
        
        for N_prec_comp={'single'}    %N_prec_comp={'single','fields'}
            
            for network_number=1:length(ntwrk_abbrev2)
                    
                colored=0;
                grd=0;
                cntr=0;
                colored_2=1;
                
                for parameter_type={'A'}
                    
                    tmp_subnet=0;
                    
                    ntwrk_sizes=cellfun('length',pos_network{:});
                    for subnet=1:length(pos_network{network_number})
                        for receiving_subnet=1:length(pos_network{network_number})
                            if subnet==receiving_subnet
                                continue;
                            end
                            
                            clear PEB;
                            load([Work_dir '/Results_paper_GSR/DCM/Basic/' name_ROI_def '/Full_model/PEB_group/PEB_' parameter_type{:} '_mean_comp_' procedure '_HCP_' cell2mat(N_prec_comp) '_' ntwrk_abbrev2{network_number} '.mat'],'PEB','DCM','GCM_full');
                            PEB1=PEB;
                            
                            %Put all values (# subject x # prmtrs) in one array
                            for subject=1:length(GCM_full)
                                PEB_group{subject}.Ep=spm_vec(GCM_full{subject}.Ep.A);
                            end
                            clear PEB_group1;
                            %Without GSR
                            tmp=0;
                            for subject=1:length(PEB_group)
                                if subject>1
                                    tmp=tmp+length(PEB_tmp2(:));
                                end
                                PEB_tmp=full(vec2mat(PEB_group{subject}.Ep,11,11)');
                                
                                clear PEB_tmp2
                                
                                PEB_tmp2=PEB_tmp(sum(ntwrk_sizes(1:receiving_subnet-1))+1:sum(ntwrk_sizes(1:receiving_subnet-1))+length(pos_network{network_number}{receiving_subnet}),sum(ntwrk_sizes(1:subnet-1))+1:sum(ntwrk_sizes(1:subnet-1))+length(pos_network{network_number}{subnet}));
                                
                                PEB_group1(tmp+1:tmp+length(PEB_tmp2(:)))=PEB_tmp2(:);
                                
                                clear PEB_tmp PEB_tmp;
                            end
                            
                            clear PEB PEB_group PEB_tmp2;
                            load([Work_dir '/Results_paper_GSR/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_group/PEB_' parameter_type{:} '_mean_comp_' procedure '_HCP_' cell2mat(N_prec_comp) '_' ntwrk_abbrev2{network_number} '.mat'],'PEB','DCM','GCM_full');
                            PEB2=PEB;
                            
                            for subject=1:length(GCM_full)
                                PEB_group{subject}.Ep=spm_vec(GCM_full{subject}.Ep.A);
                            end
                            
                            clear PEB_group2;
                            %With GSR
                            tmp=0;
                            for subject=1:length(PEB_group)
                                if subject>1
                                    tmp=tmp+length(PEB_tmp2(:));
                                end
                                
                                PEB_tmp=full(vec2mat(PEB_group{subject}.Ep,11,11)');
                                
                                PEB_tmp2=PEB_tmp(sum(ntwrk_sizes(1:receiving_subnet-1))+1:sum(ntwrk_sizes(1:receiving_subnet-1))+length(pos_network{network_number}{receiving_subnet}),sum(ntwrk_sizes(1:subnet-1))+1:sum(ntwrk_sizes(1:subnet-1))+length(pos_network{network_number}{subnet}));
                                
                                PEB_group2(tmp+1:tmp+length(PEB_tmp2(:)))=PEB_tmp2(:);
                                
                                clear PEB_tmp PEB_tmp;
                            end
                            
%                             clear PEB_tmp2;
                            
                            Y_1=spm_unvec(PEB_group1,ones(length(PEB_tmp2(:)),length(PEB_group)))';
                            Y_2=spm_unvec(PEB_group2,ones(length(PEB_tmp2(:)),length(PEB_group)))';
                            
                            Y_axis=[Y_1;Y_2];
                            %                     Y_axis=Y_axis(:);
                            
                            %X_position: Without GSR
                            X_Basic_1=[];
                            tmp=-1;
                            for prmtr1=1:size(PEB_tmp2,2)
                                for prmtr2=1:size(PEB_tmp2,1)
                                    
                                    if rem(prmtr2,sqrt(length(PEB_group{1}.Ep)))==1
                                        tmp=tmp+1;
                                    end
                                    
                                    X_Basic_1=[X_Basic_1 ones(1,length(PEB_group))'*prmtr2+size(PEB_tmp2,1)*(prmtr1-1)+tmp-0.1429];
                                    
                                    
                                end
                            end
                            
                            %X_position: With GSR
                            X_Basic_2=[];
                            tmp=-1;
                            for prmtr1=1:size(PEB_tmp2,2)
                                for prmtr2=1:size(PEB_tmp2,1)
                                    
                                    if rem(prmtr2,sqrt(length(PEB_group{1}.Ep)))==1
                                        tmp=tmp+1;
                                    end
                                    
                                    X_Basic_2=[X_Basic_2 ones(1,length(PEB_group))'*prmtr2+size(PEB_tmp2,1)*(prmtr1-1)+tmp+0.1429];
                                    
                                    
                                end
                            end
                            
                            X_axis=[X_Basic_1; X_Basic_2];
                            
                            figure('units','normalized','outerposition',[0 0 1 1]);
                            
                            for parameter=1:length(PEB_tmp2(:))
                                
                                barclr=rem(parameter-1,size(PEB_tmp2,1))+1;
                                
                                %Plot Basic
                                clr_basic=[0.2+0.6/size(PEB_tmp2,1)*barclr 0.2+0.6/size(PEB_tmp2,1)*barclr 0.2+0.6/size(PEB_tmp2,1)*barclr];
                                
                                if barclr==1
                                    clr_basic_2=[0 0 0.7]; %blue
                                elseif barclr==2
                                    clr_basic_2=[0.7 0 0]; %red
                                elseif barclr==3
                                    clr_basic_2=[0 0 0]; %black
                                elseif barclr==4
                                    clr_basic_2=[0 .7 0]; %green
                                elseif barclr==5
                                    clr_basic_2=[0.7 0 0.7]; %magenta
                                end
                                
                                
                                if colored==1
                                    plot(X_axis(1:length(PEB_group),parameter),Y_axis(1:length(PEB_group),parameter),'Color',clr_basic,'Marker','o','LineStyle','none','MarkerFaceColor',clr_basic,'MarkerEdgeColor','k');
                                elseif cntr==1
                                    plot(X_axis(1:length(PEB_group),parameter),Y_axis(1:length(PEB_group),parameter),'Color',clr_basic,'Marker','o','LineStyle','none','MarkerFaceColor',clr_basic,'MarkerEdgeColor','k','LineWidth',1.5);
                                elseif colored_2==1
                                    plot(X_axis(1:length(PEB_group),parameter),Y_axis(1:length(PEB_group),parameter),'Color',clr_basic_2,'Marker','o','LineStyle','none','MarkerFaceColor',clr_basic_2,'MarkerEdgeColor','k');
                                else
                                    plot(X_axis(1:length(PEB_group),parameter),Y_axis(1:length(PEB_group),parameter),'Color','k','Marker','o','LineStyle','none','MarkerFaceColor','k','MarkerEdgeColor','k');
                                end
                                
                                %Plot GSR
                                hold on;
                                if barclr==1
                                    clr_GSR=[0 0 1];
                                else
                                    clr_GSR=[0 0+1/(size(PEB_tmp2,1)-(size(PEB_tmp2,1)-4))*(barclr-(size(PEB_tmp2,1)-4)) 1];
                                end
                                
                                if barclr==1
                                    clr_GSR_2=[0.5 0.7 1]; %lightblue
                                elseif barclr==2
                                    clr_GSR_2=[1 0.5 0.5]; %lightred
                                elseif barclr==3
                                    clr_GSR_2=[0.5 0.5 0.5]; %grey
                                elseif barclr==4
                                    clr_GSR_2=[0.5 1 0.5]; %lightgreen
                                elseif barclr==5
                                    clr_GSR_2=[1 0.5 1]; %lightmagenta
                                end
                                
                                if colored==1
                                    plot(X_axis(length(PEB_group)+1:length(PEB_group)*2,parameter),Y_axis(length(PEB_group)+1:length(PEB_group)*2,parameter),'Color',clr_GSR,'Marker','o','LineStyle','none','MarkerFaceColor',clr_GSR,'MarkerEdgeColor','k');
                                elseif cntr==1
                                    plot(X_axis(length(PEB_group)+1:length(PEB_group)*2,parameter),Y_axis(length(PEB_group)+1:length(PEB_group)*2,parameter),'Color',clr_basic,'Marker','o','LineStyle','none','MarkerFaceColor',clr_basic,'MarkerEdgeColor',[0.65 0.65 0.65],'LineWidth',1.5);
                                elseif colored_2==1
                                    plot(X_axis(length(PEB_group)+1:length(PEB_group)*2,parameter),Y_axis(length(PEB_group)+1:length(PEB_group)*2,parameter),'Color',clr_GSR_2,'Marker','o','LineStyle','none','MarkerFaceColor',clr_GSR_2,'MarkerEdgeColor','k','LineWidth',1.5);
                                else
                                    plot(X_axis(length(PEB_group)+1:length(PEB_group)*2,parameter),Y_axis(length(PEB_group)+1:length(PEB_group)*2,parameter),'Color','k','Marker','o','LineStyle','none','MarkerFaceColor','k','MarkerEdgeColor','k');
                                end
                                hold on;
                                
                            end
                            
                            sz1=size(PEB_tmp2,1);
                            sz2=size(PEB_tmp2,2);
                            
                            %Rearrange estimates & covariances
                            clear PEB_tmp;
                            PEB_tmp=full(vec2mat(PEB1.Ep,11,11)');
                               
                            PEB_group_bas_tmp.Ep=PEB_tmp(sum(ntwrk_sizes(1:receiving_subnet-1))+1:sum(ntwrk_sizes(1:receiving_subnet-1))+length(pos_network{network_number}{receiving_subnet}),sum(ntwrk_sizes(1:subnet-1))+1:sum(ntwrk_sizes(1:subnet-1))+length(pos_network{network_number}{subnet}));
                            
                            PEB_group_bas.Ep=PEB_group_bas_tmp.Ep(:);
                            
                            clear PEB_tmp;
                            
                            PEB_tmp=full(vec2mat(diag(PEB1.Cp),11,11)');
                                
                            PEB_group_bas_tmp.Cp=PEB_tmp(sum(ntwrk_sizes(1:receiving_subnet-1))+1:sum(ntwrk_sizes(1:receiving_subnet-1))+length(pos_network{network_number}{receiving_subnet}),sum(ntwrk_sizes(1:subnet-1))+1:sum(ntwrk_sizes(1:subnet-1))+length(pos_network{network_number}{subnet}));
                            
                            PEB_group_bas.Cp=PEB_group_bas_tmp.Cp(:);
                            
                            clear PEB_tmp;
                            PEB_tmp=full(vec2mat(PEB2.Ep,11,11)');
                               
                            PEB_group_GSR_tmp.Ep=PEB_tmp(sum(ntwrk_sizes(1:receiving_subnet-1))+1:sum(ntwrk_sizes(1:receiving_subnet-1))+length(pos_network{network_number}{receiving_subnet}),sum(ntwrk_sizes(1:subnet-1))+1:sum(ntwrk_sizes(1:subnet-1))+length(pos_network{network_number}{subnet}));
                            
                            PEB_group_GSR.Ep=PEB_group_GSR_tmp.Ep(:);
                            
                            clear PEB_tmp;
                            
                            PEB_tmp=full(vec2mat(diag(PEB2.Cp),11,11)');
                                
                            PEB_group_GSR_tmp.Cp=PEB_tmp(sum(ntwrk_sizes(1:receiving_subnet-1))+1:sum(ntwrk_sizes(1:receiving_subnet-1))+length(pos_network{network_number}{receiving_subnet}),sum(ntwrk_sizes(1:subnet-1))+1:sum(ntwrk_sizes(1:subnet-1))+length(pos_network{network_number}{subnet}));
                            
                            PEB_group_GSR.Cp=PEB_group_GSR_tmp.Cp(:);
                            
                            %load random GCM for region names
                            load([Work_dir '/DatasetKirby/sub-001/ses-001/func/DCM/' procedure '/' name_ROI_def '/Full_model/GCM_full_estim_' ntwrk_abbrev2{network_number} '.mat']);
                            
                            [hh xb]=plot_ci_adapt12(full([PEB_group_bas.Ep';PEB_group_GSR.Ep']),[PEB_group_bas.Cp';PEB_group_GSR.Cp'],[],[],[],[1 1 1]*0.5,0.9,sz1,sz2,GCM{1},grd,colored,cntr,colored_2,ntwrk_sizes,subnet);
                            hold on;
                            
                            if colored==1
                                alpha 0.5
                            elseif cntr
                                alpha 0.5
                            elseif colored_2==1
                                alpha 0.7
                            end
                            %%%%%%%%%%%%%%%
                            %Create legend
                            %%%%%%%%%%%%%%%
                            clear clr_basic clr_GSR  clr_basic_2 clr_GSR_2
                            AxesHandle=findobj(gcf,'Type','axes');
                            pt1 = get(AxesHandle,{'Position'});
                            
                            clr={'r','b','g','k','y'};
                            if colored==1
                                for barclr=1:size(PEB_tmp2,1)
                                    clr_basic{barclr}=[0.2+0.6/size(PEB_tmp2,1)*barclr 0.2+0.6/size(PEB_tmp2,1)*barclr 0.2+0.6/size(PEB_tmp2,1)*barclr];
                                end
                                
                                for barclr=1:size(PEB_tmp2,1)
                                    %                         barclr=rem(connections-1,sqrt(length(hh)))+1;
                                    if barclr==1
                                        clr_GSR{barclr}=[0 0 1];
                                        %                         elseif barclr==2
                                        %                             clr_GSR{barclr}=[0 0 1];
                                    else
                                        clr_GSR{barclr}=[0 0+1/(size(PEB_tmp2,1)-(size(PEB_tmp2,1)-4))*(barclr-(size(PEB_tmp2,1)-4)) 1];
                                    end
                                end
                                
                                tmp=0;
                                for region=1:sum(ntwrk_sizes(1:receiving_subnet-1))+1:sum(ntwrk_sizes(1:receiving_subnet-1))+sz1
                                    tmp=tmp+1;
                                    names{region}=GCM{1}.xY(region).name(5:end);
                                    annotation('textbox','string','','BackgroundColor',clr_basic{tmp},'Position',[pt1{1}(3)-0.02 pt1{1}(4)+0.1-tmp*0.10/size(PEB_tmp2,1) 0.05 0.10/size(PEB_tmp2,1)],'FontWeight','bold','FontSize',12);
                                end
                                
                                tmp=0;
                                for region=1:sum(ntwrk_sizes(1:receiving_subnet-1))+1:sum(ntwrk_sizes(1:receiving_subnet-1))+sz1
                                    tmp=tmp+1;
                                    annotation('textbox','string','','BackgroundColor',clr_GSR{tmp},'Position',[pt1{1}(3)+0.03 pt1{1}(4)+0.1-tmp*0.10/size(PEB_tmp2,1) 0.05 0.10/size(PEB_tmp2,1)],'FontWeight','bold','FontSize',12);
                                end
                                
                                %Connectivity to
                                annotation('textbox','string','Target: ','Position',[pt1{1}(3)+0.05+0.03 pt1{1}(4)+0.10 0.10 0.03],'FontWeight','bold','FontSize',12,'LineStyle','none');
                                
                                %GSR and Basic above colors
                                annotation('textbox','string','w/o GSR','HorizontalAlignment','left','Position',[pt1{1}(3)-0.10+0.05+0.03 pt1{1}(4)+0.10 0.10 0.03],'FontWeight','bold','FontSize',12,'LineStyle','none');
                                annotation('textbox','string','with GSR','HorizontalAlignment','left','Position',[pt1{1}(3)-0.05+0.05+0.03 pt1{1}(4)+0.10 0.10 0.03],'FontWeight','bold','FontSize',12,'LineStyle','none');
                                
                                %strings
                                tmp=0;
                                for region=1:size(PEB_tmp2,1)
                                    tmp=tmp+1;
                                    annotation('textbox','string',names{region},'VerticalAlignment','middle','Position',[pt1{1}(3)+0.05+0.03 pt1{1}(4)+0.1-tmp*0.10/size(PEB_tmp2,1) 0.10 0.10/size(PEB_tmp2,1)],'FontWeight','bold','FontSize',12,'LineStyle','none');
                                end
                            elseif cntr==1
                                for barclr=1:size(PEB_tmp2,1)
                                    clr_basic{barclr}=[0.2+0.6/size(PEB_tmp2,1)*barclr 0.2+0.6/size(PEB_tmp2,1)*barclr 0.2+0.6/size(PEB_tmp2,1)*barclr];
                                end
                                
                                for barclr=1:size(PEB_tmp2,1)
                                    clr_GSR{barclr}=[0.2+0.6/size(PEB_tmp2,1)*barclr 0.2+0.6/size(PEB_tmp2,1)*barclr 0.2+0.6/size(PEB_tmp2,1)*barclr];
                                end
                                
                                tmp=0;
                                for region=sum(ntwrk_sizes(1:receiving_subnet-1))+1:sum(ntwrk_sizes(1:receiving_subnet-1))+sz1
                                    tmp=tmp+1;
                                    names{region}=GCM{1}.xY(region).name(5:end);
                                    annotation('textbox','string','','BackgroundColor',clr_basic{tmp},'Position',[pt1{1}(3)-0.02 pt1{1}(4)+0.1-tmp*0.10/size(PEB_tmp2,1) 0.05 0.10/size(PEB_tmp2,1)],'FontWeight','bold','FontSize',12,'EdgeColor',[0 0 0],'LineWidth',1.25);
                                end
                                
                                tmp=0;
                                for region=1:sum(ntwrk_sizes(1:receiving_subnet-1))+1:sum(ntwrk_sizes(1:receiving_subnet-1))+sz1
                                    tmp=tmp+1;
                                    annotation('textbox','string','','BackgroundColor',clr_basic{tmp},'Position',[pt1{1}(3)+0.03 pt1{1}(4)+0.1-tmp*0.10/size(PEB_tmp2,1) 0.05 0.10/size(PEB_tmp2,1)],'FontWeight','bold','FontSize',12,'EdgeColor',[0.65 0.65 0.65],'LineWidth',1.25);
                                end
                                
                                %Connectivity to
                                annotation('textbox','string','Target: ','Position',[pt1{1}(3)+0.05+0.03 pt1{1}(4)+0.10 0.10 0.03],'FontWeight','bold','FontSize',12,'LineStyle','none');
                                
                                %GSR and Basic above colors
                                annotation('textbox','string','w/o GSR','HorizontalAlignment','left','Position',[pt1{1}(3)-0.10+0.05+0.03 pt1{1}(4)+0.10 0.10 0.03],'FontWeight','bold','FontSize',12,'LineStyle','none');
                                annotation('textbox','string','with GSR','HorizontalAlignment','left','Position',[pt1{1}(3)-0.05+0.05+0.03 pt1{1}(4)+0.10 0.10 0.03],'FontWeight','bold','FontSize',12,'LineStyle','none');
                                
                                %strings
                                tmp=0;
                                for region=1:size(PEB_tmp2,1)
                                    tmp=tmp+1;
                                    annotation('textbox','string',names{region},'VerticalAlignment','middle','Position',[pt1{1}(3)+0.05+0.03 pt1{1}(4)+0.1-tmp*0.10/size(PEB_tmp2,1) 0.10 0.10/size(PEB_tmp2,1)],'FontWeight','bold','FontSize',12,'LineStyle','none');
                                end
                            elseif colored_2==1
                                for barclr=1:size(PEB_tmp2,1)
                                    if barclr==1
                                        clr_basic{barclr}=[0 0 0.7]; %blue
                                    elseif barclr==2
                                        clr_basic{barclr}=[0.7 0 0]; %red
                                    elseif barclr==3
                                        clr_basic{barclr}=[0 0 0]; %black
                                    elseif barclr==4
                                        clr_basic{barclr}=[0 .7 0]; %green
                                    elseif barclr==5
                                        clr_basic{barclr}=[0.7 0 0.7]; %magenta
                                    end
                                end
                                
                                for barclr=1:size(PEB_tmp2,1)
                                    if barclr==1
                                        clr_GSR{barclr}=[0.5 0.7 1]; %lightblue
                                    elseif barclr==2
                                        clr_GSR{barclr}=[1 0.5 0.5]; %lightred
                                    elseif barclr==3
                                        clr_GSR{barclr}=[0.5 0.5 0.5]; %grey
                                    elseif barclr==4
                                        clr_GSR{barclr}=[0.5 1 0.5]; %lightgreen
                                    elseif barclr==5
                                        clr_GSR{barclr}=[1 0.5 1]; %lightmagenta
                                    end
                                end
                                
                                %color basic blocks
                                tmp=0;
                                for region=sum(ntwrk_sizes(1:receiving_subnet-1))+1:sum(ntwrk_sizes(1:receiving_subnet-1))+sz1
                                    tmp=tmp+1;
                                    names{region-sum(ntwrk_sizes(1:receiving_subnet-1))}=GCM{1}.xY(region).name(5:end);
                                    annotation('textbox','string','','BackgroundColor',clr_basic{tmp},'Position',[pt1{1}(3)-0.02 pt1{1}(4)+0.15-tmp*0.15/size(PEB_tmp2,1) 0.075 0.15/size(PEB_tmp2,1)],'FontWeight','bold','FontSize',12,'FaceAlpha',1);
                                end
                                
                                %color GSR blocks
                                tmp=0;
                                for region=sum(ntwrk_sizes(1:receiving_subnet-1))+1:sum(ntwrk_sizes(1:receiving_subnet-1))+sz1
                                    tmp=tmp+1;
                                    annotation('textbox','string','','BackgroundColor',clr_GSR{tmp},'Position',[pt1{1}(3)+0.055 pt1{1}(4)+0.15-tmp*0.15/size(PEB_tmp2,1) 0.075 0.15/size(PEB_tmp2,1)],'FontWeight','bold','FontSize',12,'FaceAlpha',1);
                                end
                                
                                %Connectivity to
                                annotation('textbox','string','Target: ','Position',[pt1{1}(3)+0.07+0.07 pt1{1}(4)+0.16 0.10 0.03],'FontWeight','bold','FontSize',20,'LineStyle','none');
                                
                                %GSR and Basic above colors
                                annotation('textbox','string','w/o GSR','HorizontalAlignment','left','Position',[pt1{1}(3)-0.10+0.05+0.033 pt1{1}(4)+0.16 0.10 0.03],'FontWeight','bold','FontSize',20,'LineStyle','none');
                                annotation('textbox','string','with GSR','HorizontalAlignment','left','Position',[pt1{1}(3)-0.05+0.06+0.045 pt1{1}(4)+0.16 0.10 0.03],'FontWeight','bold','FontSize',20,'LineStyle','none');
                                
                                %regions
                                tmp=0;
                                for region=1:size(PEB_tmp2,1)
                                    tmp=tmp+1;
                                    annotation('textbox','string',names{region},'VerticalAlignment','middle','Position',[pt1{1}(3)+0.07+0.07 pt1{1}(4)+0.15-tmp*0.15/size(PEB_tmp2,1) 0.10 0.15/size(PEB_tmp2,1)],'FontWeight','bold','FontSize',20,'LineStyle','none');
                                end
                            end
                            
                            set(gcf,'PaperPositionMode','auto');
                            
                            clear names;
                            if colored==1
                                saveas(gcf,[direct '/' char(N_prec_comp) '_' parameter_type{:} '_from_' orig_network_abbrev2{network_number}(sum(ntwrk_sizes(1:subnet-1))+1,:) '_to_' orig_network_abbrev2{network_number}(sum(ntwrk_sizes(1:receiving_subnet-1))+1,:) '_HCP_colored.bmp']);
                            elseif cntr==1
                                saveas(gcf,[direct '/' char(N_prec_comp) '_' parameter_type{:} '_from_' orig_network_abbrev2{network_number}(sum(ntwrk_sizes(1:subnet-1))+1,:) '_to_' orig_network_abbrev2{network_number}(sum(ntwrk_sizes(1:receiving_subnet-1))+1,:) '_HCP_contour_grey.bmp']);
                            elseif colored_2==1
                                saveas(gcf,[direct '/' char(N_prec_comp) '_' parameter_type{:} '_from_' orig_network_abbrev2{network_number}(sum(ntwrk_sizes(1:subnet-1))+1,:) '_to_' orig_network_abbrev2{network_number}(sum(ntwrk_sizes(1:receiving_subnet-1))+1,:) '_HCP_colored_2.bmp']);
                            end
                            close;
                            
                            clearvars -except parameter_type N_prec_comp procedure ntwrk_name pos_network network_number ntwrk_size name_ROI_def ROI_list direct procedure subnet receiving_subnet ntwrk_abbrev2 ntwrk_sizes colored cntr colored_2 grd  orig_network_abbrev2 Work_dir SPM_dir
                        end
                    end
                end
                
               
            end
            
        end
    end
end

end