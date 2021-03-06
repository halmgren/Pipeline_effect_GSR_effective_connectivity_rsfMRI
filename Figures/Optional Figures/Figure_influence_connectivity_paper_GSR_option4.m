function Figure_influence_connectivity_paper_GSR_option4(SPM_dir,Work_dir)

all_ROI_defs={'Smith'}; %ROI definition based on Smith et al., 2009

all_procedure_names={'GSR_regr'};

for number_ROI_def=1:length(all_ROI_defs)
    
    name_ROI_def=all_ROI_defs{number_ROI_def};
    
    [ROI_list]=Define_ROIs_paper_GSR(name_ROI_def);
    
    tmp=0;
    
    for VOI_number=1:size(ROI_list,1)
        ntwrk=ROI_list{VOI_number,1}(1:3);
        
        if VOI_number>1 && strcmp(ROI_list{VOI_number,1}(1:3),ROI_list{VOI_number-1,1}(1:3))
            ntwrk_size(tmp)=ntwrk_size(tmp)+1;
            continue
            
        else
            tmp=tmp+1;
            ntwrk_size(tmp)=1;
            ntwrk_name{tmp}=ROI_list{VOI_number,1}(1:3);
        end
    end
    
    for number_procedure=1:length(all_procedure_names)
        procedure=all_procedure_names{number_procedure};
        
        direct=[Work_dir '/Figures_paper_GSR/PEB_Figures/' procedure '/' name_ROI_def '/Full_model/Figure_influence_connect'];
        mkdir(direct);
        
        for N_prec_comp={'single'}    %N_prec_comp={'single','fields','all'}
            
            for network_number=1:length(ntwrk_name)
                    
                colored=0;
                grd=0;
                cntr=0;
                colored_2=1;
                
                for parameter_type={'A'}
                    clear PEB;
                    load([Work_dir '/Results_paper_GSR/DCM/Basic/' name_ROI_def '/Full_model/PEB_group/PEB_' parameter_type{:} '_mean_comp_' procedure '_group_' cell2mat(N_prec_comp) '_' ntwrk_name{network_number} '.mat'],'PEB','DCM','PEB_group');
                    PEB1=PEB;
                    
                    tmp=0;
                    for subject=1:length(PEB_group)
                        if subject>1
                            tmp=tmp+length(PEB_group{subject}.Ep);
                        end
                        PEB_group1(tmp+1:tmp+length(PEB_group{subject}.Ep))=PEB_group{subject}.Ep;
                    end
                    
                    %put matrix in correct order for plotting
                    
                    
                    clear PEB PEB_group;
                    load([Work_dir '/Results_paper_GSR/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_group/PEB_' parameter_type{:} '_mean_comp_' procedure '_group_' cell2mat(N_prec_comp) '_' ntwrk_name{network_number} '.mat'],'PEB','DCM','PEB_group');
                    PEB2=PEB;
                    
                    tmp=0;
                    for subject=1:length(PEB_group)
                        if subject>1
                            tmp=tmp+length(PEB_group{subject}.Ep);
                        end
                        PEB_group2(tmp+1:tmp+length(PEB_group{subject}.Ep))=PEB_group{subject}.Ep;
                    end
                    
                    Y_1=spm_unvec(PEB_group1,ones(length(PEB_group{1}.Ep),length(PEB_group)))';
                    Y_2=spm_unvec(PEB_group2,ones(length(PEB_group{1}.Ep),length(PEB_group)))';
                    
                    Y_axis=[Y_1;Y_2];
%                     Y_axis=Y_axis(:);
                    
                    X_Basic_1=[];
                    tmp=-1;
                    for prmtr1=1:sqrt(length(PEB_group{subject}.Ep))
                        for prmtr2=1:sqrt(length(PEB_group{subject}.Ep))
                            
                            if rem(prmtr2,sqrt(length(PEB_group{1}.Ep)))==1
                                tmp=tmp+1;
                            end
                            
                            X_Basic_1=[X_Basic_1 ones(1,length(PEB_group))'*prmtr2+sqrt(length(PEB_group{1}.Ep))*(prmtr1-1)+tmp-0.1429];
                            
                            
                        end
                    end
                    
                    
                    X_Basic_2=[];
                    tmp=-1;
                    for prmtr1=1:sqrt(length(PEB_group{subject}.Ep))
                        for prmtr2=1:sqrt(length(PEB_group{subject}.Ep))
                            
                            if rem(prmtr2,sqrt(length(PEB_group{1}.Ep)))==1
                                tmp=tmp+1;
                            end
                            
                            
                            
                            X_Basic_2=[X_Basic_2 ones(1,length(PEB_group))'*prmtr2+sqrt(length(PEB_group{1}.Ep))*(prmtr1-1)+tmp+0.1429];
                        end
                    end
                    
                    X_axis=[X_Basic_1; X_Basic_2];
                    
                    figure('units','normalized','outerposition',[0 0 1 1]);
                    
                    for parameter=1:length(PEB_group{1}.Ep)
                        
                        barclr=rem(parameter-1,sqrt(length(PEB_group{1}.Ep)))+1;
                        
                        %Plot Basic
                        clr_basic=[0.2+0.6/sqrt(length(PEB_group{1}.Ep))*barclr 0.2+0.6/sqrt(length(PEB_group{1}.Ep))*barclr 0.2+0.6/sqrt(length(PEB_group{1}.Ep))*barclr];
                        
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
                            clr_GSR=[0 0+1/(sqrt(length(PEB_group{1}.Ep))-(sqrt(length(PEB_group{1}.Ep))-4))*(barclr-(sqrt(length(PEB_group{1}.Ep))-4)) 1];
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
                    
                    sz=sqrt(length(PEB_group{1}.Ep));
                    
                    %load random GCM for region names
                    load([Work_dir '/DatasetKirby/sub-001/ses-001/func/DCM/' procedure '/' name_ROI_def '/Full_model/GCM_full_estim_' ntwrk_name{network_number} '.mat']);
                    [hh xb]=plot_ci_adapt6(full([PEB1.Ep';PEB2.Ep']),[diag(PEB1.Cp)';diag(PEB2.Cp)'],[],[],[],[1 1 1]*0.5,0.9,sz,GCM{1},grd,colored,cntr,colored_2);
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
                        for barclr=1:length(GCM{1}.xY)
                            clr_basic{barclr}=[0.2+0.6/sqrt(length(hh))*barclr 0.2+0.6/sqrt(length(hh))*barclr 0.2+0.6/sqrt(length(hh))*barclr];
                        end
                        
                        for barclr=1:length(GCM{1}.xY)
                            %                         barclr=rem(connections-1,sqrt(length(hh)))+1;
                            if barclr==1
                                clr_GSR{barclr}=[0 0 1];
                                %                         elseif barclr==2
                                %                             clr_GSR{barclr}=[0 0 1];
                            else
                                clr_GSR{barclr}=[0 0+1/(sqrt(length(hh))-(sz-4))*(barclr-(sz-4)) 1];
                            end
                        end
                        
                        tmp=0;
                        for region=1:length(GCM{1}.xY)
                            tmp=tmp+1;
                            names{region}=GCM{1}.xY(region).name(5:end);
                            annotation('textbox','string','','BackgroundColor',clr_basic{tmp},'Position',[pt1{1}(3)-0.02 pt1{1}(4)+0.1-tmp*0.10/length(GCM{1}.xY) 0.05 0.10/length(GCM{1}.xY)],'FontWeight','bold','FontSize',12);
                        end
                        
                        tmp=0;
                        for region=1:length(GCM{1}.xY)
                            tmp=tmp+1;
                            annotation('textbox','string','','BackgroundColor',clr_GSR{tmp},'Position',[pt1{1}(3)+0.03 pt1{1}(4)+0.1-tmp*0.10/length(GCM{1}.xY) 0.05 0.10/length(GCM{1}.xY)],'FontWeight','bold','FontSize',12);
                        end
                        
                        %Connectivity to
                        annotation('textbox','string','Target: ','Position',[pt1{1}(3)+0.05+0.03 pt1{1}(4)+0.10 0.10 0.03],'FontWeight','bold','FontSize',12,'LineStyle','none');
                        
                        %GSR and Basic above colors
                        annotation('textbox','string','w/o GSR','HorizontalAlignment','left','Position',[pt1{1}(3)-0.10+0.05+0.03 pt1{1}(4)+0.10 0.10 0.03],'FontWeight','bold','FontSize',12,'LineStyle','none');
                        annotation('textbox','string','with GSR','HorizontalAlignment','left','Position',[pt1{1}(3)-0.05+0.05+0.03 pt1{1}(4)+0.10 0.10 0.03],'FontWeight','bold','FontSize',12,'LineStyle','none');
                        
                        %strings
                        tmp=0;
                        for region=1:length(GCM{1}.xY)
                            tmp=tmp+1;
                            annotation('textbox','string',names{region},'VerticalAlignment','middle','Position',[pt1{1}(3)+0.05+0.03 pt1{1}(4)+0.1-tmp*0.10/length(GCM{1}.xY) 0.10 0.10/length(GCM{1}.xY)],'FontWeight','bold','FontSize',12,'LineStyle','none');
                        end
                    elseif cntr==1
                        for barclr=1:length(GCM{1}.xY)
                            clr_basic{barclr}=[0.2+0.6/sqrt(length(hh))*barclr 0.2+0.6/sqrt(length(hh))*barclr 0.2+0.6/sqrt(length(hh))*barclr];
                        end
                        
                        for barclr=1:length(GCM{1}.xY)
                            clr_GSR{barclr}=[0.2+0.6/sqrt(length(hh))*barclr 0.2+0.6/sqrt(length(hh))*barclr 0.2+0.6/sqrt(length(hh))*barclr];
                        end
                        
                        tmp=0;
                        for region=1:length(GCM{1}.xY)
                            tmp=tmp+1;
                            names{region}=GCM{1}.xY(region).name(5:end);
                            annotation('textbox','string','','BackgroundColor',clr_basic{tmp},'Position',[pt1{1}(3)-0.02 pt1{1}(4)+0.1-tmp*0.10/length(GCM{1}.xY) 0.05 0.10/length(GCM{1}.xY)],'FontWeight','bold','FontSize',12,'EdgeColor',[0 0 0],'LineWidth',1.25);
                        end
                        
                        tmp=0;
                        for region=1:length(GCM{1}.xY)
                            tmp=tmp+1;
                            annotation('textbox','string','','BackgroundColor',clr_basic{tmp},'Position',[pt1{1}(3)+0.03 pt1{1}(4)+0.1-tmp*0.10/length(GCM{1}.xY) 0.05 0.10/length(GCM{1}.xY)],'FontWeight','bold','FontSize',12,'EdgeColor',[0.65 0.65 0.65],'LineWidth',1.25);
                        end
                        
                        %Connectivity to
                        annotation('textbox','string','Target: ','Position',[pt1{1}(3)+0.05+0.03 pt1{1}(4)+0.10 0.10 0.03],'FontWeight','bold','FontSize',12,'LineStyle','none');
                        
                        %GSR and Basic above colors
                        annotation('textbox','string','w/o GSR','HorizontalAlignment','left','Position',[pt1{1}(3)-0.10+0.05+0.03 pt1{1}(4)+0.10 0.10 0.03],'FontWeight','bold','FontSize',12,'LineStyle','none');
                        annotation('textbox','string','with GSR','HorizontalAlignment','left','Position',[pt1{1}(3)-0.05+0.05+0.03 pt1{1}(4)+0.10 0.10 0.03],'FontWeight','bold','FontSize',12,'LineStyle','none');
                        
                        %strings
                        tmp=0;
                        for region=1:length(GCM{1}.xY)
                            tmp=tmp+1;
                            annotation('textbox','string',names{region},'VerticalAlignment','middle','Position',[pt1{1}(3)+0.05+0.03 pt1{1}(4)+0.1-tmp*0.10/length(GCM{1}.xY) 0.10 0.10/length(GCM{1}.xY)],'FontWeight','bold','FontSize',12,'LineStyle','none');
                        end
                    elseif colored_2==1
                        for barclr=1:length(GCM{1}.xY)
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
                        
                        for barclr=1:length(GCM{1}.xY)
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
                        for region=1:length(GCM{1}.xY)
                            tmp=tmp+1;
                            names{region}=GCM{1}.xY(region).name(5:end);
                            annotation('textbox','string','','BackgroundColor',clr_basic{tmp},'Position',[pt1{1}(3)-0.02 pt1{1}(4)+0.15-tmp*0.15/length(GCM{1}.xY) 0.075 0.15/length(GCM{1}.xY)],'FontWeight','bold','FontSize',12,'FaceAlpha',1);
                        end
                        
                        %color GSR blocks
                        tmp=0;
                        for region=1:length(GCM{1}.xY)
                            tmp=tmp+1;
                            annotation('textbox','string','','BackgroundColor',clr_GSR{tmp},'Position',[pt1{1}(3)+0.055 pt1{1}(4)+0.15-tmp*0.15/length(GCM{1}.xY) 0.075 0.15/length(GCM{1}.xY)],'FontWeight','bold','FontSize',12,'FaceAlpha',1);
                        end
                        
                        %Connectivity to
                        annotation('textbox','string','Target: ','Position',[pt1{1}(3)+0.07+0.07 pt1{1}(4)+0.16 0.10 0.03],'FontWeight','bold','FontSize',20,'LineStyle','none');
                        
                        %GSR and Basic above colors
                        annotation('textbox','string','w/o GSR','HorizontalAlignment','left','Position',[pt1{1}(3)-0.10+0.05+0.033 pt1{1}(4)+0.16 0.10 0.03],'FontWeight','bold','FontSize',20,'LineStyle','none');
                        annotation('textbox','string','with GSR','HorizontalAlignment','left','Position',[pt1{1}(3)-0.05+0.06+0.045 pt1{1}(4)+0.16 0.10 0.03],'FontWeight','bold','FontSize',20,'LineStyle','none');
                        
                        %regions
                        tmp=0;
                        for region=1:length(GCM{1}.xY)
                            tmp=tmp+1;
                            annotation('textbox','string',names{region},'VerticalAlignment','middle','Position',[pt1{1}(3)+0.07+0.07 pt1{1}(4)+0.15-tmp*0.15/length(GCM{1}.xY) 0.10 0.15/length(GCM{1}.xY)],'FontWeight','bold','FontSize',20,'LineStyle','none');
                        end
                    end
                    
                    set(gcf,'PaperPositionMode','auto');
                    
                    clear names;
                    if colored==1
                        saveas(gcf,[direct '/' char(N_prec_comp) '_' parameter_type{:} '_' ntwrk_name{network_number} '_colored.bmp']);
                    elseif cntr==1
                        saveas(gcf,[direct '/' char(N_prec_comp) '_' parameter_type{:} '_' ntwrk_name{network_number} '_contour_grey.bmp']);
                    elseif colored_2==1
                        saveas(gcf,[direct '/' char(N_prec_comp) '_' parameter_type{:} '_' ntwrk_name{network_number} '_colored_2.bmp']);
                    end
                    
                    close;
                    
                    clearvars -except parameter_type N_prec_comp procedure ntwrk_name  network_number ntwrk_size name_ROI_def ROI_list direct procedure subnet receiving_subnet Work_dir SPM_dir
                    
                end
                
               
            end
            
        end
    end
end

end