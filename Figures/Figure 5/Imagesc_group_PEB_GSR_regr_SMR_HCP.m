function Imagesc_group_PEB_GSR_regr_SMR_HCP(SPM_dir,Work_dir)

procedure='GSR_regr';
name_ROI_def='Smith';

load([Work_dir '/Results_paper_GSR/DCM/' procedure '/' name_ROI_def '/Full_model/PEB_group/PEB_A_mean_comp_GSR_regr_HCP_single_SMR.mat'],'PEB','DCM','PEB_group');
load([Work_dir '/DatasetHCP/sub-001_summary/DCM/Basic/Smith/Full_model/GCM_SMR_full_estim.mat']);
DCM=GCM;

figure('units','centimeters','outerposition',[0 0 28 28]);

% B=vec2mat(PEB.Ep,11)';
% B=B([1,4,6,8,10,2,5,9,3,7,11],[1,4,6,8,10,2,5,9,3,7,11]);
% PEB.Ep=B(:);
% 
% for Nrow=1:121
%     A_tmp=vec2mat(PEB.Cp(:,Nrow),11)';
%     A_tmp=A_tmp([1,4,6,8,10,2,5,9,3,7,11],[1,4,6,8,10,2,5,9,3,7,11]);
%     PEB.Cp(:,Nrow)=A_tmp(:);
% end
% 
% for Ncol=1:121
%     B_tmp=vec2mat(PEB.Cp(Ncol,:),11)';
%     B_tmp=B_tmp([1,4,6,8,10,2,5,9,3,7,11],[1,4,6,8,10,2,5,9,3,7,11]);
%     PEB.Cp(Ncol,:)=B_tmp(:);
% end

% reshape(vec2mat(PEB.Cp(:,1)',11),[121,1]);
% PEB.Cp()=PEB.Cp();

ci=spm_invNcdf(1-0.05);
EP=full(vec2mat(PEB.Ep(1:9),3)');
CP=diag(PEB.Cp);
CP=full(vec2mat(CP(1:9),3)');
sgn=sign(EP-ci*sqrt(CP)).*sign(EP+ci*sqrt(CP));

A_matrix=full(vec2mat(PEB.Ep(1:9),3)');

A_matrix_raw=A_matrix;

% A_matrix=A_matrix+diag(diag(A_matrix)*NaN); %discard diagonal elements
% A_matrix(1:5,1:5)=NaN(5,5);
% A_matrix(6:8,6:8)=NaN(3,3);
% A_matrix(9:11,9:11)=NaN(3,3);
A_matrix=A_matrix-diag(diag(A_matrix))-diag(diag(exp(A_matrix)))/2; %rescale diagonal elements

h=imagesc(A_matrix);
set(h,'alphadata',~isnan(A_matrix))

colorbar;

axis square;
xlabel('\textbf{\underline{From}}','FontSize',40,'Fontweight','bold','Interpreter','latex'); ylabel('\textbf{\underline{To}}','FontSize',40,'Fontweight','bold','Interpreter','latex');
set(gca,'XAxisLocation','top','Tag','connectivity');

title_str = ['Connectivity'];

title(title_str,'FontSize',46);

% region_order=([1,4,6,8,10,2,5,9,3,7,11]);

for region=1:length({DCM{1}.xY.name})
    regions(region)={DCM{1}.xY(region).name(5:end)};
end

if size(A_matrix,1) == size(A_matrix,2) && ...
        size(A_matrix,1) == length({DCM{1}.xY.name})
    set(gca,'YTickLabel',regions,...
        'YTick',1:length({DCM{1}.xY.name}),...
        'XTickLabel',regions,'fontweight','bold','fontsize',40,'XTick',...
        1:length({DCM{1}.xY.name}),'TickLabelInterpreter', 'none');
end


%%%%%%%%%%%%%%%%%%%
%specify colormap
%%%%%%%%%%%%%%%%%%%
L=length(A_matrix)^2;
indexValue = 0;     % value for which to set a particular color

topColor = [1 0 0];         % color for maximum data value (red = [1 0 0])
indexColor = [1 1 1];       % color for indexed data value (white = [1 1 1])
bottomcolor = [0 0 1];      % color for minimum data value (blue = [0 0 1])

% Calculate where proportionally indexValue lies between minimum and
% maximum values
largest = max(max(A_matrix));
smallest = min(min(A_matrix));
index = L*abs(indexValue-smallest)/(largest-smallest);

% Create color map ranging from bottom color to index color
% Multipling number of points by 100 adds more resolution
customCMap1 = [linspace(bottomcolor(1),indexColor(1),100*index)',...
            linspace(bottomcolor(2),indexColor(2),100*index)',...
            linspace(bottomcolor(3),indexColor(3),100*index)'];

% Create color map ranging from index color to top color
% Multipling number of points by 100 adds more resolution
customCMap2 = [linspace(indexColor(1),topColor(1),100*(L-index))',...
            linspace(indexColor(2),topColor(2),100*(L-index))',...
            linspace(indexColor(3),topColor(3),100*(L-index))'];

customCMap = [customCMap1;customCMap2];  % Combine colormaps

colormap(customCMap);
c=colorbar;
% 
%%%%%%%%%%%%%%%%%%%
A_matrix_NaNs=A_matrix;
A_matrix_NaNs(sgn.*(ones(3,3)-diag([1 1 1]))==-1)=NaN;

for side1=1:3
    for side2=1:3
        if ~isnan(A_matrix_NaNs(side2,side1))||(side1==side2)
            text(side1,side2,num2str(round(A_matrix(side2,side1),2)),'HorizontalAlignment','center','Fontweight','bold','Fontsize',40);
%         elseif A_matrix(side2,side1)==Inf;
%             text(side1,side2,'','HorizontalAlignment','center','Fontweight','bold','Fontsize',18);
%         elseif side2~=side1&&isnan(A_matrix(side2,side1))&&side2<6&&side1<6
%             text(side1,side2,'','HorizontalAlignment','center','Fontweight','bold','Fontsize',18);
%         elseif side2~=side1&&isnan(A_matrix(side2,side1))&&side2>5&&side2<9&&side1>5&&side1<9
%             text(side1,side2,'','HorizontalAlignment','center','Fontweight','bold','Fontsize',18);
%         elseif side2~=side1&&isnan(A_matrix(side2,side1))&&side2>8&&side1>8
%             text(side1,side2,'','HorizontalAlignment','center','Fontweight','bold','Fontsize',18);
        elseif isnan(A_matrix_NaNs(side2,side1))
            text(side1,side2,num2str(round(A_matrix(side2,side1),2)),'HorizontalAlignment','center','Fontweight','normal','Fontsize',20);
%         elseif side2==side1&&isnan(A_matrix(side2,side1))
%             text(side1,side2,'','HorizontalAlignment','center','Fontweight','bold','Fontsize',18);
        end
    end
end

for side1=1:3
    for side2=1:3
%         if side1==side2
%             continue;
            
        if isnan(A_matrix(side1,side2))
            h.AlphaData(side1,side2)=0;
        end
    end
end

c.Limits=[-0.6 0.6];

set(gcf,'PaperPositionMode','auto');

direct=[Work_dir '/Figures_paper_GSR/PEB_Figures/GSR_regr/' name_ROI_def '/Full_model/Figure_influence_connect'];
saveas(gcf,[direct '/Single_A_SMR_matrix_GSR_regr_comp_HCP.bmp']);
close;

end