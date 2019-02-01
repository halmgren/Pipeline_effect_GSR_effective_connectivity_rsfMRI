function Mesh_DMN_paper_GSR(SPM_dir,Work_dir)

%%%%%%%%%%%%
%Adapted version of spm_dcm_graph
%%%%%%%%%%%%
direct=[Work_dir '/Figures_paper_GSR/PEB_Figures/Mesh_network_Figure'];
mkdir(direct);


m     = 4;

ROI_list=Define_ROIs_paper_GSR('Smith');

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


colopt   = {'b','g','r','c','m','y','k','w'};

% for network_number=1:length(ntwrk_size)
    
    L     = [];
%     
%     for nd=1:ntwrk_size(network_number)
%         
%         L=[L ROI_list{sum(ntwrk_size(1:network_number-1))+nd,2}'];
%         
%         name{nd}=ROI_list{sum(ntwrk_size(1:network_number-1))+nd,1}(5:end);
%         
%     end
    
tmp2=0;
    for network_number2=1:length(ntwrk_size)
    
        for nd=1:ntwrk_size(network_number2)
            tmp2=tmp2+1;
            L=[L ROI_list{sum(ntwrk_size(1:network_number2-1))+nd,2}'];
        
            name{tmp2}=ROI_list{sum(ntwrk_size(1:network_number2-1))+nd,1}(5:end);
            col(tmp2)=colopt{network_number2};
        end
    end

    
    m     = size(L,2);
    
    
    %-Render graph in anatomical space
    %==========================================================================
    h1=figure('units','normalized','outerposition',[0 0 1 1]);%subplot(1,1,1);
    ax=gca;
    
    cla(ax);
    %set(ax,'position',[0 .5 1 .5])
    options.query = [];
    options.hfig  = ancestor(ax,'figure');
    options.ParentAxes = ax;
    options.markersize = 32;
    options.meshsurf = fullfile(spm('Dir'),'canonical','iskull_2562.surf.gii');
    Pos=L(:,1);
    Orient=[];
    Var=0;
    Names=[];
    % spm_eeg_displayECD(L(:,1),[],0,[],options);
    
    %%%%%%%%%%%%%%
    %Cut from spm_eeg_displayECD
    %%%%%%%%%%%%%%
    
    options.meshsurf = fullfile(spm('Dir'),'canonical','cortex_8196.surf.gii');
    Pos=L;
    Orient=[];
    Var=8;
    Names=name;
    
    hfig       = [];
    ParentAxes = [];
    query      = [];
    handles    = [];
    tag        = '';
    try, options; catch, options = [];      end
    try, hfig        = options.hfig;        end
    try, tag         = options.tag;         end
    try, ParentAxes  = options.ParentAxes;  end
    try, query       = options.query;       end
    try, handles     = options.handles;     end
    try
        figure(hfig);
    catch
        hfig  = spm_figure('GetWin','Graphics');
        spm_figure('Clear',hfig);
        ParentAxes = axes('parent',hfig);
    end
    try
        markersize = options.markersize;
    catch
        markersize = 20;
    end
    try
        meshsurf = options.meshsurf;
    catch
        meshsurf = fullfile(spm('Dir'),'canonical','cortex_5124.surf.gii');
    end
    
    if isscalar(Var), Var = Pos*0 + Var^2;   end
    try, Pos{1};    catch, Pos = {Pos};      end
    try, Orient{1}; catch, Orient = {Orient};end
    try, Var{1};    catch, Var = {Var};      end
    
    ndip = size(Pos{1},2);
    if ~exist('Names','var') || isempty(Names)
        for i=1:ndip
            Names{i} = num2str(i);
        end
    end
    
    
%     col = ['b','g','r','y','c','m','k','w'];
    tmp = ceil(ndip./numel(col));
    col = repmat(col,1,tmp);
    pa  = get(ParentAxes,'position');
    hold on;
    if ndip > 0
        
        if isempty(query)
            opt.hfig = hfig;
            opt.ParentAxes = ParentAxes;
            opt.visible = 'on';
            pos2 = [pa(1),pa(2)+0.25*pa(4),0.03,0.5*pa(4)];
            out  = spm_eeg_render(meshsurf,opt);
            out.handles.ParentAxes.CameraPosition=[0 -100 300];
            handles.mesh = out.handles.p;
            handles.BUTTONS.transp = out.handles.transp;
            handles.hfig = out.handles.fi;
            handles.ParentAxes  = out.handles.ParentAxes;
            
            %         set(handles.mesh,...
            %             'facealpha',0.4,...
            %             'visible','on')
        end
        
        set(ParentAxes,'nextplot','add')
        
        for j=1:length(Pos)
            for i =1:ndip
                
                
                
                %INNER CIRCLE (darker)
                
                handles.hp(j,i) = plot3(handles.ParentAxes,...
                    Pos{j}(1,i),Pos{j}(2,i),Pos{j}(3,i),...
                    [col(i),'.'],...
                    'markerSize',markersize,...
                    'visible','on');
                
                [x,y,z]= ellipsoid(Pos{j}(1,i),Pos{j}(2,i),Pos{j}(3,i),...
                    6,6,6,20);
                
                %OUTER CIRCLE (lighter)
                handles.hs(j,i) = surf(handles.ParentAxes,...
                    x,y,z,...
                    'edgecolor','none',...
                    'facecolor',col(i),...
                    'facealpha',1,...
                    'visible','on');
                
                %ROI NAMES
                %
                axHidden = axes('Visible','off','hittest','on'); % Invisible axes
                
                linkprop([handles.ParentAxes axHidden],{'CameraPosition' 'CameraUpVector' 'CameraViewAngle' 'CameraTarget' 'XLim' 'CLim' 'YLim' 'ZLim' 'Position' 'DataAspectRatio'});
                
                %             handles.inv_axes=handles.ParentAxes;
                %             set(handles.inv_axes,'visible','off');
%                 if i~=2
%                     handles.ht(j,i) = text(...
%                         Pos{j}(1,i),Pos{j}(2,i),Pos{j}(3,i),...
%                         Names{i},...
%                         'Parent',handles.ParentAxes,...
%                         'visible','on','color','k','FontSize',44,'FontWeight','bold');
% % %                 else
                if Pos{j}(1,i)<-10
                    handles.ht(j,i) = text(...
                        Pos{j}(1,i)-17.5,Pos{j}(2,i)-4,Pos{j}(3,i),...
                        Names{i},...
                        'Parent',handles.ParentAxes,...
                        'visible','on','color','k','FontSize',20,'FontWeight','bold');
                elseif Pos{j}(1,i)>10
                    handles.ht(j,i) = text(...
                        Pos{j}(1,i)+17.5,Pos{j}(2,i)-4,Pos{j}(3,i),...
                        Names{i},...
                        'Parent',handles.ParentAxes,...
                        'visible','on','color','k','FontSize',20,'FontWeight','bold');
                elseif Pos{j}(2,i)>54
                    handles.ht(j,i) = text(...
                        Pos{j}(1,i),Pos{j}(2,i)+12.5,Pos{j}(3,i),...
                        Names{i},...
                        'Parent',handles.ParentAxes,...
                        'visible','on','color','k','FontSize',20,'FontWeight','bold');
                else
                    handles.ht(j,i) = text(...
                        Pos{j}(1,i),Pos{j}(2,i)-12.5,Pos{j}(3,i),...
                        Names{i},...
                        'Parent',handles.ParentAxes,...
                        'visible','on','color','k','FontSize',20,'FontWeight','bold');
                end
% %                 end

                    
                %             handles.ht(i).FontSize=15;
                handles.ht(i).HorizontalAlignment='center';
                
                set(handles.ht(i),'Parent',axHidden);
                
            end
        end
        
        try, set(handles.hp(1,:),'visible','on'); end
        try, set(handles.hq(1,:),'visible','on'); end
        try, set(handles.hs(1,:),'visible','on'); end
        try, set(handles.ht(1,:),'visible','on'); end
        
        
    end
    
    try
        clear out
        out.handles = handles;
    catch
        out = [];
    end
    
    
%     set(out.handles.ht,'FontWeight','bold');
    set(out.handles.mesh,'FaceAlpha',1/8);
    
    set(gcf,'PaperPositionMode','auto');
    saveas(gcf,[direct '/Mesh_regions.bmp']);
    close;
% end

end