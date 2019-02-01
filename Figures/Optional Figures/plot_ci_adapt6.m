function [hh xb]=plot_ci_adapt6(E,C,x,j,s,barcolor,width,sz,GCM,grd,colored,cntr,colored_2)
%adapted script from spm_plot_ci


% get axis
%--------------------------------------------------------------------------
ax = gca;

% confidence region plotting
%--------------------------------------------------------------------------
if size(E,1) == 1 && size(E,2) == 2
    E  = E';
    CR = true;
else
    CR = false;
end
    

% unpack expectations into a matrix
%--------------------------------------------------------------------------
if isstruct(E),       E = spm_vec(E);    end
if iscell(E),         E = spm_vec(E);    end

if ~exist('x','var'), x = 1:size(E,2);   end
if ~exist('j','var'), j = 1:size(E,1);   end
if ~exist('s','var'), s = '';            end

if isempty(x),        x = 1:size(E,2);   end
if isempty(j),        j = 1:size(E,1);   end


% order and length of sequence
%--------------------------------------------------------------------------
O     = E;
% E     = E(j,:);
[n,N] = size(E);

% unpack conditional covariances
%--------------------------------------------------------------------------
ci    = spm_invNcdf(1 - 0.05);               % confidence interval
gr    = 0.9;                                 % grey level
% if iscell(C)
%     
%     % try cell array of covariances (from spm_DEM amd spm_LAP)
%     %----------------------------------------------------------------------
%     try
%         for i = 1:N
%             c(:,i) = ci*sqrt(diag(C{i}(j,j)));
%         end
%     catch
%         
%         % try cell array of variances
%         %------------------------------------------------------------------
%         c = ci*sqrt(spm_unvec(spm_vec(C),O));
%         c = c(j,:);
%     end
%     
% elseif isstruct(C)
%     
%     % try structure of variances
%     %----------------------------------------------------------------------
%     c = ci*sqrt(spm_unvec(spm_vec(C),O));
%     c = c(j,:);
%     
% elseif isnumeric(C)
%     
%     % try matrix of variances
%     %----------------------------------------------------------------------
%     if all(size(C) == size(O))
%         c = ci*sqrt(C(j,:));
%     elseif all(size(C') == size(O))
%         c = ci*sqrt(C(:,j));
%     else
%         
%         % try covariance matrix
%         %------------------------------------------------------------------
%         C = diag(C);
%         c = ci*sqrt(C(j,:));
%     end
%     
% end

c=ci*sqrt(C);

% set plot parameters
%--------------------------------------------------------------------------
switch lower(get(ax,'NextPlot'))
    case 'add'
%         col   = [1 1/4 1/4];
%         width = .9;
    otherwise
%         col   = [1 3/4 3/4];
%         width = .8;
end

% plot elliptical confidence region
%--------------------------------------------------------------------------
if CR
    [x,y] = ellipsoid(E(1),E(2),1,c(1),c(2),0,32);
    fill(x(16,:)',y(16,:)',[1 1 1]*gr,'EdgeColor',[1 1 1]*.5,'Parent',ax);
    hold(ax,'on');
    plot(ax,E(1),E(2),'.','MarkerSize',16);
    hold(ax,'off'); drawnow
    return
end


% plot bar chart
%--------------------------------------------------------------------------
% if N >= 8
%     
%     % time-series plot
%     %======================================================================
%     if strcmpi(s,'exp')
%         fill([x fliplr(x)],exp([full(E + c) fliplr(full(E - c))]),...
%             [.95 .95 1],'EdgeColor',[.8 .8 1],'Parent',ax);
%         hold(ax,'on');
%         plot(x,exp(E));
%         
%     else
%         fill([x fliplr(x)],[full(E + c) fliplr(full(E - c))],...
%             [.95 .95 1],'EdgeColor',[.8 .8 1],'Parent',ax);
%         hold(ax,'on');
%         plot(ax,x,E,s);
%     end    
%     
% else
    
    
    % bar
    %======================================================================
%     if N == 1
        
        if strcmpi(s,'exp')
            
            % conditional means
            %--------------------------------------------------------------
            bar(ax,exp(E),width,'Edgecolor',[1 1 1]/2,'Facecolor',[1 1 1]*.8);
            hold(ax,'on');
            
            % conditional variances
            %--------------------------------------------------------------
            for k = 1:n
                line([k k],exp([-1 1]*c(k) + E(k)),'LineWidth',4,'Color',col,'Parent',ax);
            end
            
        else
                        
            if n > 1
                
                % conditional means
                %----------------------------------------------------------
                tmp=0;
                tmp2=0;
                for regions=1:sz
                    if regions>1
                        tmp2=tmp2+1;
                    end
                    
                    for regions2=1:sz
                        tmp=tmp+1;
%                         for method=1:2
%                             if method==1
                                hh{(regions-1)*sz+regions2}=bar(ax,[tmp+tmp2:sz^2+sz],[E(:,tmp) zeros(2,sz^2+sz-tmp-tmp2)]',width,'Edgecolor',[1 1 1]/2,'Facecolor',barcolor,'BarWidth',0.02);
%                             end
%                         end
                        hold on;
                    end
                    
                end
                hold on;
                
                %%%%%%%%%%%%%%%%%%%%%%%%%
                %lines wihtin bar graphs
                %%%%%%%%%%%%%%%%%%%%%%%%%
                if grd==1
                    for connections=1:length(hh)
                        %             h(barclr).FaceColor=barcolor*barclr/ceil(length(h)/2);
                        %                 if all(barcolor==[0.5 0.5 0.5])
                        for method=1:2
                            if method==1
                                barclr=rem(connections-1,sqrt(length(hh)))+1;
                                hh{connections}(method).FaceColor=[0.2+0.6/sqrt(length(hh))*barclr 0.2+0.6/sqrt(length(hh))*barclr 0.2+0.6/sqrt(length(hh))*barclr];
                                hh{connections}(method).EdgeColor=[0 0 0];
                            elseif method==2
                                barclr=rem(connections-1,sqrt(length(hh)))+1;
                                
                                hh{connections}(method).FaceColor=[0.2+0.6/sqrt(length(hh))*barclr 0.2+0.6/sqrt(length(hh))*barclr 0.2+0.6/sqrt(length(hh))*barclr];
                                hh{connections}(method).EdgeColor=[0 0 0];
                                
                                POS=hh{connections}(1).XData(1);
                                
                                if E(2,connections)<0
                                    
                                    vertic_loc=E(2,connections);
                                    while vertic_loc<0
                                        vertic_loc=vertic_loc+0.02;
                                        
                                        if vertic_loc>0
                                            %                                 Y=(-vertic_loc+0.02+0.0373*-(0.129+0.25))/0.0373;
                                            Y=0.25/0.02*(0-(vertic_loc-0.02))+POS+0.1429-0.25/2;
                                            line([POS+0.1429-0.25/2 Y],[vertic_loc-0.02 0],'Color','k','LineWidth',1.5);
                                            break
                                        end
                                        
                                        line([POS+0.1429-0.25/2 POS+0.1429+0.25/2],[vertic_loc-0.02 vertic_loc],'Color','k','LineWidth',1.5);
                                        
                                        
                                    end
                                end
                                
                                
                                if E(2,connections)>0
                                    
                                    vertic_loc=0;
                                    while vertic_loc<E(2,connections)
                                        
                                        vertic_loc=vertic_loc+0.02;
                                        
                                        if vertic_loc>E(2,connections)
                                            Y=0.25/0.02*(E(2,connections)-(vertic_loc-0.02))+POS+0.1429-0.25/2;
                                            line([POS+0.1429-0.25/2 Y],[vertic_loc-0.02 E(2,connections)],'Color','k','LineWidth',1.5);
                                            break
                                        end
                                        
                                        line([POS+0.1429-0.25/2 POS+0.1429+0.25/2],[vertic_loc-0.02 vertic_loc],'Color','k','LineWidth',1.5);
                                        
                                    end
                                end
                                %                     if barclr==1
                                %                         hh{connections}(method).FaceColor=[0 0 1];
                                % %                     elseif barclr==2
                                % %                         hh{connections}(method).FaceColor=[0 0 1];
                                %                     else
                                %                         hh{connections}(method).FaceColor=[0 0+1/(sqrt(length(hh))-(sz-4))*(barclr-(sz-4)) 1];
                                %                     end
                                
                                hh{connections}(method).EdgeColor=[0 0 0];
                            end
                        end
                        
                        %                 end
                        %
                        %                 if strcmp(barcolor,'b')
                        %                     hh(barclr).FaceColor=[0 0 0.1+0.9/length(h)*barclr];
                        %                 end
                    end
                elseif colored==1
                    for connections=1:length(hh)
                        %             h(barclr).FaceColor=barcolor*barclr/ceil(length(h)/2);
                        %                 if all(barcolor==[0.5 0.5 0.5])
                        for method=1:2
                            if method==1
                                barclr=rem(connections-1,sqrt(length(hh)))+1;
                                hh{connections}(method).FaceColor=[0.2+0.6/sqrt(length(hh))*barclr 0.2+0.6/sqrt(length(hh))*barclr 0.2+0.6/sqrt(length(hh))*barclr];
                                hh{connections}(method).EdgeColor=[0 0 0];
                            elseif method==2
                                
                                                    if barclr==1
                                                        hh{connections}(method).FaceColor=[0 0 1];
                                %                     elseif barclr==2
                                %                         hh{connections}(method).FaceColor=[0 0 1];
                                                    else
                                                        hh{connections}(method).FaceColor=[0 0+1/(sqrt(length(hh))-(sz-4))*(barclr-(sz-4)) 1];
                                                    end
                                
                                hh{connections}(method).EdgeColor=[0 0 0];
                            end
                        end
%                         
%                         
%                         if strcmp(barcolor,'b')
%                             hh(barclr).FaceColor=[0 0 0.1+0.9/length(h)*barclr];
%                         end
                    end
                
                    
                elseif cntr==1
                    
                    for connections=1:length(hh)
                        for method=1:2
                            if method==1
                                barclr=rem(connections-1,sqrt(length(hh)))+1;
                                hh{connections}(method).FaceColor=[0.2+0.6/sqrt(length(hh))*barclr 0.2+0.6/sqrt(length(hh))*barclr 0.2+0.6/sqrt(length(hh))*barclr];
                                hh{connections}(method).EdgeColor=[0 0 0];
                                hh{connections}(method).LineWidth=1.25;
                            elseif method==2
                                
                               barclr=rem(connections-1,sqrt(length(hh)))+1;
                                hh{connections}(method).FaceColor=[0.2+0.6/sqrt(length(hh))*barclr 0.2+0.6/sqrt(length(hh))*barclr 0.2+0.6/sqrt(length(hh))*barclr];
                                hh{connections}(method).EdgeColor=[0.65 0.65 0.65];
                                hh{connections}(method).LineWidth=1.25;
                            end
                        end
                    end
                elseif colored_2==1
                    for connections=1:length(hh)
                        for method=1:2
                            barclr=rem(connections-1,sqrt(length(hh)))+1;
                            if method==1
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
                                
                                hh{connections}(method).FaceColor=clr_basic_2;
                                hh{connections}(method).EdgeColor=[0 0 0];
                                hh{connections}(method).LineWidth=2;
                            elseif method==2
                                
                                if barclr==1
                                    clr_GSR_2=[0.5 0.7 1]; %lightblue
                                elseif barclr==2
                                    clr_GSR_2=[1 0.5 0.5]; %lightred
                                elseif barclr==3
                                    clr_GSR_2=[0.5 0.5 0.5]; %grey
                                elseif barclr==4
                                    clr_GSR_2=[0.5 1 0.5]; %green
                                elseif barclr==5
                                    clr_GSR_2=[1 0.5 1]; %magenta
                                end
                                
                                hh{connections}(method).FaceColor=clr_GSR_2;
                                hh{connections}(method).EdgeColor=[0 0 0];
                                hh{connections}(method).LineWidth=2;
                            end
                        end
%                         
%                         
%                         if strcmp(barcolor,'b')
%                             hh(barclr).FaceColor=[0 0 0.1+0.9/length(h)*barclr];
%                         end
                    end
                end
                r=refline(0,0);
                set(r,'LineWidth',3,'Color','k');
                
                
                %%%%%%%%%%%%%%%%%%%%%%%
                %ERRORBARS
                %%%%%%%%%%%%%%%%%%%%%%%
                hold on;
                %             hold(ax,'off');
                for connections=1:length(hh)
                    xb = bsxfun(@plus, hh{connections}(1).XData, [hh{connections}.XOffset]');
                    %             ax1 = axes('Visible','off','hittest','on');
                    %             linkprop([ax,ax1],{'XLim','YLim','Position'});
                    %                 E=spm_vec(E');
                    for method=1:2
                        
%                         if connections>1
                            if method==1
%                                  line([xb(method,1)-0.1429 xb(method,1)-0.1429],[-1 1]*c(method,connections)' + E(method,connections)','LineWidth',4,'Color',[1 0 0],'Parent',ax,'Visible','on');
                                    line([xb(method,1) xb(method,1)],[-1 1]*c(method,connections)' + E(method,connections)','LineWidth',4,'Color',[1 0 0],'Parent',ax,'Visible','on');
                            elseif method==2
%                                 line([xb(method,1)+0.1429 xb(method,1)+0.1429],[-1 1]*c(method,connections)' + E(method,connections)','LineWidth',4,'Color',[1 0 0],'Parent',ax,'Visible','on');
                                    line([xb(method,1) xb(method,1)],[-1 1]*c(method,connections)' + E(method,connections)','LineWidth',4,'Color',[1 0 0],'Parent',ax,'Visible','on');
                                %                 set(handles.ax1,'Parent',ax1)
                                %                 set(handles.ax1,'visible','on');
                            end
                        
%                         else
%                            line([xb(method,1) xb(method,1)],[-1 1]*c(method,connections)' + E(method,connections)','LineWidth',4,'Color',[1 2/4 2/4],'Parent',ax,'Visible','on'); 
                            
%                         end
                    end
                    
                    
                end
                
                hold(ax,'on');
                
                
                
                % conditional means
                %----------------------------------------------------------
%                 bar(ax,E,'Edgecolor',[1 1 1]/2,'Facecolor','b');
%                 hold(ax,'on');
                
            end
            
            % conditional variances
            %--------------------------------------------------------------
        
            
            
            
        end
        
        
        
        
        Xl=xlabel('Connectivity From');
        Yl=ylabel('Posterior Expectation (Hz)');
        
        xlim([0 sz^2+sz]);
        
        ax.XTick=[];
        
        for region=1:length(GCM.xY)
            ax.XTick(region)=[ceil(sz/2)+(region-1)*(sz+1)];
            ax.XTickLabel{region}=GCM.xY(region).name(5:end);
        end
        
        ax.XAxis.FontSize=30;
        ax.XAxis.FontWeight='bold';
        
        ax.YAxis.FontSize=36;
        ax.XAxis.FontWeight='bold';
        Xl.FontSize=30;
        Xl.FontWeight='bold';
        Yl.FontSize=30;
        Yl.FontWeight='bold';
        
        box(ax,'off');
        
        YL=ylim;
        
        for region=1:length(GCM.xY)-1
            line([(sz+1)*region (sz+1)*region],[YL(1) YL(2)],'Color',[0.7 0.7 0.7],'LineStyle','--','LineWidth',2);
        end
%         set(ax,'XLim',[0 n + 1]);
        
%     else
        
        
%         % conditional means
%         %------------------------------------------------------------------
%         h = bar(ax,E); hold(ax,'on');
%         
%         % conditional variances
%         %------------------------------------------------------------------
%         for m = 1:N
%             x = mean(get(get(h(m),'Children'),'Xdata'));
%             for k = 1:n
%                 line([x(k) x(k)],[-1 1]*c(k,m) + E(k,m),'LineWidth',4,'Color',col,'Parent',ax);
%             end
%         end
%     end
%     
% end

drawnow
