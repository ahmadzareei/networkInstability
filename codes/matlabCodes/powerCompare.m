%plot different m and n's network simulation results
set(0,'defaultAxesFontSize',10)
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
% set(gca,'color','none')
fs = 12;
figFolder = fullfile("E:","TempCode","MatlabFlow","Figs","MNCompare2");
a = 0.2;
T = 6001;
% dir = fullfile("E:","TempCode","MatlabFlow","matData","erosion","DelaunayNet",strcat("200by100T",num2str(T)));
dir = fullfile("E:","TempCode","MatlabFlow","matData","erosion","DelaunayNet",strcat("50by50T",num2str(T),"d0.2"));
subdir = fullfile(dir,strcat('N',num2str(1,'%.1f')),strcat('M',num2str(0,'%.1f')),strcat('a',num2str(a,'%.1f')));
ST = fullfile(subdir,'ST1.mat');
posData = fullfile(subdir,strcat('configArrayS',num2str(1),'.mat'));
load(ST);
load(posData);
A = graph(s,t);

Ns = [1.5 2 2.5 3 3.5 4 4.5];
% Ns = flip(Ns);
Ms = [0 0.25 0.5 0.75 1.0 1.25 1.5 1.75 2];
% Ns = [1.5 5];
% Ms = [0 2];
ns = length(Ns);
ms = length(Ms);
width = 1/ms;
height = 1/ns;
% fig = figure;
% axis off
% box off
% hold on 
% scatter( -1, -1, 0.01);
% xlim([0,2]);
% ylim([1.5,5]);
% ha = tight_subplot(ns,ms,[.0 -.02],[0.08,0.0],[0.14,0.04]);
for i = 1:1
    for j = 1:1
        close all 
        figure('Position',[0,0,300,300])
        set(gcf, 'color', 'none');
        set(gca, 'color', 'none');
%         ax = axes('Position',[width*(j-1) height*(i-1) width height]);
%         ax.PositionConstraint = 'innerposition';
%         axes(ha((i-1)*ms+j));
        k = ms*(i-1)+j;
        subdir = fullfile(dir,strcat('N',num2str(Ns(i),'%.1f')),strcat('M',num2str(Ms(j),'%.1f')),strcat('a',num2str(a,'%.1f')));
        timeData = fullfile(subdir,strcat('matLargeDataS',num2str(1),'.mat'));
        result = isfile(timeData);
        if ~result
            subdir = fullfile(dir,strcat('N',num2str(Ns(i),'%.1f')),strcat('M',num2str(Ms(j),'%.2f')),strcat('a',num2str(a,'%.1f')));
            timeData = fullfile(subdir,strcat('matLargeDataS',num2str(1),'.mat'));
        end
        load(timeData);
        Rave_t = mean(transpose(R_t));
%         [Rmax,index] = max(Rave_t);
        index = find(Rave_t > 15);
        if isempty(index)
            [Ns(i) Ms(j)]
%             Rave_t(end)
%             continue
            [Rmax,index] = max(Rave_t);
            Rmax
        else
            index = index(1);
        end
        index = 1;
%         axes(ha(k));
        h = plot(A,'XData',posArray(:,1),'YData',posArray(:,2)); %
        h.NodeLabel='';
        h.Marker='none';

%         axis([0 40 -5 40]);
        

        edgeWidth = log(1+R_t(index,:)/mean(R_t(index,:)));
        h.LineWidth = edgeWidth;
        h.EdgeCData= abs(WeightP_t(index,:))/mean(abs(WeightP_t(index,:)));
        caxis([0,10]);
        h.EdgeAlpha = 1;
        J = customcolormap_preset('white-blue-red');
        colormap(J);
        
%         set(colorbar,'position',[0.88 .17 .03 .7])
        axis equal
        axis([5 45 5 45]);
        box off
        axis off
        ax = gca;
%         exportgraphics(gca,fullfile(figFolder,strcat('G_N_',num2str(Ns(i)),'_M_',num2str(Ms(j)),'initial.emf')),'BackgroundColor','none');
%         saveas(gcf,fullfile(figFolder,strcat('G_N_',num2str(Ns(i)),'_M_',num2str(Ms(j)),'_test.pdf'))) % ,'epsc');
%         if Ms(j)>((1+Ns(i))/4)
%             edgeWidth = log(1+R_t(index,:)/mean(R_t(index,:)));
%             h.LineWidth = edgeWidth;
% %             h.LineWidth = 3;
% %             title(strcat('N',num2str(Ns(i),'%.1f'),'M',num2str(Ms(j),'%.2f')),'Color','blue','FontSize',5);
%         else
%             h.LineWidth = 1;
% %             title(strcat('N',num2str(Ns(i),'%.1f'),'M',num2str(Ms(j),'%.2f')),'Color','red','FontSize',5);
%         end
        
%         colorbar;
    end
end

ax = axes;
J = customcolormap_preset('white-blue-red');
c = colormap(J);
colorbar;
caxis([0,10]);
ax.Visible = 'off';

% han=axes(fig,'visible','off'); 
% han.XLim = [0,2];
% han.YLim =[1.5,5] ;
% han.Title.Visible='on';
% han.XLabel.Visible='on';
% han.YLabel.Visible='on';
% ylabel(han,'$m$','Interpreter','latex');
% xlabel(han,'$n$','Interpreter','latex');
% % title(han,'2D Erosion');