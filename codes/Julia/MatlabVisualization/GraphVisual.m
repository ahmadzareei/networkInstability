%visulazation of 3d networks

set(0,'defaultAxesFontSize',10)
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
fs = 12;
figFolder = fullfile("E:","TempCode","MatlabFlow","Figs","Voro3d");
a = 0.02;
T = 12001;
% dir = fullfile("E:","TempCode","MatlabFlow","matData","erosion","DelaunayNet",strcat("200by100T",num2str(T)));
dir = fullfile("E:","TempCode","MatlabFlow","matData","erosion","VoronoiNet",strcat("20by12by12T",num2str(T),"d0.20"));
for N=2:5
    close all


    subdir = fullfile(dir,strcat('N',num2str(N,'%.1f')),strcat('a',num2str(a,'%.2f')));
    timeData = fullfile(subdir,strcat('matLargeDataS',num2str(1),'.mat'));
    orderData = fullfile(subdir,'A_Matlab_Data.mat');
    result = isfile(timeData);
    if ~result
        subdir = fullfile(dir,strcat('N',num2str(N,'%.2f')),strcat('a',num2str(a,'%.1f')));
        orderData = fullfile(subdir,'A_Matlab_Data.mat');
        timeData = fullfile(subdir,strcat('matLargeDataS',num2str(1),'.mat'));
    end
    timeData = fullfile(subdir,strcat('matLargeDataS',num2str(1),'.mat'));
    posData = fullfile(dir,strcat('configArrayS',num2str(1),'.mat'));
    graphData = fullfile(dir,strcat('ST',num2str(1),'.mat'));

    load(graphData);
    load(orderData);
    load(timeData);
    load(posData);

    A = graph(s,t);
    h = plot(A,'XData',posArray(:,1),'YData',posArray(:,2),'ZData',posArray(:,3)); %
    h.NodeLabel='';
    h.Marker='none';

    axis equal;
    % axis([5 15 0 12 0 12]);

    Rave_t = mean(transpose(R_t));
    index = find(Rave_t > 15);
    index = index(1);
    edgeWidth = log(1+R_t(index,:)/mean(R_t(index,:)));
    edgeWidth = R_t(index,:)/mean(R_t(index,:));
    h.LineWidth = edgeWidth;
    % h.LineWidth = 1;
    h.EdgeCData= log(abs(WeightP_t(index,:))+1);
    caxis([0,15])
    h.EdgeAlpha = 0.3;
    J = customcolormap_preset("white-blue-red");
    colormap(J);
    % colormap(bluewhitered(256));
    % colorbar;
    exportgraphics(gca,fullfile(figFolder,strcat('G_N_',num2str(N),'.pdf')),'BackgroundColor','none')
end

% 
% hold on;
% C_t = R_t.^4;
% E_t = WeightP_t.^2./C_t;
% E_t_s = E_t./(Outflow_t.^2);
% plot(mean(transpose(R_t)),sum(transpose(E_t_s)));