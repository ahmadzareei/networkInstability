%plot change of average order parameters vs different m and n's

set(0,'defaultAxesFontSize',10)
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
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
Odata = [];
Rave_data = [];
% ha = tight_subplot(ns,ms,[.0 -.02],[0.08,0.0],[0.14,0.04]);
seeds = 10;
for seed = 1:seeds
    Odata_ = [];
    Rave_data_ = [];
    for i = 1:ns
        orderData = [];
        Rave_d = [];
        for j = 1:ms
    %         close all 
    %         figure('Position',[0,0,300,300])
    %         set(gcf, 'color', 'none');
    %         set(gca, 'color', 'none');
    %         ax = axes('Position',[width*(j-1) height*(i-1) width height]);
    %         ax.PositionConstraint = 'innerposition';
    %         axes(ha((i-1)*ms+j));
            k = ms*(i-1)+j;
            subdir = fullfile(dir,strcat('N',num2str(Ns(i),'%.1f')),strcat('M',num2str(Ms(j),'%.1f')),strcat('a',num2str(a,'%.1f')));
            timeData = fullfile(subdir,strcat('matLargeDataS',num2str(seed),'.mat'));
            result = isfile(timeData);
            if ~result
                subdir = fullfile(dir,strcat('N',num2str(Ns(i),'%.1f')),strcat('M',num2str(Ms(j),'%.2f')),strcat('a',num2str(a,'%.1f')));
                timeData = fullfile(subdir,strcat('matLargeDataS',num2str(seed),'.mat'));
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
%                 Rave_t(index)
            end
            Ne = size(WeightP_t,2);
            Qs = transpose(WeightP_t);
            orderParas = 1/(Ne-1)*(Ne - sum(Qs.^2).^2./sum(Qs.^4));
            orderData = [orderData orderParas(index)-orderParas(1)];
            Rave_d = [Rave_d Rave_t(index)];
    %         axes(ha(k));
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
        Odata_ = [Odata_; orderData];
        Rave_data_ = [Rave_data_;Rave_d];
    end
    if seed == 1
        Odata = Odata_;
        Rave_data = Rave_data_;
    else
        Odata = Odata + Odata_;
        Rave_data = Rave_data + Rave_data_;
    end
end
Odata = Odata./seeds;
Rave_data = Rave_data/seeds;
figure;
h=heatmap(Ms,Ns,Odata,'Colormap',customcolormap_preset('red-white-blue'),'ColorLimits',[-0.3 0.3]);
figure;
h2=heatmap(Ms,Ns,Rave_data,'Colormap',customcolormap_preset('red-white-blue'));

