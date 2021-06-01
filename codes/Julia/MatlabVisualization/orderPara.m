%plot order parameters vs average radius
set(0,'defaultAxesFontSize',10)
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
fs = 12;
figFolder = fullfile("E:","TempCode","MatlabFlow","Figs","OrderPara","2d");
a = 0.2;
T = 8001;
% dir = fullfile("E:","TempCode","MatlabFlow","matData","erosion","DelaunayNet",strcat("200by100T",num2str(T)));
dir = fullfile("E:","TempCode","MatlabFlow","matData","erosion","DelaunayNet",strcat("100by50T",num2str(T),"d0.2"));
R_c = 0;
Ns = [];
Rs = [];
Os = [];
indices = [1 2 3 4];
Os1 = [];
Os0 = [];
Os2 = [];
Os3 = [];
Os4 = [];
for N = 1:0.2:6
    subdir = fullfile(dir,strcat('N',num2str(N,'%.1f')),strcat('a',num2str(a,'%.1f')));
    orderData = fullfile(subdir,'A_Matlab_Data.mat');
    result = isfile(orderData);
    if ~result
        subdir = fullfile(dir,strcat('N',num2str(N,'%.2f')),strcat('a',num2str(a,'%.1f')));
        orderData = fullfile(subdir,'A_Matlab_Data.mat');
    end
    timeData = fullfile(subdir,strcat('matLargeDataS',num2str(1),'.mat'));
    load(orderData);
    load(timeData);
    Rave_t = mean(transpose(R_t));
    for i = 2:1:10
        timeData = fullfile(subdir,strcat('matLargeDataS',num2str(i),'.mat'));
        load(orderData);
        load(timeData);
        Rave_t = Rave_t + mean(transpose(R_t));
    end
    Rave_t = Rave_t/10;
%     Rave_t(1)
    indices = find(Rave_t>7 & Rave_t<25);
    Qs = transpose(WeightP_t);
    Ne = size(WeightP_t,2);
    orderParas = 1/(Ne-1)*(Ne - sum(Qs.^2).^2./sum(Qs.^4));
    
    [ d, ix1 ] = min( abs( Rave_t - 3*7.5 ) );
    Os1 = [Os1 orderParas(ix1)];
    
    [ d, ix0 ] = min( abs( Rave_t - 7.5 ) );
    Os0 = [Os0 orderParas(ix0)];
    [ d, ix2 ] = min( abs( Rave_t - 1.1*7.5 ) );
    Os2 = [Os2 orderParas(ix2)];
    [ d, ix3 ] = min( abs( Rave_t - 1.3*7.5 ) );
    Os3 = [Os3 orderParas(ix3)];
    [ d, ix4] = min( abs( Rave_t - 1.5*7.5 ) );
    Os4 = [Os4 orderParas(ix4)];
%     indices = indices(indices~=ix1);
%     indices = indices(indices~=ix0);
    indices = indices(indices~=ix2);
    indices = indices(indices~=ix3);
    indices = indices(indices~=ix4);
%     orderParas = Data(1,1,:);
%     orderParas = orderParas(:,:);
    Rs = [Rs Rave_t(indices)/7.5];
    Os = [Os orderParas(indices)];
    Ns = [Ns N*ones(size(indices))];
end

% T = 601;
% dir = fullfile("E:","TempCode","MatlabFlow","matData","erosion","DelaunayNet",strcat("100by50T",num2str(T),"d0.2"));
% 

% for N = 1:0.2:6
%     subdir = fullfile(dir,strcat('N',num2str(N,'%.1f')),strcat('a',num2str(a,'%.1f')));
%     timeData = fullfile(subdir,strcat('matLargeDataS',num2str(1),'.mat'));
% 
%     result = isfile(timeData);
%     if ~result
%         subdir = fullfile(dir,strcat('N',num2str(N,'%.2f')),strcat('a',num2str(a,'%.1f')));
% 
%         timeData = fullfile(subdir,strcat('matLargeDataS',num2str(1),'.mat'));
%     end
%     load(timeData);
%     Rave_t = mean(transpose(R_t));
%     indices = find(Rave_t>0.1 & Rave_t<25);
%     Qs = transpose(WeightP_t);
%     Ne = size(WeightP_t,2);
%     orderParas = 1/(Ne-1)*(Ne - sum(Qs.^2).^2./sum(Qs.^4));

%     Rs = [Rs Rave_t(indices)/7.5];
%     Os = [Os orderParas(indices)];
%     Ns = [Ns N*ones(size(indices))];
% end

sc = 40;
% fig = figure;
figure('position',[100,100,600*1.6,600]);
s = scatter(Ns,Os,sc,Rs,'filled','MarkerFaceAlpha',0.7);
s.Annotation.LegendInformation.IconDisplayStyle = 'off';
xlabel("$n$",'Interpreter','latex','FontSize',20);
ylabel("$\mathcal{O}$",'Interpreter','latex','FontSize',20);
title("2D Erosion");
J = customcolormap_preset("pasteljet");
colormap(J);
cb = colorbar;
cb.Label.String = "$\langle r \rangle/r_0$ ";
cb.Label.Interpreter = 'latex';
set(cb,'TickLabelInterpreter','latex')
m = length(J);
cmin = min(Rs);
cmax = max(Rs);
Cindex1 = fix((3-cmin)/(cmax-cmin)*m)+1; 
Cindex2 = fix((1.1-cmin)/(cmax-cmin)*m)+1;
Cindex3 = fix((1.3-cmin)/(cmax-cmin)*m)+1;
Cindex4 = fix((1.5-cmin)/(cmax-cmin)*m)+1;
color0 = J(1,:);
color1 = J(Cindex1,:);
color2 = J(Cindex2,:);
color3 = J(Cindex3,:);
color4 = J(Cindex4,:);

figure('position',[100,100,600*1.6,600]);

plot(1:0.2:6,Os0,'Marker','.','Color',color0,'LineWidth',2,'MarkerSize',10);
hold on 
plot(1:0.2:6,Os2,'Marker','x','Color',color2,'LineWidth',2,'MarkerSize',10);
plot(1:0.2:6,Os3,'Marker','^','Color',color3,'LineWidth',2,'MarkerSize',10);
plot(1:0.2:6,Os4,'Marker','d','Color',color4,'LineWidth',2,'MarkerSize',10);
plot(1:0.2:6,Os1,'Marker','v','Color',color1,'LineWidth',2,'MarkerSize',10);
xlabel("$n$",'Interpreter','latex','FontSize',20);
ylabel("$\mathcal{O}$",'Interpreter','latex','FontSize',20);
legend('$r_0$','$1.1r_0$','$1.3r_0$','$1.5r_0$','$3r_0$','Interpreter','latex','FontSize',15)