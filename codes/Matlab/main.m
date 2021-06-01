
%% Clear all and close all
clear all
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Setting Diameter Distribution
%%%%%%%%%%% Uniform
dDist = 'uniform';
dVar.Rmax = 14;
dVar.Rmin = 1;
%%%%%%%%%%%% Constant
% dDist = 'const';
% dVar.D0 = 1;
%%%%%%%%%%%%% Gaussian
% dDist = 'gaussian'
% dVar.sigma = 3.0;
% dVar.mu = 7.0;
% dVar.smallest = 1e-2;
%%%%%%%%%%%%% Log Normal
% dDist = 'lognormal'
% dVar.sigma = .30;
% dVar.mu = 3.0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% L distribution
lDist = 'const';
lVar.L0 = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Size & Grid Type
Lengths.Ly = 50;
Lengths.Lx = 100;
grid_type = 'random'; % 'diamond'; % 'struc'; % '3D'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% P boundary condition
PBCs.left = 10; % Presure at left
PBCs.right = 1; % Pressure at right
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Running Initializer
[A,maxflow,order,points,Boundary] = initializer(dDist,dVar, lDist,lVar,Lengths, PBCs,grid_type);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot  result
plot_graph(A,points,A.Edges.Weight/maxflow*Lengths.Ly)


%%

function [A,maxflow,order,points,Boundary] = initializer(dDist,dVar, lDist,lVar,Lengths, P_BCs,grid_type)


    %% Defining the grid
    disp('defining the grid....')

    
    [A, points,n_tot,Boundary] = make_graph(Lengths,grid_type);
    % [A, points,n_tot,Boundary] = make_graph(Lx,Ly,'struc');
    % plot_graph(A,points,'');
    %%%%%%%% Definding the boundary nodes (This is for a square grid)
    left_nodes = Boundary.left;
    right_nodes = Boundary.right;
    bottom_nodes = Boundary.bottom;
    top_nodes = Boundary.top;
    middle_nodes = Boundary.middle;


    %% Initializing the pressure at left and right boundary - Diameter of edges, and length of each edge
    disp('initializing the parameters...')

    A = set_diam(A,dVar,dDist);
    A = set_leng(A,lVar,lDist);



    %% LHS*P = RHS -- finding LHS matrix, and RHS matrix
    disp('setting up the matrix ...')
    [LHS,RHS] = set_LHS_RHS(A,n_tot);

    disp('Setting up the boundary condition ...');
    [LHS,RHS] = set_dP(LHS,RHS,Boundary, P_BCs.left,P_BCs.right);

    disp('solving for the result ...');
    P = full(LHS\RHS); % solving for P


    % calculating flux in each node and maximum flow
    [A,maxflow,order] = post_process(A,P,Boundary);


    % plot_graph(A,points,A.Edges.Weight);

    % save('initializer.mat');
    
end    