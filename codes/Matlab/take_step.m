

function [A,maxflow,order] = take_step(A,N,maxDeltaD,type,PBCs,Boundary)
    % take a step in cloggin or erosion
    % type can be type='clog' or 'erode'
    % A is the graph
    % N is the power for the change
    % maxDeltaD is the maximum change in diameters
    % PBCs are the boundary conditions
    % Boundary are the boundaries
    
    deltaD = A.Edges.Weight./A.Edges.Diameter.^N;
    alpha = maxDeltaD/max(deltaD);
    
    if strcmp(type,'clog')
        alpha = -1*alpha;
    elseif strcmp(type,'erode')
        alpha = 1*alpha;
    end
    
    A.Edges.Diameter = A.Edges.Diameter + alpha*deltaD;
    
    n_tot =  max(max(A.Edges.EndNodes));
    [LHS,RHS] = set_LHS_RHS(A,n_tot);
    [LHS,RHS] = set_dP(LHS,RHS,Boundary, PBCs.left,PBCs.right);
    P = full(LHS\RHS); % solving for P
    [A,maxflow,order] = post_process(A,P,Boundary);

end