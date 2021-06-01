

function [A, points,n_tot,Boundary] = make_graph(Lengths,type)
% Creates a graph of size Lx X Ly 
% returns graph A, with Edges 
% each node is at location x = points(:,1) and y = points(:,2);

    if strcmp(type,'diamond')
        Ly = Lengths.Ly;
        Lx = Lengths.Lx;
        m = Ly;
        n = Lx;    

        n_tot = n*(2*m+1)+m;

        pointX = repmat([1:2:2*m,0:2:2*m],[1,n+1]);
        pointX = pointX(1:n_tot);
        pointY = [];
        for j = 0:2*n
            if mod(j,2) == 0
                pointY = [pointY,j*ones(1,m)];
            else
                pointY = [pointY,j*ones(1,m+1)];        
            end
        end

        points = [pointX;pointY];


        s = [];
        t = [];

        for j = 1:n_tot - m
            if mod(j,2*m+1)<m+1 && mod(j,2*m+1)>0  
                s = [s,j,j];
                t = [t,j+m,j+m+1];
            elseif mod(j,2*m+1) == m+1
                s = [s,j];
                t = [t,j+m+1];
            elseif mod(j,2*m+1) == 0
                s = [s,j];
                t = [t,j+m];
            else
                s = [s,j,j];
                t = [t,j+m,j+m+1];
            end
        end

        A = graph(s,t);
        points = flipud(points);
        numEdge = size(A.Edges,1); % which is basically m*(n-1) + n*(m-1) for a square graph

        % Definding the boundary nodes (This is for a square grid)
        Boundary.left = 1:1:m; ; % left nodes  m:(2*m+1):n_tot
        Boundary.right = n_tot-m+1:1:n_tot; % right nodes

        Boundary.bottom = Ly+1:(2*Ly+1):n_tot;
        Boundary.top = 2*Ly+1:(2*Ly+1):n_tot;
        all_nodes = 1:n_tot;
        Boundary.middle =  all_nodes(~ismember(all_nodes,[Boundary.top,Boundary.bottom,Boundary.left,Boundary.right]));
    elseif strcmp(type,'struc')
        Lx = Lengths.Lx;
        Ly = Lengths.Ly;
        
        x = 1:1:Lx; m = length(x);
        y = 1:1:Ly; n = length(y);
        points = [repmat(x,1,n); reshape(repmat(y',1,m)',m*n,1)'];
        % create the adjency matrix (adj) for grid connections
        % it will be used for creating the graph of conections
        r= n;c=m;                          % Get the matrix size
        diagVec1 = repmat([ones(c-1, 1); 0], r, 1);  % Make the first diagonal vector
        %   (for horizontal connections)
        diagVec1 = diagVec1(1:end-1);                % Remove the last value
        diagVec2 = ones(c*(r-1), 1);                 % Make the second diagonal vector
        %   (for vertical connections)
        adj = diag(diagVec1, 1)+diag(diagVec2, c);   % Add the diagonals to a zero matrix
        adj = adj+adj.';                             % Add the matrix to a transposed copy of
        % Creating the graph
        A = graph(adj);
        numEdge = size(A.Edges,1); % which is basically m*(n-1) + n*(m-1) for a square graph
        n_tot = m*n;
        % Definding the boundary nodes (This is for a square grid)
        Boundary.left = 1:m:m*n; % left nodes
        Boundary.right = m:m:m*n; % right nodes
        Boundary.bottom = 1:1:m;
        Boundary.top = (n*m-m+1):1:n*m;
    elseif strcmp(type,'3D')
        Lx = Lengths.Lx;
        Ly = Lengths.Ly;
        Lz = Lengths.Lz;

        x = 1:1:Lx; m = length(x);
        y = 1:1:Ly; n = length(y);
        z = 1:1:Lz; o = length(z);


        %     1 -----------2 ------------3 ---------- ..... -------m-1 --------------  m
        %    m+1 ---------m+2 ----------m+3 ----------..... ------ 2m-1 ------------- 2m
        %   2m+1 --------2m+2 ---------2m+3 ----------..... ------ 3m-1 ------------- 3m
        %  ....
        % (n-1)m+1 ---- (n-1)m+2 ----(n-1)m+3 ----------..... -----nm-1 --------------nm

        % the second row starts from nm+1 and goes to 2nm
        % and this goes on until the last row

        pointsXY = [repmat(x,1,n); reshape(repmat(y',1,m)',m*n,1)']; % this is the xy plane points
        pointsXY = repmat(pointsXY,1,o);
        points = [pointsXY;reshape(repmat(z',1,m*n)',m*n*o,1)'];

        % plotting the points if you like
        % scatter3(points(1,:),points(2,:), points(3,:))

        % create the adjency matrix (adj) for grid connections
        % it will be used for creating the graph of conections
        r= n;c=m;                          % Get the matrix size
        diagVec1 = repmat([ones(c-1, 1); 0], r, 1);  % Make the first diagonal vector
        %   (for horizontal connections)
        diagVec1 = diagVec1(1:end-1);                % Remove the last value
        diagVec2 = ones(c*(r-1), 1);                 % Make the second diagonal vector
        %   (for vertical connections)
        adj = spdiags([0;diagVec1], 1,m*n,m*n)+spdiags([zeros(m*n-length(diagVec2),1);diagVec2], c,m*n,m*n);   % Add the diagonals to a zero matrix
        adj = adj+adj.';                             % Add the matrix to a transposed copy of
        % Creating the graph

        adjA =sparse([]);
        for i = 1:o
            adjA = blkdiag(adjA,adj);
        end

        adjA = adjA + spdiags([zeros(m*n,1); ones(m*n*(o-1),1)],m*n,m*n*o,m*n*o) + spdiags([ones(m*n*(o-1),1);zeros(m*n,1)],-m*n,m*n*o,m*n*o);

        A = graph(adjA);

        %  plot(A,'XData',points(1,:),'YData',points(2,:),'ZData', points(3,:));

        numEdge = size(A.Edges,1); % which is basically m*(n-1) + n*(m-1) for a square graph
        n_tot = m*n*o;
        % Definding the boundary nodes (This is for a square grid)
        Boundary.left = 1:m:m*n*o; % left nodes
        Boundary.right = m:m:m*n*o; % right nodes
    else strcmp(type,'random')
        % Lengths.Lx = 20;
        % Lengths.Ly = 10;

        Lx = Lengths.Lx;
        Ly = Lengths.Ly;


        pointX = 1 + Lx*rand(1,Lx*Ly); % [0,Lx]
        pointY = 1 + Ly*rand(1,Lx*Ly); % [0,Ly]

        pointX = [(Lx+2)*rand(1,Lx), (Lx+2)*ones(1,Ly),(Lx+2)*rand(1,Lx),zeros(1,Ly), 0, Lx+2, Lx+2,0, pointX];
        pointY = [zeros(1,Lx), (Ly+2)*rand(1,Ly),(Ly+2)*ones(1,Lx),(Ly+2)*rand(1,Ly), 0, 0, Ly+2, Ly+2, pointY];
        points = [pointX;pointY];

        n_tot = length(pointX);
        Boundary.bottom = 1:Lx;
        Boundary.right = (Lx+1):(Lx+Ly);
        Boundary.top = (Lx+Ly+1):(2*Lx+Ly);
        Boundary.left = (2*Lx+Ly+1):(2*Lx+2*Ly);
        all_nodes = 1:n_tot;
        Boundary.middle =  all_nodes(~ismember(all_nodes,[Boundary.top,Boundary.bottom,Boundary.left,Boundary.right]));


        delauny_triangulation= delaunay(pointX,pointY);
        s = [delauny_triangulation(:,1);delauny_triangulation(:,2);delauny_triangulation(:,3)]; % first points of the edges
        t = [delauny_triangulation(:,2);delauny_triangulation(:,3);delauny_triangulation(:,1)]; % second points of the edges
        ST = [s,t];
        ST = sort(ST,2);
        ST = unique(ST,'rows');

        A = graph(ST(:,1),ST(:,2));
        % scatter(pointX,pointY,'.')
        % triplot(dt,pointX,pointY)

    end

end
