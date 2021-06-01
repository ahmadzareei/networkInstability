


function A = set_diam(A,vars,type)
   numEdge = size(A.Edges.EndNodes,1);       
    if strcmp(type,'uniform')
       A.Edges.Diameter = (vars.Rmax-vars.Rmin)*rand(numEdge,1) + vars.Rmin;
    elseif strcmp(type,'uniform_setAve')
        % sigma = vars.d_ave*sqrt(3);
        sigma = vars.d_ave*1.75;
        A.Edges.Diameter = vars.d_ave + sigma*(rand(numEdge,1)-0.5);
    elseif strcmp(type,'const');
       A.Edges.Diameter = vars.D0*ones(numEdge,1);
    elseif strcmp(type,'gaussian')
       A.Edges.Diameter = vars.mu + vars.sigma*randn(numEdge,1); % The diameter
       % making sure no diameter is below zero
       A.Edges.Diameter = A.Edges.Diameter.*(A.Edges.Diameter>0) + (A.Edges.Diameter<=0)*vars.smallest;       
    elseif strcmp(type,'lognormal')
        pd = makedist('Lognormal','mu',vars.mu,'sigma',vars.sigma);
        A.Edges.Diameter = random(pd,numEdge,1); % The diameter
    end
end
