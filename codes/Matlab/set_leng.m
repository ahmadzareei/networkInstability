


function A = set_leng(A,vars,type)

   numEdge = size(A.Edges.EndNodes,1);       
    if strcmp(type,'uniform')
       A.Edges.Length = (vars.Lmax-vars.Lmin)*rand(numEdge,1) + vars.Lmin;
    elseif strcmp(type,'uniform_setAve')
        sigma = vars.l_ave*sqrt(3);
        A.Edges.Length = vars.l_ave + sigma*(rand(numEdge,1)-0.5);
    elseif strcmp(type,'const');
       A.Edges.Length = vars.L0*ones(numEdge,1);
    elseif strcmp(type,'gaussian')
       A.Edges.Length = vars.mu + vars.sigma*randn(numEdge,1); % The diameter
       % making sure no diameter is below zero
       A.Edges.Length = A.Edges.Length.*(A.Edges.Length>0) + (A.Edges.Length<=0)*vars.smallest;       
    elseif strcmp(type,'lognormal')
        pd = makedist('Lognormal','mu',vars.mu,'sigma',vars.sigma);
        A.Edges.Length = random(pd,numEdge,1); % The diameter
    end

end
