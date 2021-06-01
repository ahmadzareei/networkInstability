
function [A,maxflow,order] = post_process(A,P,Boundary)
    
    A.Edges.Weight = (P(A.Edges.EndNodes(:,1)) - P(A.Edges.EndNodes(:,2))).*(A.Edges.Diameter).^4./A.Edges.Length;
    maxflow = sum(A.Edges.Weight(ismember(A.Edges.EndNodes(:,1),Boundary.left)));
    N = size(A.Edges.EndNodes,1);
    order = 1/(N-1)*(N - sum(A.Edges.Weight.^2).^2/sum(A.Edges.Weight.^4));    

end