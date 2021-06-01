
function [LHS,RHS] = set_LHS_RHS(A,n_tot)
    LHS = sparse(n_tot,n_tot);
    for i = 1:n_tot
        % We know that the flux is Q =  d^4 / (128 mu ) ( P_neigh - P_node)
        % and we set the total flux going into a node should addup to zero
        [eid,nid] = outedges(A,i); % find the outgoing edge from node i, the edgeID is stored in eid, and the nodeid is stored in nid
        LHS(i,i) = -sum((A.Edges.Diameter(eid)).^4./(A.Edges.Length(eid)));
        for k = 1:length(nid)
            LHS(i,nid(k)) = (A.Edges.Diameter(eid(k)))^4/A.Edges.Length(eid(k));
        end
    end
    RHS = sparse(n_tot,1);
end





