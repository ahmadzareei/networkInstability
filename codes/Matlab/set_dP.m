
function [LHS,RHS] = set_dP(LHS,RHS,Boundary, P_left,P_right)
    % We know that the Pressure on the left boundary nodes should be P_left
    % we remove the equation corresponding to flux equation for those nodes
    % and replace it with P = P_left
    for i = Boundary.left
        LHS(i,:) = 0;
        LHS(i,i) = 1;
        RHS(i,1) = P_left;
    end
    % Similar procedure as above for right bouPlotPlotndary nodes
    for j = Boundary.right
        LHS(j,:) = 0;
        LHS(j,j) = 1;
        RHS(j,1) = P_right;
    end
end