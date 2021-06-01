function [C,Q,R,Cn,Qn,Rn,Cb,Qb,Rb] = make_CQR(A,n_tot,Boundary)

    Nn = n_tot;
    Ne = size(A.Edges.EndNodes,1);
    C = zeros(Ne,Nn);
    C((1:Ne)' + (A.Edges.EndNodes(:,2)-1)*Ne)=1;
    C((1:Ne)' + (A.Edges.EndNodes(:,1)-1)*Ne)=-1;
    Q = C';
    R = diag(A.Edges.Diameter.^4./A.Edges.Length);

    Cn = C;
    Cn(:,[Boundary.left,Boundary.right])= [];
    Qn = Cn';
    Rn = R;
    
    Cb = C;
    Cb = C(:,[Boundary.left,Boundary.right]);
    Qb = Cb';
    Rb = R;   
    
end