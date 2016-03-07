function [F,G]=OraclePG(qc,ind)
    B=[-AdI*AdC;eye(n-md,n-md)];
    q=q0+B*qc;
    select ind
    case 2 then
        F=1/3*sum(q.*r.*q.*abs(q))+sum(pr.*(Ar*q));
    case 3 then
        G=(1/3*...
        sum(B.*repmat(r.*q.*abs(q)...
        ,1,n-md),1)...
        +2/3*sum(B.*repmat((2*sign(q)).*r.*q.*q,1,n-md),1)...
        +sum(repmat(pr,1,n-md).*(Ar*B),1))';
    case 4 then
        F=1/3*sum(q.*r.*q.*abs(q))+sum(pr.*(Ar*q));
        G=(1/3*sum(B.*repmat(r.*q.*abs(q),1,n-md),1)...
        +2/3*sum(B.*repmat((2*sign(q)).*r.*q.*q,1,n-md),1)...
        +sum(repmat(pr,1,n-md).*(Ar*B),1))';

    end
endfunction