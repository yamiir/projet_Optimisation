function [F,G,H]=OraclePG(qc,ind)
    q=q0+B*qc;
    F=0
    G=0
    H=0
    select ind
    case 2 then
        F=1/3*sum(q.*r.*q.*abs(q))+pr'*(Ar*q);
    case 3 then
        G=(1/3*(B'*(r.*q.*abs(q)))...
        +2/3*(B'*((sign(q)).*r.*q.*q))...
        +(Ar*B)'*pr);
    case 4 then
        F=1/3*sum(q.*r.*q.*abs(q))+pr'*(Ar*q);
        G=(B'*(r.*q.*abs(q)))+(Ar*B)'*pr;
        //+2/3*(B'*((sign(q)).*r.*q.*q))...
    case 5 then
        H=2*(B'*(diag(sign(q).*r.*q))*B);
       // +4/3*(B'*(sign(q).*r.*q)';
    case 6 then
        G=(B'*(r.*q.*abs(q)))+(Ar*B)'*pr;
        H=2*(B'*(diag(sign(q).*r.*q))*B);
     case 7 then
        F=1/3*sum(q.*r.*q.*abs(q))+pr'*(Ar*q);
        G=(B'*(r.*q.*abs(q)))+(Ar*B)'*pr;
        H=2*(B'*(diag(sign(q).*r.*q))*B);
    end
endfunction
