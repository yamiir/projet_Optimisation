function [F,G,H]=OraclePG(qc,ind)
    q=q0+B*qc;
    select ind
    case 2 then
        F=1/3 * (q' * r.*q.*abs(q)) + pr' * (Ar*q);
    case 3 then
        G=(B'*(r.*q.*abs(q)))+(Ar*B)'*pr;
    case 4 then
        F=1/3 * (q' * r.*q.*abs(q)) + pr' * (Ar*q);
        G=(B'*(r.*q.*abs(q)))+(Ar*B)'*pr;
        //+2/3*(B'*((sign(q)).*r.*q.*q))...
    case 5 then
        H=2*(B'*(repmat(sign(q).*r.*q,1,n-md).*B));
       // +4/3*(B'*(sign(q).*r.*q)';
    case 6 then
        G=(B'*(r.*q.*abs(q)))+(Ar*B)'*pr;
        H=2*(B'*(repmat(sign(q).*r.*q,1,n-md).*B));
     case 7 then
        F=1/3 * (q' * r.*q.*abs(q)) + pr' * (Ar*q);
        G=(B'*(r.*q.*abs(q)))+(Ar*B)'*pr;
        H=2*(B'*(repmat(sign(q).*r.*q,1,n-md).*B));
    end
endfunction
