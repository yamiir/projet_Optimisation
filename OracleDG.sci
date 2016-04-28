function [F,G,H]=OracleDG(lambda,ind)
    q=get_Q(lambda)
    F=0
    G=0
    H=0
    T=-Ad'*lambda-Ar'*pr
    q_lambda=diag(T./(2*(T./r).^0.5))*(-Ad')
    select ind
    case 2 then
        F=1/3*sum(q.*q.*abs(q).r)+pr'*(Ar*q)+lambda'*(Ad*q-fd)
    case 3 then
        G=q_lambda*(r.*q.*abs(q)+Ar'*pr+Ad'*lambda)+Ad*q-fd

    case 4 then
        F=1/3*sum(q.*q.*abs(q).r)+pr'*(Ar*q)+lambda'*(Ad*q-fd)
        G=q_lambda*(r.*q.*abs(q)+Ar'*pr+Ad'*lambda)+Ad*q-fd
    case 5 then
        H=
    case 6 then
        G=
        H=
    case 7 then
        F=
        G=
        H=
    end
endfunction
