function [F,G,H]=OracleDG(lambda,ind)
    q=get_Q(lambda)
    //disp("q=",q)
    F=0
    G=0
    H=0
    T=-Ad'*lambda-Ar'*pr
    q_lambda=diag(1 ./(2*sqrt(r.*abs(T))))*(-Ad')
    select ind
    case 2 then
        F=-(1/3*(q'*(q.*abs(q).*r))+pr'*(Ar*q)+lambda'*(Ad*q-fd))
    case 3 then
        G=-(Ad*q-fd)

    case 4 then
        F=-(1/3*(q'*(q.*abs(q).*r))+pr'*(Ar*q)+lambda'*(Ad*q-fd))
        G=-(Ad*q-fd)
    case 5 then
        H=-Ad*q_lambda
    case 6 then
        G=-(Ad*q-fd)
        H=-Ad*q_lambda
    case 7 then
        F=-(1/3*(q'*(q.*abs(q).*r))+pr'*(Ar*q)+lambda'*(Ad*q-fd))
        G=-(Ad*q-fd)
        H=-Ad*q_lambda
    end
endfunction
