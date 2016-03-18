function [fopt,xopt,gopt]=Gradient_Conjuge(OraclePG,xini)
    alphai=0.001
    G=OraclePG(x,3)
    d=-G
    alpha=Wolfe(alphai,x,d,OraclePG)
    iter=5000
    x=xini
    for k=1:iter
        Beta=d'*G
        d=-OraclePG(x,3)+
        
    end
