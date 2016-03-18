function [PasOp]=Pas_Cauchy(OraclePG,xini,D)
    x1=xini
    g1=OraclePG(x1,3)
    G1=g'*D
    iter=500
    alpha=0.01
    for k=1:iter
        x1=x1-alpha*G1
        g1=OraclePG(x1,3)
        G1=g'*D
    end
