function [fopt,xopt,gopt]=Gradient_Conjuge(OraclePG,xini)
    alphai=0.99
    x0=xini
    [F0,G0]=OraclePG(x0,3)
    d0=-G0
    alpha=Wolfe(alphai,x0,d0,OraclePG)
    x1=x0+alpha*d0
    dk0=d0
    tol=0.00001
    iter=5000
    logG=[]
    logP=[]
    Cout=[]
    for k=1:iter

        [F1,G1]=OraclePG(x1,3)
        Beta=(G1-G0)'*G1/(norm(G0)^2)
        dk1=-G1+Beta*dk0
        //disp("G1=",G1)
        if norm(G1) <= tol then
            break
        end
        alpha=Wolfe(alphai,x1,dk1,OraclePG)
        //disp('alpha=',alpha)
        x2=x1+alpha*dk1
        x1=x2
        x0=x1
        dk0=dk1
        F0=F1
        G0=G1  
        //disp('G1=',G1)
        logG = [ logG ; log10(norm(G1)) ];
        logP = [ logP ; log10(alpha) ];
        Cout = [ Cout ; F1 ];
    end
    fopt=F1
    xopt=x2
    gopt=G1
    tcpu = timer();

    cvge = ['Iteration         : ' string(k);...
    'Temps CPU         : ' string(tcpu);...
    'Critere optimal   : ' string(fopt);...
    'Norme du gradient : ' string(norm(gopt))];
    disp('Fin de la methode de gradient a pas fixe')
    disp(cvge)
    // - visualisation de la convergence

    Visualg(logG,logP,Cout);

endfunction
