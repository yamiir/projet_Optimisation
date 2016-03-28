function [fopt,xopt,gopt]=BFGS(OraclePG,xini)
    x0=xini
    [F0,G0,H0]=OraclePG(xini,7)
    W0=inv(H0)
    d0=-W0*G0
    alphai=1
    alpha=Wolfe(alphai,xini,d0,OraclePG)
    x1=x0+alpha*d0
    iter=5000
    tol=0.00001
    logG=[]
    logP=[]
    Cout=[]
    for k=1:iter
        [F1,G1]=OraclePG(x1,3)
        
        if norm(G1)<tol then
            break
        end
        
        thetaU=x1-x0
        thetaG=G1-G0;
        //disp('G1=',thetaG)
        Dim=size(thetaU);
        I=eye(Dim(1),Dim(1));
        W1=(I-(thetaU*thetaG')/(thetaG'*thetaU))*W0*(I-(thetaG*thetaU')/(thetaG'*thetaU))...
        +(thetaU*thetaU')/(thetaG'*thetaG);
        d1=-W1*G1;
        alpha=Wolfe(alphai,x1,d1,OraclePG);
        //disp('alpha=',alpha)
        x2=x1+alpha*d1;

        F0=F1;
        G0=G1;
        x0=x1;
        x1=x2;
        W0=W1;
      logG = [ logG ; log10(norm(G1)) ];
      logP = [ logP ; log10(alpha) ];
      Cout = [ Cout ; F1 ];
    end
    fopt=F1
    gopt=G1
    xopt=x2
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
