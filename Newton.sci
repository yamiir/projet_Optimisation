function [fopt,xopt,gopt]=Newton(OraclePG,xini)
x0=xini;
iter=5000;
tol=0.0001;
alphai=1;
    logG=[];
    logP=[];
    Cout=[];
for k=1:iter
    [F0,G0,H0]=OraclePG(x0,7)
    if norm(G0)<tol then
        break;
    end
    
    d=-inv(H0)*G0;
    alpha=Wolfe(alphai,x0,d,OraclePG);
    x1=x0+alpha*d;
      logG = [ logG ; log10(norm(G0)) ];
      logP = [ logP ; log10(alpha) ];
      Cout = [ Cout ; F0 ];
      x0=x1
end
    fopt=F0;
    gopt=G0;
    xopt=x1;
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
