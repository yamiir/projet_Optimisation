function [fopt,xopt,gopt]=BFGS(OraclePG,xini)
    titre = "Parametres du gradient a pas fixe";
    labels = ["Nombre maximal d''iterations";...
    "Valeur du pas de gradient";...
    "Seuil de convergence sur ||G||"];
    typ = list("vec",1,"vec",1,"vec",1);
    default = ["500";"1";"0.000001"];
    [ok,iter,alphai,tol] = getvalue(titre,labels,typ,default);

    // ----------------------------
    // Initialisation des variables
    // ----------------------------
    logG = [];
    logP = [];
    Cout = [];

    timer();

    // -------------------------
    // Boucle sur les iterations
    // -------------------------

    dim = length(xini);
    x = xini;
    delta_u=0;
    delta_G=0;
    H = OraclePG(xini,5);
    disp(H);
    Wk=inv(H);
    kstar = iter;
    for k = 1:iter
        [Fn,Gn] = OraclePG(x,4);
        if norm(Gn) <= tol then
            kstar = k;
            break
        end
        Dk=-Wk*Gn;
        alpha = Wolfe(alphai,x,Dk,OraclePG);
        delta_u = alpha * Dk;
        x=x+delta_u;
        [Fn_1,Gn_1] = OraclePG(x,4);
        delta_G = Gn_1 - Gn;

        uG = delta_u * delta_G';
        denom = delta_G'*delta_u;
        Wk_1=(eye(dim,dim)-uG/denom) * Wk * (eye(dim,dim) - uG'/denom) + delta_u*delta_u' / denom;

        Wk=Wk_1;    
        logG = [ logG ; log10(norm(G)) ];
        logP = [ logP ; log10(alpha) ];
        Cout = [ Cout ; F ];
    end


    // ---------------------------
    // Resultats de l'optimisation
    // ---------------------------

    fopt = F;
    xopt = x;
    gopt = G;

    tcpu = timer();

    cvge = ['Iteration         : ' string(kstar);...
    'Temps CPU         : ' string(tcpu);...
    'Critere optimal   : ' string(fopt);...
    'Norme du gradient : ' string(norm(gopt))];
    disp('Fin de la methode de gradient a pas fixe')
    disp(cvge)

    // - visualisation de la convergence

    Visualg(logG,logP,Cout);



endfunction
