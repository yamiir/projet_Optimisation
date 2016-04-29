///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  MONITEUR D'ENCHAINEMENT POUR LE CALCUL DE L'EQUILIBRE D'UN RESEAU D'EAU  //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

// --------------------------------------
// Dimensionnement de l'espace de travail
// --------------------------------------

stacksize(10000000);

// ------------------------------------------
// Fonctions fournies dans le cadre du projet
// ------------------------------------------

// Donnees du problemes

exec('Probleme_R.sce');
exec('Structures_R.sce');

// Affichage des resultats

exec('Visualg.sci');

// Verification  des resultats

exec('HydrauliqueP.sci');
exec('HydrauliqueD.sci');
exec('Verification.sci');

// ------------------------------------------
// Fonctions a ecrire dans le cadre du projet
// ------------------------------------------

// ---> Charger les fonctions  associees a l'oracle du probleme,
//      aux algorithmes d'optimisation et de recherche lineaire.
//
// Exemple : la fonction "optim" de Scilab
//
exec('OraclePG.sci');
exec('Gradient_F.sci');
exec('Wolfe_Skel.sci');
exec('Gradient_Conjuge.sci');
exec('BFGS.sci');
exec('OracleDG.sci')
exec('get_Q.sci')
exec('Newton.sci');
titrgr = "Fonction optim de Scilab sur le probleme primal";

// -----> A completer...
// -----> A completer...
// -----> A completer...

// ------------------------------
// Initialisation de l'algorithme
// ------------------------------

// La dimension (n-md) est celle du probleme primal

xini_P = 0.1 * rand(n-md,1);

// La dimension (md) est celle du probleme dual

xini_D = 0.1 * rand(md,1);
// ----------------------------
// Minimisation proprement dite
// ----------------------------

// Exemple : la fonction "optim" de Scilab
methode=8
select methode
    case 1 then
        [fopt,xopt,gopt] = Gradient_F(OraclePG,xini_P);
    case 2 then
        [fopt,xopt,gopt] = Gradient_Conjuge(OraclePG,xini_P);
    case 3 then
        [fopt,xopt,gopt] = Newton(OraclePG,xini_P);
    case 4 then
        [fopt,xopt,gopt] = BFGS(OraclePG,xini_P);
    case 5 then
        [fopt,xopt,gopt] = Gradient_F(OracleDG,xini_D);
    case 6 then
        [fopt,xopt,gopt] = Gradient_Conjuge(OracleDG,xini_D);
    case 7 then
        [fopt,xopt,gopt] = Newton(OracleDG,xini_D);
    case 8 then
        [fopt,xopt,gopt] = BFGS(OracleDG,xini_D);

end
// -----> A completer...

// --------------------------
// Verification des resultats
// --------------------------
//if methode==1 then
//[q,z,f,p] = HydrauliqueP(xopt);
//else
[q,z,f,p] = HydrauliqueD(xopt);
//end
Verification(q,z,f,p);

//
