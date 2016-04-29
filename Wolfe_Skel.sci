function [alphan,ok]=Wolfe(alpha,x,D,Oracle)


//////////////////////////////////////////////////////////////
//                                                          //
//   RECHERCHE LINEAIRE SUIVANT LES CONDITIONS DE WOLFE     //
//                                                          //
//                                                          //
//  Arguments en entree                                     //
//  -------------------                                     //
//    alpha  : valeur initiale du pas                       //
//    x      : valeur initiale des variables                //
//    D      : direction de descente                        //
//    Oracle : nom de la fonction Oracle                    //
//                                                          //
//  Arguments en sortie                                     //
//  -------------------                                     //
//    alphan : valeur du pas apres recherche lineaire       //
//    ok     : indicateur de reussite de la recherche       //
//             = 1 : conditions de Wolfe verifiees          //
//             = 2 : indistinguabilite des iteres           //
//                                                          //
//                                                          //
//    omega1 : coefficient pour la 1-ere condition de Wolfe //
//    omega2 : coefficient pour la 2-eme condition de Wolfe //
//                                                          //
//////////////////////////////////////////////////////////////


// -------------------------------------
// Coefficients de la recherche lineaire
// -------------------------------------

   omega1 = 0.1;
   omega2 = 0.9;

   alphamin = 0.0;
   alphamax = %inf;

   ok = 0;
   dltx = 0.00000001;

// ---------------------------------
// Algorithme de Fletcher-Lemarechal
// ---------------------------------

   // Appel de l'oracle au point initial
   
   ind = 4;
   [F,G] = Oracle(x,ind);

   // Initialisation de l'algorithme

   alphan = alpha;
   xn = x;

   // Boucle de calcul du pas
   //
   // xn represente le point pour la valeur courante du pas,
   // xp represente le point pour la valeur precedente du pas.

   while ok == 0
      
      xp = xn;
      xn = x + (alphan*D);

      // Calcul des conditions de Wolfe
      
      //Condition 1 :
      cond1 = 0;
      Phi_alpha = Oracle(xn,2);
      //disp("xn",xn)
      Phi_max = F + omega1*alphan*(G'*D);
      //disp("G=",D)
      //disp("phi_alpha=",Phi_alpha)
      //disp("phi_max=",Phi_max)
      if Phi_alpha <= Phi_max then
          cond1 = 1;
      end
      
      //Condition 2 :
      cond2 = 0;
      [F2,G2] = Oracle(xn,4);
      max2 = G2'*D;
      min2 = omega2*(G'*D)
      if min2 <= max2 then
          cond2 = 1;
      end
      
      // -----> A completer...
      // -----> A completer...

      // Test de la valeur de alphan :
      // - si les deux conditions de Wolfe sont verifiees,
      //   faire ok = 1 : on sort alors de la boucle while
      // - sinon, modifier la valeur de alphan : on reboucle.
      // -----> A completer...
      // -----> A completer...
      
      if cond1 == 0 then
          alphamax = alphan;
          alphan = (alphamin + alphamax)/2;
          
      elseif cond2 == 0 then
          alphamin = alphan;
          if alphamax == %inf then
              alphan = 2*alphamin;
          else
              alphan = (alphamin + alphamax)/2;
          end
          
      else
          ok = 1;
          break
      end

      
      // -----> A completer...
      // -----> A completer...
      // Test d'indistinguabilite

      if norm(xn-xp) < dltx then
        ok = 2;
      end

   end

endfunction
