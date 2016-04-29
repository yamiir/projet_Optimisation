function Q=get_Q(lambda)
    T=(-Ar'*pr-Ad'*lambda)
    //disp("T=",T)
    Q=sqrt((abs(T)./r)).*sign(T)
endfunction
