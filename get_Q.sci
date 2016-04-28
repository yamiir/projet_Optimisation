function Q=get_Q(lambda)
    T=(-Ar'*pr-Ad'*lambda)
    Q=T.^(0.5).*sign(T)
endfunction
