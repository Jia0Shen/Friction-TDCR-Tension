function K = getKmat(Kb,Kt,n)

    K = diag(repmat([Kb,Kb,Kt],1,n));
    
end