function [P,F,T,theta,beta,lamb] = stateExtract_t(tdcr, x, z)
    m = tdcr.m; 
    n = tdcr.n;
    
    P = x(1:3*m*n);
    F = x(3*m*n+1:6*m*n);
    T = x(6*m*n+1:7*m*n);
    theta = x(7*m*n+1:8*m*n);
    
    beta = z(1:2*m*n);
    lamb = z(2*m*n+1:3*m*n);

end