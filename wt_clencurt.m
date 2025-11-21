function [x,w] = wt_clencurt(N)
% (N+1) quadrature points: n=0,1,...,N
    
 theta = pi*(0:N)/N; theta = fliplr(theta); x = cos(theta);       % abscissae
 
 dvec = linspace(1,1,N+1); dvec(1)=0.5; dvec(end)=0.5;
 bvec = linspace(0,0,N+1); bvec(1) = pi;
 avec = (pi/2)*linspace(1,1,N+1); avec(1) = pi;
 w = (pi/(N))*(bvec(1)/avec(1))*dvec;    
end