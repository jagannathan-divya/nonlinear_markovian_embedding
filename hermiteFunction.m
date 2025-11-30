% returns a column vector of first n-modes hermite functions at x

function hvec = hermiteFunction(N, x)
    hvec = linspace(0,0,N); % row vector
    hvec(1) = pi^(-.25)*exp(-0.5*x^2);
    hvec(2) = sqrt(2)*pi^(-.25).*x.*exp(-0.5*x^2);
    iter = 2;
    while iter<N
            n = iter-1;
            a = sqrt(2/(n+1)); b = sqrt((n)/(n+1));
            hvec(iter+1) = a*x*hvec(iter) - b*hvec(iter-1);   
            iter = iter + 1;
    end
    hvec = hvec'; %column vector
end

   
    
