% generate N matrix in the evolution equation for the coefficient vector
% this is a lower bidiagonal matrix

function [Nmatrix] = Nmat(N)
    n = 0:1:(N-1);
    leadDiag = sqrt(0.5*(2*n+1));
    downleadDiag = sqrt(n(1:end-1)+1);
    A0 = diag(leadDiag);
    A1 = diag(downleadDiag,-1);
    Nmatrix = A0 + A1;
end