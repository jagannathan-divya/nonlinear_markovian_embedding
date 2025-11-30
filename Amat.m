% generate A1, A2 matrix in the evolution equation for the coefficient vector
% A1 is a upper bidiagonal matrix, and A2 is lower bidiagonal.

function [Amat1, Amat2] = Amat(N)
    n = 0:1:(N-1);
    leadDiag = -((-1).^n).*sqrt(0.5*(2*n+1));
    upleadDiag = ((-1).^(n(1:end-1))).*sqrt(n(1:end-1)+1);
    downleadDiag =  -upleadDiag;
    A0 = diag(leadDiag);
    A1 = diag(downleadDiag,-1);
    A2 = diag(upleadDiag,1);
    Amat1 = A0 + A2;
    Amat2 = A0 + A1;
end