function [L, N, A1, A2] = operatorMatrices(Nherm)
    % matrix operators in the forcing term
    N = Nmat(floor(Nherm/2));
    [A1, A2] = Amat(floor(Nherm/2));
    % matrix propagator
    L1 = Lmat(floor(Nherm/2), 2); L2 = Lmat(floor(Nherm/2), 1); 
    L = vertcat(horzcat(L1, 0*L1), horzcat(0*L2, L2));   
end
