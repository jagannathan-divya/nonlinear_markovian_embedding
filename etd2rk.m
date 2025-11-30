function [expLh1, expLh2, expLh, M1, M2, N21, nodes, aij, b] = etd2rk(L, dt, opt)

   if opt==1
        % Heun's
       nodes = [0, 1];
       aij = zeros(2,2); aij(2,1) = 1;
       b = linspace(0,0,2);
       b(1) = 0.5; b(2) = 0.5;
   elseif opt==2
       % Raltson
       nodes = [0, 2/3];
       aij = zeros(2,2); aij(2,1) = 2/3;
       b = linspace(0,0,2);
       b(1) = 1/4; b(2) = 3/4;
   else
       % Midpoint
       nodes = [0, 1/2];
       aij = zeros(2,2); aij(2,1) = 1/2;
       b = linspace(0,0,2);
       b(1) = 0; b(2) = 1;
   end

   % Computing the exponential of matrix
   [V,D] = eig(dt*L);
   expLh = (V*expm(D))/(V);
   %expLh = real(expLh);

   [V1,D1] = eig(nodes(1)*dt*L);
   expLh1 = (V1*expm(D1))/(V1);
   %expLh1 = real(expLh1);

   
   [V2,D2] = eig(nodes(2)*dt*L);
   expLh2 = (V2*expm(D2))/(V2);
   %expLh2 = real(expLh2);

   N21 = (1/dt)*(L\(expLh2 - eye(size(L))));
 
   L2 = L*L;
   M2 = (1/dt^2)*(L2\(expLh - (eye(size(L)) + dt*L)))/nodes(2);
   M1 = (1/dt)*(L\(expLh - eye(size(L)))) - M2;
  
   % checkvec = (diag(D)==0);
   % check = 0;
   % if (any(checkvec))
   %     check=1; 
   % end
end