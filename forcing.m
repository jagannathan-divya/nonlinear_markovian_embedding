function f = forcing(u, interface, vel, t, bta, N, A1, A2, Nherm, whichEx)

    % Need to construct the Nmat, Fmat, and the Hermite Function vector
    x = u(1:floor(Nherm/2)); y = u(floor(Nherm/2)+1:end);
    
    psiatl = hermiteFunction(Nherm, interface);
    psiat2l = hermiteFunction(Nherm, 2*interface);
    psiat0 = hermiteFunction(Nherm, 0);

    somen = 0:1:(floor(Nherm/2)-1); alter = (-1).^somen;
    C = diag(alter);

   if whichEx==1
    dfdt = exp(t);
   else
       dfdt = 0;
   end
    %dfdt = (pi/4)*cos((pi/2)*t);
    % A = 1; mu = 0.1; t0 = 5;
    % dfdt = -2*A*mu*exp(-mu*(t-t0)^2)*(t-t0);

    f1 = -vel*N*y + (bta*vel*sqrt(2*pi))*(A2*psiat2l(2:2:end)) - (2*dfdt*sqrt(2*pi))*C*psiatl(1:2:end);
    f2 = vel*N'*x + (bta*vel*sqrt(2*pi))*(A1*(psiat0(1:2:end) - psiat2l(1:2:end))) - (2*dfdt*sqrt(2*pi))*C*psiatl(2:2:end);

    f = vertcat(f1, f2);
end
