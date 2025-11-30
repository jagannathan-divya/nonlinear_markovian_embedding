function [x, v, hr, hi] = swd_evolve(oldstate,hr0,hi0,params,k,w,nstage)
    c1 = params(1); c2 = params(2); dt = params(3); 
    x0 = oldstate(1); v0 = oldstate(2);
    if nstage == 1
        % Forward Euler
        hint0 = dot(w,hr0);
        x = x0 + dt*v0;
        v = v0 + dt*(-v0+c1*hint0);
        hr = hr0 + dt*(-c2*hr0-v0*k.*hi0);
        hi = hi0 + dt*(-c2*hi0-k/pi + v0*k.*hr0);
    elseif nstage == 2
        % Heun's method
        a21 = 2./3.;
        b = [0.25,.75];
        % constructing the first stage
        hints1 = dot(w,hr0);
        xs1 = x0;
        vs1 = v0;
        hrs1 = hr0;
        his1 = hi0;
        % constructing the second stage
        hints2 = dot(w,hrs1);
        xs2 = x0 + a21*dt*vs1;
        vs2 = v0 + a21*dt*(-vs1+c1*hints1);
        hrs2 = hr0 + a21*dt*(-c2*hrs1-vs1*k.*his1);
        his2 = hi0 + a21*dt*(-c2*his1-k/pi + vs1*k.*hrs1);
        % constructing the new step using the above two stages
        x = x0 + dt*(b(1)*vs1+b(2)*vs2);
        v = v0 + dt*(b(1)*(-vs1+c1*hints1)+b(2)*(-vs2+c1*hints2));
        hr = hr0 + dt*(b(1)*(-c2*hrs1-vs1*k.*his1)+b(2)*(-c2*hrs2-vs2*k.*his2));
        hi = hi0 + dt*(b(1)*(-c2*his1-k/pi + vs1*k.*hrs1)+b(2)*(-c2*his2-k/pi + vs2*k.*hrs2));
    else
        disp("not a valid option.")
    end
end
