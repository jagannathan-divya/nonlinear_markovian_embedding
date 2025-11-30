% Solving the Markovian system for the Superwalking droplet problem with
% the Bessel-J1 kernel
% Written by: Divya Jagannathan at ICTS on 27 June 2023

clear all
clc
close all

% physical system parameters:
c1 = 1.5; c2 = 0.1; x0 = 1; v0 = 1;

% scheme parameters
Ndiv = 31; [k,w] = wt_clencurt(Ndiv); Nquad = size(k,2);                    % Nquad is the number of quadrature points used to compute the integral in 3.3b

% time-integration parameters
tfinal = 2000; dt = 1/2^8; nsteps = floor(tfinal/dt)+1; stage = 2;          % Toggle stage between 1 and 2 for Euler forward and Heun's method resp.

params = [c1,c2,dt];

% state declaration
x = linspace(0,0,nsteps); v = linspace(0,0,nsteps);                         % original state vars
hr = linspace(0,0,Nquad); hi = linspace(0,0,Nquad);                         % addt. vars in extended state; hr and hi are real and imag parts of the H-function in the draft
time = linspace(0,0,nsteps);

% initial conditions
x(1) = x0; v(1) = v0; 

% ----------------------------%
% TIME-INTEGRATION BEGINS HERE

istep = 2;

while istep<=nsteps
   oldstate = [x(istep-1),v(istep-1)];
   [x(istep),v(istep),hr,hi] = swd_evolve(oldstate,hr,hi,params,k,w,stage);
   time(istep) = time(istep-1) + dt;
   istep = istep + 1;
end

% TIME-INTEGRATION ENDS HERE 
% ----------------------------%


% PLOTS
hop = 1;

% Figure 1 for position and velocity versus time
figure(1)
p1 = plot(time(1:hop:end),x(1:hop:end),'b-','LineWidth',2,'MarkerSize',2); hold on
p2 = plot(time(1:hop:end),v(1:hop:end),'r-','LineWidth',2,'MarkerSize',2);
xlabel('time', 'FontSize', 18);
legend([p1,p2],{'position $x_d$', 'velocity $\dot{x}_d$'},'Location','south', 'Interpreter','latex', 'FontSize',18);
set(gca, 'FontSize', 16);

% Figure 2 for v, vdot
figure(2)
vdot = diff(v)./diff(time);
plot(v(2:hop:end),vdot(1:hop:end),'k-','MarkerSize',2,'MarkerFaceColor','k','LineWidth',0.5); hold on
ylabel('$\ddot{x}_d$','Interpreter','latex', 'FontSize', 18); xlabel('$\dot{x}_d$','Interpreter','latex', 'FontSize',18);
set(gca, 'FontSize', 16);

figure(3)
trunc_step = ceil(nsteps/2);
p2 = plot(time(trunc_step:hop:end),v(trunc_step:hop:end),'r-','LineWidth',1,'MarkerSize',2);
xlabel('time');
ylim([-5,5]); pbaspect([10, 3, 1]);
%legend([p2],{'velocity'},'Location','north');
set(gca, 'FontSize', 16);
