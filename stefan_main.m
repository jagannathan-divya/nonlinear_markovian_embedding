% Solve H-equation using ETDRK & AB2AM2 for Stefan problem
% Written by Divya Jagannathan on 15 Sept 2023

%clear all
clc

alp = 0.620063; whichEx = 2;

% Physical and iterative parameters
bta = 1; t0 = .1; 

if whichEx==1
    init_len = t0; % Ex 1
else
    init_len = 2*alp*sqrt(t0) ; % Ex 2
end

Tfinal = 10; dt = 0.001 ; nsteps = ceil((Tfinal - t0)/dt);

% Keep Nherm EVEN.
Nherm = 500; opt = 1;


% Matrices given number of quadrature points, Nherm. 
[L, N, A1, A2] = operatorMatrices(Nherm);

% Scheme-specific matrices
[expL1, expL2, expL, M1, M2, N21, nodes, aij, b] = etd2rk(-L, dt, opt);

% This one's to reconstruct the velocity from the history integral 
somen = 0:1:(floor(Nherm/2)-1); alter = (-1).^somen; 

evenHermite = hermiteFunction(Nherm, 0);

basis4y = sqrt(2*pi)*alter'.*evenHermite(1:2:end);

% Initialisation
t = linspace(0, 0, nsteps); t(1) = t0;

interface = linspace(0, 0, nsteps); vel = linspace(0, 0, nsteps);

u = linspace(0, 0, Nherm); u = u';

% Initial condition
interface(1) = init_len; 
if whichEx==1
    vel(1)= init_int(init_len,t0,t0);
else
    vel(1) = init_int_2(init_len, init_len, t0, t0);
end

monebybtapi = -1/(bta*pi);

% ITERATION BEGINS HERE

istep = 1; 

while istep <= nsteps
    % stage 1
    t1 = t(istep) + nodes(1)*dt;
    u1 = expL1*u;
    l1 = interface(istep);
    if whichEx==1
        v1 = init_int(l1,t1,t0) + monebybtapi*dot(u1(1:floor(Nherm/2)), basis4y);
    else
        v1 = init_int_2(l1,init_len,t1,t0) + monebybtapi*dot(u1(1:floor(Nherm/2)), basis4y);
    end
    f1 = forcing(u1, l1, v1, t1, bta, N, A1, A2, Nherm, whichEx);

    % stage 2
    t2 = t(istep) + nodes(2)*dt;
    u2 = expL2*u+ dt*(N21*f1);
    l2 = interface(istep) + dt*(aij(2,1)*v1);
    if whichEx==1
        v2 = init_int(l2,t2,t0) + monebybtapi*dot(u2(1:floor(Nherm/2)), basis4y);
    else
        v2 = init_int_2(l2,init_len,t2,t0) + monebybtapi*dot(u2(1:floor(Nherm/2)), basis4y);
    end
    f2 = forcing(u2, l2, v2, t2, bta, N, A1, A2, Nherm, whichEx);

    % next step 
    t(istep + 1) = t(istep) + dt;
    u = expL*u + dt*(M1*f1 + M2*f2);
    interface(istep + 1) = interface(istep) + dt*(b(1)*v1 + b(2)*v2);
    if whichEx==1
        vel(istep + 1) = init_int(interface(istep+1),  t(istep+1), t0) + monebybtapi*dot(u(1:floor(Nherm/2)), basis4y);
    else
        vel(istep + 1) = init_int_2(interface(istep+1), init_len, t(istep+1),t0) + monebybtapi*dot(u(1:floor(Nherm/2)), basis4y);
    end
    % if (mod(istep,500)==0)
    %     plot(u(1:floor(Nherm/2)),'k-o','MarkerSize',2); hold on; plot(basis4y,'r-o','MarkerSize',2); ylim([-100,100]); title(t(istep)); pause(0.1); hold off
    % end
    istep = istep + 1;
end
%pause(1)
hop = 1;

%figure
tref = linspace(t0,Tfinal,100);

if whichEx==1
    subplot(1,2,1)
    plot(t(1:hop:end),vel(1:hop:end),'b-o',t(1:hop:end), interface(1:hop:end),'r-o','MarkerSize',3); grid on; xlim([0,Tfinal+0.5]); ylim([0,max(interface)+1]); pbaspect([2,1,1]); hold on;
    plot(tref, tref,'k', t, ones(size(t)),'k','LineWidth',1); hold on
else
    subplot(1,2,1)
    plot(t(1:hop:end),vel(1:hop:end),'b-o',t(1:hop:end), interface(1:hop:end),'r-o','MarkerSize',3); grid on; xlim([0,Tfinal+0.5]); ylim([0,max(interface)+1]); pbaspect([2,1,1]); hold on;
    plot(tref,2*alp*sqrt(tref),'k', tref,alp./sqrt(tref),'k','LineWidth',1); hold on
end

ylabel('$v(t), x(t)$', 'Interpreter','latex','FontSize',14);
xlabel('$t$','Interpreter','latex','FontSize',14);
%plot(t, vel, 'k-o', t, cos(t), 'b', 'MarkerSize',2, 'LineWidth',1); grid on; axis square;


if whichEx==1
    subplot(1,2,2)
    semilogy(t, abs(vel-1),'b-o', t, abs(interface-t),'r-o','MarkerSize',2); grid on; axis square;
else
    subplot(1,2,2)
    %semilogy(t,abs(interface-2*alp*sqrt(t)),'r-o','MarkerSize',2); grid on; axis square;
    %semilogy(t, abs(vel-alp./sqrt(t)),'m-o','MarkerSize',2); grid on; axis square;
    semilogy(t, abs(vel-alp./sqrt(t)),'b-o',t,abs(interface-2*alp*sqrt(t)),'r-o','MarkerSize',2); grid on; axis square;
end
ylim([1e-6,1e1]);
xlabel('$t$', 'Interpreter','latex','FontSize',14);
ylabel('$|v(t)-\alpha/\sqrt{t}|, |x(t)-2\alpha \sqrt{t}|$','Interpreter','latex','FontSize',14);
hold on
%pbaspect([6,4,1]);
