%--------------------------------------------------------------------------
%----------- fractional DE D^2 x(t) + c D^alpha x(t) + (delta + epsilon*cos(omega t))x(t) + d*x(t - 2*pi) = 0 using RK4
%----------- code by balaji adireddi
%----------- Indian Institute of Technology Hyderabad
clc
clear all %#ok
set(groot,'defaultAxesTickLabelInterpreter','latex')
set(groot,'defaulttextInterpreter','latex')
set(groot,'defaultLegendInterpreter','latex')
global c delta epsi d alpha A_mat B_mat C_mat m omega tau CoeMat N %#ok
% Parameters
alpha = 0.5;  % Fractional order
tau = 2*pi;   % Delay of 2*pi
dt = 0.001;        % Time step
tf = 100;         % Final time
t = 0:dt:tf+tau;      % Time vectorSimpson_1_3

c = 0.1;          % Coefficient of fractional derivative
delta = 1;      % Delta coefficient
epsi = 0.4;       % Epsilon coefficient
omega = 1;        % Omega (frequency)
d = -0.04;       % d coefficient for delayed term

% Initial conditions
x = zeros(1, length(t));      % Displacement
x_dot = zeros(1, length(t));  % Velocity
D_alpha = zeros(1, length(t));% Fractional derivative

tau_pos = ceil(tau/dt); % Delay in terms of position index

for i = 1:tau_pos
    x(i+1) = 1;   
    x_dot(i+1) = 0; 
end

% Runge-Kutta 4th order (RK4) Method
for ii = tau_pos:length(t)-1
    k1_x = x_dot(ii);
    k1_x_dot = -c*D_alpha(ii)-(delta+epsi*cos(omega*t(ii)))*x(ii)-d*x(ii+1-tau_pos);

    k2_x = x_dot(ii)+0.5*dt*k1_x_dot;
    k2_x_dot = -c*D_alpha(ii)-(delta+epsi*cos(omega*(t(ii)+0.5*dt)))*(x(ii)+0.5*dt*k1_x)- d*x(ii +1- tau_pos);

    k3_x = x_dot(ii)+0.5*dt*k2_x_dot;
    k3_x_dot = -c*D_alpha(ii)-(delta+epsi*cos(omega*(t(ii)+0.5*dt)))*(x(ii)+0.5*dt*k2_x)-d*x(ii+1-tau_pos);

    k4_x = x_dot(ii)+dt*k3_x_dot;
    k4_x_dot = -c*D_alpha(ii)-(delta+epsi*cos(omega*(t(ii)+ dt)))*(x(ii)+dt*k3_x)-d*x(ii +1- tau_pos);

    % Update x and x_dot using RK4
    x(ii+1) = x(ii)+(dt/6)*(k1_x+2*k2_x+2*k3_x+k4_x);
    x_dot(ii+1) = x_dot(ii)+(dt/6)*(k1_x_dot+2*k2_x_dot+2*k3_x_dot+k4_x_dot);

    % Compute fractional derivative D^alpha
    sum_D_alpha = 0;
    for jj = 1:ii
        sum_D_alpha = sum_D_alpha+((x_dot(ii+1-jj)+x_dot(ii+1-(jj-1)))/2)/(((t(jj)/2+t(jj+1)/2)^alpha))*dt;
    end
    D_alpha(ii+1) = 1/gamma(1-alpha)*sum_D_alpha;
end
toc
%%
N=14;
m=7;
load(sprintf('Coeff_mat_N%d.mat',N));% load coefficient matrix for Galerkin
[A_mat, B_mat, C_mat] = frac_sys_mat(alpha, m);  % Assuming this function is defined elsewhere
tspan = [0 tf];
init=zeros(N+m+1,1);
init(1)=1;
lags = 1;
% Ode solution
optODE = odeset('AbsTol',1e-6,'RelTol',1e-6);
optDDE = ddeset('AbsTol',1e-6,'RelTol',1e-6);
tic
sol_Ode=ode15s(@OdeFun,tspan,init,optODE);
toc
% Solve using dde23
tic;
sol_Dde = dde23(@dde_Fun, 1, @history, tspan, optDDE);
toc;

%% Plot results from dde23 and RK4 and ode
figure(7)
plot(t(1:end)-tau,x(1:end),'-b','LineWidth',2);
hold on
% plot(sol_Dde.x*tau, sol_Dde.y(1,:), '--b', 'LineWidth', 2);
plot(sol_Ode.x*tau,sol_Ode.y(1,:),'--r','LineWidth', 2)
xlabel('$t$', 'Interpreter', 'latex');
ylabel('$x(t)$', 'Interpreter', 'latex', 'Rotation', 0);
legend('Direct numerical integration','Proposed method','Orientation','vertical','Location','best','Interpreter','latex', 'FontSize', 16);
box on
grid on;
xlim([0 100])
% xticks([0 20 40 60 80 100])
% xticks([0 25 50 75 100])
set(gca,'fontsize',20)

%% Ode function
function dy  = OdeFun(t,y)
global c delta epsi d CoeMat N omega m  A_mat B_mat C_mat tau %#ok
dy = zeros(N+m+1,1);
dy(1:N,1) = CoeMat*[y(1:N+1,1)];
dy(N+1,1)=-c*tau^2*C_mat'*y(N+2:N+m+1,1)-tau^2*(delta+epsi*cos(omega*tau*t))*y(1,1)-d*tau^2*(y(1,1)+y(2,1));
dy(N+2:N+m+1,1)=-tau*(A_mat\B_mat)*y(N+2:N+m+1,1)+(A_mat\C_mat)*y(N+1,1);
end

function du = dde_Fun(t,u,Z)
global c delta epsi d A_mat B_mat C_mat m omega tau; %#ok
du = zeros(m+2,1);
du(1) = u(2);
du(2) = -c*tau^2*C_mat'*u(3:m+2)-tau^2*(delta +epsi*cos(omega*tau*t))*u(1)-d*tau^2*Z(1);
du(3:m+2) = -tau*(A_mat\B_mat)*u(3:m+2)+(A_mat\C_mat)*u(2);
end

% History function for initial conditions
function s = history(~)
global A_mat B_mat C_mat m tau; %#ok
s = zeros(m+2, 1);
s(1) = 1;   % Initial condition for x(t)
s(2) = 0;   % Initial condition for x'(t)
s(3:m+2) =0; % Initial conditions for fractional derivatives
end
