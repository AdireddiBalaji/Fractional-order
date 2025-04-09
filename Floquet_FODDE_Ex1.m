%--------------------------------------------------------------------------
%-----------  D^2 x(t) + c D^(1/2)x(t)+(delta+epsilon*cos(omega t))x(t)+k1 x(t-tau)=0 using Ode15s floquet
%----------- code by balaji adireddi
%----------- Indian Institute of Technology Hyderabad
clc
clear all %#ok
global c delta epsi k1 CoeMat N omega A_mat B_mat C_mat Td m%#ok

set(groot,'defaultAxesTickLabelInterpreter','latex')
set(groot,'defaulttextInterpreter','latex')
set(groot,'defaultLegendInterpreter','latex')

m=7;   % Number of shape functions for first Galerkin approximation
N=14;  % Number of shape functions in second Galerkin approximation
Td=2*pi; % time delay in the system
% parameters used 
c=0.1;k1=-0.04;omega=1;
Delta=linspace(-0.2,1.5,200); 
Epsi=linspace(0,1.6,200);
alpha=0.5; % farctional order
load(sprintf('Coeff_mat_N%d.mat',N));% load coefficient matrix for 2nd Galerkin appriximation from Maple code
[A_mat, B_mat, C_mat] =frac_sys_mat(alpha,m); % Matrix from 1st Galerkin approximation

%%
red= [1 0 0]; %color specification for the stability chart
figure(11)
hold on
set(groot,'DefaultAxesTickLabelInterpreter','latex');
set(gca,'FontSize',20);
set(get(gca,'YLabel'),'Rotation',0)
xlabel('$\delta$','Interpreter','latex', 'FontSize', 25)
ylabel('$\epsilon$','Interpreter','latex','FontSize', 25)
axis([Delta(1) Delta(end) Epsi(1) Epsi(end)]);

tspan=[0 (2*pi)/Td]; % integration time period 
options = odeset('RelTol',1e-6,'AbsTol',1e-6); % Tolerances
lambda = cell(length(Delta),length(Epsi)); % pre allocation for eigen values
for i=length(Delta):-1:1
    for j=1:length(Epsi)
        delta = Delta(i);
        epsi=Epsi(j);
        init=eye(N+m+1);
        for k=1:N+m+1
            [t,y]=ode15s(@OdeFun,tspan,init(:,k),options);
            M(:,k)=y(end,:)';
        end
        L=max(abs(eig(M))); % find maximum absolute eigen values
        lambda{i,j}=L;
        save('floquet_Frac_alpha0pt5_Delay_Damped','lambda_Re_Im','lambda','i','j','Epsi','c','Delta','k1','N','m'); % save data
        
        if lambda{i,j} <1+1e-6 % stability criterion
            figure(11)
            plot(Delta(i),Epsi(j),'.','MarkerEdge',red,'MarkerFace',red);
       end
    end
end

%% Ode function
function dy  = OdeFun(t,y)
global c delta epsi k1 CoeMat N omega m  A_mat B_mat C_mat Td %#ok
dy = zeros(N+m+1,1);
dy(1:N,1) = CoeMat*[y(1:N+1,1)];
dy(N+1,1)=-c*(Td^2)*C_mat'*y(N+2:N+m+1,1)-Td^2*(delta+epsi*cos(omega*Td*t))*y(1,1)-k1*Td^2*(y(1,1)+y(2,1));
dy(N+2:N+m+1,1)=-Td*(A_mat\B_mat)*y(N+2:N+m+1,1)+(A_mat\C_mat)*y(N+1,1);
end

function [A, B, c] = frac_sys_mat(x,N)
char_freq = 1; decs = 4;
    
    % Comments on parameter values:
    % N is the number of elements, and the size of the approximating matrices;
    % char_freq (characteristic frequency) and decs (width in decades) roughly
    % determine the frequency range over which the approximation is sought.
    % Large values of decs (larger than, say, 8) may cause trouble within finite
    % precision arithmetic in some analyses of some problems. 
    % Note that the role of these parameters is merely to define beta1 and beta2 below.

    yc = char_freq^x;
    beta1 = log10(yc) - 0.4 * decs;
    beta2 = beta1 + 0.6 * decs;
    d = logspace(beta1, beta2, N - 1);
    p = d.^2 ./ (1 + d.^2); 
    p = [0, p, 1];

    A = zeros(N); 
    B = A; 
    c = zeros(N, 1);
    
    for i = 1:N-1
        a = p(i); b = p(i+1); a0 = -a / (b - a);
        a1 = 1 / (b - a);
        
        A(i, i+1) = x * (a0 * (1 - a0) * mb(x - 1, -x - 1, a, b) ...
            + (1 - 2 * a0) * a1 * mb(x, -x - 1, a, b) ...
            - a1^2 * mb(x + 1, -x - 1, a, b));
        
        A(i+1, i+1) = A(i+1, i+1) + x * (a0^2 * mb(x - 1, -x - 1, a, b) ...
            + 2 * a0 * a1 * mb(x, -x - 1, a, b) ...
            + a1^2 * mb(x + 1, -x - 1, a, b));
        
        A(i, i) = A(i, i) + x * ((1 - a0)^2 * mb(x - 1, -x - 1, a, b) ...
            - 2 * (1 - a0) * a1 * mb(x, -x - 1, a, b) ...
            + a1^2 * mb(x + 1, -x - 1, a, b));

        B(i, i+1) = x * (a0 * (1 - a0) * mb(x, -x - 2, a, b) ...
            + (1 - 2 * a0) * a1 * mb(x + 1, -x - 2, a, b) ...
            - a1^2 * mb(x + 2, -x - 2, a, b));

        B(i+1, i) = B(i, i+1); 
        A(i+1, i) = A(i, i+1);

        B(i+1, i+1) = B(i+1, i+1) + x * (a0^2 * mb(x, -x - 2, a, b) ...
            + 2 * a0 * a1 * mb(x + 1, -x - 2, a, b) ...
            + a1^2 * mb(x + 2, -x - 2, a, b));

        B(i, i) = B(i, i) + x * ((1 - a0)^2 * mb(x, -x - 2, a, b) ...
            - 2 * (1 - a0) * a1 * mb(x + 1, -x - 2, a, b) ...
            + a1^2 * mb(x + 2, -x - 2, a, b));

        c(i+1) = c(i+1) + x * (a0 * mb(x - 1, -x - 1, a, b) ...
            + a1 * mb(x, -x - 1, a, b));

        c(i) = c(i) + x * ((1 - a0) * mb(x - 1, -x - 1, a, b) ...
            - a1 * mb(x, -x - 1, a, b));
    end

    a = p(N); 
    b = 1; 
    a0 = 1 / (b - a);
    
    A(N, N) = A(N, N) + a0^2 * x * mb(x - 1, -x + 1, a, b);
    B(N, N) = B(N, N) + a0^2 * x * mb(x, -x, a, b);
    c(N) = c(N) + a0 * x * mb(x - 1, -x, a, b);

    c = c / sqrt(gamma(1 - x) * gamma(1 + x));
end

function s = mb(m, n, a, b)
    if (n > -1) && (m > -1)
        bb = beta(1 + m, 1 + n);
        s = (betainc(b, 1 + m, 1 + n) - betainc(a, 1 + m, 1 + n)) * bb;
    elseif (n <= -1) && (m > 0)
        s = (a^m * (1 - a)^(n + 1) - b^m * (1 - b)^(n + 1)) / (n + 1) ...
            + m / (n + 1) * mb(m - 1, n + 1, a, b);
    else
        s = mb(m, n + 1, a, b) + mb(m + 1, n, a, b);
    end
end
