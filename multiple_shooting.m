% clear vars from from workspace and close all figure windows
clear
close all

% load casadi
import casadi.*

% problem constants
T=2;
a=1;
b=-2.694528;
c=-1.155356;

% grid size
N=50;
% length of a control interval
dt = T/N; 

%opti class
opti=Opti();

% decision variables
X=opti.variable(2,N+1);
U=opti.variable(1,N);

% intial condition
X0=[0;0];
opti.subject_to(X(:,1)==X0);

% RK 45 integration
for k=1:N 
    x=X(:,k);
    k1 = ode_fun(x,         U(:,k));
    k2 = ode_fun(x+dt/2*k1, U(:,k));
    k3 = ode_fun(x+dt/2*k2, U(:,k));
    k4 = ode_fun(x+dt*k3,   U(:,k));
    x_next = x + dt/6*(k1+2*k2+2*k3+k4);
% Euler forward    
%   x_next=x+dt*(ode_fun(x,U(:,k)));
% continuity constraints
    opti.subject_to(x_next-X(:,k+1)==0);
end

% Objective function
obj=dt*sum(U.^2);
opti.minimize(obj)

% boundary constriant
opti.subject_to(a*X(1,end)+b*X(2,end)-c==0)

%IPOPT solver
opti.solver('ipopt')
sol=opti.solve();

% Compare the numerical and analytical solution
t=linspace(0,T,N+1);
X=sol.value(X);
U=sol.value(U);
[Xa,Ua]=analytical_solution(t);
plot(Xa(1,:),Xa(2,:))
hold on
scatter(X(1,:),X(2,:))
grid
xlabel('$x_1$')
ylabel('$x_2$')
legend('analytical','single shooting','Location','northwest')

%control input
figure
plot(t,Ua)
hold on
stairs(t,[U NaN])
xlabel('time [s]')
ylabel('f')
legend('analytical','single shooting','Location','northwest')

% Sparsity
figure
Lag=opti.f+ opti.lam_g'*opti.g;
H=hessian(Lag,opti.x);
spy(H) 
legend('Hessian sparsity') 
figure 
Jac=jacobian(opti.g,opti.x);
spy(Jac)
legend('Jacobian sparsity') 
