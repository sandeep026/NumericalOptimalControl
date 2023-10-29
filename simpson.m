% clear workspace and close existing figure windows
clear
close all

% load casadi package
import casadi.*

% problem parameters
T=2;
a=1;
b=-2.694528;
c=-1.155356;
% grid size
N=50;
% length of a control interval
dt = T/N; 

% opti class
opti=Opti();

% decision vectors
X=opti.variable(2,N+1);
U=opti.variable(1,N+1);

% initial condition
X0=[0;0];
opti.subject_to(X(:,1)-X0==0)

% Collocation constraint 
for k=1:N 
Xmid=0.5*(X(:,k)+X(:,k+1))+dt/8*(ode_fun(X(:,k),U(:,k))-ode_fun(X(:,k+1),U(:,k+1)));
Umid=0.5*(U(:,k)+U(:,k+1));
%simpson quadrature
opti.subject_to( X(:,k+1)-X(:,k)-dt/6*(ode_fun(X(:,k),U(:,k))+4*ode_fun(Xmid,Umid)+ode_fun(X(:,k+1),U(:,k+1)))==0)
end

% objective function - Lagrange term
s=0;
for k=1:N
Xmid=0.5*(X(:,k)+X(:,k+1))+dt/8*(ode_fun(X(:,k),U(:,k))-ode_fun(X(:,k+1),U(:,k+1)));
Umid=0.5*(U(:,k)+U(:,k+1));
s=s+dt/6*(obj(X(:,k),U(k))+obj(X(k+1),U(k+1))+4*obj(Xmid,Umid));
end
obj=s;
opti.minimize(obj)

% boundary condition
opti.subject_to(a*X(1,end)+b*X(2,end)-c==0)

% NLP solver
opti.solver('ipopt')
sol=opti.solve();

% Compare the numerical and analytical solution
t=linspace(0,T,N+1);
X=sol.value(X);
U=sol.value(U);
[Xa,Ua,Ja]=analytical_solution(t);
plot(Xa(1,:),Xa(2,:))
hold on
scatter(X(1,:),X(2,:))
grid
xlabel('$x_1$')
ylabel('$x_2$')
legend('analytical','simpson','Location','northwest')
title(['Objective error=' num2str(Ja-sol.value(obj))])
 
%Control input 
figure
plot(t,Ua)
hold on
plot(t,U)
xlabel('time [s]')
ylabel('f')
grid
legend('analytical','simpson','Location','northwest')

% Sparsity pattern
figure
Lag=opti.f+ opti.lam_g'*opti.g;
H=hessian(Lag,opti.x);
spy(H)
grid 
legend('Hessian sparsity') 
figure 
Jac=jacobian(opti.g,opti.x);
spy(Jac)
grid
legend('Jacobian sparsity') 
