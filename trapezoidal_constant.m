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
%grid size
N=50;
% length of a control interval
dt = T/N;

% opti class
opti=Opti();

% decision vectors
X=opti.variable(2,N+1);
U=opti.variable(1,N);

% initial condition
X0=[0;0];
opti.subject_to(X(:,1)-X0==0)


% Collocation constraint 
for k=1:N 
opti.subject_to(X(:,k+1)-X(:,k)-0.5*dt*(ode_fun(X(:,k),U(:,k))+ode_fun(X(:,k+1),U(:,k)))==0)
end

% objective function - Lagrange term
s=0;
for k=1:N
s=s+0.5*dt*(U(:,k)'*U(:,k)+U(:,k)'*U(:,k));
end
obj=s;
opti.minimize(obj)

% boundary condition
opti.subject_to(a*X(1,end)+b*X(2,end)-c==0)

% NLP solver
opti.solver('ipopt')
sol=opti.solve();

% Compare the numerical and analytical solution
t=linspace(0,2,N+1);
X=sol.value(X);
U=sol.value(U);
[Xa,Ua,Ja]=analytical_solution(t);
plot(Xa(1,:),Xa(2,:))
hold on
scatter(X(1,:),X(2,:))
grid
xlabel('$x_1$')
ylabel('$x_2$')
legend('analytical','Trapezoidal constant control','Location','northwest')
title(['Objective error=' num2str(Ja-sol.value(obj))])

%Control input 
figure
plot(t,Ua)
hold on
stairs(t,[U NaN])
xlabel('time [s]')
ylabel('f')
legend('analytical','Trapezoidal constant control','Location','northwest')

% Sparsity pattern
figure
Lag=opti.f+ opti.lam_g'*opti.g;
H=hessian(Lag,opti.x);
spy(H) 
legend('Hessian sparsity') 
figure 
Jac=jacobian(opti.g,opti.x);
spy(Jac)
legend('Jacobian sparsity') 
