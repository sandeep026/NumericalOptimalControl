% clear vars from from workspace and close all figure windows
clear
close all

% load casadi
import casadi.*

% parameters
T=2;
a=1;
b=-2.694528;
c=-1.155356;
%grid size
N=80;

%opti class
opti=Opti();

%decision variables
W=opti.variable(3,N+1);
X=W(1:2,:)
U=W(end,:)


%boundary condition 
X0=[0;0];
opti.subject_to(X(:,1)-X0==0)

%LGL grid points and quadrature weights
[tau,wi]=legslb(N+1);
tau=tau';
wi=wi';
% differentiation matrix
D=legslbdiff(N+1,tau);

%Collocate dynamics
rhs=(ode_fun(X,U)')';
lhs=2/T*(D*X')';
error=lhs-rhs;
opti.subject_to(error(:)==0)

%objective function
obj=T/2*sum(wi.*U.^2);
opti.minimize(obj)

%boundary condition
opti.subject_to(a*X(1,end)+b*X(2,end)-c==0)

%NLP solver
opti.solver('ipopt')
sol=opti.solve();

%Compare numerical and analytical solution
t=T/2*(tau+1);
X=sol.value(X);
U=sol.value(U);
[Xa,Ua]=analytical_solution(t);
plot(Xa(1,:),Xa(2,:))
hold on
scatter(X(1,:),X(2,:))
grid
xlabel('$x_1$')
ylabel('$x_2$')
legend('analytical','LGL pseudospectral','Location','northwest')

%control input
figure
plot(t,[Ua],'-')
hold on
plot(t,[U ],'--o')
legend('analytical','LGL pseudospectral','Location','northwest')
grid
xlabel('$t~[s]$')
ylabel('$u$')

%sparsity
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



