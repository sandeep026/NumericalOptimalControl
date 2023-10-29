function [xd] = ode_fun(x,u)
x1=x(1,:);
x2=x(2,:);
xd=[x2;-x2+u];
end

