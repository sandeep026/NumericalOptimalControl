function [x,u,J] = analytical_solution(t)
u=[exp(t)/4-1/2];
x=[-3/8*exp(-t)+exp(t)/8-1/2*t+1/4;3/8*exp(-t)+1/8*exp(t)-1/2];
J=1/32*(8*2-8*exp(2)+exp(2*2))-1/32*(8*0-8*exp(0)+exp(2*0));
end

