function [t,x] = LV_SSA(x, c1, c2, c3, T)
% Generate a single trajectory of the stocastic Lotka-Volterra model
% Inputs:
%   x - initial condition.
%   c1, c2, c2 - rate parameters
%   T - duration of trajectory
% Output:
%   t - sequence of reaction times 
%   x - sequence of system configurations, after the rection.

dx = [1, -1, 0;
    0, 1, -1];    

t=0;
j=0;
while t(end) < T
    j=j+1;
    rates = [c1*x(end,1), c2*x(end,1)*x(end,2), c3*x(end,2)];
    escape_rate = sum(rates);
    
    % Sample next event from an exponential distribution
    t(j+1) = t(end) - log(rand)/escape_rate;
    % (alternatively: t(end)=t(end)+exprnd(escape_rate);
    
    % Sample the reaction:
    i = find((escape_rate*rand)<=cumsum(rates), 1);
    x(j+1,:) = x(j,:) + dx(:, i)';
end


