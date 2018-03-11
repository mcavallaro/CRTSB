function x = LV_SSA_stepper(x, c1, c2, c3, dt)
% Evolves a single Lotka-Volterra model trajectory for a certain time
% Inputs:
%   x - initial condition.
%   c1, c2, c2 - rate parameters
%   dt - duration of trajectory
% Output:
%   x - system configurations after time dt.

dx = [1, -1, 0;
      0, 1, -1];

t = 0;
while t < dt   
    rates = [c1*x(1), c2*x(1)*x(2), c3*x(2)];
    escape_rate = sum(rates);
    
    % Sample next event from an exponential distribution
    t = t - log(rand)/escape_rate;
    % (alternatively: t=t+exprnd(escape_rate);
    
    % Sample the reaction:
    i = find((escape_rate*rand)<=cumsum(rates), 1);
    x = x + dx(:, i)';
end


