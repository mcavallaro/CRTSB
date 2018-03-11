function dx = LV_derivative(x, c1, c2, c3)
% Function to be used with the ODE solver.
% Inputs:
%   x - Independent variable.
%   c1, c2, c2 - rate parameters
% Output:
%   dx - First derivative: the rate of change of the populations

  dx = [0; 0];

  dx(1) = c1*x(1)-c2*x(1)*x(2);
  dx(2) = c2*x(1)*x(2)-c3*x(2);