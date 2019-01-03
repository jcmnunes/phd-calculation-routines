function [tout, yout] = rk4arg(FunFcn, tspan, y0, nstep, varargin)
% RK4ARG  Runge-Kutta 4th order method to integrate ODE's.
%   The function was adaptaded from [1] and slightly modified
%   to accept varargin. 
%
%   References
%
%   [1] D. Arnold and J. C. Polking, Ordinary differential 
%   equations using MATLAB, Prentice Hall, 2003

% Initialization
t0 = tspan(1);
tfinal = tspan(2);
pm = sign(tfinal - t0); % Which way are we computing?
ssize = (tfinal - t0)/nstep;
if ssize < 0, ssize = -ssize; end
h = pm*ssize;
t = t0;
y = y0(:);
% We need to compute the number of steps.
dt = abs(tfinal - t0);
N = floor(dt/ssize) + 1;
if (N-1)*ssize < dt
N = N + 1;
end
% Initialize the output.
tout = zeros(N,1);
tout(1) = t;
yout = zeros(N,size(y,1));
yout(1,:) = y.';
k = 1;
% The main loop
while (k < N)
if pm*(t + h - tfinal) > 0
h = tfinal - t;
tout(k+1) = tfinal;
else
tout(k+1) = t0 +k*h;
end
k = k + 1;
% Compute the slopes
s1 = feval(FunFcn, t, y,varargin{:}); s1 = s1(:);
s2 = feval(FunFcn, t + h/2, y + h*s1/2, varargin{:}); s2 = s2(:);
s3 = feval(FunFcn, t + h/2, y + h*s2/2, varargin{:}); s3 = s3(:);
s4 = feval(FunFcn, t + h, y + h*s3, varargin{:}); s4 = s4(:);
y = y + h*(s1 + 2*s2 + 2*s3 +s4)/6;
t = tout(k);
yout(k,:) = y.';
end