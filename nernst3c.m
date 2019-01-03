function df = nernst3c(~,f,Cp1,Cp2,Cp3,z1,z2,z3,F,R,...
    T,Jv,D1,D2,D3)
% NERNST3C  System of ODE's defined by the extended
% Nersnt-Planck equation (3-component system).
%   For further details on the extended Nernst-Planck
%   equation see [1].
%
%   To see the ODE system type:
%     >> type nernst3c
% 
%   References
% 
%   [1] Morao et al, Ultrafiltration of supercoiled plasmid
%       DNA: Modelling and application, Journal of Membrane
%       Science, 378 (2011) 280 - 289.
%
%   see also nernst4c

% Author: J. Nunes <josenunes34@gmail.com>
% Date: 2013

df=zeros(4,1);

df(1)=-z1*f(1)*F*(Jv*z1*(f(1)-Cp1)/D1+Jv*z2*(f(2)-Cp2)/D2...
    +Jv*z3*(f(3)-Cp3)/D3)/(F*(z1^2*f(1)+z2^2*f(2)...
    +z3^2*f(3))/(R*T))/(R*T)+Jv*(f(1)-Cp1)/D1;

df(2)=-z2*f(2)*F*(Jv*z1*(f(1)-Cp1)/D1+Jv*z2*(f(2)-Cp2)/D2...
    +Jv*z3*(f(3)-Cp3)/D3)/(F*(z1^2*f(1)+z2^2*f(2)...
    +z3^2*f(3))/(R*T))/(R*T)+Jv*(f(2)-Cp2)/D2;

df(3)=-z3*f(3)*F*(Jv*z1*(f(1)-Cp1)/D1+Jv*z2*(f(2)-Cp2)/D2...
    +Jv*z3*(f(3)-Cp3)/D3)/(F*(z1^2*f(1)+z2^2*f(2)...
    +z3^2*f(3))/(R*T))/(R*T)+Jv*(f(3)-Cp3)/D3;

df(4)=(Jv*z1*(f(1)-Cp1)/D1+Jv*z2*(f(2)-Cp2)/D2...
    +Jv*z3*(f(3)-Cp3)/D3)/(F*(z1^2*f(1)+z2^2*f(2)...
    +z3^2*f(3))/(R*T));