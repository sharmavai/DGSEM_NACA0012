function Qg = farfield_bc(~,~,~,Mach,alpha_deg,gamma)
% Robust baseline: prescribed freestream ghost state. Validate reflections before using for unsteady studies.
a=alpha_deg*pi/180; rho=1; p=1/gamma; u=Mach*cos(a); v=Mach*sin(a); E=p/(gamma-1)+.5*rho*(u^2+v^2);
Qg=[rho;rho*u;rho*v;E];
end
