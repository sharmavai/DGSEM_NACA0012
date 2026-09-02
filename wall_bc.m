function Qg = wall_bc(Qi,gamma)
% Euler slip-wall ghost state. No-slip requires a verified viscous DG formulation.
Qi=Qi(:); rho=Qi(1); if rho<=0,error('Non-positive density.');end
u=Qi(2)/rho; v=Qi(3)/rho; E=Qi(4); p=(gamma-1)*(E-.5*rho*(u^2+v^2)); if p<=0,error('Non-positive pressure.');end
Qg=Qi; % Normal reflection must be performed by caller with normal; retained only for compatibility.
end
