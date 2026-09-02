function dt = compute_timestep(Q,mesh,CFL,CFLv,gamma,mu,N)
assert(N==mesh.N,'N must match mesh.N.'); p=N+1; dt=inf;
for e=1:mesh.n_elem
 rho=squeeze(Q(e,:,:,1)); rhou=squeeze(Q(e,:,:,2)); rhov=squeeze(Q(e,:,:,3)); E=squeeze(Q(e,:,:,4));
 if any(rho(:)<=0), error('Non-positive density in element %d.',e); end
 u=rhou./rho; v=rhov./rho; pr=(gamma-1)*(E-.5*rho.*(u.^2+v.^2));
 if any(pr(:)<=0), error('Non-positive pressure in element %d.',e); end
 c=sqrt(gamma*pr./rho);
 ax=abs(mesh.xi_x(:,:,e).*u+mesh.xi_y(:,:,e).*v)+c.*sqrt(mesh.xi_x(:,:,e).^2+mesh.xi_y(:,:,e).^2);
 ae=abs(mesh.eta_x(:,:,e).*u+mesh.eta_y(:,:,e).*v)+c.*sqrt(mesh.eta_x(:,:,e).^2+mesh.eta_y(:,:,e).^2);
 lam=max(ax(:)+ae(:)); dtInv=CFL/(p^2*max(lam,eps));
 if mu>0
  nu=max(mu./rho,[],'all'); g2=max(mesh.xi_x(:,:,e).^2+mesh.xi_y(:,:,e).^2+mesh.eta_x(:,:,e).^2+mesh.eta_y(:,:,e).^2,[],'all');
  dtVis=CFLv/(p^4*max(nu*g2,eps));
 else, dtVis=inf; end
 dt=min(dt,min(dtInv,dtVis));
end
if ~isfinite(dt)||dt<=0,error('Invalid timestep.');end
end
