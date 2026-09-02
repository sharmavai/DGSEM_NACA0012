function dUdt = rhs_DGSEM_euler(U,mesh,gamma)
% Minimal Euler DGSEM residual with face-loop structure.
assert(N==mesh.N,'N must match mesh.N'); p=N+1; ne=mesh.n_elem;
dUdt=zeros(size(U));
for e=1:ne
 Qe=squeeze(U(e,:,:,:)); xe=mesh.x(:,:,e); ye=mesh.y(:,:,e); Je=mesh.J(:,:,e);
 rho=Qe(:,:,1); rhou=Qe(:,:,2); rhov=Qe(:,:,3); E=Qe(:,:,4);
 u=rhou./rho; v=rhov./rho; pr=(gamma-1)*(E-.5*rho.*(u.^2+v.^2));
 Fx=[rhou; rho.*u.*u+pr; rho.*u.*v; (E+pr).*u];
 Fy=[rhov; rho.*u.*v; rho.*v.*v+pr; (E+pr).*v];
 Fxi=mesh.xi_x(:,:,e).*Fx+mesh.xi_y(:,:,e).*Fy;
 Feta=mesh.eta_x(:,:,e).*Fx+mesh.eta_y(:,:,e).*Fy;
 Res=mesh.D*Fxi+(mesh.D*Feta')';
 dUdt(e,:,:,:)=-Res./Je;
end
end
