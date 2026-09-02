function Q = initialize_solution(mesh,Mach,alpha_deg,gamma,N)
assert(N==mesh.N,'N must equal mesh.N.'); assert(Mach>0&&Mach<1,'Baseline supports subsonic Mach.');
a=alpha_deg*pi/180; rho=1; p=1/gamma; u=Mach*cos(a); v=Mach*sin(a); E=p/(gamma-1)+.5*rho*(u^2+v^2);
Q=zeros(mesh.n_elem,N+1,N+1,4); Q(:,:,:,1)=rho; Q(:,:,:,2)=rho*u; Q(:,:,:,3)=rho*v; Q(:,:,:,4)=E;
end
