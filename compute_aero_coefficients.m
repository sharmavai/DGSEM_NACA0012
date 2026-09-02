function [CL,CD,CM] = compute_aero_coefficients(Q,mesh,Mach,alpha_deg,gamma,~,N)
% Pressure contribution only. Do not report this CD as total viscous drag.
assert(N==mesh.N); q=.5*Mach^2; a=alpha_deg*pi/180; Fx=0;Fy=0;M=0; D=mesh.D;w=mesh.w_lgl;
for k=1:numel(mesh.wall_elems)
 e=mesh.wall_elems(k); xe=mesh.x(:,:,e);ye=mesh.y(:,:,e); tx=D*xe(:,1);ty=D*ye(:,1); ds=sqrt(tx.^2+ty.^2).*w;
 nx=-ty./max(ds./w,eps); ny=tx./max(ds./w,eps); % outward normal from fluid sign must be verified per mesh
 for i=1:N+1
  rho=Q(e,i,1,1); u=Q(e,i,1,2)/rho;v=Q(e,i,1,3)/rho;E=Q(e,i,1,4);p=(gamma-1)*(E-.5*rho*(u^2+v^2));
  if rho<=0||p<=0,error('Invalid wall state.');end
  fx=-p*nx(i)*ds(i); fy=-p*ny(i)*ds(i); Fx=Fx+fx;Fy=Fy+fy;M=M+(xe(i,1)-.25)*fy-ye(i,1)*fx;
 end
end
CL=(-Fx*sin(a)+Fy*cos(a))/q; CD=(Fx*cos(a)+Fy*sin(a))/q; CM=M/q;
end
