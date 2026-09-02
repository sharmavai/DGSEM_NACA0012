function [Cp,xs,ys] = compute_pressure_coefficient(Q,mesh,gamma,Mach,N)
assert(N==mesh.N); pInf=1/gamma;q=.5*Mach^2; xs=[];ys=[];Cp=[];
for k=1:numel(mesh.wall_elems)
 e=mesh.wall_elems(k); ids=1:N+1; if k>1,ids=2:N+1;end
 for i=ids
  rho=Q(e,i,1,1);u=Q(e,i,1,2)/rho;v=Q(e,i,1,3)/rho;E=Q(e,i,1,4);p=(gamma-1)*(E-.5*rho*(u^2+v^2));
  if rho<=0||p<=0,error('Invalid surface state.');end
  xs(end+1,1)=mesh.x(i,1,e);ys(end+1,1)=mesh.y(i,1,e);Cp(end+1,1)=(p-pInf)/q;
 end
end
end
