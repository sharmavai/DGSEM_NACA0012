function [xi,w]=LGL_quadrature(N)
Np1=N+1; xi=zeros(Np1,1); xi(1)=-1; xi(end)=1;
for j=1:N-1
 x=-cos(j*pi/N);
 for k=1:100, [~,~,dP,d2P]=leg_derivs(x,N); dx=-dP/d2P; x=x+dx; if abs(dx)<1e-15,break,end,end
 xi(j+1)=x;
end
w=zeros(Np1,1);
for i=1:Np1, [~,Pn,~,~]=leg_derivs(xi(i),N); w(i)=2/(N*Np1*Pn^2); end
end
function [Pnm1,Pn,dPn,d2Pn]=leg_derivs(x,N)
if N==0, Pnm1=0;Pn=1;dPn=0;d2Pn=0;return,end
if N==1, Pnm1=1;Pn=x;dPn=1;d2Pn=0;return,end
Pkm2=1;Pkm1=x;dPkm2=0;dPkm1=1;d2Pkm2=0;d2Pkm1=0;
for k=2:N
 Pk=((2*k-1)*x*Pkm1-(k-1)*Pkm2)/k;
 dPk=((2*k-1)*(Pkm1+x*dPkm1)-(k-1)*dPkm2)/k;
 d2Pk=((2*k-1)*(2*dPkm1+x*d2Pkm1)-(k-1)*d2Pkm2)/k;
 Pkm2=Pkm1;Pkm1=Pk;dPkm2=dPkm1;dPkm1=dPk;d2Pkm2=d2Pkm1;d2Pkm1=d2Pk;
end
Pnm1=Pkm2;Pn=Pkm1;dPn=dPkm1;d2Pn=d2Pkm1;
end
