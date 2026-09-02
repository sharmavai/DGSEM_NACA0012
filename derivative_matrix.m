function [D,M,B]=derivative_matrix(N)
Np1=N+1; [xi,w]=LGL_quadrature(N);
[~,PN,~,~]=arrayfun(@(x) leg_derivs(x,N),xi,'UniformOutput',false);
PN=cell2mat(PN);
D=zeros(Np1);
for i=1:Np1, for j=1:Np1, if i~=j, D(i,j)=PN(i)/(PN(j)*(xi(i)-xi(j))); end, end, end
for i=1:Np1, D(i,i)=-sum(D(i,:)); end
D(1,1)=-N*(N+1)/4; D(end,end)=N*(N+1)/4;
M=diag(w); B=diag([-1;zeros(N-1,1);1]);
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
