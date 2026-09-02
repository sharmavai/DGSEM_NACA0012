function E=compute_total_energy(Q,mesh)
w=mesh.w_lgl; E=0;
for e=1:mesh.n_elem, for j=1:length(w), for i=1:length(w)
 E=E+w(i)*w(j)*mesh.J(i,j,e)*Q(e,i,j,4);
end,end,end
end
