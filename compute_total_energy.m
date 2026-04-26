% =========================================================================
%  compute_total_energy.m
%  Total energy in domain — conservation check.
%
%  FIX: Previous code did sum(U(:,:,:,4),'all') without Jacobian or
%       quadrature weighting — not a meaningful integral.
%       Corrected to: E = sum_{e,i,j}  w_i * w_j * J(i,j,e) * rhoE(i,j,e)
% =========================================================================

function E_total = compute_total_energy(Q, mesh)
% Q:    (n_elem, Np1, Np1, 4)
% mesh: struct with J (Np1,Np1,n_elem) and w_lgl (Np1,1)

w = mesh.w_lgl;
n_elem = mesh.n_elem;
Np1 = length(w);

E_total = 0;
for e = 1:n_elem
    for jj = 1:Np1
        for ii = 1:Np1
            E_total = E_total + w(ii) * w(jj) * mesh.J(ii,jj,e) * Q(e,ii,jj,4);
        end
    end
end
end
