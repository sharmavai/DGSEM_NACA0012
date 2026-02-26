function E_total = compute_total_energy(U, mesh)
% COMPUTE_TOTAL_ENERGY  Total energy in domain (conservation check).

E_total = sum(U(:,:,:,4), 'all');
end
