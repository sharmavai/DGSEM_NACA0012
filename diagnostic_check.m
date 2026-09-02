clear;clc;close all;
% Stage-1 diagnostic: mesh and state validation only; not a Navier-Stokes validation.
Re=1e4; Mach=.3; alpha=5; gamma=1.4; N=3; mu=1/Re;
mesh=generate_mesh_NACA0012(48,16,20,30,5/Re,N); mesh.Mach=Mach;mesh.alpha_deg=alpha;
Q=initialize_solution(mesh,Mach,alpha,gamma,N);
dt=compute_timestep(Q,mesh,.05,.02,gamma,mu,N);
fprintf('elements=%d, min(J)=%.3e, min(metric length)=%.3e, dt=%.3e\n',mesh.n_elem,min(mesh.J(:)),min(mesh.h_min),dt);
assert(all(isfinite(Q(:)))&&min(mesh.J(:))>0);
fprintf('Stage-1 mesh/state checks passed. Do not run airfoil NS production case until Euler and viscous verification tests are added.\n');
