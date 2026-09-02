function Qnew = rk4_step(Q,dt,mesh,gamma,mu,kcond)
% Classical RK4. For Euler, pass mu=0.
k1=rhs_DGSEM_euler(Q,mesh,gamma);
k2=rhs_DGSEM_euler(Q+0.5*dt*k1,mesh,gamma);
k3=rhs_DGSEM_euler(Q+0.5*dt*k2,mesh,gamma);
k4=rhs_DGSEM_euler(Q+dt*k3,mesh,gamma);
Qnew=Q+(dt/6.0)*(k1+2*k2+2*k3+k4);
end
