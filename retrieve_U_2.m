% Retrieve the physical displacements from the wave amplitudes
function    u_vector = retrieve_U_2(N, lb, phi_p, phi_n, qp, qn)
% The function u_vector = retrieve_U(N, lb, phi_p, phi_n, qp, qn)
% takes as inputs the number of unit-cells (N), the propagation constants
% (lb) and wave vectors (phi) and the wave amplitudes computed previously,
% to return the physical displacements using the Bloch wave decomposition
% method. One has to be careful of the positions n=0 and n=N, since the
% matrix powers operations may result in M^0 matrices.
% The resulting u_vector is in fact a matrix: u(DOF_index, n_index), where 
% DOF_index is the DOF number on the UC edge (uL), while n_index is the
% unit-cell position from 0 (i.e., left of the 1st UC) to N (right of the
% last UC).

    m = length(qp);
    u_vector = zeros(m,N+1);

    u_vector(:,1)= phi_p*qp + phi_n*(lb.^N)*qn;
    u_vector(:,N+1)= phi_p*(lb.^N)*qp + phi_n*qn;


    %u_vector = [phi_p*qp + phi_n*(lb.^N)*qn];

    for n=1:N-1
        u_vector(:,n+1) = phi_p*(lb.^n)*qp + phi_n*(lb.^(N-n))*qn;
        %u_vector = [u_vector,phi_p*(lb.^n)*qp + phi_n*(lb.^(N-n))*qn];

    end
