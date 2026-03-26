function [qp, qn] = wave_Amplitudes_2(DLL, DLR, DRL, DRR, N, lb, phi_p, phi_n, BC, type)
% Function [qp, qn]=wave_Amplitudes(DLL, DLR, DRL, DRR, N, lb, phi_p, phi_n, BC, type)
% Determines, for each 'type' of boundary conditions, ex. type='U-U' or 'U-F', 'F-U', 'F-F',
% the wave amplitudes q+, q-, solutions of the problem H*[q+ ; q-] = BC,
% where BC is the imposed boundary conditions vector. BC can be of the form
% [U0 ; UN], [U0 ; FN],  [F0 ; UN] or [F0 ; FN], depending on what is imposed.
% H is a matrix of size 2*m x 2*m, and BC is a vector of size 2*m.
% It returns the two wave amplitudes qp and qn, each of size m.
m = length(lb);
switch type
  case 'U-U'
        H = [phi_p,phi_n*lb^N ; phi_p*lb^N, phi_n];

  case 'F-F'
        H = [DLL*phi_p + DLR*phi_p*lb,  DLL*phi_n*lb^N + DLR*phi_n*lb^(N-1) ...
            ; DRR*phi_p*lb^N + DRL*phi_p*lb^(N-1), DRR*phi_n + DRL*phi_n*lb];

  case 'U-F'
        H = [ phi_p, phi_n*lb^N ...
            ; DRR*phi_p*lb^N + DRL*phi_p*lb^(N-1), DRR*phi_n + DRL*phi_n*lb];

  case 'F-U'
        H = [DLL*phi_p + DLR*phi_p*lb,  DLL*phi_n*lb^N + DLR*phi_n*lb^(N-1) ...
            ; phi_p*lb^N, phi_n];
end

Qv = H\BC;
qp = Qv(1:m); 
qn = Qv(m+1:2*m);