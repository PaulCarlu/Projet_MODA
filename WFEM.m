function sol = WFEM(K,M,freqrange,uL,uR,ui,dx)
%  wavebasis = WFEM(K,M,freqrange,uL,uR,ui,dx) is a structures which contains :
% 1D matrix : f(frequency_index) the frequencies of the solutions
% 2D matrix : lbpos(wave_index, frequency_index), the propagation constants
% 3D matrix : phipos(DOF_index, wave_index, frequency_index), the eigenvectors
% optional 2D matrix : k_pos(wave_index, frequency_index), which contains
% the wavenumbers of a few mainly propagating solution (for easy plot of
% the dispersion curves)

for i=1:length(freqrange)  % frfr=nbfreq
    % On commence à 2 car si trop proche de 0, une erreur survient
    om=sparse(2*pi*freqrange(i));
    D=K - om^2*M;

    %% Dynamic condensation
    [DLL, DLR, DRL, DRR] = condensation_dyn(D, uL, ui, uR);
    % Use as : [DLL, DLR, DRL, DRR, vL, vR] = condensation_dyn(D, uL, ui, uR)
    % D is the original dynamic stiffness matrix, uL, ui, uR the left, inner
    % and right DOFs.
    % This function create the condensed submatrices DLL, DLR, DRL, DRR, and
    % ensures that the new left and right DOFs are uL=1:n and uR=(n+1):(2*n).
    % Resulting matrices are in 'full' format.


    %% Dispersion relation resolution
    [PHI,Lambdas] = solve_palindrome_2(DRL,DLL+DRR,DLR);
    % 'solve_palindrome' solves the quadratic eigenvalue problem (EVP):
    % [A0*1/lambda + A1 + A2*lambda] PHI = 0
    % It takes square matrices A of dimension n, and uses a linearization
    % of the EVP to retrieve [PHI,Lambdas],
    % where PHI=[phi_1, ..., phi_2n] is of size n x 2*n, and
    % where Lambdas is a vector of size 2*n


    %% Sorting and normalizing the wave basis
    [lbpos, lbneg, phipos, phineg] = wavesorting(Lambdas,PHI);
    % function 'wavesorting' collects the PHI and Lambdas variables as defined
    % produced by solve_palindrome, and does three operations:
    % First, it separates the positive-going and negative-going propagation
    % constants and vectors (associated with abs(lambda)< or > 1).
    % Second, it sorts the two sets (+ and - waves) in a similar way, i.e.,
    % from the most propagating to the most evanescent waves.
    % Third, it normalizes each vector phi to the unity.
    % lbpos and lbneg are vectors of length n
    % phipos, phineg are matrices of size n x n.


    %% Test the palindomic nature of the solutions:
    % test=abs([lbpos 1./lbneg]); % propagation constants (should be =)
    % test=abs([phipos ; phineg]); % wave vectors (same norms)

    %% Fill the wavebasis structure at each frequency
    sol.lbpos(:,i)=lbpos;
    sol.lbneg(:,i)=lbneg;
    sol.phipos(:,:,i)=phipos;
    sol.phineg(:,:,i)=phineg;
    sol.k_pos(:,i)=1i*log(lbpos(1:min(2,length(lbpos))))/dx; % extracts main propagating wavenumbers for easy plot
    sol.f=freqrange(i);
end

end









