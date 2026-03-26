
%% Loading a pre-existing unit-cell model :
clear; clc

load UC_model2D_2x2.mat

%load UC_model2D_10x10.mat

freqrange = linspace(5,50000,300); % Given frequency range for these two models


%% Parametrization of the unit-cell
clear; clc

% Geometrical parameters
lx=0.1; % Unit-cell length
ly=0.1; % Unit-cell height
r=0.035; % Inclusion's radius
geom=struct('lx',lx,'ly',ly,'r',r);
geom

% Material parameters
E=210e9; % Young modulus
eta=0.01; % Damping loss factor
nu=0.3; % Poisson coefficient
rho=7800; % Material density
t=0.1; % Plate width
mater=struct('E',E,'nu',nu,'rho',rho,'eta',eta,'t',t);

plotopt=1; % Option to visualize unit-cell

model = generate_model(geom,mater,plotopt);

% Possible frequency parameters for this model :
freqrange = linspace(1,10000,200);


%% Test & view the free eigenmodes of the unit-cell

[PHI,omeg]=eigs(model.K,model.M,20,'sm');
f=sqrt(diag(omeg))/2/pi; phi=PHI(:,8);

visualize_mode(model.mesh,phi,'animate', true)


%% Extract simplified variables and compute waves

K=model.K;
M=model.M;
uL=model.dofs.uL;
ui=model.dofs.ui;
uR=model.dofs.uR;
d=model.geom.lx;

nf=length(freqrange);

% Computing waves from dispersion relation
wavebasis = WFEM(K,M,freqrange,uL,uR,ui,d);

plot(freqrange,wavebasis.k_pos,'k.') % ploting the results Re(k)

%figure
plot(freqrange,abs(wavebasis.lbpos),'k--') % ploting the results Re(k)

%% Forced response
m=length(uL)

N=20; % Number of unit-cells in the finite structure
U = zeros(m,N+1,nf);

% Define the boundary conditions (imposed displacements or forces)
type='F-F'; % Case where forces are applied on both edges of the waveguide
F0 = zeros(m,1);
FN = 1e6*ones(m,1);
BC = [F0 ; FN]; % vector of size 2*m containing either U0 or F0 and UN or FN

for i=1:nf

    % Creating D matrices
    om=sparse(2*pi*freqrange(i));
    D = K - om^2*M;

    % Dynamic condensation (see details in WFEM_MODA function)
    %[DLL, DLR, DRL, DRR, vL, vR] = condensation_dyn(D, uL, ui, uR);
    [DLL, DLR, DRL, DRR,vL, vR] = condensation_dyn_2(D, uL, ui, uR);

    % Defining wave matrices for each frequency
    phi_p = wavebasis.phipos(:,:,i);
    phi_n = wavebasis.phineg(:,:,i);
   % lb = sparsediag(wavebasis.lbpos(:,i));
% alternatively :
    lb = sparse(diag(wavebasis.lbpos(:,i)));

    % Compute the wave amplitudes for the selected boundary conditions
    [qp, qn] = wave_Amplitudes(DLL, DLR, DRL, DRR, N, lb, phi_p, phi_n, BC, type);
% Function [qp, qn]=wave_Amplitudes(DLL, DLR, DRL, DRR, N, lb, phi_p, phi_n, BC, type)
% Determines, for each 'type' of boundary conditions, ex. type='U-U' or 'U-F', 'F-U', 'F-F',
% the wave amplitudes q+, q-, solutions of the problem H*[q+ ; q-] = BC,
% where BC is the imposed boundary conditions vector. BC can be of the form
% [U0 ; UN], [U0 ; FN],  [F0 ; UN] or [F0 ; FN], depending on what is imposed.
% H is a matrix of size 2*m x 2*m, and BC is a vector of size 2*m.
% It returns the two wave amplitudes qp and qn, each of size m.


    % Retrieve the physical displacements from the wave amplitudes
    u_vector = retrieve_U(N, lb, phi_p, phi_n, qp, qn);
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

    U(:,:,i) = u_vector; % collects the displacement vector at ferquency index i.

end

figure(1); plot(real(U(2,:,end)),'k-')

figure(2); semilogy(freqrange,squeeze(abs(U(end,N+1,:))),'k-')

%% Create global matrices for finite periodic structure (case 'F-F')

[Kg, Mg] = assemble_structure(model,N);
% Use: [Kg, Mg] = assemble_structure(UC,N)
% This function takes a unit-cell model as defined above, and returns, for
% a given integer N, the global Stiffness and Mass matrices resulting from
% an assembly of N unit-cells.

tot_dofs=length(Mg);

nui=length([uL ui]);
Fg=zeros(tot_dofs,1);
Fg(end-m+1:end)=FN; % Apply non-zero forces on right edge of the structure
Ug = zeros(tot_dofs,nf);
for i=1:nf
    % Creating D matrices
    om=sparse(2*pi*freqrange(i));
    Dg= Kg - om^2*Mg;
    % Solve the linear harmonic system
    Ug(:,i)=Dg\Fg;
end

figure(1); hold on; plot(real(Ug(2:nui:end,end)),'r--')

figure(2); hold on; semilogy(freqrange,squeeze(abs(Ug(end,:))),'r--')


