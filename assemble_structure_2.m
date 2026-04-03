%% Create global matrices for finite periodic structure (case 'F-F')

function [Kg, Mg] = assemble_structure_2(model,N)
% Use: [Kg, Mg] = assemble_structure(UC,N)
% This function takes a unit-cell model as defined above, and returns, for
% a given integer N, the global Stiffness and Mass matrices resulting from
% an assembly of N unit-cells.

K=model.K;
M=model.M;
uL=model.dofs.uL;
ui=model.dofs.ui;
uR=model.dofs.uR;

dim_K=length(K);
nuL = length(uL);
nui = length(ui);

dim_cell = 2*nuL+nui;
dim_sys = N*dim_cell - (N-1)*nuL;
Kg = sparse(dim_sys,dim_sys);
Mg = sparse(dim_sys,dim_sys);


for i=0:N-1
    I = (nuL+nui)*i+1;
    I_p1 = I+dim_cell-1;
    ind = I:I_p1;
    length(ind);
    Kg(ind,ind) = Kg(ind,ind) + K;
    Mg(ind,ind) = Mg(ind,ind) + M;
end

