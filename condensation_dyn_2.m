
%% Dynamic condensation

[DLL, DLR, DRL, DRR] = condensation_dyn(D, uL, ui, uR);
    % Use as : [DLL, DLR, DRL, DRR, vL, vR] = condensation_dyn(D, uL, ui, uR)
    % D is the original dynamic stiffness matrix, uL, ui, uR the left, inner
    % and right DOFs.
    % This function create the condensed submatrices DLL, DLR, DRL, DRR, and
    % ensures that the new left and right DOFs are uL=1:n and uR=(n+1):(2*n).
    % Resulting matrices are in 'full' format.


% Condensation dynamique dans le cas Fi = 0

m = length(uL);

% Sous matrices élémentaires de D original
DLL=D(uL,uL);
DLi = D(uL,ui);
DLR = D(uL,uR);
DiL = D(ui,uL);
Dii = D(ui,ui);
DiR = D(ui,uR);
DRL=D(uR,uL);
DRi=D(uR,ui);
DRR=D(uR,uR);

% Sous matrices élémentaires de D condensée
DLL = DLL-DLi*inv(Dii)*DiL;
DLR = -DLi*inv(Dii)*DiR;
DRL = -DRi*inv(Dii)*DiL;
DRR = DRR-DRi*inv(Dii)*DiR;

uL = 1:m;
uR = (m+1):(2*m);





