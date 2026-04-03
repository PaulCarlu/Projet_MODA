function [PHI,Lambdas] = solve_palindrome_2(A0,A1,A2)
% 'solve_palindrome' solves the quadratic eigenvalue problem (EVP):
% [A0*1/lambda + A1 + A2*lambda] PHI = 0
% It takes square matrices A of dimension n, and uses a linearization
% of the EVP to retrieve [PHI,Lambdas], 
% where PHI=[phi_1, ..., phi_2n] is of size n x 2*n, and
% where Lambdas is a vector of size 2*n

% Variables
[n,m] = size(A0);
matnulle = zeros(n,n);
matiden = eye(n,n);

% Ordre de grandeur de D
ordreGrandeur = norm(A2,'inf')/length(A2);

% Résolution
N = [-A0,matnulle;matnulle,ordreGrandeur*matiden];
L = [A1,A2;ordreGrandeur*matiden,matnulle];
[PSI,LAMBDA] = eig(L,N);

% Attribution
PHI = PSI(1:n,:);
Lambdas = diag(LAMBDA);