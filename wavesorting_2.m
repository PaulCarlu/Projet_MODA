function [lbpos, lbneg, phipos, phineg] = wavesorting_2(Lambdas,PHI)
% function 'wavesorting' collects the PHI and Lambdas variables as defined
% produced by solve_palindrome, and does three operations:
% First, it separates the positive-going and negative-going propagation
% constants and vectors (associated with abs(lambda)< or > 1).
% Second, it sorts the two sets (+ and - waves) in a similar way, i.e.,
% from the most propagating to the most evanescent waves.
% Third, it normalizes each vector phi to the unity.
% lbpos and lbneg are vectors of length n
% phipos, phineg are matrices of size n x n.
n = length(Lambdas); n=n/2;
lbpos = zeros(n,1);
lbneg = zeros(n,1);
phipos = zeros(n,n);
phineg = zeros(n,n);

% Tri de 0 à Inf pour récupérer les lambdas des ondes progressives
abs_Lambdas = abs(Lambdas);
[abs_Lambdas,I]=sort(abs_Lambdas);

lbpos = Lambdas(I(1:n));
phipos = PHI(:,I(1:n));
phineg = PHI(:,I(n+1:2*n));

% Tri dans le bon sens de lbpos
[lbpos,I]=sort(lbpos,'descend');
phipos = phipos(:,I);
phineg = phineg(:,I);

% Comme seulement les lambdas des ondes progressives vont etre utilisés
% on peut reconstruire lbneg à partir de lbpos
lbneg = 1./lbpos;

% Normalisation des phi
for i=1:n
    phipos(:,i) = phipos(:,i)/norm(phipos(:,i));
    phineg(:,i) = phineg(:,i)/norm(phineg(:,i));
end