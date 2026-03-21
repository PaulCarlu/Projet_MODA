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

abs_Lambdas = abs(Lambdas);
tol = 1e-6; % Tolérance pour la recherche de 1/lambda dans abs_Lambdas à partir de lambda

for i=1:n
  [~,wheremin] = min(abs_Lambdas);
  Lambdasmin = Lambdas(wheremin);
  lbpos(i,1) = Lambdasmin;
  phipos(:,i) = PHI(:,wheremin);

  wheremax = find(abs(Lambdas - 1/Lambdasmin) < tol);
  % Test pour savoir si wheremax contient plusieurs indices associés à
  % 1/Lambdasmin pour compléter la paire (Lambdasmin,1/Lambdasmin)
  if (length(wheremax) > 1)
    [~,wheremax] = min(abs(Lambdasmin - 1/Lambdas(wheremax)));
  end
  % ATTENTION : wheremax peut être vide avec 'find'
  % Si c'est le cas, réduire tol pour accepter plus de valeurs
  lbneg(n-i+1,1) = Lambdas(wheremax);
  phineg(:,n-i+1) = PHI(:,wheremax);

  abs_Lambdas(wheremin) = Inf; % Changement dans abs_Lambdas pour pouvoir rechercher denouveau le minimum
end
