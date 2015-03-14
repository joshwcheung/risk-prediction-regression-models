%m SNPs, n individuals
%p = vector with m minor allele frequencies (MAF)
%G = n x m matrix with m SNPs, n individuals
%B = vector with m betas, 1% causal
%E = vector with n epsilons
%Y = vector with n phenotypes

%Change m, n, varG
m = 1000;
n = 100000;
varg = 0.7;
varE = 1 - varg;
numCausal = m * 0.01;
causalSNPs = randperm(m, numCausal);

%p ~ U(0.05, 0.5)

p = 0.05 + (0.5 - 0.05)*rand(m,1);
dlmwrite('p.csv', p);

%Training Set

%G ~ B(n, p)

G = NaN(n, m);
for i = 1:n;
    for j = 1:m;
        G(i, j) = binornd(2, p(j));
        %disp([i, j]);
    end
end
dlmwrite('TrainingG.csv', G);

%B ~ N(0, varg/numCausal) for causal SNPs, beta = 0 otherwise

B = zeros(m, 1);
for i = 1:numCausal;
    B(causalSNPs(i)) = normrnd(0, sqrt(varg/numCausal));
end
dlmwrite('TrainingB.csv', B);

%E ~ N(0, varE)

E = normrnd(0, sqrt(varE), n, 1);
dlmwrite('TrainingE.csv', E);

%Y_i = sum(B_j * G_ij) + E_i

Y = NaN(n, 1);
for i = 1:n;
    Gi = (G(i, :)).';
    Y(i) = sum(Gi .* B) + E(i);
end
dlmwrite('TrainingY.csv', Y);

%Validation Set

G = NaN(n, m);
for i = 1:n;
    for j = 1:m;
        G(i, j) = binornd(2, p(j));
        %disp([i, j]);
    end
end
dlmwrite('ValidationG.csv', G);

B = zeros(m, 1);
for i = 1:numCausal;
    B(causalSNPs(i)) = normrnd(0, sqrt(varg/numCausal));
end
dlmwrite('ValidationB.csv', B);

E = normrnd(0, sqrt(varE), n, 1);
dlmwrite('ValidationE.csv', E);

Y = NaN(n, 1);
for i = 1:n;
    Gi = (G(i, :)).';
    Y(i) = sum(Gi .* B) + E(i);
end
dlmwrite('ValidationY.csv', Y);