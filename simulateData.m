function simulateData(m, n, varG)
%simulateData simulates genotype and phenotype data
%   stores p, G, B, E, and Y in .csv files for training and validation sets
%   p: vector of m minor allele frequencies (MAFs)
%   G: n x m matrix of genotypes
%   B: vector of m betas, the SNP effect sizes, where 1% are causal
%   E: vector of n epsilons, the environmental (non-genetic) contribution
%   Y: vector of n phenotypes
%Args:
%   m: the number of SNPs
%   n: the number of individuals
%   varG: the variance in phenotype due to genetics

varE = 1 - varG;
%Assume 1% of SNPs are causal
numCausal = m * 0.01;
causalSNPs = randperm(m, numCausal);

%p ~ U(0.05, 0.5)
p = 0.05 + (0.5 - 0.05)*rand(m,1);
dlmwrite('p.csv', p);

trainingG = simulateG(m, n, p, 'Training');
trainingB = simulateB(m, numCausal, causalSNPs, varG, 'Training');
trainingE = simulateE(n, varE, 'Training');
calculateY(n, trainingG, trainingB, trainingE, 'Training');

validationG = simulateG(m, n, p, 'Validation');
validationB = simulateB(m, numCausal, causalSNPs, varG, 'Validation');
validationE = simulateE(n, varE, 'Validation');
calculateY(n, validationG, validationB, validationE, 'Validation');
end

function G = simulateG(m, n, p, set)
%simulateG simulates genotypes
%   G ~ B(n, p)
%Args:
%   m: the number of SNPs
%   n: the number of individuals
%   p: vector of m minor allele frequencies (MAFs)
%   set: string describing set ('Training' or 'Validation')
%Returns:
%   G: n x m matrix of genotypes

G = NaN(n, m);
for i = 1:n;
    for j = 1:m;
        G(i, j) = binornd(2, p(j));
    end
end
dlmwrite(strcat(set, 'G.csv'), G);
end

function B = simulateB(m, numCausal, causalSNPs, varG, set)
%simulateB simulates betas
%   B ~ N(0, varg/numCausal) for causal SNPs, beta = 0 otherwise
%Args:
%   m: the number of SNPs
%   numCausal: the number of causal SNPs
%   causalSNPs: vector of positions of causal SNPs
%   varG: the variance in phenotype due to genetics
%   set: string describing set ('Training' or 'Validation')
%Returns:
%   B: vector of m betas, the SNP effect sizes, where 1% are causal

B = zeros(m, 1);
for i = 1:numCausal;
    B(causalSNPs(i)) = normrnd(0, sqrt(varG/numCausal));
end
dlmwrite(strcat(set, 'B.csv'), B);
end

function E = simulateE(n, varE, set)
%simulateE simulates epsilons
%   E ~ N(0, varE)
%Args:
%   n: the number of individuals
%   varE: the variance due to environment
%   set: string describing set ('Training' or 'Validation')
%Returns:
%   E: vector of n epsilons, the environmental (non-genetic) contribution

E = normrnd(0, sqrt(varE), n, 1);
dlmwrite(strcat(set, 'E.csv'), E);
end

function calculateY(n, G, B, E, set)
%calculateY calculates phenotypes
%   Y_i = sum(B_j * G_ij) + E_i
%Args:
%   n: the number of individuals
%   G: n x m matrix of genotypes
%   B: vector of m betas, the SNP effect sizes, where 1% are causal
%   E: vector of n epsilons, the environmental (non-genetic) contribution
%   set: string describing set ('Training' or 'Validation')

Y = NaN(n, 1);
for i = 1:n;
    Gi = (G(i, :)).';
    Y(i) = sum(Gi .* B) + E(i);
end
dlmwrite(strcat(set, 'Y.csv'), Y);
end