function estimateBetas()
%estimateBetas uses regression models to calculate beta hat values
%   reads 'TrainingG.csv' and 'TrainingY.csv' to calculate beta hat values 
%   using Ordinary Least Squares (OLS), Ridge, and Lasso; also produces a
%   modified OLS output

G = csvread('TrainingG.csv');
Y = csvread('TrainingY.csv');

%OLS
%BHat = (inv(G.'G))G.'Y
%(inv(G.'G))G.' = pinv(G), the Moore-Penrose pseudoinverse of G
%BHat = pinv(G) * Y

BHatOLS = pinv(G) * Y;
dlmwrite('BHatOLS.csv', BHatOLS);

%Modified OLS
%1% of SNPs are causal
%Include SNPs with highest effect sizes (2% = 20 SNPs)

absBHat = abs(BHatOLS);
[~, sortingIndeces] = sort(absBHat, 'descend');
maxIndeces = sortingIndeces(1:20);
BHatModifiedOLS = zeros(1000, 1);
for i = 1:length(maxIndeces);
    BHatModifiedOLS(maxIndeces(i)) = BHatOLS(maxIndeces(i));
end
dlmwrite('BHatModifiedOLS.csv', BHatModifiedOLS);

%Ridge

lambda = 0.002;
BHatRidge = ridge(Y, G, lambda);
dlmwrite('BHatRidge.csv', BHatRidge);

%Lasso

[beta, ~] = lasso(G, Y);
%Change this depending on desired lambda
BHatLasso = beta(:, 93);
dlmwrite('BHatLasso.csv', BHatLasso);
end