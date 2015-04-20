function validateModel()
%validateModel produces R^2 and MSE calculations for regression models

BHatOLS = csvread('BHatOLS.csv');
BHatModifiedOLS = csvread('BHatModifiedOLS.csv');
BHatRidge = csvread('BHatRidge.csv');
BHatLasso = csvread('BHatLasso.csv');
G = csvread('ValidationG.csv');
Y = csvread('ValidationY.csv');

YHatOLS = calculateYHat(G, BHatOLS, 'OLS');
YHatModifiedOLS = calculateYHat(G, BHatModifiedOLS, 'ModifiedOLS');
YHatRidge = calculateYHat(G, BHatRidge, 'Ridge');
YHatLasso = calculateYHat(G, BHatLasso, 'Lasso');

RSquaredOLS = calculateRSquared(Y, YHatOLS);
RSquaredModifiedOLS = calculateRSquared(Y, YHatModifiedOLS);
RSquaredRidge = calculateRSquared(Y, YHatRidge);
RSquaredLasso = calculateRSquared(Y, YHatLasso);
MSEOLS = calculateMSE(Y, YHatOLS);
MSEModifiedOLS = calculateMSE(Y, YHatModifiedOLS);
MSERidge = calculateMSE(Y, YHatRidge);
MSELasso = calculateMSE(Y, YHatLasso);

printRSquaredAndMSE(RSquaredOLS, MSEOLS, 'OLS');
printRSquaredAndMSE(RSquaredModifiedOLS, MSEModifiedOLS, 'ModifiedOLS');
printRSquaredAndMSE(RSquaredRidge, MSERidge, 'Ridge');
printRSquaredAndMSE(RSquaredLasso, MSELasso, 'Lasso');
end

function YHat = calculateYHat(G, BHat, model)
%calculateYHat calculates YHat
%Args:
%   G: n x m matrix of genotypes
%   BHat: vector of m beta hats
%   model: string describing model ('OLS', 'Ridge', 'Lasso', ...)
%Returns:
%   YHat: vector of n predicted phenotypes

YHat = G * BHat;
dlmwrite(strcat('YHat', model, '.csv'), YHat);
end

function RSquared = calculateRSquared(Y, YHat)
%calculateRSquared calculates R^2 between Y and YHat
%Args:
%   Y: vector of n true phenotypes
%   YHat: vector of n predicted phenotypes
%Returns:
%   RSquared: R^2

YResid = Y - YHat;
SSResid = sum(YResid.^2);
SStotal = (length(Y) - 1) * var(Y);
RSquared = 1 - SSResid/SStotal;
end

function MSE = calculateMSE(Y, YHat)
%calculateMSE calculates mean squared error (MSE) of YHat
%Args:
%   Y: vector of n true phenotypes
%   YHat: vector of n predicted phenotypes
%Returns:
%   MSE: mean squared error

YResid = Y - YHat;
SSResid = sum(YResid.^2);
MSE = SSResid/length(Y);
end

function printRSquaredAndMSE(RSquared, MSE, model)
%printRSquaredAndMSE prints RSquared and MSE
%Args:
%   RSquared: R^2
%   MSE: mean squared error
%   model: model: string describing model ('OLS', 'Ridge', 'Lasso', ...)

fprintf('%s RSquared = %f\n', model, RSquared);
fprintf('%s MSE = %f\n', model, MSE);
end