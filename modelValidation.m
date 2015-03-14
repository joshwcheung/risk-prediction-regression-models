BHatOLS = csvread('BHatOLS.csv');
BHatModifiedOLS = csvread('BHatModifiedOLS.csv');
BHatRidge = csvread('BHatRidge.csv');
BHatLasso = csvread('BHatLasso.csv');
G = csvread('ValidationG.csv');
Y = csvread('ValidationY.csv');

YHatOLS = G * BHatOLS;
dlmwrite('YHatOLS.csv', YHatOLS);

YHatModifiedOLS = G * BHatModifiedOLS;
dlmwrite('YHatModifiedOLS.csv', YHatModifiedOLS);

YHatRidge = G * BHatRidge;
dlmwrite('YHatRidge.csv', YHatRidge);

YHatLasso = G * BHatLasso;
dlmwrite('YHatLasso.csv', YHatLasso);

%R^2 and MSE

%OLS

Yresid = Y - YHatOLS;
SSresid = sum(Yresid.^2);
SStotal = (length(Y) - 1) * var(Y);
Rsquared = 1 - SSresid/SStotal;
MSE = SSresid/length(Y);
disp('OLS Rsquared = ');
disp(Rsquared);
disp('OLS MSE = ');
disp(MSE);

%Modified OLS

Yresid = Y - YHatModifiedOLS;
SSresid = sum(Yresid.^2);
SStotal = (length(Y) - 1) * var(Y);
Rsquared = 1 - SSresid/SStotal;
MSE = SSresid/length(Y);
disp('Modified OLS Rsquared = ');
disp(Rsquared);
disp('Modified OLS MSE = ');
disp(MSE);

%Ridge

Yresid = Y - YHatRidge;
SSresid = sum(Yresid.^2);
SStotal = (length(Y) - 1) * var(Y);
Rsquared = 1 - SSresid/SStotal;
MSE = SSresid/length(Y);
disp('Ridge Rsquared = ');
disp(Rsquared);
disp('Ridge MSE = ');
disp(MSE);

%Lasso

Yresid = Y - YHatLasso;
SSresid = sum(Yresid.^2);
SStotal = (length(Y) - 1) * var(Y);
Rsquared = 1 - SSresid/SStotal;
MSE = SSresid/length(Y);
disp('Lasso Rsquared = ');
disp(Rsquared);
disp('Lasso MSE = ');
disp(MSE);