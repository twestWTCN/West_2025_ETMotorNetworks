function [X,Y] = switchVariables(X,Y)
Xtmp = X; Ytmp = Y;
X = Ytmp;
Y = Xtmp;