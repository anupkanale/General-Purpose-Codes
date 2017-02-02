% Matlab program to calculate the coefficients for the least squares fit
% for a line and plot

close all; clear all;

%Take input values of x and y
x = [1,2,3.5,4.1,5];
y = [6.5,7,8,8,10];
N = length(x);

%Calculate the coefficients a and b
sumX = sum(x);
sumX2 = sum(x.^2);
sumY = sum(y);
sumXY = sum(x.*y);

a = (N*sumXY - sumX*sumY)/(N*sumX2 - sumX^2);
b = (sumY - sumX*a)/N;
yFit = a*x+b; % Fitted y

% Plots
plot(x,yFit,'-r', x,y, 'ob', 'LineWidth', 2);