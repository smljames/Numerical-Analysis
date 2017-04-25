y3 = importdata('../hw08/output/f3.txt'); % yl: using Lagrange to calculate y
y5 = importdata('../hw08/output/f5.txt');
y7 = importdata('../hw08/output/f7.txt');
y13 = importdata('../hw08/output/f13.txt');
y21 = importdata('../hw08/output/f21.txt');

x = 475:775;
y = importdata('../hw08/code/f301.dat');
y = y.data(:, 2);
%{
figure
plot(x, y3);
hold;
%}
%{
plot(x, y5);
hold;
%}
%{
plot(x, y7);
hold;
%}
%{
plot(x, y13);
hold;
%}
%
plot(x, y21);
hold;
%}
plot(x, y);
xlabel('x');
ylabel('y');
title('Lagrange Interpolation');