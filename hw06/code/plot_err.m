clc
close all
err1 = importdata('../output/err1.txt');
err2 = importdata('../output/err2.txt');
err3 = importdata('../output/err3.txt');
err4 = importdata('../output/err4.txt');
lambda1 = importdata('../output/lambda1.txt');
lambda2 = importdata('../output/lambda2.txt');
lambda3 = importdata('../output/lambda3.txt');
lambda4 = importdata('../output/lambda4.txt');

iter = (0:length(err1)-1);

figure
fig = loglog(iter, err1);
title('err1 v.s. iteration');
%saveas(fig, 'err1_iter.png');

figure
fig = loglog(iter, err2);
title('err2 v.s. iteration');
%saveas(fig, 'err2_iter.png');

figure
fig = loglog(iter, err3);
title('err3 v.s. iteration');
%saveas(fig, 'err3_iter.png');

figure
fig = loglog(iter, err4);
title('err4 v.s. iteration');
%saveas(fig, 'err4_iter.png');

figure
fig = plot(iter, lambda1);
title('lambda v.s. iteration');
%saveas(fig, 'lambda_iter.png');

%{
figure
plot(iter, lambda2);

figure
plot(iter, lambda3);

figure
plot(iter, lambda4);
%}