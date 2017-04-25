A = importdata('../hw07/code/m4.dat');
A = A(2:end,1);
A_size = size(A);
side = sqrt(A_size(1));
A = reshape(A, [side,side]);
e = eig(A)