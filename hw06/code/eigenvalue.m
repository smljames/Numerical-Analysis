A = importdata('../exam/net.txt');
A(1,1)=1;
%A = A(2:end-14);
%A = reshape(A, [14, 14]);
%A = A';
%A = A*A;
e = eig(A);

max(e)