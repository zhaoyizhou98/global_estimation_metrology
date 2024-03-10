function X = mysylvester(A,B,C)
%         Solve the equation A*X + X*B = C 
A = (A + A') / 2;
B = (B + B') / 2;
[QA, dA] = eig(A, 'vector');
[QB, dB] = eig(B, 'vector');

CC = QA'*C*QB;
tmp = dA + dB';
idx = find(abs(tmp) <= 0.00001);
tmp(idx) = inf;
X = CC ./ tmp;
X = QA*X*QB';

end