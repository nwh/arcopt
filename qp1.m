function [f g] = qp1(A,b,c,x)
Ax = A*x;
f = 0.5*x'*Ax + b'*x + c;
g = Ax + b;
