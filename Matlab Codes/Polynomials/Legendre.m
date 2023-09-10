%% parameter field
L=1;
dx=0.01;
%% variable field
x=-L:dx:L;
x=x';
A=[x.^0,x.^1,x.^2,x.^3,x.^4];
N=size(x,1);
%% computation field
[Q,R]=qr(A,0);
Q1=Q(N,:);
P=Q*diag(1./Q1);
%% display field
plot(x,P);