%%% update(34)
function out_xi = updateXi(M,N,TJ,xi,tau,c,pc1)
G = [pc1/M, ones(1,N), zeros(1,N)];
a = G*TJ - 1/M*ones(1,M) + (tau - xi)/c;
bar_a = sum(a)/M;
out_xi = xi +c*M*bar_a/(M+1); %%(34)