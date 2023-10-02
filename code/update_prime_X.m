%%%this function update variables x. (29)
function X_out = update_prime_X(pc1,c,N,M ,tau,xi,TJ)
G = [pc1/M, ones(1,N), zeros(1,N)];
a = G*TJ - 1/M*ones(1,M) + (tau - xi)/c;
bar_a = sum(a)/M;
X = a-bar_a*M/(M+1); %(29)
X_out = X';