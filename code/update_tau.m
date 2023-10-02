%%%update (35)
function  out_tau = update_tau(tau, c, M,N,pc1,TJ,X)
G = [pc1/M, ones(1,N), zeros(1,N)];
out_tau = tau + c*(G*TJ -X'- 1/M*ones(1,M));