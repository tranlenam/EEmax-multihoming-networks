function [EE_store] = runADMM(C_M,Pmax,rmin,B,M,N,pc,po,pr)
c = 70; %70 %120 %50 %250  %%penalty parameter of ADMM
%% initial value
Va_in = 0.05; 
u_in = Va_in*kron(ones(1,M),zeros(N,1));
v_in = Va_in*kron(ones(N,1),B/N);
pc1 = sum(sum(pc))+sum(po)+sum(pr);
tau = zeros(1,M);
zetai = zeros(N,1+2*M);
zetaj = zeros(1+2*N,M);
xi = 0;
ZJ = [Va_in*ones(1,M); u_in ;v_in];
X = zeros(M,1);
%%
IterationNum = 50;
EE_store = zeros(1,IterationNum); %%contain EE of each iteration
for iterationIndex=1:IterationNum
    [TJ,TI,~] = update_prime_T(C_M,pc1,c,Pmax,B,rmin,N,M ,tau,zetai,zetaj,ZJ,X);
    X = update_prime_X(pc1,c,N,M ,tau,xi,TJ);
    [ZJ,ZI] = update_prime_Z(c,N,M ,zetai,zetaj,TJ,TI);
    xi = updateXi(M,N,TJ,xi,tau,c,pc1);
    tau = update_tau(tau, c, M,N,pc1,TJ,X);
    zetaj = update_zetaj(zetaj, c,TJ,ZJ);
    zetai = update_zetai(zetai,TI,ZI,c);
    EE = cal_EE(C_M, TJ, ZJ,pc1,N); %calculate EE each iteration. Can use "out_EE" alternatively.
    EE_store(iterationIndex+1)= EE;
end
