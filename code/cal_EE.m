function EE = cal_EE(C_M, TJ, ZJ,pc1,N)
P = TJ(2:N+1,:)/ZJ(1,1);
W = TJ(N+2:end,:)/ZJ(1,1);
D = sum(sum(W.*log2(1+(C_M.*P)./W)));
EE=D/(sum(sum(P))+pc1);