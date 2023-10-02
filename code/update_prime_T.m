%%% This function solve  (21) and (22). Author Quang-Doanh Vu
% Because the variables are decoupable, we solve all problems of user and
% BS at the same time without loss of distributed property. The reader can
% modify this fuction such as each user and BS solve one problem as written
% in the paper. Reader can also solve this problem by YALMIP with solver
% SNOPT or fmincom (use function "kullbackleibler" instead of "rel_entr").

function [TJ,TI,out_EE] = update_prime_T(C_M,pc1,c,Pmax,B,rmin,N,M ,tau,zetai,zetaj,ZJ,X)
chan = reshape(C_M,[],1);
G = [eye(M)*pc1/M, kron(eye(M),ones(1,N)), kron(eye(M),zeros(1,N))];
Aj = [zetaj(1,:)'; reshape(zetaj(2:N+1,:),[],1);  reshape(zetaj(N+2:end,:),[],1)];
Ai = reshape(zetai,[],1);
AA = kron(ones(1,M),eye(N));
BB = kron(eye(M),ones(1,N));
VA = ZJ(1,1); %vartheta
U = reshape(ZJ(2:N+1,:),[],1);
V = reshape(ZJ(N+2:end,:),[],1);

cvx_solver sedumi
cvx_begin quiet
variables P(N*M,1)  P1(N*M,1) W(N*M,1) W1(N*M,1) VJ(M,1) VI(N,1)   %power, bandwidth, extra. P,W,VJ for BS ; P1 W1 and VI for users
minimize (sum(rel_entr(W,chan.*P+W)/log(2))  + sum(tau*(G*[VJ;P;W]-X-ones(M,1)/M) + c/2*(G*[VJ;P;W]-X-ones(M,1)/M).^2 )+ ...
    Aj'*([VJ;P;W] - [VA*ones(M,1);U;V]) + c/2*sum(([VJ;P;W] - [VA*ones(M,1);U;V]).^2)...
    + Ai'*([VI;P1;W1] - [VA*ones(N,1);U;V])+ c/2*sum(([VI;P1;W1] - [VA*ones(N,1);U;V]).^2)  )
subject to
BB*W - VJ.*B' <= 0
AA*P1 - VI.*Pmax'<=0
AA*(-rel_entr(W1, chan.*P1+W1)/log(2)) - VI.*rmin' >=0  %data rate constraints at user
P>=0
P1>=0
W>=0
W1>=0
VI>=0
VJ>=0
cvx_end
TJ = [VJ'; reshape(P,N,M); reshape(W,N,M)]; %BS output
TI = [VI,reshape(P1,N,M),reshape(W1,N,M) ]; %user output
out_EE = sum( W.*log2(1+(chan.*P)./W));