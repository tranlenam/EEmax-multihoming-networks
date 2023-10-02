% Reader can change the penalty paremater c in function "runADMM". 
% Some parameters c require more iterations to converge, 
% you can change "iterationNum" in file "runADMM" if it requires.
% The convex problems in this example is solved by CVX with inner solver SEDUMI
clear 
clc
%% setup the system, do not change these paramters 
M = 2;   %number of subsystem
N = 3;   %number of user
pc1 = 5; %dBm
pc = ones(N,M)*10^(pc1/10);%mW
po_dB = 10; %dBm
po = ones(1,N)*10^(po_dB/10);%mW active circuit power
pr_dB = 10; %dBm
pr = ones(1,M)*10^(pr_dB/10);%mW active circuit power
B = [2.4 1.2 ];%bandwidth RAT allocate to multi homing service Mhz
rmin = [2, 1.4, 1.1];%QoS data rate
P_peak = 10^(24/10); % maximum transmit power 24dBm
Pmax = P_peak*ones(1,N); %mW maximum power transmission, same for all user 
%%
%  example random channel (effective channel)
C_M =[0.134898423125817,6.18586310294603;4.59823100370149,106.806386079692;0.320769348148690,44.7400514027452;];

EE_central = 0.2341; %% centralized result obtained by CVX

EE_sequence = runADMM(C_M,Pmax,rmin,B,M,N,pc,po,pr);


plot(1:size(EE_sequence,2),EE_sequence)
hold on
plot(1:size(EE_sequence,2),ones(1,size(EE_sequence,2))*EE_central,'--r')
legend('Centralized','Proposed decentralized solution')
saveas(gcf,'../results/Fig3.png')