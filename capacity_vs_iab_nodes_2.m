clc;
close all;
clear all;
NOI2=1000;
d=40;
d_m=100;
Pm=1;
N_dB=-74;
N=10^(-3)*10^(-7.4);
f=28*10^9;
W=10^9 ;
n=4;
c=3*10^8;
lamda=c/f;
d_ref=1.5;
k0=(lamda/4*pi*d_ref)^2;

node_dens=[];
cap_npc_arr=[];
cap_pc_arr=[];

for n_it=2:6
    cap_pc=0;
    cap_npc=0;
    for j=1:NOI2
        capacity_vs_iab_nodes_1;
        cap_npc=cap_npc + log2(1+sinr_npc);
        cap_pc=cap_pc + log2(1+sinr_pc);
    end
    cap_npc_arr=[cap_npc_arr  cap_npc/NOI2];
    cap_pc_arr=[cap_pc_arr  cap_pc/NOI2];
    node_dens=[node_dens n_it];
end

plot(node_dens,cap_npc_arr,'*--');
hold on;
plot(node_dens,cap_pc_arr,'<--');
xlabel('Node Density (/km^2)');
ylabel('Capacity (bps/Hz)');
legend("With out power control","With power control")
grid on;
