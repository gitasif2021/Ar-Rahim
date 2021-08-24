% Call this file from capacity_vs_iab_nodes_2.m file
d=40;
d_m=100;
Pm=1;
N_dB=-74;
N=10^(-3)*10^(-7.4);
f=28*10^9;
W=10^9 ;
n=2;
c=3*10^8;
lamda=c/f;
d_ref=1.5;
k0=(lamda/4*pi*d_ref)^2;

count_node=0;
x_node=[];
y_node=[];
x_typical=0;
y_typical=0;
r_u=unifrnd(0,d);
th_u=unifrnd(0,2*pi);
x_typ_u=r_u*cos(th_u);
y_typ_u=r_u*sin(th_u);
r_d=unifrnd(2*d,d_m);
th_d=unifrnd(0,2*pi);
x_d=r_d*cos(th_d);
y_d=r_d*sin(th_d);
r_d_u=unifrnd(0,d);
th_d_u=unifrnd(0,2*pi);
x_d_u=r_u*cos(th_u);
y_d_u=r_u*sin(th_u);

while 1
    r=unifrnd(2*d,d_m);
    th=unifrnd(0,2*pi);
    x=r*cos(th);
    y=r*sin(th);
    r1=sqrt((x_d-x)^2+(y_d-y)^2);
    if r1>=80 && r1<=100
        x_node=[x_node x];
        y_node=[y_node y];
        count_node=count_node+1;
        if count_node==(n_it-1)
            break;
        end
    end
    
end

r_nu=unifrnd(0,d,1,n_it-1);
t_nu=unifrnd(0,2*pi,1,n_it-1);
x_nu=x_node+r_nu.*cos(t_nu);
y_nu=y_node+r_nu.*sin(t_nu);

t1=x_node-x_typical;
t2=y_node-y_typical;
t1=t1.^2;
t2=t2.^2;
t3=t1+t2;
t3=sqrt(t3);
d_d_n1=sqrt((x_typical-x_d)^2+(y_typical-y_d)^2);
t3=[d_d_n1 t3];
h_sq_npc_l1=exprnd(1,1,n_it);
t3=t3.^-n;
t3=k0*t3;
g_npc_l1=h_sq_npc_l1.*t3;

des_npc_l1=Pm*g_npc_l1(1);
sum_temp=0;
for v=2:n_it
    sum_temp=sum_temp+g_npc_l1(v);
end
int_npc_l1=Pm*sum_temp+N;
sinr_npc_l1=des_npc_l1/int_npc_l1;

t1=x_node-x_typ_u;
t2=y_node-y_typ_u;
t1=t1.^2;
t2=t2.^2;
t3=t1+t2;
t3=sqrt(t3);
d_d_typu=sqrt((x_d-x_typ_u)^2+(y_d-y_typ_u)^2);
t3=[d_d_typu t3];
d_typn_u=sqrt((x_typical-x_typ_u)^2+(y_typical-y_typ_u)^2);
t3=[d_typn_u t3];
h_sq_npc_l2=exprnd(1,1,n_it+1);
t3=t3.^-n;
t3=k0*t3;
g_npc_l2=h_sq_npc_l2.*t3;

des_npc_l2=Pm*g_npc_l2(1);
sum_temp=0;
for v=2:(n_it+1)
    sum_temp=sum_temp+g_npc_l2(v);
end
int_npc_l2=Pm*sum_temp+N;
sinr_npc_l2=des_npc_l2/int_npc_l2;

sinr_npc=(1/sinr_npc_l1 + 1/sinr_npc_l2)^-1;

t1=x_node-x_d_u;
t2=y_node-y_d_u;
t1=t1.^2;
t2=t2.^2;
t3=t1+t2;
t3=sqrt(t3);
d_n1_du=sqrt((x_typical-x_d_u)^2+(y_typical-y_d_u)^2);
t3=[d_n1_du t3];
d_d_u=sqrt((x_d-x_d_u)^2+(y_d-y_d_u)^2);
t3=[d_d_u t3];
h_sq_temp1=exprnd(1,1,n_it+1);
t3=t3.^-n;
t3=k0*t3;
g_temp1=h_sq_npc_l2.*t3;
des_temp1=Pm*g_temp1(1);
sum_temp=0;
for v=2:(n_it+1)
    sum_temp=sum_temp+g_temp1(v);
end
int_temp1=Pm*sum_temp+N;
sinr_temp1_donor=des_temp1/int_temp1;

sinr_th=1;
P_donor_dbm=30;
P_donor=10^(-3)*10^(P_donor_dbm/10);

if (sinr_temp1_donor<=sinr_th)
    P_donor=Pm;
else
    
    while 1
        
        P_donor_dbm=P_donor_dbm-1;
        P_donor=10^(-3)*10^(P_donor_dbm/10);
        des_temp1=P_donor*g_temp1(1);
        int_temp1=Pm*sum_temp+N;
        sinr_temp1_donor=des_temp1/int_temp1;
        
        if(sinr_temp1_donor<= sinr_th)
            P_donor_dbm=P_donor_dbm+1;
            P_donor=10^(-3)*10^(P_donor_dbm/10);
            break;
        end
    end
end

P_controlled_n=[];
for s=1:n_it-1
    d_temp_ar1=[];
    temp1=sqrt((x_node(s)-x_nu(s))^2+(y_node(s)-y_nu(s))^2);
    temp2=sqrt((x_d-x_nu(s))^2+(y_d-y_nu(s))^2);
    temp3=sqrt((x_typical-x_nu(s))^2+(y_typical-y_nu(s))^2);
    d_temp_ar1=[d_temp_ar1 temp1 temp2 temp3];
    for s1=1:n_it-1
        if s1~=s
            temp1=sqrt((x_node(s1)-x_nu(s))^2+(y_node(s1)-y_nu(s))^2);
            d_temp_ar1=[d_temp_ar1 temp1];
        end
    end
    h_sq_temp_ar1=exprnd(1,1,n_it+1);
    d_temp_ar1=d_temp_ar1.^-n;
    d_temp_ar1=k0*d_temp_ar1;
    g_temp_ar1=h_sq_temp_ar1.*d_temp_ar1;
    des_temp_ar1=Pm*g_temp_ar1(1);
    sum_temp=0;
    for v=2:(n_it+1)
        sum_temp=sum_temp+g_temp_ar1(v);
    end
    int_temp_ar1=Pm*sum_temp+N;
    sinr_temp_ar1=des_temp_ar1/int_temp_ar1;
    
    sinr_th=1;
    P_node_dbm=30;
    P_node=10^(-3)*10^(P_node_dbm/10);
    
    if (sinr_temp_ar1<=sinr_th)
        P_node=Pm;
    else
        
        while 1
            
            P_node_dbm=P_node_dbm-1;
            P_node=10^(-3)*10^(P_node_dbm/10);
            des_temp_ar1=P_node*g_temp_ar1(1);
            int_temp_ar1=Pm*sum_temp+N;
            sinr_temp_ar1=des_temp_ar1/int_temp_ar1;
            
            if(sinr_temp_ar1<= sinr_th)
                P_node_dbm=P_node_dbm+1;
                P_node=10^(-3)*10^(P_node_dbm/10);
                break;
            end
        end
    end
    P_controlled_n=[P_controlled_n P_node];
end

des_pc_l1=Pm*g_npc_l1(1);
sum_temp=0;
for it1=1:(n_it-1)
    sum_temp=sum_temp+g_npc_l1(it1+1)*P_controlled_n(it1);
end
int_pc_l1=sum_temp+N;
sinr_pc_l1=des_pc_l1/int_pc_l1;

des_pc_l2=Pm*g_npc_l2(1);
sum_temp=0;
for it2=1:(n_it-1)
    sum_temp=sum_temp+g_npc_l2(it2+2)*P_controlled_n(it2);
end
sum_temp=sum_temp+P_donor*g_npc_l2(2);
int_pc_l2=sum_temp+N;
sinr_pc_l2=des_pc_l2/int_pc_l2;

sinr_pc=(1/sinr_pc_l1 + 1/sinr_pc_l2)^-1;
sinr_npc_dB=10*log10(sinr_npc);
sinr_pc_dB=10*log10(sinr_pc);
