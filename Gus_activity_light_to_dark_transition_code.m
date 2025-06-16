clear all;

%If you run this code you will get the figure only for WT

%timespan

    
domain = [0 200];

xmax = 148;
xmin = 120;

color = ["#000000" "#D95319" "#4DBEEE" "#EDB120"];

m = 1; %Color code counter



%parameters values

p_b = 10.0;
m_b = 1;
k_r = 0.232;
d_e = 27.2;
m_c =  1;
p_cl = 1.00;
p_cd = 112;
d_c  = 1.79;
p_p  = 1;
p_pe = 0.332; 
d_p = 4.91;
k_pc = 34.3;
d_pb = 0.313;
p_g = 0.009;
k_g = 0.113;
p_gp = 2.93;
p_ge = 0.465;
p_gb = 10.7;
p_hc = 0.180;
p_f = 450; 
p_fp = 22;
d_f = 0.0009;
k_0 = 10;
p_0 = 25;
m_f = 1;
d_ec = 0.01;

p_th = [0.2 0.04 0.2 1.3];

%WT parameters

m_p = 1;
p_t1 = 0.2;



%initial condition 

tc1 = 0;
tc2 = 0;           
tc3 = 0;               
tc4 = 0;
tc5 = 0;
tc6 = 0;

tc = [tc1 tc2 tc3 tc4 tc5 tc6];

[IVSOL, DVSOL] = ode45(@(t,dp) Model_18_b(t,dp,p_b,m_b,k_r,d_e,m_c,p_cl,p_cd,d_c...
    ,m_p,p_p,p_pe,d_p,k_pc,d_pb,p_g,k_g,p_gp,p_ge,p_gb,p_f,p_fp,d_f,k_0,p_0,p_t1,m_f, ...
    d_ec),domain,tc);  
                               



figure(1)

%subplot (2,3,6);
hold on;
plot(IVSOL, DVSOL(:,6),'R','LineWidth',10,'Color',color(m))
%xlabel('Time(Hour)','Fontsize',14)
%ylabel('F','Fontsize',14)



%xlim([1 24])
%xticks([0 4 8 12 16 20 24])
%yticks([5000 9000 12000])
xlim([112 144.2]);
xticks([112 120 128 136 144])
set(gca,'xticklabel',{'-8','0','8','16','24'})
ylim([1500 15000]);
yticks([0 7500 15000])
set(gca,'yticklabel',{'0','7500','15000'})
hold on;


