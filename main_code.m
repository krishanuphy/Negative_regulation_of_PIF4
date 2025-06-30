clear all;
a1 = [];
i = 1;
d1 = [];


    
domain = [0 120];

xmax = 120;
xmin = 0;

color = ["#000000" "#FF0000" "#0072BD" "#EDB120"];

%parameter Values

p_b = 10.0;       %Light-induced activation rate of phyB
m_b = 1;
k_r = 0.232;      %Deactivation rate of phyB in dark
d_e = 27.2;       %Decay rate of ELF3
m_c =  1;
p_cl = 1.00;      %Production rate of COP1 in light
p_cd = 112;       %Production rate of COP1 in dark
d_c  = 1.79;      %Decay rate of COP1
p_p  = 1;         %Production rate of PIF4
p_pe = 0.332;     %Intensity of ELF3’s inhibition of PIF4 production
d_p = 4.91;       %Decay rate of PIF4
k_pc = 34.3;      %Intensity of COP1’s inhibition of PIF4 degradation
d_pb = 0.313;     %Inhibition rate of PIF4 by phyB
p_g = 0.009;      %Basal rate of hypocotyl growth
k_g = 0.113;      %Conversion factor between PIF4-targeted gene expression and hypocotyl growth
p_gp = 2.93;      %
p_ge = 0.465;     %Intensity of ELF3’s inhibition of hypocotyl growth
p_gb = 10.7;      %Intensity of phyB’s inhibition of hypocotyl growth
p_f = 450;        %Production rate of GUS 
p_fp = 22;        %Intensity of PIF4’s inhibition of GUS expression
d_f = 0.0009;     %Decay rate of GUS
k_0 = 10;         %Basal rate of GUS production
p_0 = 25;         %Intensity of autoinhibition for PIF4 synthesis
k_01 = k_0;       
d_ec = 0.01;      %Inhibition rate of ELF3 by COP1
 
p_th = [0.2 0.04 0.2 1.3];       %Threshold Value for four genotype









%Experimental Data

exp_D = [4 8 12 16];                         %Day legth
exp_hyp = [7.3538 4.3948 2.68285 2.0395];    %Hypocotyle length for wildtype(WT)
exp_dev = [0.668411125 0.406120101 0.277211163 0.190259795]; %standard devation of the experimental measurement of hypocotyl for WT 

exp_hyp_mut = [5.532076923 2.7465 2.177263158 1.208666667];  %Hypocotyl legth of pif4-101
exp_hyp_mut_dev  = [0.74006649 0.361545546 0.372886912 0.184124628]; %standard devation of the experimental measurement of hypocotyl for pif4-101 

exp_hyp_oe1 = [12.2035 7.22295 5.8414 3.5197];  %Hypocotyl legth of PIF4-OE1
exp_hyp_oe1_dev = [0.365437409 0.785811647 0.455545179 0.302029469];%standard devation of the experimental measurement of hypocotyl for PIF4-OE1  

exp_hyp_oe2 = [13.16542857 11.64195 9.8363 6.335]; %Hypocotyl legth of PIF4-OE2
exp_hyp_oe2_dev = [1.291944162 0.564585745 0.522554367 0.544857312];  %standard devation of the experimental measurement of hypocotyl for PIF4-OE2

m=0;
for m_p = [1 0.65 4 5.9]

    p_t1 = p_th(i)
    i = i+1;
    j=0;
for D = 1:24

j = j+1;
%initial condition 

tc1 = 0;   %phyB initial condition
tc2 = 0;   %ELF3 initial condition        
tc3 = 0;   %COP1 initial condition        
tc4 = 0;   %PIF4 initial condition
tc5 = 0;   %Hypocotyl growth initial condition
tc6 = 0;   %GUS activity initial condition

tc = [tc1 tc2 tc3 tc4 tc5 tc6];

[IVSOL, DVSOL] = ode45(@(t,dp) hyp_function(t,dp,p_b,m_b,k_r,d_e,m_c,p_cl,p_cd,d_c...
    ,m_p,p_p,p_pe,d_p,k_pc,d_pb,p_g,k_g,p_gp,p_ge,p_gb,p_f,p_fp,d_f,k_0,p_0,D,p_t1,k_01,...
    d_ec),domain,tc);    
                               




[Row,Col] = size(DVSOL);

g_end = DVSOL(Row,5);

g_end_data(j) = g_end;
day_length(j) = D;



end
m = m+1;
plot(day_length,g_end_data,'-','Color',color(m),'LineWidth',5)
hold on;


end


e = errorbar(exp_D,exp_hyp,exp_dev,"o",'LineWidth',1);
e.Marker = 'o';
e.MarkerSize = 5;
e.Color = "#000000";
e.CapSize = 10;

e = errorbar(exp_D,exp_hyp_mut,exp_hyp_mut_dev,"o",'LineWidth',1);
e.Marker = 'o';
e.MarkerSize = 5;
e.Color = "#FF0000";
e.CapSize = 10;

e = errorbar(exp_D,exp_hyp_oe1,exp_hyp_oe1_dev,"o",'LineWidth',1);
e.Marker = 'o';
e.MarkerSize = 5;
e.Color = "#0072BD";
e.CapSize = 10;

e = errorbar(exp_D,exp_hyp_oe2,exp_hyp_oe2_dev,"o",'LineWidth',1);
e.Marker = 'o';
e.MarkerSize = 5;
e.Color = "#EDB120";
e.CapSize = 10;

xlabel('Day Length(Hour)','Fontsize',14)
ylabel('Growth','Fontsize',14)
xlim([1 24])

txt = ['Hypocotyl Length vs Day Length Graph,Model Code No - 18a, Mutp - 0.65,1,4,5.9 and' ...
    '  and P_s_e_l_f = 25'];
sgtitle(txt);





