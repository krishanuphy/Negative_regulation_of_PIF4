clear all;

%timespan


D = 8;  %day length
   
domain = [0 48];


color = ['g' 'b' 'k' 'red'];


%parameters values;  already explained in main code(hypocotyl fitting code)
p_b = 10.0;
k_r = 0.232;
p_cl = 1.00;
p_pe = 0.332; 
p_e1 = 108;
p_e2 = 39.8;
p_t1 = 0.2;
m_b = 1;
d_e = 27.2;
m_c =  1;
p_cd = 112;
d_c  = 1.79;
p_p  = 1;
d_p = 4.91;
k_pc = 34.3;
d_pb = 0.313; 
p_g = 0.009;
k_g = 0.113;
p_gp = 2.93;
p_ge = 0.465;
p_gb = 10.7;
p_fp = 22;
d_f = 0.3;
d_ec = 0.01;
p_f = 450;
k_0 = 10;
%Negative Feedback parameters

p_max = 70; %MAX feedback
p_min = 0.1; %Min feedback
p_max_end = 80;
D_1 = 2;
D_2 = 6;

%Gus coefficient parameters
p_max1 = 30;    %
p_min1 = 0.5;    %Min
p_max1_end = 10;
D_3 = 2;
D_4 = 6;
m_p = 1; 

m = 3;    %Color change plot 





%initial condition 

tc1 = 0;
tc2 = 0;           
tc3 = 0;               
tc4 = 0;
tc5 = 0;
tc6 = 0;

tc = [tc1 tc2 tc3 tc4 tc5 tc6];

[IVSOL, DVSOL] = ode45(@(t,dp) model(t,dp,p_b,m_b,k_r,d_e,m_c,p_cl,p_cd,d_c...
    ,m_p,p_p,p_pe,d_p,k_pc,d_pb,p_g,k_g,p_gp,p_ge,p_gb,p_f,d_f,k_0,D,p_t1, ...
    d_ec,p_max,p_min,D_1,p_e1,p_e2,D_2,D_3,D_4,p_max1,p_min1,p_max_end,p_max1_end),domain,tc);  
                               

figure(1)




   

hold on;
plot(IVSOL-24, DVSOL(:,4)/0.629,'R','LineWidth',3,'Color',color(m))   %0.629 isthe normalization factor
%plot(IVSOL-24, DVSOL(:,axis),'R','LineWidth',2,'Color',color(m))
xlabel('Time(Hour)','Fontsize',14)
ylabel('PIF4 Protein','Fontsize',14)
%set(gca,'Xscale','log')
xlim([0 24])
xticks([0 8 16 24]);
%xlim([xmin xmax]);
hold off;





%27 degree pif4 code,parameter is going to change




p_b = 5;
k_r = 0.251;  
p_cl = 2.37; 
p_pe = 0.028; 
p_e1 = 127;  
p_e2 = 7.29;

p_t1 = 0.2;


m_b = 1;
d_e = 27.2;
m_c =  1;
p_cd = 112;
d_c  = 1.79;
p_p  = 1;
d_p = 4.91;
k_pc = 34.3;
d_pb = 0.313; 
p_g = 0.009;
k_g = 0.113;
p_gp = 2.93;
p_ge = 0.465;
p_gb = 10.7;
p_fp = 22;
d_f = 0.3;

p_f = 450; 
k_0 =  100; 

m_p = 1;
d_ec = 0.01;
%Negative Feedback parameters

p_max = 20; 
p_min = 0.1; 
p_max_end = 20;
D_1 = 2;
D_2 = 6;

%Gus coefficient parameters
p_max1 = 5;
p_min1 = 0;    %Min
p_max1_end = 5;
D_3 = 2;
D_4 = 6;

m = 4;    %Color change plot 


%initial condition 

tc1 = 0;
tc2 = 0;           
tc3 = 0;               
tc4 = 0;
tc5 = 0;
tc6 = 0;

tc = [tc1 tc2 tc3 tc4 tc5 tc6];

[IVSOL, DVSOL] = ode45(@(t,dp) model(t,dp,p_b,m_b,k_r,d_e,m_c,p_cl,p_cd,d_c...
    ,m_p,p_p,p_pe,d_p,k_pc,d_pb,p_g,k_g,p_gp,p_ge,p_gb,p_f,d_f,k_0,D,p_t1, ...
    d_ec,p_max,p_min,D_1,p_e1,p_e2,D_2,D_3,D_4,p_max1,p_min1,p_max_end,p_max1_end),domain,tc);  
                               

figure(1)




   
hold on;
plot(IVSOL-24, DVSOL(:,4)/0.629,'R','LineWidth',3,'Color',[0.92, 0, 0])   %0.629 isthe normalization factor
xlabel('Time(Hour)','Fontsize',14)
ylabel('PIF4 Protein','Fontsize',14)
%set(gca,'Xscale','log')
xlim([0 24])
xticks([0 8 16 24]);
xlabel('Time(hour)', 'FontSize', 20); % X-axis label with font size 14
ylabel('PIF4 Protein', 'FontSize', 20); % Y-axis label with font size 14
%xticks([0 4 8 12 16 20 24])


pbaspect([1 1.5 1]);


% main equatins function

function [diff] = model(t,dp,p_b,m_b,k_r,d_e,m_c,p_cl,p_cd,d_c...
    ,m_p,p_p,p_pe,d_p,k_pc,d_pb,p_g,k_g,p_gp,p_ge,p_gb,p_f,d_f,k_0,D,p_t1, ...
    d_ec,p_max,p_min,D_1,p_e1,p_e2,D_2,D_3,D_4,p_max1,p_min1,p_max_end,p_max1_end) 


B = dp(1);     
E = dp(2);
C = dp(3);
P = dp(4);
G = dp(5);
F = dp(6);


%write the nested function

          

    function p_s2 = calc_ps2(P,E,p_0)

        p_s2 = m_p*(p_p/(1+p_pe*E+p_0*P));

    end

    function p_s1 = calc_ps1(P,p_1)

        p_s1 = k_0 + (p_f)/(1+p_1*P) ;

    end
%Defining equations

m_e = 1;
t0 = mod(t,24);
t1 = t0 - D;
t2 = t0 - 24;
k0 = 5;


if t0 >= D
    L = 0;
end

if t0 < D
    L = 1;
end

if 0 < D < 24

p_e = m_e*p_e1 - p_e2*(-1+(2/(1+exp(-k0*t0)))-(2/(1+exp(-k0*t1)))+...
     (2/(1+exp(-k0*t2))));
end

if D == 0

    p_e = m_e*p_e1 + p_e2;

end

if D == 24

    p_e = m_e*p_e1 - p_e2;

end

%PIF4 protein first half
if t0 <= D_1  
    p_0 = -((p_max-p_min)/D_1)*t0+p_max ;
    %p_s2 = calc_ps2(P,E,p_0); 
end
%PIF4 second part 
if t0 > D_1 & t0 <= D_2
    p_0 = p_min;
    %p_s2 = calc_ps2(P,E,p_0); 
end
%PIF4 third part
if  t0 > D_2 
    t1 = t0-D_2;
    %t0 = t0-D_2;
    p_0 = ((p_max_end-p_min)/(24-D_2))*t1+p_min ;
    %p_s2 = calc_ps2(P,E,p_0);  
end

%GUS first half
if  t0 <= D_3 
    p_1 = -((p_max1-p_min1)/D_3)*t0+p_max ;
    %p_s1 = calc_ps1(P,p_1); 
end
%GUS second half
if t0 > D_3 & t0 <= D_4
    p_1 = p_min1;
    %p_s1 = calc_ps1(P,p_1);
end
%GUS third part
if t0 > D_4 
    %t0 = t0-D_4;
    t2 = t0-D_4;
    p_1 = ((p_max1_end-p_min1)/(24-D_4))*t2+p_min ;
    %p_s1 = calc_ps1(P,p_1);
end

if P > p_t1
    p_s2 = calc_ps2(P,E,p_0);
    p_s1 = calc_ps1(P,p_1);
end
if P <= p_t1
    p_s2 = m_p*(p_p/(1+p_pe*E));  %%%%% activator
    p_s1 = k_0 ; %%%% activator
     % p_s2 = m_p*(p_p/(1+p_pe*E));
     % p_s1 = (p_f) ; 
end
 



   
diff = [p_b*L*(m_b - B) - k_r*B
       p_e - d_ec*E*C - d_e*E
       m_c*(p_cl*L + p_cd*(1-L))- d_c*C
       p_s2 - (d_p/(1+k_pc*C))*P - d_pb*B*P
       p_g + k_g *((p_gp*P)/(1+p_ge*E+p_gb*B))
       p_s1 - d_f*F];

  
end