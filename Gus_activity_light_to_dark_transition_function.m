%Plant Biology

% p_self included here only
function [diff] = model(t,dp,p_b,m_b,k_r,d_e,m_c,p_cl,p_cd,d_c...
    ,m_p,p_p,p_pe,d_p,k_pc,d_pb,p_g,k_g,p_gp,p_ge,p_gb,p_f,p_fp,d_f,k_0,p_0,p_t1,m_f, ...
    d_ec) 


B = dp(1);     
E = dp(2);
C = dp(3);
P = dp(4);
G = dp(5);
F = dp(6);

%Defining equations

m_e = 1;
p_e1 = 108;
p_e2 = 39.8;



if t <= 120

    L = 1;
    D = 24;
    t0 = mod(t,24);
    t1 = t0 - D;
    t2 = t0 - 24;
    k0 = 5;
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


if P >= p_t1
    
    p_s2 = m_p*(p_p/(1+p_pe*E+p_0*P));
    
    p_s1 = k_0 + (p_f)/(1+p_fp*P) ;   
end

if P < p_t1
    p_s2 = m_p*(p_p/(1+p_pe*E));  %%%%% activator
    
    p_s1 = m_f*k_0 ; %%%% activator
end

end

if t > 120

    L = 0;
    D = 0;
    t0 = mod(t,24);
    t1 = t0 - D;
    t2 = t0 - 24;
    k0 = 5;

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


if P >= p_t1
    
    p_s2 = m_p*(p_p/(1+p_pe*E+p_0*P));
    
    p_s1 = k_0 + (p_f)/(1+p_fp*P) ;   
end

if P < p_t1
    p_s2 = m_p*(p_p/(1+p_pe*E));  %%%%% activator
    
    p_s1 = m_f*k_0 ; %%%% activator
end

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


if P >= p_t1
    
    p_s2 = m_p*(p_p/(1+p_pe*E+p_0*P));
    
    p_s1 = k_0 + (p_f)/(1+p_fp*P) ;   
end

if P < p_t1
    p_s2 = m_p*(p_p/(1+p_pe*E));  %%%%% activator
    
    p_s1 = m_f*k_0 ; %%%% activator
end
 



   
diff = [p_b*L*(m_b - B) - k_r*B 
       p_e - d_ec*E*C - d_e*E
       m_c*(p_cl*L + p_cd*(1-L))- d_c*C
       p_s2 - (d_p/(1+k_pc*C))*P - d_pb*B*P
       p_g + k_g *((p_gp*P)/(1+p_ge*E+p_gb*B))
       p_s1 - d_f*F];

  
end