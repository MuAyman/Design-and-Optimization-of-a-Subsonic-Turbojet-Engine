clc;
%% Initial Values and Spacs
AR = 8; %aspect ratio
b =  280; %wing area in m2
S = b^2/AR; %wing span
gammac = 1.4;   gammat = 1.3;
cpc = 1004.5;   cpt = 1148.86;
Rc = 287;       Rt = cpt*(gammat-1)/gammat;
pi_d = 0.96;    pi_b = 0.96;    pi_n = 0.96;
eta_c = 0.88;   eta_b = 0.99;   eta_t = 0.9;    eta_m = 0.98;
hpr = 4.28e7;
W_max = 141700;
W_empty = 64600;
W_payload = 21150;
W_f = W_max - W_empty - W_payload;
W_endClimb = W_max - 0.1*W_f;   % 10% fuel comsumed
W_endCruise = W_max - 0.9*W_f;   % 90% fuel comsumed
% take-off(start mid), climb(start mid), & cruise(start mid end)
W = [W_max 140301.25 138902.5 137503.75 W_endClimb 113725 W_endCruise]*9.81; % N
H = [0 3 8.7]*10^3; % m
RC = [6.35 10]; % rate of climb
V = [241 216];

%% Evaluating M4
To3 = 288 - 6.5e-3*H(1);
To5 = 288 - 6.5e-3*(2*H(2));
M3 = V(1)/sqrt(gammac*Rc*To3);
M5 = V(2)/sqrt(gammac*Rc*To5);
M4 = 0.5*(M3+M5); % M mid climb 
M = [0.286  M4 0.89];

%% Thrust Required
theta = asin(RC./V); % rad
theta_avg = 0.5*(theta(1)+theta(2)); 
To = 288 - 6.5e-3*H;    % K
Po = (101.325*(1-(2.257e-5)*H).^5.256)*1000; % Pa
pho = Po./(0.2869*To);  % kg/m3
V = M.*sqrt(gammac*Rc*To); % update the V vector with the mid points speeds

% cruise
a3 = 0.5*pho(3)*V(3)^2*b;
CL(3) = W(6)/a3;
CD(3) = (0.01 + 0.065*CL(3)^2);
F_r(3) = 1/4*(a3*CD(3)); % thrust required for one engine 

%% Finding the min Tt4 & Pi_c
% Assumptions
Tt4_3 = 660:0.2:1100;

for o = 1:length(Tt4_3)

% cruise @ 6
% ram effect
tau_r_3 = 1 + (gammac - 1)/2 * M(3)^2;
pi_r_3 = tau_r_3^(gammac/(gammac-1));
Tt2_3 = tau_r_3 * To(3);
Pt2_3 = pi_r_3 * Po(3);

% comb.
pi_c_3 = 8:0.01:30;
tau_c_3 = (pi_c_3.^((gammac-1)/gammac)-1)./eta_c + 1;

% c.c.
tau_lambda_3 = cpt*Tt4_3(1,o)/(cpc*To(3));
f_3 = (tau_lambda_3-tau_r_3.*tau_c_3)./(((eta_b*hpr)/(cpc*To(3)))-tau_lambda_3);

% turb.
tau_t_3 = 1 - tau_r_3.*(tau_c_3-1)./((eta_m.*(1+f_3))*tau_lambda_3);
pi_t_3 = (1-(1-tau_t_3)./eta_t).^(gammat/(gammat-1));

% nozzle
pt9p9_cr_3 = (0.5*gammat+0.5)^(gammat/(gammat-1));
pt9p0_3 = pi_r_3.*pi_d.*pi_c_3.*pi_b.*pi_t_3.*pi_n;
% Asuuming full expansion
p9p0_3 = pt9p0_3./pt9p9_cr_3;
M9_3 = 1;
T9_3 = Tt4_3(1,o).*tau_t_3/(0.5+0.5*gammat);
V9_3 = sqrt(gammat*Rt*T9_3);
V9e_3 = V9_3 + V9_3./(gammat.*M9_3^2).*(1-1./p9p0_3).*(1+f_3);

% preformance
Vo_3 = M(3) .* sqrt(gammac*Rc*To(3));  
Fmo = (1+f_3).*V9e_3 - Vo_3;
SFC = f_3./Fmo;

plot(pi_c_3,SFC);
xlabel('Pi_c'); ylabel('SFC'); title('Pi_C vs SFC');

%% geting mo3 & A9 
[MinSFC, best] = min(SFC);
pi_cBest = pi_c_3(best);
ct3 = MinSFC;
mo_3 = F_r(3)/Fmo(best);
A9 = mo_3*(1+f_3)/(pho(3)*V(3));
R(3) = 2*sqrt(2/(pho(3)*S))*1/ct3*CL(3)^0.5/CD(3)*(W(5)^0.5 - W(7)^0.5); 

%% off-design procedure
F_a(1) = 0; F_r(1) = 1; F_a(2) = 0; F_r(2) = 1;
mo_ref_req = mo_3;
while ((F_a(1) < F_r(1)) || (F_a(2) < F_r(2)))
    mo_3 = mo_ref_req;
    Tt2_ref = Tt2_3;
    Tt4_ref = Tt4_3(1,o);
    tau_c_ref = tau_c_3(best);
    pi_c_ref = pi_cBest;
    Pt2_ref = Po(1)*pi_r_3;
    mo_ref = mo_3;

    %take off
    % ram effect
    tau_r_off_1 = 1+((gammac-1)/2)*M(1)^2;
    pi_r_off_1 = (1+((gammac-1)/2)*M(1)^2)^(1.4/0.4);
    Tt2_off_1 = tau_r_off_1*To(1);
    Tt4_off_1 = Tt4_3(1,o); % same as Tt4_3
    sqrt_Tt2_rel_1 = sqrt(Tt2_off_1)/sqrt(Tt2_ref);
    sqrt_Tt4_rel_1 = sqrt(Tt4_off_1)/sqrt(Tt4_ref);
    
    syms tau_c_off_1
    eqn2 = ((tau_c_off_1-1)/(tau_c_ref-1))==((Tt4_off_1/Tt2_off_1)/(Tt4_ref/Tt2_ref));
    tau_c_off_1 = double(solve(eqn2,tau_c_off_1));
    pi_c_off_1=tau_c_off_1^(gammac/(gammac-1));
    pi_c_rel_1=pi_c_off_1/pi_c_ref;
    
    syms mo_rel_1
    Pt2_off_1=pi_r_off_1*Po(1);
    Pt2_rel_1=Pt2_off_1/Pt2_ref;
    eqn1=mo_rel_1*sqrt_Tt2_rel_1/Pt2_rel_1==pi_c_rel_1/(sqrt_Tt4_rel_1/sqrt_Tt2_rel_1);
    mo_rel_1=double(solve(eqn1,mo_rel_1));
    mo_off_1=mo_rel_1*mo_ref;
    
    % available thrust
    M9_off_1=1;
    Tt9_off_1=tau_t_3(best)*Tt4_off_1;
    T9_off_1=Tt9_off_1/(1+((gammat-1)/2));
    V9_off_1=sqrt(gammat*Rt*T9_off_1);
    Pt9P9_cr_off = (0.5*gammat+0.5)^(gammat/(gammat-1));
    Pt9P0_off_1=pi_r_off_1*pi_d*pi_c_off_1*pi_b*pi_t_3(best)*pi_n;
    P0P9_off_1=Pt9P9_cr_off/Pt9P0_off_1;
    F_a(1)=mo_off_1*((1+f_3(best))*V9_off_1-V(1)+((1+f_3(best))*(V9_off_1/(gammat*(M9_off_1^2)))*(1-(P0P9_off_1))));
    ct1 = f_3(best)/(F_a(1)/mo_off_1);
    
    % required thrust
    a1 = 0.5*pho(1)*V(1)^2*b;
    CL(1) = W(2)/a1;
    CD(1) = (0.01 + 0.065*CL(1)^2);
    F_r(1) = 1/4*(a1*CD(1)); % thrust required for one engine
    R(1) = 2*sqrt(2/(pho(1)*S))*1/ct1*CL(1)^0.5/CD(1)*(W(1)^0.5 - W(3)^0.5);
    
    if F_a(1) < F_r(1)
        mo_off_r_1 = F_r(1)/F_a(1)*mo_off_1;
        mo_ref_req = mo_off_r_1/mo_rel_1;
    end
    
    % climb
    % ram effect
    tau_r_off_2 = 1+((gammac-1)/2)*M(2)^2;
    pi_r_off_2 = (1+((gammac-1)/2)*M(2)^2)^(1.4/0.4);
    Tt2_off_2 = tau_r_off_2*To(2);
    Tt4_off_2 = Tt4_3(1,o);
    sqrt_Tt2_rel_2 = sqrt(Tt2_off_2)/sqrt(Tt2_ref);
    sqrt_Tt4_rel_2 = sqrt(Tt4_off_2)/sqrt(Tt4_ref);

    syms tau_c_off_2
    eqn2 = ((tau_c_off_2-1)/(tau_c_ref-1))==(Tt4_off_2/Tt2_off_2)/(Tt4_ref/Tt2_ref);
    tau_c_off_2 = double(solve(eqn2,tau_c_off_2));
    pi_c_off_2 = tau_c_off_2^(gammac/(gammac-1));
    pi_c_rel_2 = pi_c_off_2/pi_c_ref;

    syms mo_rel_2
    Pt2_off_2 = pi_r_off_2*Po(2);
    Pt2_rel_2 = Pt2_off_2/Pt2_ref;
    eqn1 = mo_rel_2*sqrt_Tt2_rel_2/Pt2_rel_2==pi_c_rel_2/(sqrt_Tt4_rel_2/sqrt_Tt2_rel_2);
    mo_rel_2 = double(solve(eqn1,mo_rel_2));
    mo_off_2 = mo_rel_2*mo_ref;
    
    % available thrust
    M9_off_2 = 1;
    Tt9_off_2 = tau_t_3(best)*Tt4_off_2;
    T9_off_2 = Tt9_off_2/(1+((gammat-1)/2));
    V9_off_2 = sqrt(gammat*Rt*T9_off_2);
    Pt9P0_off_2 = pi_r_off_2*pi_d*pi_c_off_2*pi_b*pi_t_3(best)*pi_n;
    P0P9_off_2 = (1/Pt9P0_off_2)*Pt9P9_cr_off;
    F_a(2) = mo_off_2*((1+f_3(best))*V9_off_2-V(2)+((1+f_3(best))*(V9_off_2/(gammat*(M9_off_2^2)))*(1-(P0P9_off_2))));
    ct2 = f_3(best)/(F_a(2)/mo_off_2);
    
    % required thrust
    a2 = 0.5*pho(2)*V(2)^2*b;
    CL(2) = W(4)*cos(theta_avg)/a2;
    CD(2) = (0.01 + 0.065*CL(2)^2);
    F_r(2) = 1/4*(a2*CD(2)+W(4)*sin(theta_avg));
    R(2) = 2*sqrt(2/(pho(2)*S))*1/ct2*CL(2)^0.5/CD(2)*(W(3)^0.5 - W(5)^0.5);
        
    if F_a(2) < F_r(2)
        mo_off_r_2 = F_r(2)/F_a(2)*mo_off_2;
        mo_ref_req = mo_off_r_2/mo_rel_2;
    end

end


end








