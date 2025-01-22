%% Simulation of SG formation in 3D space
rng('shuffle')

% parameter
kT = 1.0;
C = 0.1;
eta = 1;

zeta_40s = 1;
zeta_eIF3 = 1;
zeta_43s = 1;
zeta_RNA = 1;
zeta_60s = 1;
zeta_80s1RNA = 1;
zeta_80s2RNA = 1;
zeta_80s3RNA = 1;
zeta_80s4RNA = 1;
zeta_80s5RNA = 1;
zeta_80s = 1;
zeta_G3BP = 1;

% size
L = 31;
L3 = L^3;
dL = 1;

% total concentration
phi_40s_tot = 17.58;
phi_60s_tot = 15.92;
phi_RNA_tot = 15.65*0.9/5;
phi_eIF3_tot = 1.2;
phi_G3BP_tot = 10;

% parameters before and after stress on
k_40s_to_43s_ini = 0.07;
k_80sRNA_form_ini = 5;
k_80sRNA_break_ini = 0.012;
k_80s_break_ini = 0.002;
k_40s_bind_60s_ini = 0.008;
k_60s_break_ini = 0;

k_40s_to_43s_stress = k_40s_to_43s_ini;
k_80sRNA_form_stress = 4e-4;
k_80sRNA_break_stress = k_80sRNA_break_ini;
k_80s_break_stress = k_80s_break_ini;
k_40s_bind_60s_stress = 0.02;
k_60s_break_stress = 0;

% parameters of chemical reaction for initial condition
k_40s_to_43s = k_40s_to_43s_ini;
k_80sRNA_form = k_80sRNA_form_ini;
k_80sRNA_break = k_80sRNA_break_ini;
k_80s_break = k_80s_break_ini;
k_40s_bind_60s = k_40s_bind_60s_ini;
k_60s_break = k_60s_break_ini;

% lower limit for chemical reaction
phi_chem_thresh = 1e-2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parameters for interaction energy
lambda_mu_penalty = 0.012;
lambda_mu_penalty_RNA = 0.025;
lambda_mu_penalty_G3BP = 0.025;

% generate chi
chi = zeros(11,11); % 40s=1, 43s=2, RNA=3, 60s=4, 80s1-5RNA=5-9, 80s=10, G3BP=11 from 1 to 11
chi(3,3) = 0;
chi(11,11) = -0.3;
chi(11,3) = -0.7;       chi(3,11) = chi(11,3);
chi(1,11) = -0.039;     chi(11,1) = chi(1,11);
chi(2,11) = -0.09;     chi(11,2) = chi(2,11);
chi(4,11) = 0.138;      chi(11,4) = chi(4,11);
chi(5,11) = 0.0;      chi(11,5) = chi(5,11);
chi(6,11) = 0.0;      chi(11,6) = chi(6,11);
chi(7,11) = 0.0;      chi(11,7) = chi(7,11);
chi(8,11) = 0.0;      chi(11,8) = chi(8,11);
chi(9,11) = 0.0;      chi(11,9) = chi(9,11);
chi(10,11) = 0.217;     chi(11,10) = chi(10,11);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initial condition
phi_43s = 6e-1*ones(L,L,L); % PIC
phi_eIF3 = phi_eIF3_tot*ones(L,L,L)-phi_43s;
phi_80s1RNA = 2e-1*ones(L,L,L);
phi_80s2RNA = 2e-1*ones(L,L,L);
phi_80s3RNA = 2e-1*ones(L,L,L);
phi_80s4RNA = 2e-1*ones(L,L,L);
phi_80s5RNA = 1.5*ones(L,L,L);
phi_80s = 5e-1*ones(L,L,L);

randL_40s = 0.4 + 0.2*rand(L,L,L);
phi_40s = (phi_40s_tot*L3-sum(1*phi_80s1RNA(:))-sum(2*phi_80s2RNA(:))-sum(3*phi_80s3RNA(:)) ...
    -sum(4*phi_80s4RNA(:))-sum(5*phi_80s5RNA(:))-sum(phi_80s(:))-sum(phi_43s(:))).*randL_40s./sum(randL_40s(:));

randL_60s = 0.4 + 0.2*rand(L,L,L);
phi_60s = (phi_60s_tot*L3-sum(1*phi_80s1RNA(:))-sum(2*phi_80s2RNA(:))-sum(3*phi_80s3RNA(:)) ...
    -sum(4*phi_80s4RNA(:))-sum(5*phi_80s5RNA(:))-sum(phi_80s(:))).*randL_60s./sum(randL_60s(:));

randL_RNA = 0.4 + 0.2*rand(L,L,L);
phi_RNA = (phi_RNA_tot*L3-sum(phi_80s1RNA(:))-sum(phi_80s2RNA(:))-sum(phi_80s3RNA(:)) ...
    -sum(phi_80s4RNA(:))-sum(phi_80s5RNA(:))).*randL_RNA./sum(randL_RNA(:));

phi_G3BP = (phi_G3BP_tot-0.1) + 0.2*rand(L,L,L);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% simulation info
t_total = 5000;
t_stress_on = 500;
t_stress_duration = 9500;
dt = 0.01;
dt_rec = 10;

N_rec = round(t_total/dt_rec);
t_rec = zeros(1,N_rec);
phi_40s_rec = zeros(L3,N_rec);
phi_eIF3_rec = zeros(L3,N_rec);
phi_43s_rec = zeros(L3,N_rec);
phi_RNA_rec = zeros(L3,N_rec);
phi_60s_rec = zeros(L3,N_rec);
phi_80s1RNA_rec = zeros(L3,N_rec);
phi_80s2RNA_rec = zeros(L3,N_rec);
phi_80s3RNA_rec = zeros(L3,N_rec);
phi_80s4RNA_rec = zeros(L3,N_rec);
phi_80s5RNA_rec = zeros(L3,N_rec);
phi_80s_rec = zeros(L3,N_rec);
phi_G3BP_rec = zeros(L3,N_rec);

zetaphi_coef = ones(L,L,L);

t = 0;
t_next = 0;
count = 1;
tic;
while t <= t_total
    if t >= t_next
        disp(t);

        t_rec(count) = t;
        phi_40s_rec(:,count) = phi_40s(:);
        phi_eIF3_rec(:,count) = phi_eIF3(:);
        phi_43s_rec(:,count) = phi_43s(:);
        phi_RNA_rec(:,count) = phi_RNA(:);
        phi_60s_rec(:,count) = phi_60s(:);
        phi_80s1RNA_rec(:,count) = phi_80s1RNA(:);
        phi_80s2RNA_rec(:,count) = phi_80s2RNA(:);
        phi_80s3RNA_rec(:,count) = phi_80s3RNA(:);
        phi_80s4RNA_rec(:,count) = phi_80s4RNA(:);
        phi_80s5RNA_rec(:,count) = phi_80s5RNA(:);
        phi_80s_rec(:,count) = phi_80s(:);
        phi_G3BP_rec(:,count) = phi_G3BP(:);

        count = count + 1;
        t_next = t_next + dt_rec;

        k_80sRNA_break = 0.09/10 + 0.06/10.*rand(L,L,L);

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if t >= t_stress_on & t < t_stress_on+t_stress_duration
        k_40s_to_43s = k_40s_to_43s_stress;
        k_80sRNA_form = k_80sRNA_form_stress;
        % k_80sRNA_break = k_80sRNA_break_stress;
        k_80s_break = k_80s_break_stress;
        k_40s_bind_60s = k_40s_bind_60s_stress;
        k_60s_break = k_60s_break_stress;
        dt = 0.005;
    end

    if t >= t_stress_on+t_stress_duration
        k_40s_to_43s = k_40s_to_43s_ini;
        k_80sRNA_form = k_80sRNA_form_ini;
        k_80sRNA_break = k_80sRNA_break_ini;
        k_80s_break = k_80s_break_ini;
        k_40s_bind_60s = k_40s_bind_60s_ini;
        k_60s_break = k_60s_break_ini;
        dt = 0.01;
    end

    if t >= t_stress_on+t_stress_duration+500
        dt = 0.01;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % chemical reaction fluxes
    J_40s_to_43s = k_40s_to_43s.*phi_40s.*(phi_40s>phi_chem_thresh).*phi_eIF3.*(phi_eIF3>phi_chem_thresh);
    J_80s1RNA_form = k_80sRNA_form.*phi_43s.*(phi_43s>phi_chem_thresh).*phi_60s.*(phi_60s>phi_chem_thresh)...
        .*phi_RNA.*(phi_RNA>phi_chem_thresh);
    J_80s1RNA_break = k_80sRNA_break.*phi_80s1RNA.*(phi_80s1RNA>5*phi_chem_thresh);
    J_80s2RNA_form = k_80sRNA_form.*phi_43s.*(phi_43s>phi_chem_thresh).*phi_60s.*(phi_60s>phi_chem_thresh)...
        .*phi_80s1RNA.*(phi_80s1RNA>phi_chem_thresh);
    J_80s2RNA_break = k_80sRNA_break.*phi_80s2RNA.*(phi_80s2RNA>5*phi_chem_thresh/2);
    J_80s3RNA_form = k_80sRNA_form.*phi_43s.*(phi_43s>phi_chem_thresh).*phi_60s.*(phi_60s>phi_chem_thresh)...
        .*phi_80s2RNA.*(phi_80s2RNA>phi_chem_thresh);
    J_80s3RNA_break = k_80sRNA_break.*phi_80s3RNA.*(phi_80s3RNA>5*phi_chem_thresh/3);
    J_80s4RNA_form = k_80sRNA_form.*phi_43s.*(phi_43s>phi_chem_thresh).*phi_60s.*(phi_60s>phi_chem_thresh)...
        .*phi_80s3RNA.*(phi_80s3RNA>phi_chem_thresh);
    J_80s4RNA_break = k_80sRNA_break.*phi_80s4RNA.*(phi_80s4RNA>5*phi_chem_thresh/4);
    J_80s5RNA_form = k_80sRNA_form.*phi_43s.*(phi_43s>phi_chem_thresh).*phi_60s.*(phi_60s>phi_chem_thresh)...
        .*phi_80s4RNA.*(phi_80s4RNA>phi_chem_thresh);
    J_80s5RNA_break = k_80sRNA_break.*phi_80s5RNA.*(phi_80s5RNA>5*phi_chem_thresh/5);
    J_40s_bind_60s = k_40s_bind_60s.*phi_40s.*(phi_40s>phi_chem_thresh).*phi_60s.*(phi_60s>phi_chem_thresh);
    J_80s_break = k_80s_break.*phi_80s.*(phi_80s>phi_chem_thresh);
    J_60s_break = k_60s_break.*(phi_60s>0.135);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % interactions between components
    idx = 1;
    U_40s = chi(idx,1).*phi_40s + chi(idx,2).*phi_43s + chi(idx,3).*phi_RNA + chi(idx,4).*phi_60s + ...
        chi(idx,5).*phi_80s1RNA + chi(idx,6).*phi_80s2RNA + chi(idx,7).*phi_80s3RNA + chi(idx,8).*phi_80s4RNA ...
        + chi(idx,9).*phi_80s5RNA + chi(idx,10).*phi_80s + chi(idx,11).*phi_G3BP;
    [v_x_temp,v_y_temp,v_z_temp] = delta(-(log(phi_40s) + lambda_mu_penalty*phi_40s.^2 - C*delta2(phi_40s,dL) + U_40s),dL);
    v_x = kT./(zeta_40s.*zetaphi_coef).*v_x_temp;
    v_y = kT./(zeta_40s.*zetaphi_coef).*v_y_temp;
    v_z = kT./(zeta_40s.*zetaphi_coef).*v_z_temp;
    J = - J_40s_to_43s + J_80s1RNA_break + J_80s2RNA_break + J_80s3RNA_break + J_80s4RNA_break + J_80s5RNA_break ...
        + J_80s_break - J_40s_bind_60s;
    dphi_40s = dt*(-div(phi_40s.*v_x,phi_40s.*v_y,phi_40s.*v_z,dL) + J);

    [v_x_temp,v_y_temp,v_z_temp] = delta(-(log(phi_eIF3) + lambda_mu_penalty*phi_eIF3.^2 - C*delta2(phi_eIF3,dL)),dL);
    v_x = kT./(zeta_eIF3.*zetaphi_coef).*v_x_temp;
    v_y = kT./(zeta_eIF3.*zetaphi_coef).*v_y_temp;
    v_z = kT./(zeta_eIF3.*zetaphi_coef).*v_z_temp;
    J = - J_40s_to_43s + J_80s1RNA_form + J_80s2RNA_form + J_80s3RNA_form + J_80s4RNA_form + J_80s5RNA_form;
    dphi_eIF3 = dt*(-div(phi_eIF3.*v_x,phi_eIF3.*v_y,phi_eIF3.*v_z,dL) + J);
    
    idx = 2;
    U_43s = chi(idx,1).*phi_40s + chi(idx,2).*phi_43s + chi(idx,3).*phi_RNA + chi(idx,4).*phi_60s + ...
        chi(idx,5).*phi_80s1RNA + chi(idx,6).*phi_80s2RNA + chi(idx,7).*phi_80s3RNA + chi(idx,8).*phi_80s4RNA ...
        + chi(idx,9).*phi_80s5RNA + chi(idx,10).*phi_80s + chi(idx,11).*phi_G3BP;
    [v_x_temp,v_y_temp,v_z_temp] = delta(-(log(phi_43s) + lambda_mu_penalty*phi_43s.^2 - C*delta2(phi_43s,dL) + U_43s),dL);
    v_x = kT./(zeta_43s.*zetaphi_coef).*v_x_temp;
    v_y = kT./(zeta_43s.*zetaphi_coef).*v_y_temp;
    v_z = kT./(zeta_43s.*zetaphi_coef).*v_z_temp;
    J = + J_40s_to_43s - J_80s1RNA_form - J_80s2RNA_form - J_80s3RNA_form - J_80s4RNA_form - J_80s5RNA_form;
    dphi_43s = dt*(-div(phi_43s.*v_x,phi_43s.*v_y,phi_43s.*v_z,dL) + J);
    
    idx = 3;
    U_RNA = chi(idx,1).*phi_40s + chi(idx,2).*phi_43s + chi(idx,3).*phi_RNA + chi(idx,4).*phi_60s + ...
        chi(idx,5).*phi_80s1RNA + chi(idx,6).*phi_80s2RNA + chi(idx,7).*phi_80s3RNA + chi(idx,8).*phi_80s4RNA ...
        + chi(idx,9).*phi_80s5RNA + chi(idx,10).*phi_80s + chi(idx,11).*phi_G3BP;
    [v_x_temp,v_y_temp,v_z_temp] = delta(-(log(phi_RNA) + lambda_mu_penalty_RNA*phi_RNA.^2 - C*delta2(phi_RNA,dL) + U_RNA),dL);
    v_x = kT./(zeta_RNA.*zetaphi_coef).*v_x_temp;
    v_y = kT./(zeta_RNA.*zetaphi_coef).*v_y_temp;
    v_z = kT./(zeta_RNA.*zetaphi_coef).*v_z_temp;
    J = - J_80s1RNA_form + J_80s1RNA_break;
    dphi_RNA = dt*(-div(phi_RNA.*v_x,phi_RNA.*v_y,phi_RNA.*v_z,dL) + J);
    
    idx = 4;
    U_60s = chi(idx,1).*phi_40s + chi(idx,2).*phi_43s + chi(idx,3).*phi_RNA + chi(idx,4).*phi_60s + ...
        chi(idx,5).*phi_80s1RNA + chi(idx,6).*phi_80s2RNA + chi(idx,7).*phi_80s3RNA + chi(idx,8).*phi_80s4RNA ...
        + chi(idx,9).*phi_80s5RNA + chi(idx,10).*phi_80s + chi(idx,11).*phi_G3BP;
    [v_x_temp,v_y_temp,v_z_temp] = delta(-(log(phi_60s) + lambda_mu_penalty*phi_60s.^2 - C*delta2(phi_60s,dL) + U_60s),dL);
    v_x = kT./(zeta_60s.*zetaphi_coef).*v_x_temp;
    v_y = kT./(zeta_60s.*zetaphi_coef).*v_y_temp;
    v_z = kT./(zeta_60s.*zetaphi_coef).*v_z_temp;
    J = - J_80s1RNA_form - J_80s2RNA_form - J_80s3RNA_form - J_80s4RNA_form - J_80s5RNA_form ...
        + J_80s1RNA_break + J_80s2RNA_break + J_80s3RNA_break + J_80s4RNA_break + J_80s5RNA_break + J_80s_break ...
        - J_40s_bind_60s - J_60s_break;
    dphi_60s = dt*(-div(phi_60s.*v_x,phi_60s.*v_y,phi_60s.*v_z,dL) + J);
    
    idx = 5;
    U_80s1RNA = chi(idx,1).*phi_40s + chi(idx,2).*phi_43s + chi(idx,3).*phi_RNA + chi(idx,4).*phi_60s + ...
        chi(idx,5).*phi_80s1RNA + chi(idx,6).*phi_80s2RNA + chi(idx,7).*phi_80s3RNA + chi(idx,8).*phi_80s4RNA ...
        + chi(idx,9).*phi_80s5RNA + chi(idx,10).*phi_80s + chi(idx,11).*phi_G3BP;
    [v_x_temp,v_y_temp,v_z_temp] = delta(-(log(phi_80s1RNA) + lambda_mu_penalty*phi_80s1RNA.^2 - C*delta2(phi_80s1RNA,dL) + U_80s1RNA),dL);
    v_x = kT./(zeta_80s1RNA.*zetaphi_coef).*v_x_temp;
    v_y = kT./(zeta_80s1RNA.*zetaphi_coef).*v_y_temp;
    v_z = kT./(zeta_80s1RNA.*zetaphi_coef).*v_z_temp;
    J = + J_80s1RNA_form - J_80s1RNA_break - J_80s2RNA_form + J_80s2RNA_break;
    dphi_80s1RNA = dt*(-div(phi_80s1RNA.*v_x,phi_80s1RNA.*v_y,phi_80s1RNA.*v_z,dL) + J);

    idx = 6;
    U_80s2RNA = chi(idx,1).*phi_40s + chi(idx,2).*phi_43s + chi(idx,3).*phi_RNA + chi(idx,4).*phi_60s + ...
        chi(idx,5).*phi_80s1RNA + chi(idx,6).*phi_80s2RNA + chi(idx,7).*phi_80s3RNA + chi(idx,8).*phi_80s4RNA ...
        + chi(idx,9).*phi_80s5RNA + chi(idx,10).*phi_80s + chi(idx,11).*phi_G3BP;
    [v_x_temp,v_y_temp,v_z_temp] = delta(-(log(phi_80s2RNA) + lambda_mu_penalty*phi_80s2RNA.^2 - C*delta2(phi_80s2RNA,dL) + U_80s2RNA),dL);
    v_x = kT./(zeta_80s2RNA.*zetaphi_coef).*v_x_temp;
    v_y = kT./(zeta_80s2RNA.*zetaphi_coef).*v_y_temp;
    v_z = kT./(zeta_80s2RNA.*zetaphi_coef).*v_z_temp;
    J = + J_80s2RNA_form - J_80s2RNA_break - J_80s3RNA_form + J_80s3RNA_break;
    dphi_80s2RNA = dt*(-div(phi_80s2RNA.*v_x,phi_80s2RNA.*v_y,phi_80s2RNA.*v_z,dL) + J);

    idx = 7;
    U_80s3RNA = chi(idx,1).*phi_40s + chi(idx,2).*phi_43s + chi(idx,3).*phi_RNA + chi(idx,4).*phi_60s + ...
        chi(idx,5).*phi_80s1RNA + chi(idx,6).*phi_80s2RNA + chi(idx,7).*phi_80s3RNA + chi(idx,8).*phi_80s4RNA ...
        + chi(idx,9).*phi_80s5RNA + chi(idx,10).*phi_80s + chi(idx,11).*phi_G3BP;
    [v_x_temp,v_y_temp,v_z_temp] = delta(-(log(phi_80s3RNA) + lambda_mu_penalty*phi_80s3RNA.^2 - C*delta2(phi_80s3RNA,dL) + U_80s3RNA),dL);
    v_x = kT./(zeta_80s3RNA.*zetaphi_coef).*v_x_temp;
    v_y = kT./(zeta_80s3RNA.*zetaphi_coef).*v_y_temp;
    v_z = kT./(zeta_80s3RNA.*zetaphi_coef).*v_z_temp;
    J = + J_80s3RNA_form - J_80s3RNA_break - J_80s4RNA_form + J_80s4RNA_break;
    dphi_80s3RNA = dt*(-div(phi_80s3RNA.*v_x,phi_80s3RNA.*v_y,phi_80s3RNA.*v_z,dL) + J);

    idx = 8;
    U_80s4RNA = chi(idx,1).*phi_40s + chi(idx,2).*phi_43s + chi(idx,3).*phi_RNA + chi(idx,4).*phi_60s + ...
        chi(idx,5).*phi_80s1RNA + chi(idx,6).*phi_80s2RNA + chi(idx,7).*phi_80s3RNA + chi(idx,8).*phi_80s4RNA ...
        + chi(idx,9).*phi_80s5RNA + chi(idx,10).*phi_80s + chi(idx,11).*phi_G3BP;
    [v_x_temp,v_y_temp,v_z_temp] = delta(-(log(phi_80s4RNA) + lambda_mu_penalty*phi_80s4RNA.^2 - C*delta2(phi_80s4RNA,dL) + U_80s4RNA),dL);
    v_x = kT./(zeta_80s4RNA.*zetaphi_coef).*v_x_temp;
    v_y = kT./(zeta_80s4RNA.*zetaphi_coef).*v_y_temp;
    v_z = kT./(zeta_80s4RNA.*zetaphi_coef).*v_z_temp;
    J = + J_80s4RNA_form - J_80s4RNA_break - J_80s5RNA_form + J_80s5RNA_break;
    dphi_80s4RNA = dt*(-div(phi_80s4RNA.*v_x,phi_80s4RNA.*v_y,phi_80s4RNA.*v_z,dL) + J);

    idx = 9;
    U_80s5RNA = chi(idx,1).*phi_40s + chi(idx,2).*phi_43s + chi(idx,3).*phi_RNA + chi(idx,4).*phi_60s + ...
        chi(idx,5).*phi_80s1RNA + chi(idx,6).*phi_80s2RNA + chi(idx,7).*phi_80s3RNA + chi(idx,8).*phi_80s4RNA ...
        + chi(idx,9).*phi_80s5RNA + chi(idx,10).*phi_80s + chi(idx,11).*phi_G3BP;
    [v_x_temp,v_y_temp,v_z_temp] = delta(-(log(phi_80s5RNA) + lambda_mu_penalty*phi_80s5RNA.^2 - C*delta2(phi_80s5RNA,dL) + U_80s5RNA),dL);
    v_x = kT./(zeta_80s5RNA.*zetaphi_coef).*v_x_temp;
    v_y = kT./(zeta_80s5RNA.*zetaphi_coef).*v_y_temp;
    v_z = kT./(zeta_80s5RNA.*zetaphi_coef).*v_z_temp;
    J = + J_80s5RNA_form - J_80s5RNA_break;
    dphi_80s5RNA = dt*(-div(phi_80s5RNA.*v_x,phi_80s5RNA.*v_y,phi_80s5RNA.*v_z,dL) + J);

    idx = 10;
    U_80s = chi(idx,1).*phi_40s + chi(idx,2).*phi_43s + chi(idx,3).*phi_RNA + chi(idx,4).*phi_60s + ...
        chi(idx,5).*phi_80s1RNA + chi(idx,6).*phi_80s2RNA + chi(idx,7).*phi_80s3RNA + chi(idx,8).*phi_80s4RNA ...
        + chi(idx,9).*phi_80s5RNA + chi(idx,10).*phi_80s + chi(idx,11).*phi_G3BP;
    [v_x_temp,v_y_temp,v_z_temp] = delta(-(log(phi_80s) + lambda_mu_penalty*phi_80s.^2 - C*delta2(phi_80s,dL) + U_80s),dL);
    v_x = kT./(zeta_80s.*zetaphi_coef).*v_x_temp;
    v_y = kT./(zeta_80s.*zetaphi_coef).*v_y_temp;
    v_z = kT./(zeta_80s.*zetaphi_coef).*v_z_temp;
    J = + J_40s_bind_60s - J_80s_break;
    dphi_80s = dt*(-div(phi_80s.*v_x,phi_80s.*v_y,phi_80s.*v_z,dL) + J);

    idx = 11;
    U_G3BP = chi(idx,1).*phi_40s + chi(idx,2).*phi_43s + chi(idx,3).*phi_RNA + chi(idx,4).*phi_60s + ...
        chi(idx,5).*phi_80s1RNA + chi(idx,6).*phi_80s2RNA + chi(idx,7).*phi_80s3RNA + chi(idx,8).*phi_80s4RNA ...
        + chi(idx,9).*phi_80s5RNA + chi(idx,10).*phi_80s + chi(idx,11).*phi_G3BP;
    [v_x_temp,v_y_temp,v_z_temp] = delta(-(log(phi_G3BP) + lambda_mu_penalty_G3BP*phi_G3BP.^2 - C*delta2(phi_G3BP,dL) + U_G3BP),dL);
    v_x = kT./(zeta_G3BP.*zetaphi_coef).*v_x_temp;
    v_y = kT./(zeta_G3BP.*zetaphi_coef).*v_y_temp;
    v_z = kT./(zeta_G3BP.*zetaphi_coef).*v_z_temp;
    J = 0;
    dphi_G3BP = dt*(-div(phi_G3BP.*v_x,phi_G3BP.*v_y,phi_G3BP.*v_z,dL) + J);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    t = t + dt;

    phi_40s = phi_40s + dphi_40s;
    phi_eIF3 = phi_eIF3 + dphi_eIF3;
    phi_43s = phi_43s + dphi_43s;
    phi_RNA = phi_RNA + dphi_RNA;
    phi_60s = phi_60s + dphi_60s;
    phi_80s1RNA = phi_80s1RNA + dphi_80s1RNA;
    phi_80s2RNA = phi_80s2RNA + dphi_80s2RNA;
    phi_80s3RNA = phi_80s3RNA + dphi_80s3RNA;
    phi_80s4RNA = phi_80s4RNA + dphi_80s4RNA;
    phi_80s5RNA = phi_80s5RNA + dphi_80s5RNA;
    phi_80s = phi_80s + dphi_80s;
    phi_G3BP = phi_G3BP + dphi_G3BP;

    phi_check = phi_40s + phi_43s + phi_RNA + phi_60s + phi_80s1RNA + phi_80s2RNA + phi_80s3RNA ...
        + phi_80s4RNA + phi_80s5RNA + phi_80s + phi_G3BP;
    if max(abs(imag(phi_check(:)))) > 0 || ~(max(abs(phi_check(:))) > 0)
        disp('Result Not Real');
        break;
    end
end
toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data Processing

%% Plot 3D results
for i = 10:10:N_rec
    phi_plot = reshape(phi_G3BP_rec(:,i),L,L,L);
    h = slice(phi_plot, [], [], 1:size(phi_plot,3));colorbar;
    set(h, 'EdgeColor','none', 'FaceColor','interp');
    alpha(.15);
    hold on;
    title(['t=',num2str(i*dt_rec)])
    hold off;
    pause(0.05)
end


%% Calculate total amounts of ribosomal subunits
phi_80splot_rec = phi_80s_rec+phi_80s1RNA_rec+2*phi_80s2RNA_rec+3*phi_80s3RNA_rec+4*phi_80s4RNA_rec+5*phi_80s5RNA_rec;

i = 50;
phi_idx_1 = 1:L;
phi_idx_2 = 1:L;
phi_idx_3 = 1:L;
phi_40s_plot = circshift(reshape(phi_40s_rec(:,i),L,L,L),[0,0,6]);
phi_40s_plot = phi_40s_plot(phi_idx_1,phi_idx_2,phi_idx_3);
phi_43s_plot = circshift(reshape(phi_43s_rec(:,i),L,L,L),[0,0,6]);
phi_43s_plot = phi_43s_plot(phi_idx_1,phi_idx_2,phi_idx_3);
phi_60s_plot = circshift(reshape(phi_60s_rec(:,i),L,L,L),[0,0,6]);
phi_60s_plot = phi_60s_plot(phi_idx_1,phi_idx_2,phi_idx_3);
phi_80s_plot = circshift(reshape(phi_80splot_rec(:,i),L,L,L),[0,0,6]);
phi_80s_plot = phi_80s_plot(phi_idx_1,phi_idx_2,phi_idx_3);
phi_RNA_ref = circshift(reshape(phi_RNA_rec(:,i),L,L,L),[0,0,6]);
phi_RNA_ref = phi_RNA_ref(phi_idx_1,phi_idx_2,phi_idx_3);
phi_RNA_plot = phi_RNA_ref;
[Lx,Ly,Lz] = size(phi_RNA_ref);
[pos_SG,N_SG] = find_droplet_3D(phi_RNA_ref,5,Lx,Ly,Lz);
ball_num_mat = zeros(4,N_SG+1);
% inside
for SG_i = 1:N_SG
    ball_num_mat(1,SG_i) = sum(phi_40s_plot(pos_SG{SG_i}));
    ball_num_mat(2,SG_i) = sum(phi_43s_plot(pos_SG{SG_i}));
    ball_num_mat(3,SG_i) = sum(phi_60s_plot(pos_SG{SG_i}));
    ball_num_mat(4,SG_i) = sum(phi_80s_plot(pos_SG{SG_i}));
    % ball_num_mat(5,SG_i) = sum(phi_RNA_plot(pos_SG{SG_i}));
end
% outside
pos_Cyto = 1:length(phi_40s_plot(:));
for SG_i = 1:N_SG
    pos_Cyto = setdiff(pos_Cyto,pos_SG{SG_i});
end
ball_num_mat(1,N_SG+1) = sum(phi_40s_plot(pos_Cyto));
ball_num_mat(2,N_SG+1) = sum(phi_43s_plot(pos_Cyto));
ball_num_mat(3,N_SG+1) = sum(phi_60s_plot(pos_Cyto));
ball_num_mat(4,N_SG+1) = sum(phi_80s_plot(pos_Cyto));
% ball_num_mat(5,N_SG+1) = sum(phi_RNA_plot(pos_Cyto));

ball_num_mat = round(ball_num_mat.*100); % .*100 for real number
particle_num_before = ball_num_mat;


i = 500;
phi_idx_1 = 1:L;
phi_idx_2 = 1:L;
phi_idx_3 = 1:L;
phi_40s_plot = circshift(reshape(phi_40s_rec(:,i),L,L,L),[0,0,6]);
phi_40s_plot = phi_40s_plot(phi_idx_1,phi_idx_2,phi_idx_3);
phi_43s_plot = circshift(reshape(phi_43s_rec(:,i),L,L,L),[0,0,6]);
phi_43s_plot = phi_43s_plot(phi_idx_1,phi_idx_2,phi_idx_3);
phi_60s_plot = circshift(reshape(phi_60s_rec(:,i),L,L,L),[0,0,6]);
phi_60s_plot = phi_60s_plot(phi_idx_1,phi_idx_2,phi_idx_3);
phi_80s_plot = circshift(reshape(phi_80splot_rec(:,i),L,L,L),[0,0,6]);
phi_80s_plot = phi_80s_plot(phi_idx_1,phi_idx_2,phi_idx_3);
phi_RNA_ref = circshift(reshape(phi_RNA_rec(:,i),L,L,L),[0,0,6]);
phi_RNA_ref = phi_RNA_ref(phi_idx_1,phi_idx_2,phi_idx_3);
phi_RNA_plot = phi_RNA_ref;
[Lx,Ly,Lz] = size(phi_RNA_ref);
[pos_SG,N_SG] = find_droplet_3D(phi_RNA_ref,5,Lx,Ly,Lz);
ball_num_mat = zeros(4,N_SG+1);
% inside
for SG_i = 1:N_SG
    ball_num_mat(1,SG_i) = sum(phi_40s_plot(pos_SG{SG_i}));
    ball_num_mat(2,SG_i) = sum(phi_43s_plot(pos_SG{SG_i}));
    ball_num_mat(3,SG_i) = sum(phi_60s_plot(pos_SG{SG_i}));
    ball_num_mat(4,SG_i) = sum(phi_80s_plot(pos_SG{SG_i}));
    % ball_num_mat(5,SG_i) = sum(phi_RNA_plot(pos_SG{SG_i}));
end
% outside
pos_Cyto = 1:length(phi_40s_plot(:));
for SG_i = 1:N_SG
    pos_Cyto = setdiff(pos_Cyto,pos_SG{SG_i});
end
ball_num_mat(1,N_SG+1) = sum(phi_40s_plot(pos_Cyto));
ball_num_mat(2,N_SG+1) = sum(phi_43s_plot(pos_Cyto));
ball_num_mat(3,N_SG+1) = sum(phi_60s_plot(pos_Cyto));
ball_num_mat(4,N_SG+1) = sum(phi_80s_plot(pos_Cyto));
% ball_num_mat(5,N_SG+1) = sum(phi_RNA_plot(pos_Cyto));

ball_num_mat = round(ball_num_mat.*100); % .*100 for real number
particle_num_after = [sum(ball_num_mat(:,1:end-1),2),ball_num_mat(:,end)];
particle_num_sum = [particle_num_before,particle_num_after];

particle_num_sum


%% BarPlot Control and inside and outside concentrations after SG formation
load("ExpData.mat");

phi_RNA_in_thresh = 10;
phi_RNA_out_thresh = 1.5;

data_40s_bar = [data_40s_exp(1),mean(phi_40s_rec(:,50));data_40s_exp(2),mean(phi_40s(phi_RNA(:)<phi_RNA_out_thresh));data_40s_exp(3),mean(phi_40s(phi_RNA(:)>phi_RNA_in_thresh))];
error_40s_bar = [error_40s_exp(1),std(phi_40s_rec(:,50));error_40s_exp(2),std(phi_40s(phi_RNA(:)<phi_RNA_out_thresh));error_40s_exp(3),std(phi_40s(phi_RNA(:)>phi_RNA_in_thresh))];
figure;
hb = bar(data_40s_bar); % get the bar handles
hold on;
pause(1);
for k = 1:size(data_40s_bar,2)
% get x positions per group
xpos = hb(k).XData + hb(k).XOffset;
% draw errorbar
errorbar(xpos, data_40s_bar(:,k), error_40s_bar(:,k), 'LineStyle', 'none', ...
    'Color', 'k', 'LineWidth', 1);
end
set(gca,'xticklabel',{'Control'; 'Cyto';'SG'});
set(gca,'FontSize',15)
ylabel('Concentration/\muM')
title('40S')
legend('Experiment','Simulation')

%%%%%%%%%%

data_43s_bar = [data_43s_exp(1),mean(phi_43s_rec(:,50));data_43s_exp(2),mean(phi_43s(phi_RNA(:)<phi_RNA_out_thresh));data_43s_exp(3),mean(phi_43s(phi_RNA(:)>phi_RNA_in_thresh))];
error_43s_bar = [error_43s_exp(1),std(phi_43s_rec(:,50));error_43s_exp(2),std(phi_43s(phi_RNA(:)<phi_RNA_out_thresh));error_43s_exp(3),std(phi_43s(phi_RNA(:)>phi_RNA_in_thresh))];
figure;
hb = bar(data_43s_bar); % get the bar handles
hold on;
pause(1);
for k = 1:size(data_43s_bar,2)
% get x positions per group
xpos = hb(k).XData + hb(k).XOffset;
% draw errorbar
errorbar(xpos, data_43s_bar(:,k), error_43s_bar(:,k), 'LineStyle', 'none', ...
    'Color', 'k', 'LineWidth', 1);
end
set(gca,'xticklabel',{'Control'; 'Cyto';'SG'});
set(gca,'FontSize',15)
ylabel('Concentration/\muM')
title('PIC')
legend('Experiment','Simulation')

%%%%%%%%%%

data_60s_bar = [data_60s_exp(1),mean(phi_60s_rec(:,50));data_60s_exp(2),mean(phi_60s(phi_RNA(:)<phi_RNA_out_thresh));data_60s_exp(3),mean(phi_60s(phi_RNA(:)>phi_RNA_in_thresh))];
error_60s_bar = [error_60s_exp(1),std(phi_60s_rec(:,50));error_60s_exp(2),std(phi_60s(phi_RNA(:)<phi_RNA_out_thresh));error_60s_exp(3),std(phi_60s(phi_RNA(:)>phi_RNA_in_thresh))];
figure;
hb = bar(data_60s_bar); % get the bar handles
hold on;
pause(1);
for k = 1:size(data_60s_bar,2)
% get x positions per group
xpos = hb(k).XData + hb(k).XOffset;
% draw errorbar
errorbar(xpos, data_60s_bar(:,k), error_60s_bar(:,k), 'LineStyle', 'none', ...
    'Color', 'k', 'LineWidth', 1);
end
set(gca,'xticklabel',{'Control'; 'Cyto';'SG'});
set(gca,'FontSize',15)
ylabel('Concentration/\muM')
title('60S')
legend('Experiment','Simulation')

%%%%%%%%%%

phi_80s_plot = phi_80s+phi_80s1RNA+2*phi_80s2RNA+3*phi_80s3RNA+4*phi_80s4RNA+5*phi_80s5RNA;
phi_80splot_rec = phi_80s_rec+phi_80s1RNA_rec+2*phi_80s2RNA_rec+3*phi_80s3RNA_rec+4*phi_80s4RNA_rec+5*phi_80s5RNA_rec;

data_80s_bar = [data_80s_exp(1),mean(phi_80splot_rec(:,50));data_80s_exp(2),mean(phi_80s_plot(phi_RNA(:)<phi_RNA_out_thresh));data_80s_exp(3),mean(phi_80s_plot(phi_RNA(:)>phi_RNA_in_thresh))];
error_80s_bar = [error_80s_exp(1),std(phi_80splot_rec(:,50));error_80s_exp(2),std(phi_80s_plot(phi_RNA(:)<phi_RNA_out_thresh));error_80s_exp(3),std(phi_80s_plot(phi_RNA(:)>phi_RNA_in_thresh))];
figure;
hb = bar(data_80s_bar); % get the bar handles
hold on;
pause(1);
for k = 1:size(data_80s_bar,2)
% get x positions per group
xpos = hb(k).XData + hb(k).XOffset;
% draw errorbar
errorbar(xpos, data_80s_bar(:,k), error_80s_bar(:,k), 'LineStyle', 'none', ...
    'Color', 'k', 'LineWidth', 1);
end
set(gca,'xticklabel',{'Control'; 'Cyto';'SG'});
set(gca,'FontSize',15)
ylabel('Concentration/\muM')
title('80S + 80S-mRNA')
legend('Experiment','Simulation')

%%%%%%%%%%

phi_80s_ratio_before = sum(phi_80s_rec(:,50))/sum(phi_80splot_rec(:,50));
disp(phi_80s_ratio_before);
phi_80s_ratio_after = sum(phi_80s_rec(:,end))/sum(phi_80splot_rec(:,end));
disp(1-phi_80s_ratio_after);




%% Demanded function
function divA = div(Ax,Ay,Az,d)
    Axp1 = circshift(Ax,[0,-1,0]);
    Axm1 = circshift(Ax,[0,1,0]);
    Axp2 = circshift(Ax,[0,-2,0]);
    Axm2 = circshift(Ax,[0,2,0]);
    dAx = (-Axp2+Axm2+8.*Axp1-8.*Axm1)/(12*d);
    Ayp1 = circshift(Ay,[-1,0,0]);
    Aym1 = circshift(Ay,[1,0,0]);
    Ayp2 = circshift(Ay,[-2,0,0]);
    Aym2 = circshift(Ay,[2,0,0]);
    dAy = (-Ayp2+Aym2+8.*Ayp1-8.*Aym1)/(12*d);
    Azp1 = circshift(Az,[0,0,-1]);
    Azm1 = circshift(Az,[0,0,1]);
    Azp2 = circshift(Az,[0,0,-2]);
    Azm2 = circshift(Az,[0,0,2]);
    dAz = (-Azp2+Azm2+8.*Azp1-8.*Azm1)/(12*d);
    divA = dAx + dAy + dAz;
end
function [dAx,dAy,dAz] = delta(A,d)
    Axp1 = circshift(A,[0,-1,0]);
    Axm1 = circshift(A,[0,1,0]);
    Axp2 = circshift(A,[0,-2,0]);
    Axm2 = circshift(A,[0,2,0]);
    dAx = (-Axp2+Axm2+8.*Axp1-8.*Axm1)/(12*d);
    Ayp1 = circshift(A,[-1,0,0]);
    Aym1 = circshift(A,[1,0,0]);
    Ayp2 = circshift(A,[-2,0,0]);
    Aym2 = circshift(A,[2,0,0]);
    dAy = (-Ayp2+Aym2+8.*Ayp1-8.*Aym1)/(12*d);
    Azp1 = circshift(A,[0,0,-1]);
    Azm1 = circshift(A,[0,0,1]);
    Azp2 = circshift(A,[0,0,-2]);
    Azm2 = circshift(A,[0,0,2]);
    dAz = (-Azp2+Azm2+8.*Azp1-8.*Azm1)/(12*d);
end
function d2A = delta2(A,d)
    Axp1 = circshift(A,[0,-1,0]);
    Axm1 = circshift(A,[0,1,0]);
    Ayp1 = circshift(A,[-1,0,0]);
    Aym1 = circshift(A,[1,0,0]);
    Azp1 = circshift(A,[0,0,-1]);
    Azm1 = circshift(A,[0,0,1]);
    d2A = (Axp1+Axm1+Ayp1+Aym1+Azp1+Azm1-6*A)/(d^2);
end
