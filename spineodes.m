function xprime = spineodes(t,x,p)
% SPINEODES: system of ODEs for dynamic spine model.

ampa_rise = p(1); % s
ampa_decay = p(2); % s
nmda_rise = p(3); % s
nmda_decay = p(4); % s
esyn = p(5); % mV; reversal potential of synapse
Cm = p(6); % uF/cm^2
Ra = p(7); % Ohm.cm
gbar_leak = p(8); % S/m^2 = pS/um^2
e_leak = p(9); % mV

% R-channel parameters
e_CaR = p(10); % mV; Ohmic eCa is higher, but this reduces GHK error (Oertner)
gbar_CaRspine = p(11); % pS/um^2
gbar_CaRdend = p(12); % pS/um^2
dend_area = p(13); % um^2

gbar_ampa = p(15); % S
gbar_nmda = p(16); % S
alpha = p(17); % gbar_ampa/headrad ratio
beta_sp = p(18);
beta_neck = p(19);
beta_d = p(20);
Ca0 = p(21);
Dca = p(22);
dend_vol = p(23);
diam = p(24); % m
L = p(25); % m
tau_mR = p(26); % s
tau_hR = p(27); % s
VmR_half = p(28); % mV
kmR = p(29); % 1/mV
VhR_half = p(30); % mV
khR = p(31); % 1/mV
kf = p(32); % /M/s
kb = p(33); % /s
Bt_sp = p(34); % M
Bt_neck = p(35); % M
Bt_d = p(36); % M
eta_ltp = p(37); % rate of potentiation
ltp_midCa = p(38);
ltp_slope = p(39);
eta_ltd = p(40); % rate of depression
ltd_midCa = p(41);
ltd_slope = p(42);
Istim = p(43); % stimulus current to dendrite

% T-channel parameters
tau_mT_0 = p(44);
tau_mT_A = p(45);
tau_mT_k = p(46);
tau_hT_0 = p(47);
tau_hT_A = p(48);
tau_hT_k = p(49);
VmT_half = p(50);
kmT = p(51);
VhT_half = p(52);
khT = p(53);
gbar_CaT_spine = p(54);
gbar_CaT_dend = p(55);

% L-channel parameters
tau_mL_0 = p(56);
tau_mL_A = p(57);
tau_mLVhalf = p(58);
tau_mL_k = p(59);
tau_hL_0 = p(60);
tau_hL_A = p(61);
tau_hLVhalf = p(62);
tau_hL_k = p(63);
VmL_half = p(64);
kmL = p(65);
VhL_half = p(66);
khL = p(67);
gbar_CaL_spine = p(68);
gbar_CaL_dend = p(69);

% N-channel parameters
alpha_m_N_A = p(70);
alpha_m_N_Vhalf = p(71); % mV; activation forward rate Vhalf
alpha_m_N_k = p(72); % mV; activation forward rate slope
beta_m_N_A = p(73); % s; activation backward rate prefactor
beta_m_N_Vhalf = p(74); % mV; activation backward rate threshold
beta_m_N_k = p(75); % mV; activation backward rate exponent
alpha_h_N_A = p(76); % s; activation forward rate prefactor
alpha_h_N_Vhalf = p(77); % mV; activation forward rate Vhalf
alpha_h_N_k = p(78); % mV; activation forward rate slope
beta_h_N_A = p(79); % s; activation backward rate prefactor
beta_h_N_Vhalf = p(80); % mV; activation backward rate threshold
beta_h_N_k = p(81); % mV; activation backward rate exponent
gbar_CaN_spine = p(82); % S/m^2  = pS/um^2
gbar_CaN_dend = p(83); % S/m^2 = pS/um^2



Cdend = dend_area*1E-12*Cm*1E-2;

neck_area = L*diam*pi;
neck_vol = L*pi*(diam)^2/4;

tpeakampa = (ampa_rise*ampa_decay/(ampa_decay-ampa_rise))*log(ampa_decay/ampa_rise);
fampa = 1/(-exp(-tpeakampa/ampa_rise) + exp(-tpeakampa/ampa_decay));
tpeaknmda = (nmda_rise*nmda_decay/(nmda_decay-nmda_rise))*log(nmda_decay/nmda_rise);
fnmda = 1/(-exp(-tpeaknmda/nmda_rise) + exp(-tpeaknmda/nmda_decay));


% List of ODEs to be solved (one for each of 31 state variables).
xprime = [-x(1)/ampa_decay;
-x(2)/ampa_rise;
-x(3)/nmda_decay;
-x(4)/nmda_rise;
1*( -x(31)*fampa*(x(1)-x(2))*(x(5)-esyn) - 1*(100e-12)*fnmda*(x(3)-x(4))*(1/(1+exp(-0.062*x(5))*1/3.57))*(x(5)-esyn) - gbar_leak*4*pi*x(30)^2*(x(5)-e_leak) - (x(5)-x(6))/(2*Ra*4*L/(pi*diam^2)) - gbar_CaRspine*x(11)^3*x(12)*(x(5)-e_CaR)*(4*pi*x(30)^2) - gbar_CaT_spine*x(15)*x(16)*(x(5)-e_CaR)*(4*pi*x(30)^2) - gbar_CaL_spine*x(19)^2*x(20)*(x(5)-e_CaR)*(4*pi*x(30)^2) - gbar_CaN_spine*x(23)^2*x(24)*(x(5)-e_CaR)*(4*pi*x(30)^2) ) / (4*pi*x(30)^2*Cm*1E-2);
1*( - (x(6)-x(5))/(2*Ra*4*L/(pi*diam^2)) - (x(6)-x(7))/(2*Ra*4*L/(pi*diam^2)) - gbar_leak*(L*diam*pi)*(x(6)-e_leak) )/ (L*pi*diam*Cm*1E-2);
1*( 1e3*Istim + -gbar_leak*dend_area*1E-12*(x(7)-e_leak) - (x(7)-x(6))/(2*Ra*4*L/(pi*diam^2)) - gbar_CaRdend*x(12)^3*x(13)*(x(7)-e_CaR)*dend_area*1e-12 - gbar_CaT_dend*x(15)*x(16)*(x(7)-e_CaR)*dend_area*1e-12 - gbar_CaL_dend*x(19)^2*x(20)*(x(7)-e_CaR)*dend_area*1e-12 - gbar_CaN_dend*x(23)^2*x(24)*(x(7)-e_CaR)*dend_area*1e-12)/Cdend;
(-1E-3*(1E-3*0.1*(100e-12)*fnmda*(x(3)-x(4))*(1/(1+exp(-0.062*x(5))*1/3.57))*(x(5)-e_CaR))/(2*96489*(4/3)*pi*x(30)^3) - beta_sp*(x(8)-Ca0)*(3/x(30)) + kb*(Bt_sp-x(27)) - kf*x(27)*x(8) - (Dca*(x(8)-x(9))*(pi*diam^2/4))/((L/2)*(4/3)*pi*x(30)^3) - (1e-3*1E-3*gbar_CaRspine*x(10)^3*x(11)*(x(5)-e_CaR)*(4*pi*x(30)^2) + 1e-3*1E-3*gbar_CaT_spine*x(15)*x(16)*(x(5)-e_CaR)*(4*pi*x(30)^2) + 1e-3*1E-3*gbar_CaL_spine*x(19)^2*x(20)*(x(5)-e_CaR)*(4*pi*x(30)^2) + 1e-3*1E-3*gbar_CaN_spine*x(23)^2*x(24)*(x(5)-e_CaR)*(4*pi*x(30)^2) )/(2*96489*(4/3)*pi*(x(30)^3)));
(-beta_neck*(x(9)-Ca0)*(neck_area/neck_vol) + kb*(Bt_neck-x(28)) - kf*x(28)*x(9) - (Dca*(x(9)-x(8))*(pi*diam^2/4))/((L/2)*neck_vol) - (Dca*(x(9)-x(10))*(pi*diam^2/4))/((L/2)*neck_vol));
(-beta_d*(x(10)-Ca0)*(dend_area/(dend_vol*1e-6)) + kb*(Bt_d-x(29)) - kf*x(29)*x(10) - (Dca*(x(10)-x(9))*(pi*diam^2/4))/((L/2)*dend_vol*1e-18) - ( 1e-3*1E-3*gbar_CaRdend*x(13)^3*x(14)*(x(7)-e_CaR)*dend_area + 1e-3*1E-3*gbar_CaT_dend*x(17)*x(18)*(x(7)-e_CaR)*dend_area + 1e-3*1E-3*gbar_CaL_dend*x(21)*x(22)*(x(7)-e_CaR)*dend_area + 1e-3*1E-3*gbar_CaN_dend*x(25)*x(26)*(x(7)-e_CaR)*dend_area )/(2*96489*dend_vol*1e-6));
(-x(11) + 1/(1+exp(-(x(5)-VmR_half)/kmR)))/tau_mR;
(-x(12) + 1/(1+exp(-(x(5)-VhR_half)/khR)))/tau_hR;
(-x(13) + 1/(1+exp(-(x(7)-VmR_half)/kmR)))/tau_mR;
(-x(14) + 1/(1+exp(-(x(7)-VhR_half)/khR)))/tau_hR;
(-x(15) + (1 /(1+exp(-(x(5)-VmT_half)/kmT))))/ ( tau_mT_0 + (tau_mT_A * exp(-x(5)/tau_mT_k)) );
(-x(16) + (1 /(1+exp(-(x(5)-VhT_half)/khT))))/ ( tau_hT_0 + (tau_hT_A * exp(-x(5)/tau_hT_k)) );
(-x(17) + (1 /(1+exp(-(x(7)-VmT_half)/kmT))))/ ( tau_mT_0 + (tau_mT_A * exp(-x(7)/tau_mT_k)) );
(-x(18) + (1 /(1+exp(-(x(7)-VhT_half)/khT))))/ ( tau_hT_0 + (tau_hT_A * exp(-x(7)/tau_hT_k)) );
(-x(19) + (1 /(1+exp(-(x(5)-VmL_half)/kmL))))/(tau_mL_0 + tau_mL_A /(1+exp(-(x(5)-tau_mLVhalf)/tau_mL_k)));
(-x(20) + (1 /(1+exp(-(x(5)-VhL_half)/khL))))/(tau_hL_0 + tau_hL_A /(1+exp(-(x(5)-tau_hLVhalf)/tau_hL_k)));
(-x(21) + (1 /(1+exp(-(x(7)-VmL_half)/kmL))))/(tau_mL_0 + tau_mL_A /(1+exp(-(x(7)-tau_mLVhalf)/tau_mL_k)));
(-x(22) + (1 /(1+exp(-(x(7)-VhL_half)/khL))))/(tau_hL_0 + tau_hL_A /(1+exp(-(x(7)-tau_hLVhalf)/tau_hL_k)));
(alpha_m_N_A*(x(5)-alpha_m_N_Vhalf)/(1-exp(-(x(5)-alpha_m_N_Vhalf)/alpha_m_N_k))) * (1-x(23)) - beta_m_N_A*exp(-(x(5)-beta_m_N_Vhalf)/beta_m_N_k) * x(23);
(alpha_h_N_A*exp(-(x(5)-alpha_h_N_Vhalf)/alpha_h_N_k))*(1-x(24)) - (beta_h_N_A/(1+exp(-(x(5)-beta_h_N_Vhalf)/beta_h_N_k)))*x(24);
(alpha_m_N_A*(x(7)-alpha_m_N_Vhalf)/(1-exp(-(x(7)-alpha_m_N_Vhalf)/alpha_m_N_k))) * (1-x(25)) - beta_m_N_A*exp(-(x(7)-beta_m_N_Vhalf)/beta_m_N_k) * x(25);
(alpha_h_N_A*exp(-(x(7)-alpha_h_N_Vhalf)/alpha_h_N_k))*(1-x(26)) - (beta_h_N_A/(1+exp(-(x(7)-beta_h_N_Vhalf)/beta_h_N_k)))*x(26);
kb*(Bt_sp-x(27)) - kf*x(27)*x(8);
kb*(Bt_neck-x(28)) - kf*x(28)*x(9);
kb*(Bt_d-x(29)) - kf*x(29)*x(10)
((3/(4*pi))^(1/3)*(1/3)*(4*pi*(x(30)^3)/3)^(-2/3))*alpha*(eta_ltp/(1+exp(-(x(8)-ltp_midCa)/ltp_slope)) - 0*eta_ltd/(1+exp(-(x(8)-ltd_midCa)/ltd_slope)));
eta_ltp/(1+exp(-(x(8)-ltp_midCa)/ltp_slope)) - 0*eta_ltd/(1+exp(-(x(8)-ltd_midCa)/ltd_slope))];
end



















