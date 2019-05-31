% spinemodel_fig6.m
%
% Biophysical dendritic spine model which was used to generate Figure 6 of
% O'Donnell et al, J Neurosci (2011). Please send any comments or questions
% to cian@salk.edu.
%
% - Detailed biophysical model of dendritic spine as metaplastic device. 
% - Three functional compartments: spinehead, dendrite, connected by spine neck.
% - [Ca+] tracked in spinehead. 
% Influx from NMDA receptors and CaV channels. 
% Outflux from diffusion to spine neck, internal buffering and membrane pumps 
% (last two lumped together). 
% Amplitude and timecourse of Ca transients depend on spine geometery and 
% synapse properties.
% - Synaptic plasticity depends on [Ca].
% - includes neck capacitance, and multiple CaV channels



%%%%%%%%%%%%%%%%%%%%%
% SETUP
%%%%%%%%%%%%%%%%%%%%%

% Run parameters
dt = 0.1e-3; % s
duration = 1500e-3; % s; 8000s for figure
starttime = cputime; % to calculate total CPU simulation time later

headvolvec = 1e-18*[0.01:0.01:0.3]'; % List of potential spine head volumes
headradvec = (3.*headvolvec./(4*pi)).^(1/3); % corresponding list of spine head radii

% Electrical parameters
Cm = 0.75; % uF/cm^2; Membrane capacitance
Ra = 200E-2; % Ohm.m; Axial resistivity
gbar_leak = 2; % S/m^2 = pS/um^2; Membrane leak conductance density
e_leak = -70; % mV; Leak conductance reversal potential
e_CaR = 10; % mV; Nernst eCa is higher (+50mV), but this approximation reduces GHK error (Koch book, 1999)

Vsp0 = -50; % Initial spine head voltage
Vd0 = Vsp0; % Initial dendritic voltage

% Spine
L = 0.5e-6; % m; neck length
diam = 1*0.1e-6; % m; neck diameter
Rneck = Ra*L*4/(pi*diam^2); % resulting neck resistance
%headrad = 1.5299e-07; % m; Specify initial spine head radius
headrad = headradvec(2); % Choose initial spine head radius from list
min_headrad = 0.1e-6; % Lower limit on spine head radius
Chead = 4*pi*headrad^2*1E4*Cm*1E-6; % Spine head capacitance
glhead = gbar_leak*4*pi*headrad^2; % Head leak conductance
alpha = 2e-10; % Ratio of  headvol/gbar_ampa (assume always proportional)

% Dendrite
dend_area = 0.5E4; % um^2; gives input resistance of 200MOhm if gleak = 1pS/um^2
dend_diam = 2; % um
dend_length = dend_area/(pi*dend_diam); % um
dend_vol = pi*dend_diam^2*dend_length/4; % um^3
Cdend = dend_area*1E-8*Cm*1E-6; % dendritic capacitance
Rdend = dend_area/gbar_leak;  % Dendritic input resistance

Istim = 200e-12; % A; stimulus current to dendrite (with default parameters Istim give V=-70, and every 100pA depolarises by 10mV.

% Synapse
ampa_rise = 0.18e-3; % s; Rise time for synaptic AMPAR conductance
ampa_decay = 1.8e-3; % s; Decay time for synaptic AMPAR conductance
nmda_rise = 2e-3; % s; Rise time for synaptic NMDAR conductance
nmda_decay = 89e-3; % s; Decay time for synaptic NMDAR conductance

% set ampa and nmda g separately
%gbar_ampa = 0.6E-9; % S
gbar_nmda = 0.090E-9; % S
% set approx like Noguchi et al (2005)
gbar_ampa = (4/3)*pi*headrad^3/alpha; % S
%gbar_nmda = 0.035E-9 + gbar_ampa/8; % S; undercompensating; % S
esyn = 0; % mV; reversal potential of synapse

% Ca parameters
Ca0 = 0.05E-6; % M; resting Ca concentration
Casp0 = Ca0;
Cad0 = Ca0;
Dca = 2.2E-10; % m^2/s; Ca diffusion constant
beta_sp = 1.3E-4; % /s; spine extrusion rate
beta_neck = beta_sp; % /s; neck extrusion rate; assumed same as spine
beta_d = 4*beta_sp; % /s; dendrite extrusion rate. 4x spine rate (to reproduce Sabatini et al, 2002)
kf = 100E6; % /M/s; buffer forward (binding) rate
kb = 1*5E-6*kf; % /s; buffer backward (unbinding) rate
Bt_sp = 100E-6; % M; total buffer concentration in spine
Bt_neck = 100E-6; % M; assumed same as spine head concentration
Bt_d = 5*Bt_sp; % total buffer in dendrite (Sabatini et al, 2002)

% R-type Ca current
gbar_CaRspine = 75; % S/m^2  = pS/um^2
gbar_CaRdend = 10; % S/m^2 = pS/um^2
tau_mR = 0.0036; % s; activation time constant
tau_hR = 0.2; % s % inactivation time constant
VmR_half = -14; % mV; activation half-voltage
kmR = 6.7; % 1/mV; activation slope
VhR_half = -65; % mV; inactivation half-voltage
khR = -11.8; % 1/mV; inactivation slope
msp_R0 = 1/(1+exp(-(Vsp0-VmR_half)/kmR));
md_R0 = 1/(1+exp(-(Vd0-VmR_half)/kmR));
hsp_R0 = 1/(1+exp(-(Vsp0-VhR_half)/khR));
hd_R0 = 1/(1+exp(-(Vd0-VhR_half)/khR));

% T-type Ca current (CaV 3.1) from Traboulsie et al J Physiol 2007
gbar_CaT_spine = 75; % S/m^2  = pS/um^2
gbar_CaT_dend = 10; % S/m^2 = pS/um^2
tau_mT_0 = -0.855809*1e-3; % ms; zero activation time constant
tau_mT_A = 1.493527*1e-3; % ms; prefactor on activation function
tau_mT_k = 27.414182; % 1/mv; scale on activation exponential
tau_hT_0 = 9.987873*1e-3; % ms % inactivation time constant
tau_hT_A = 0.002883*1e-3; % ms; prefactor on inactivation function
tau_hT_k = 5.598574; % 1/mv; scale on inactivation exponential
VmT_half = -42.921064; % mV; activation half-voltage
kmT = 5.163208; % 1/mV; activation slope
VhT_half = -72.907420; % mV; inactivation half-voltage
khT = -4.575763; % 1/mV; inactivation slope
msp_T0 = 1 /(1+exp(-(Vsp0-VmT_half)/kmT));
md_T0 = 1 /(1+exp(-(Vd0-VmT_half)/kmT));
hsp_T0 = 1 /(1+exp(-(Vsp0-VhT_half)/khT));
hd_T0 = 1 /(1+exp(-(Vd0-VhT_half)/khT));

% L-type Ca current (CaV 1.2/1.3) from Carlin et al E J Neuro, 2000
gbar_CaL_spine = 75; % S/m^2  = pS/um^2
gbar_CaL_dend = 10; % S/m^2 = pS/um^2
tau_mL_0 = 5*1e-3; % ms; zero activation time constant
tau_mL_A = 20*1e-3; % ms; prefactor on activation function
tau_mLVhalf = -15; % mV; half voltage for activation tau curve
tau_mL_k = 5; % 1/mv; scale on activation exponential
tau_hL_0 = 20*1e-3; % ms % inactivation time constant
tau_hL_A = 50*1e-3; % ms; prefactor on activation function
tau_hLVhalf = -30; % mV; half voltage for inactivation tau curve
tau_hL_k = 7; % mv; scale on activation exponential
VmL_half = -20; % mV; activation half-voltage
kmL = 6; % mV; activation slope
VhL_half = -70; % mV; inactivation half-voltage
khL = -6.4; % 1/mV; inactivation slope
msp_L0 = 1 /(1+exp(-(Vsp0-VmL_half)/kmL));
md_L0 = 1 /(1+exp(-(Vd0-VmL_half)/kmL));
hsp_L0 = 1 /(1+exp(-(Vsp0-VhL_half)/khL));
hd_L0 = 1 /(1+exp(-(Vd0-VhL_half)/khL));

% N-type Ca current (CaV 2.2) from Huang et al, Neuroscience, 1998
gbar_CaN_spine = 75; % S/m^2  = pS/um^2
gbar_CaN_dend = 10; % S/m^2 = pS/um^2
alpha_m_N_A = 0.1*1e-3; % s; activation forward rate prefactor
alpha_m_N_Vhalf = 20; % mV; activation forward rate Vhalf
alpha_m_N_k = 10; % mV; activation forward rate slope
beta_m_N_A = 0.4*1e-3; % s; activation backward rate prefactor
beta_m_N_Vhalf = -25; % mV; activation backward rate threshold
beta_m_N_k = 18; % mV; activation backward rate exponent
alpha_h_N_A = 0.01*1e-3; % s; activation forward rate prefactor
alpha_h_N_Vhalf = -50; % mV; activation forward rate Vhalf
alpha_h_N_k = 10; % mV; activation forward rate slope
beta_h_N_A = 0.1*1e-3; % s; activation backward rate prefactor
beta_h_N_Vhalf = -17; % mV; activation backward rate threshold
beta_h_N_k = 17; % mV; activation backward rate exponent

% Initialize CaV channel state variables to steady-state values
alpha_m_0_sp = (alpha_m_N_A*(Vsp0-alpha_m_N_Vhalf)/(1-exp(-(Vsp0-alpha_m_N_Vhalf)/alpha_m_N_k)));
beta_m_0_sp = beta_m_N_A*exp(-(Vsp0-beta_m_N_Vhalf)/beta_m_N_k);
msp_N0 = alpha_m_0_sp/(alpha_m_0_sp+beta_m_0_sp);
alpha_h_0_sp = alpha_h_N_A*exp(-(Vsp0-alpha_h_N_Vhalf)/alpha_h_N_k);
beta_h_0_sp = beta_h_N_A/(1+exp(-(Vsp0-beta_h_N_Vhalf)/beta_h_N_k));
hsp_N0 = alpha_h_0_sp/(alpha_h_0_sp+beta_h_0_sp);
alpha_m_0_d = (alpha_m_N_A*(Vd0-alpha_m_N_Vhalf)/(1-exp(-(Vd0-alpha_m_N_Vhalf)/alpha_m_N_k)));
beta_m_0_d = beta_m_N_A*exp(-(Vd0-beta_m_N_Vhalf)/beta_m_N_k);
md_N0 = alpha_m_0_d/(alpha_m_0_d+beta_m_0_d);
alpha_h_0_d = alpha_h_N_A*exp(-(Vd0-alpha_h_N_Vhalf)/alpha_h_N_k);
beta_h_0_d = beta_h_N_A/(1+exp(-(Vd0-beta_h_N_Vhalf)/beta_h_N_k));
hd_N0 = alpha_h_0_d/(alpha_h_0_d+beta_h_0_d);

% Plasticity parameters
eta_ltp = 1e-9; % rate of potentiation
ltp_midCa = 5.5e-6; % Potentiation curve offset
ltp_slope = 0.2e-6; % Potentiation curve slope

eta_ltd = eta_ltp/2; % rate of depression
ltd_midCa = 4e-6; % Depression curve offset
ltd_slope = 0.2e-6; % Depression curve slope


%%%%%%%%%%%%%%%
% STIMULATION
%%%%%%%%%%%%%%%

delay = 50e-3; % s; Wait time at beginning of simulation

% Regular stimulation

% rate = 10; % Hz
% isi = 1/rate; % s
% spiketimes = (delay:isi:duration)'; % s
% spiketimes = [spiketimes; duration];
% spiketimes = (floor(spiketimes)/dt)*dt;

% Poisson stimulation 

% rate = 10; % Hz
% num = rate*duration/800;
% isi = -(1000/rate)*log(rand(num,1));
% spiketimes = cumsum(isi);
% valid = find(spiketimes<duration);
% spiketimes = spiketimes(1:valid(end));
% spiketimes = (floor(spiketimes)/dt)*dt;

% Bursts stimulation

IBI = 3; % s, interburst interval
ISIrate = 100; % Hz; intraburst rate
ISI = 1/ISIrate; % s, intra-burst interval
numBursts = 4;
spikesPerBurst = 100;
spiketimes = zeros(numBursts*spikesPerBurst,1);
% IBI counted between bursts
% for  i = 0:numBursts-1;
%     for j = 1:spikesPerBurst;
%         spiketimes(i*spikesPerBurst+j) = i*IBI + i*ISI*(spikesPerBurst-1) + (j-1)*ISI;
%     end
% end
%IBI counted from start of one burst to start of next
for  i = 0:numBursts-1;
    for j = 1:spikesPerBurst;
        spiketimes(i*spikesPerBurst+j) = i*IBI + (j-1)*ISI;
    end
end
spiketimes = spiketimes + delay;
spiketimes = spiketimes(spiketimes<duration);
% spiketimes = (floor(spiketimes)/dt)*dt;



%%%%%%%%%%%
% RUN
%%%%%%%%%%%

% Write parameter array for ODE solver
p=zeros(83,1);
p(1) = ampa_rise; % ms
p(2) = ampa_decay; % ms
p(3) = nmda_rise; % ms
p(4) = nmda_decay; % ms
p(5) = esyn; % mV; reversal potential of synapse
p(6) = Cm; % F
p(7) = Ra; % Ohm.cm
p(8) = gbar_leak; % pS/um^2
p(9) = e_leak; % mV
p(10) = e_CaR; % mV; Ohmic eCa is higher, but this reduces GHK error (Oertner)
p(11) = gbar_CaRspine; % pS/um^2
p(12) = gbar_CaRdend; % pS/um^2
p(13) = dend_area; % um^2

p(15) = gbar_ampa; % S
p(16) = gbar_nmda; % S
p(17) = alpha; % HeadVol/gbar_ampa ratio
p(18) = beta_sp;
p(19) = beta_neck;
p(20) = beta_d;
p(21) = Ca0;
p(22) = Dca;
p(23) = dend_vol;
p(24) = diam;
p(25) = L;
p(26) = tau_mR; % s
p(27) = tau_hR; % s
p(28) = VmR_half; % mV
p(29) = kmR; % 1/mV
p(30) = VhR_half; % mV
p(31) = khR; % 1/mV
p(32) = kf; % /M/s
p(33) = kb; % /s
p(34) = Bt_sp; % M
p(35) = Bt_neck; % M
p(36) = Bt_d; % M
p(37) = eta_ltp; % rate of potentiation
p(38) = ltp_midCa;
p(39) = ltp_slope;
p(40) = eta_ltd; % rate of depression
p(41) = ltd_midCa;
p(42) = ltd_slope;
p(43) = Istim;

p(44) = tau_mT_0;
p(45) = tau_mT_A;
p(46) = tau_mT_k;
p(47) = tau_hT_0;
p(48) = tau_hT_A;
p(49) = tau_hT_k;
p(50) = VmT_half;
p(51) = kmT;
p(52) = VhT_half;
p(53) = khT;
p(54) = gbar_CaT_spine;
p(55) = gbar_CaT_dend;

p(56) = tau_mL_0;
p(57) = tau_mL_A;
p(58) = tau_mLVhalf;
p(59) = tau_mL_k;
p(60) = tau_hL_0;
p(61) = tau_hL_A;
p(62) = tau_hLVhalf;
p(63) = tau_hL_k;
p(64) = VmL_half;
p(65) = kmL;
p(66) = VhL_half;
p(67) = khL;
p(68) = gbar_CaL_spine;
p(69) = gbar_CaL_dend;

p(70) = alpha_m_N_A;
p(71) = alpha_m_N_Vhalf; % mV; activation forward rate Vhalf
p(72) = alpha_m_N_k; % mV; activation forward rate slope
p(73) = beta_m_N_A; % s; activation backward rate prefactor
p(74) = beta_m_N_Vhalf; % mV; activation backward rate threshold
p(75) = beta_m_N_k; % mV; activation backward rate exponent
p(76) = alpha_h_N_A; % s; activation forward rate prefactor
p(77) = alpha_h_N_Vhalf; % mV; activation forward rate Vhalf
p(78) = alpha_h_N_k; % mV; activation forward rate slope
p(79) = beta_h_N_A; % s; activation backward rate prefactor
p(80) = beta_h_N_Vhalf; % mV; activation backward rate threshold
p(81) = beta_h_N_k; % mV; activation backward rate exponent
p(82) = gbar_CaN_spine; % S/m^2  = pS/um^2
p(83) = gbar_CaN_dend; % S/m^2 = pS/um^2

% Initialise
x0 = [0 0 0 0 Vsp0 Vsp0 Vd0 1.05*Ca0 1.05*Ca0 1.05*Ca0 msp_R0 hsp_R0 md_R0 hd_R0 msp_T0 hsp_T0 md_T0 hd_T0 msp_L0 hsp_L0 md_L0 hd_L0 msp_N0 hsp_N0 md_N0 hd_N0 0.99*Bt_sp 0.9833*Bt_neck 0.9833*Bt_d headrad gbar_ampa];
track_t = (0:dt:duration-dt)'.*1000;
track_x = [];

% Run (simulation loops over intervals between synaptic stimulation times)
tspan = [0:dt:spiketimes(1)-dt];
[t,x] = ode23s(@spineodes,tspan,x0,[],p);
track_x = [track_x; x];

for i = 1:length(spiketimes)-1;
    tspan = [spiketimes(i):dt:spiketimes(i+1)-dt];
    xnew = x(end,1:4)+ 0.5.*[1-x(end,1) 1-x(end,1) 1-x(end,3) 1-x(end,3)];
    x0 = [xnew x(end,5:end)];
    [t,x] = ode23s(@spineodes,tspan,x0,[],p);
    
  if(x(end,30)<min_headrad)
        x(end,30:31) = [min_headrad 4*pi*min_headrad^3/(3*alpha)];
  end
  track_x = [track_x; x];
end

% complete simualtion after last spike event
if(spiketimes(end)<duration)
    tspan = [spiketimes(end):dt:duration];
    xnew = x(end,1:4)+ 0.5.*[1-x(end,1) 1-x(end,1) 1-x(end,3) 1-x(end,3)];
    x0 = [xnew x(end,5:end)];
    [t,x] = ode23s(@spineodes,tspan,x0,[],p);
    if(x(end,30)<min_headrad)
        x(end,30:31) = [min_headrad 4*pi*min_headrad^3/(3*alpha)];
    end
    track_x = [track_x; x];
end

track_t = (0:dt:(length(track_x)-1)*dt)';


%************
% PLOT FIGURES
%************


% figure(1)
% plot(track_t,track_x(:,1),'b',track_t,track_x(:,2),'r',track_t,track_x(:,3),'b.',track_t,track_x(:,4),'r.')
% xlabel('Time (s)')
% ylabel('Conductance')

figure(2)
plot(track_t,track_x(:,5),'b',track_t,track_x(:,6),'g',track_t,track_x(:,7),'r')
xlabel('Time (s)')
ylabel('Voltage (mV)')
legend('Spine head','Spine neck','Dendrite')

figure(3)
plot(track_t,1E6*track_x(:,8),'b',track_t,1E6*track_x(:,10),'r',track_t,1E6*track_x(:,9),'g')
xlabel('Time (s)')
ylabel('[Ca^2^+] (\muM)')
legend('Spine head','Spine neck','Dendrite')
%ylim([0 0.1])

% Buffered-to-free calcium ratio (20 in spines according to Sabatini)
% figure(4)
% plot(track_t,(Bt_sp-track_x(:,27))./(track_x(:,8)),'g')
% xlabel('Time (s)')
% hold on

% figure(5)
% plot(track_t,1e18*(4/3)*pi*(track_x(:,30).^3),'b')
% xlabel('Time (s)')
% ylabel('Head volume (\mum^3)')
% 
% figure(6)
% plot(track_t,1e9*track_x(:,31),'b')
% xlabel('Time (s)')
% ylabel('AMPAR conductance (nS)')

% Plot plasticity curve
figure(7)
caconvec = 1e-6*[0:0.1:12];
plasvec = zeros(length(caconvec),1);
for i = 1:length(caconvec)
plasvec(i) = eta_ltp/(1+exp(-(caconvec(i)-ltp_midCa)/ltp_slope)) - eta_ltd/(1+exp(-(caconvec(i)-ltd_midCa)/ltd_slope));
end
plot(1e6*caconvec,plasvec)
hold on
plot(1e6*caconvec,0)
xlabel('[Ca] (µM)')
ylabel('Synaptic change')



runtime = cputime - starttime

% figure(2); hold on
figure(3); hold on
% figure(5); hold on
%figure(6); hold on