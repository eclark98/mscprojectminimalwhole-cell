%%%
% 
% Code adapted from Weisse et al., "A mechanistic link between 
% cellular trade-offs, gene expression and growth", PNAS, 2015
% and a previous masters student code from their project 
% "The impact of cellular trade-offs on synthetic genetic 
% circuits across cell populations", MSc Research Report, 2017
%
%
%

% Loops the solver over a range of induction level wg values to run 
% repressillator model


clear all

wg= logspace(0,4,100); % create vector of induction level values

for i=1:length(wg) % runs loop of each induction level value
%%
% parameters
thetar= 426.8693338968694;
k_cm= 0.005990373118888;
s0= 10000;
gmax= 1260.0;
cl= 0;
thetax= 4.379733394834643;
Kt= 1.0e3;
M= 1.0e8;
we= 4.139172187824451;
Km= 1.0e3;
vm= 5800.0;
nx= 300.0;
Kq= 1.522190403737490e+05;
Kp= 180.1378030928276;
vt= 726.0;
wr= 929.9678874564831;
wq= 948.9349882947897;
wp= 0.0;
nq= 4;
nr= 7549.0;
ns= 0.5;
% Extra repressilator parameters
dg= log(2)/4;
dmg= log(2)/2;
Kg= 100;
hg= 2;
wg= 614;

% Define vector of parameters
parameters= [thetar k_cm s0 gmax cl thetax Kt M we Km vm nx Kq Kp vt wr wq wp nq nr ns dg dmg Kg hg 0 0];

% define rate constants
b= 0;
dm= 0.1;
kb= 1;
ku= 1.0;
f= cl*k_cm;
rates= [b dm kb ku f];

% define initial conditions
% rm? means conc. of complex of ribosome and mrna species m?
% mr,mt,mm,mq means conc. of mrna of either r,t,m or q proteins respectively
% r,et,em,q,p means conc. of r,et,em,q or p proteins respectively

rmr_0= 0;
em_0= 0;
rmp_0= 0;
rmq_0= 0;
rmt_0= 0;
et_0= 0;
rmm_0= 0;
zmm_0= 0;
zmr_0= 0;
zmp_0= 0;
zmq_0= 0;
zmt_0= 0;
mt_0= 0;
mm_0= 0;
q_0= 0;
p_0= 0;
si_0= 0;
mq_0= 0;
mp_0= 0;
mr_0= 0;
r_0= 10.0;
a_0= 1000.0;
% Define Vector for inital conditions
init= [rmr_0 em_0 rmp_0 rmq_0 rmt_0 et_0 rmm_0 zmm_0 zmr_0 zmp_0 zmq_0 zmt_0 mt_0 mm_0 q_0 p_0 si_0 mq_0 mp_0 mr_0 r_0 a_0 0 0 0 0 0 0 0 0 0 0 0 0];


 % call solver routine from script of odes
tf= 1e7;
[t,y]= ode15s(@(t,y) cellmodel_repress_odes(t, y, rates, parameters), [0 tf], init);

% Defining which values in array y are assigned to which variable 
rmr= y(:,1); em= y(:,2); rmp= y(:,3); rmq= y(:,4); rmt= y(:,5);
et= y(:,6); rmm= y(:,7); zmm= y(:,8); zmr= y(:,9); zmp= y(:,10); zmq= y(:,11); zmt= y(:,12); mt= y(:,13);
mm= y(:,14); q= y(:,15); p= y(:,16); si= y(:,17); mq= y(:,18); mp= y(:,19); mr= y(:,20); r= y(:,21); a= y(:,22);

% define initial conditions for GFP simulation as ss of the native model
rmr_0= rmr(end); em_0= em(end); rmp_0= rmp(end); rmq_0= rmq(end); rmt_0= rmt(end);et_0= et(end); rmm_0= rmm(end); zmm_0= zmm(end); zmr_0= zmr(end);
zmp_0= zmp(end); zmq_0= zmq(end); zmt_0= zmt(end); mt_0= mt(end); mm_0= mm(end); q_0= q(end); p_0= p(end); si_0= si(end); mq_0= mq(end); mp_0= mp(end);
mr_0= mr(end); r_0= r(end); a_0= a(end);

% randomize initial conds of the GFP species, mean equal to the
% Em protein from native model
meanGFP     = em(end);
meanmGFP    = mm(end);
meanrmGFP   = rmm(end);

mg_0    = meanmGFP + 0.3*meanmGFP.*randn(1,1);
rmg_0   = meanrmGFP + 0.3*meanrmGFP.*randn(1,1);
g_0     = meanGFP + 0.3*meanGFP.*randn(1,1);

init= [rmr_0 em_0 rmp_0 rmq_0 rmt_0 et_0 rmm_0 zmm_0 zmr_0 zmp_0 zmq_0 zmt_0 mt_0 mm_0 q_0 p_0 si_0 mq_0 mp_0 mr_0 r_0 a_0 0 0 0 0 0 0 0 0 0 rmg_0 mg_0 g_0];

% set numg0 (GFP gene induction strength)
%parameters(end)= wg(i);

[t,y]= ode15s(@(t,y) cellmodel_repress_odes(t, y, rates, parameters), [0 tf], init);

% Defining which values in array y are assigned to which variable 
rmr= y(:,1); em= y(:,2); rmp= y(:,3); rmq= y(:,4); rmt= y(:,5);
et= y(:,6); rmm= y(:,7); zmm= y(:,8); zmr= y(:,9); zmp= y(:,10); zmq= y(:,11); zmt= y(:,12); mt= y(:,13);
mm= y(:,14); q= y(:,15); p= y(:,16); si= y(:,17); mq= y(:,18); mp= y(:,19); mr= y(:,20); r= y(:,21); a= y(:,22); rmg= y(:,32); mg= y(:,33); g= y(:,34);

% define initial conditions for GFP simulation as ss of the native model
rmr_0= rmr(end); em_0= em(end); rmp_0= rmp(end); rmq_0= rmq(end); rmt_0= rmt(end);et_0= et(end); rmm_0= rmm(end); zmm_0= zmm(end); zmr_0= zmr(end);
zmp_0= zmp(end); zmq_0= zmq(end); zmt_0= zmt(end); mt_0= mt(end); mm_0= mm(end); q_0= q(end); p_0= p(end); si_0= si(end); mq_0= mq(end); mp_0= mp(end);
mr_0= mr(end); r_0= r(end); a_0= a(end);

% randomize initial conds of the repressilator species, mean equal to the
% gfp steady state from gfp model
meanGFP     = g(end);
meanmGFP    = mg(end);
meanrmGFP   = rmg(end);

%randomize init conds
initmg    = meanmGFP + 0.3*meanmGFP.*randn(1,3);
initrmg   = meanrmGFP + 0.3*meanrmGFP.*randn(1,3);
initg     = meanGFP + 0.3*meanGFP.*randn(1,3);
mg1_0= initmg(1);
mg2_0= initmg(2);
mg3_0= initmg(3);
rmg1_0= initrmg(1);
rmg2_0= initrmg(2);
rmg3_0= initrmg(3);
g1_0= initg(1);
g2_0= initg(2);
g3_0= initg(3);

init= [rmr_0 em_0 rmp_0 rmq_0 rmt_0 et_0 rmm_0 zmm_0 zmr_0 zmp_0 zmq_0 zmt_0 mt_0 mm_0 q_0 p_0 si_0 mq_0 mp_0 mr_0 r_0 a_0 g1_0 g2_0 g3_0 rmg1_0 rmg2_0 rmg3_0 mg1_0 mg2_0 mg3_0 0 0 0];

parameters(end)= 0;
%parameters(end-1)= wg(i);

tf= 5e3;
[t,y]= ode15s(@(t,y) cellmodel_repress_odes(t, y, rates, parameters), [0 tf], init);

    % Equations from paper
    Kgamma= gmax/Kp;
    gamma= gmax*y(end,22)/(Kgamma + y(end,22));%Translation elongation rate
    Gam(i)= gamma;
    % total translation rate
	ttrate= (y(end,1) + y(end,3) + y(end,4) + y(end,5) + y(end,7) + y(end,26) + y(end,27) + y(end,28))*gamma;
    lam(i)= ttrate/M;       % Growth/dilation rate
    
    % translation rates of different protein classes
    r_rate(i)= y(end,1)*gamma;
    q_rate(i)= y(end,4)*gamma;
    e_rate(i)= (y(end,5)+y(end,7))*gamma;
    g_rate(i)= (y(end,26) + y(end,27) + y(end,28))*gamma;
    
    % Find peak oscillations 
    [pks,loc]= findpeaks(y(:,23));
    amp(i)= pks(end);
    period(i)= t(loc(end))-t(loc(end-1));

end

% Plot of translation rates of different protein classes against induction
% figure(1)
% semilogx(wg,r_rate,'LineWidth',2)
% hold on
% semilogx(wg,q_rate,'LineWidth',2)
% semilogx(wg,e_rate,'LineWidth',2)
% semilogx(wg,g_rate,'LineWidth',2)
% legend('ribosome rate','housekeeping rate','enzyme rate','repress rate')
% xlabel('repressilator induction')
% ylabel('synthesis rate (aa/min)')
% 
% Plot of induction against growth rate
% figure(2)
% semilogx(wg,lam, 'LineWidth',2)
% xlabel('repressilator induction')
% ylabel('growth rate (1/min)')
% title('Growth rate with Repressilator')


% Plot of repressilator induction with am-plitude and period of oscillation
figure(3)
yyaxis left
semilogx(wg,amp,'LineWidth',2)
xlabel('repressilator induction')
ylabel('repressilator protein amplitude (mols/cell)')
yyaxis right
semilogx(wg,period,'LineWidth',2)
ylim([40,110])
ylabel('period (min)')


% Plot of concentrations
figure(4);
% plot intracellular nutrients
semilogy(t(end),si,'LineWidth',2);
hold on
% plot energy
semilogy(t(end),a,'LineWidth',2);
hold on
% plot free ribosomes
semilogy(t(end),r,'LineWidth',2);
hold on
% plot housekeeping enzymes 
semilogy(t(end),q,'LineWidth',2);
hold on
% plot free transporters
semilogy(t(end),et,'LineWidth',2);
hold on
% plot free metabolic enzymes
semilogy(t(end),em,'LineWidth',2);
 grid off
legend('Intracellular Nutrients','Energy', 'Ribosomes', 'q', 'transporter enzymes', 'metabolic enzymes')
xlabel('Time (minutes)')
ylabel('Concentration (M)')
%xlim([1e0 1e10]);
hold off

