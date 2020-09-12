%% 
% Supplementary material 
%
% Extended From Weisse et al., "A mechanistic link between cellular trade-offs, 
% gene expression and growth", PNAS, 2015
%
% This file initializes parameters and calls the solver routine with the
% inclusion of dynamic external nutrient and population growth. Plots 
% concentrations to steady-state
%%
% parameters
thetar= 427;
k_cm= 0.00599;
s0= 1.0e4;
gmax= 1260;
cl= 0;
thetax= 4.38;
Kt= 1000;
M= 1.0e8;
we= 4.14;
Km= 1000;
vm= 5800;
nx= 300;
Kq= 152219;
vt= 726;
wr= 930;
wq= 949;
nq= 300;
nr= 7549;

ns= 100;
Kgamma= 3e8;

% Alternative parameters
% ns= 0.5;
% Kgamma= 7;

% Define vector of parameters 
parameters= [thetar k_cm s0 gmax cl thetax Kt M we Km vm nx Kq vt wr wq nq nr ns Kgamma];

% define rate constants
b= 0;
dm= 0.1;
kb= 1;
ku= 1.0;
f= cl*k_cm;
ds= 0.01;
dn= 0.01;
kin= 0;
rates= [b dm kb ku f ds dn kin];

% define initial conditions
% rm? means conc. of complex of ribosome and mrna species m?
% mr,mt,mm,mq means conc. of mrna of either r,t,m or q proteins respectively
% r,et,em,q,p means conc. of r,et,em,q or p proteins respectively

% %%% Steady-state ICs
% rmr_0= 6.056e+03;
% %em_0= 1;
% em_0= 466.7;
% rmq_0= 2.8745e+03;
% rmt_0= 12.5401;
% et_0= 466.7;
% %et_0= 1;
% rmm_0= 12.5401;
% mt_0= 20.4599;
% mm_0= 20.4599;
% q_0= 1.0699e+05;
% %si_0= 128.4008;
% si_0= 0;
% mq_0= 4.68996e+03;
% mr_0= 3.339e+03;
% r_0= 2.0994;
% a_0= 3.8545e+08;

% %zero ICs
rmr_0= 0;
em_0= 0;
rmq_0= 0;
rmt_0= 0;
et_0= 0;
rmm_0= 0;
mt_0= 0;
mm_0= 0;
q_0= 0;
si_0= 0;
mq_0= 0;
mr_0= 0;
r_0= 10;
a_0= 1000;

% Additional variables, number of cells and external nutrient 
n_0= 0;
s0_0=1.0e10;


init= [rmr_0 em_0 rmq_0 rmt_0 et_0 rmm_0 mt_0 mm_0 q_0 si_0 mq_0 mr_0 r_0 a_0 n_0 s0_0];

% call solver routine - runs odes that include odes for external nutrient
% and number of cells
t0= 0;
tf= 1e4;
[t,y]= ode15s(@(t,y) cellmodel_odes_external(t, y, rates, parameters), [t0 tf], init);

rmr= y(:,1);
em= y(:,2);
rmq= y(:,3);
rmt= y(:,4);
et= y(:,5);
rmm= y(:,6);
mt= y(:,7);
mm= y(:,8);
q= y(:,9);
si= y(:,10);
mq= y(:,11);
mr= y(:,12);
r= y(:,13);
a= y(:,14);

n= y(:,15);
s0= y(:,16);


% %%%% Plot concentrations to steady-state
figure(1);
% plot intracellular nutrients
grid on
loglog(t,si,'LineWidth',2);
hold on
% plot energy
loglog(t,a,'LineWidth',2);
hold on
% plot free ribosomes
loglog(t,r,'LineWidth',2);
hold on
% plot housekeeping enzymes 
loglog(t,q,'LineWidth',2);
hold on
% plot free transporters
loglog(t,et,'LineWidth',2);
hold on
% plot free metabolic enzymes
loglog(t,em,'LineWidth',2);
% plot free metabolic enzymes
loglog(t,s0,'LineWidth',2);
 grid off
legend('Intracellular Nutrients','Energy', 'Ribosomes', 'q protein', 'Transporter enzymes', 'Metabolic enzymes', 'External Nutrient')
xlabel('Time (minutes)')
ylabel('Concentration (molecules)')
xlim([1e0 1e10]);
hold off
