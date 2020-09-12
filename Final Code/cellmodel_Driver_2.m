%% 
% Supplementary material 
%
% Adapted From Weisse et al., "A mechanistic link between cellular trade-offs, 
% gene expression and growth", PNAS, 2015
%
% This file initializes parameters and calls the solver routine with two
% strains of cells in the same environment. Used to plot protein 
% concentrations and internal nutrients of both both strains 

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

% ns= 0.5;
% Kgamma= 7;

ns=100;
Kgamma=3.0e08;

% vector of parameters 
parameters= [thetar k_cm s0 gmax cl thetax Kt M we Km vm nx Kq vt wr wq nq nr ns Kgamma];

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

% Strain 1
rmr1_0= 0;
em1_0= 0;
rmq1_0= 0;
rmt1_0= 0;
et1_0= 0;
rmm1_0= 0;
mt1_0= 0;
mm1_0= 0;
q1_0= 0;
si1_0= 0;
mq1_0= 0;
mr1_0= 0;
r1_0= 10;
a1_0= 1000;

% Strain 2
rmr2_0= 0;
em2_0= 0;
rmq2_0= 0;
rmt2_0= 0;
et2_0= 0;
rmm2_0= 0;
mt2_0= 0;
mm2_0= 0;
q2_0= 0;
si2_0= 0;
mq2_0= 0;
mr2_0= 0;
r2_0= 10.0;
a2_0= 1000.0;

% Steady-state ICs
% rmr2_0= 1.7602e03;
% em2_0= 97.6248;
% rmq2_0= 648.2996;
% rmt2_0= 2.8282;
% et2_0= 97.6248;
% rmm2_0= 2.8282;
% mt2_0= 22.7986;
% mm2_0= 22.7986;
% q2_0= 2.2378e04;
% si2_0= 128.3942;
% mq2_0= 5.2260e03;
% mr2_0= 4.6482e03;
% r2_0= 0.4439;
% a2_0= 1.2752e08;
% n2_0=0;

% Define vector of intial conditions using values above
init= [rmr1_0 em1_0 rmq1_0 rmt1_0 et1_0 rmm1_0 mt1_0 mm1_0 q1_0 si1_0 mq1_0 mr1_0 r1_0 a1_0 rmr2_0 em2_0 rmq2_0 rmt2_0 et2_0 rmm2_0 mt2_0 mm2_0 q2_0 si2_0 mq2_0 mr2_0 r2_0 a2_0];

% call solver routine - if it runs "edit" version, this is the competitive cells,
% with one strain more dominant then the other 
t0= 0;
tf= 1e9;
[t,y]= ode15s(@(t,y) cellmodel_odes_2_edit(t, y, rates, parameters), [t0 tf], init);
% [t,y]= ode15s(@(t,y) cellmodel_odes_2(t, y, rates, parameters), [t0 tf], init);

% Define which columns in array y represent which variable corresponding to
% their ODE
rmr1= y(:,1);
em1= y(:,2);
rmq1= y(:,3);
rmt1= y(:,4);
et1= y(:,5);
rmm1= y(:,6);
mt1= y(:,7);
mm1= y(:,8);
q1= y(:,9);
si1= y(:,10);
mq1= y(:,11);
mr1= y(:,12);
r1= y(:,13);
a1= y(:,14);

rmr2= y(:,15);
em2= y(:,16);
rmq2= y(:,17);
rmt2= y(:,18);
et2= y(:,19);
rmm2= y(:,20);
mt2= y(:,21);
mm2= y(:,22);
q2= y(:,23);
si2= y(:,24);
mq2= y(:,25);
mr2= y(:,26);
r2= y(:,27);
a2= y(:,28);


% %%%%%% Plot Concentrations
figure(1);
% plot intracellular nutrients
grid on
loglog(t,si1,'LineWidth',1);
hold on
loglog(t,si2,'LineWidth',1);
hold on
% plot energy
loglog(t,a1,'LineWidth',1);
hold on
loglog(t,a2,'LineWidth',1);
hold on
% plot free ribosomes
loglog(t,r1,'LineWidth',1);
hold on
loglog(t,r2,'LineWidth',1);
hold on
% plot housekeeping enzymes 
loglog(t,q1,'LineWidth',1);
hold on
loglog(t,q2,'LineWidth',1);
hold on
% plot free transporters
loglog(t,et1,'LineWidth',1);
hold on
loglog(t,et2,'LineWidth',1);
hold on
% plot free metabolic enzymes
loglog(t,em1,'LineWidth',1);
hold on
loglog(t,em2,'LineWidth',1);
 grid off
legend('Intracellular Nutrients Strain 1','Intracellular Nutrients Strain 2','Energy Strain 1','Energy Strain 2', 'Ribosomes Strain 1','Ribosomes Strain 2', 'q protein Strain 1','q protein Strain 2', 'Transporter enzymes Strain 1','Transporter enzymes Strain 2', 'Metabolic enzymes Strain 1','Metabolic enzymes Strain 2')
xlabel('Time (minutes)')
ylabel('Concentration (molecules)')
xlim([1e0 1e10]);
hold off


% %%%%%% Plot intracellular nutrients 
figure(2)
plot(t,si1,'LineWidth',1.5)
hold on
plot(t,si2,'LineWidth',1.5)
xlabel('Time (minutes)')
ylabel('Internal Nutrient Concentration (molecules)')
legend('Internalised nutrient of strain 1','Internalised nutrient of strain 2')
xlim([0 1e3])