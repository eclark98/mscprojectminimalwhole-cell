%% 
% Supplementary material 
%
% Adapted From Weisse et al., "A mechanistic link between cellular trade-offs, 
% gene expression and growth", PNAS, 2015
%
% This file initializes parameters and calls the solver routine with
% external nutrient and population growth, with two competitive strains and
% creates several plots

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
Kgamma=3.0e8;
% define vectors of parameters
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
% define vectors of intial conditions
rates= [b dm kb ku f ds dn kin];

% define initial conditions - different sets dependent on what to plot
% rm? means conc. of complex of ribosome and mrna species m?
% mr,mt,mm,mq means conc. of mrna of either r,t,m or q proteins respectively
% r,et,em,q,p means conc. of r,et,em,q or p proteins respectively


%%%%% SET 1 - zeros inital conditions
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
r1_0= 10.0;
a1_0= 1000.0;
n1_0=0;
% %%%Strain 2
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
n2_0=0;

%%%%% SET 2 - same steady-state
% %%%STRAIN 1
% rmr1_0= 6.056e+03;
% em1_0= 466.7;
% rmq1_0= 2.8745e+03;
% rmt1_0= 12.5401;
% et1_0= 466.7;
% rmm1_0= 12.5401;
% mt1_0= 20.4599;
% mm1_0= 20.4599;
% q1_0= 1.0699e+05;
% si1_0= 0;
% mq1_0= 4.68996e+03;
% mr1_0= 3.339e+03;
% r1_0= 2.0994;
% a1_0= 3.8545e+08;
% n1_0=1;
% %%% STRAIN 2
% rmr2_0= 6.056e+03;
% em2_0= 466.7;
% rmq2_0= 2.8745e+03;
% rmt2_0= 12.5401;
% et2_0= 466.7;
% rmm2_0= 12.5401;
% mt2_0= 20.4599;
% mm2_0= 20.4599;
% q2_0= 1.0699e+05;
% si2_0= 0;
% mq2_0= 4.68996e+03;
% mr2_0= 3.339e+03;
% r2_0= 2.0994;
% a2_0= 3.8545e+08;
% n2_0=1;


% %%%%% SET 3 - same steady-state
% %Strain 1
% rmr1_0= 5.7706e03;
% em1_0= 489.8298;
% rmq1_0= 2.9442e03;
% rmt1_0= 12.88440;
% et1_0= 489.8298;
% rmm1_0= 12.8440;
% mt1_0= 12.9705;
% mm1_0= 12.9705;
% q1_0= 1.1128e05;
% si1_0= 143.0744;
% mq1_0= 2.9732e03;
% mr1_0= 1.4014e03;
% r1_0= 5.2583;
% a1_0= 2.2287e08;
% n1_0=10;
% %Strain 2
% rmr2_0= 5.7706e03;
% em2_0= 489.8298;
% rmq2_0= 2.9442e03;
% rmt2_0= 12.88440;
% et2_0= 489.8298;
% rmm2_0= 12.8440;
% mt2_0= 12.9705;
% mm2_0= 12.9705;
% q2_0= 1.1128e05;
% si2_0= 143.0744;
% mq2_0= 2.9732e03;
% mr2_0= 1.4014e03;
% r2_0= 5.2583;
% a2_0= 2.2287e08;
% n2_0=10;

% %%%% SET 4 - competitive strains steady-state
% % Strain 1
% rmr1_0= 6.0562e03;
% em1_0= 466.7413;
% rmq1_0= 2.8745e03;
% rmt1_0= 12.5401;
% et1_0= 466.7413;
% rmm1_0= 12.5401;
% mt1_0= 20.4599;
% mm1_0= 20.4599;
% q1_0= 1.0699e05;
% si1_0= 128.4008;
% mq1_0= 4.6900e03;
% mr1_0= 3.3386e03;
% r1_0= 2.0994;
% a1_0= 3.8545e08;
% n1_0=1.4753e39;
% %STRAIN 2
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

s0_0=10e10;

% Define inital condition array
init= [rmr1_0 em1_0 rmq1_0 rmt1_0 et1_0 rmm1_0 mt1_0 mm1_0 q1_0 si1_0 mq1_0 mr1_0 r1_0 a1_0 n1_0 rmr2_0 em2_0 rmq2_0 rmt2_0 et2_0 rmm2_0 mt2_0 mm2_0 q2_0 si2_0 mq2_0 mr2_0 r2_0 a2_0 n2_0 s0_0];

% call solver routine 
t0= 0;
tf= 5e4;
% call solver routine - if it runs "edit" version, this is the competitive cells,
% with one strain more dominant then the other 
[t,y]= ode15s(@(t,y) cellmodel_odes_external_2_edit(t, y, rates, parameters), [t0 tf], init);

% defines which column in y array is which variable

% Strain 1
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
n1= y(:,15);

% Strain 2
rmr2= y(:,16);
em2= y(:,17);
rmq2= y(:,18);
rmt2= y(:,19);
et2= y(:,20);
rmm2= y(:,21);
mt2= y(:,22);
mm2= y(:,23);
q2= y(:,24);
si2= y(:,25);
mq2= y(:,26);
mr2= y(:,27);
r2= y(:,28);
a2= y(:,29);
n2= y(:,30);

s0= y(:,31);

    % %%%%%%%%% STRAIN 1
    % Translation elongation rate
    gamma1= gmax.*a1./(Kgamma + a1);
    % Total translation rate = nx*(gamma/nx*rmx) = gamma*rmx
	ttrate1= (rmq1 + rmr1 + rmt1 + rmm1).*gamma1;
    % Dilution rate = growth rate
	lam1= ttrate1/M;
	fr1= nr*(r1 + rmr1 + rmt1 + rmm1 + rmq1) / ( nr*(r1 + rmr1 + rmt1 + rmm1 + rmq1) + nx * (q1 + et1 + em1));
	% Rate of import of nutrients
    nucat1= em1.*vm.*si1./(Km + si1);
    % Rate of metabolism into energy
    nuimp1= et1.*vt.*s0./(Kt + s0);
    
    
    % %%%%%%%% STRAIN 2 - multiplied by 2
    % Translation elongation rate
    gamma2= 2*gmax.*a2./(Kgamma + a2);
    % Total translation rate = nx*(gamma/nx*rmx) = gamma*rmx
	ttrate2= 2*(rmq2 + rmr2 + rmt2 + rmm2).*gamma2;
    % Dilution rate = growth rate
	lam2= 2*ttrate2/M;
	fr2= 2*nr*(r2 + rmr2 + rmt2 + rmm2 + rmq2) / ( nr*(r2 + rmr2 + rmt2 + rmm2 + rmq2) + nx * (q2 + et2 + em2));
	% Rate of import of nutrients
    nucat2= 2*em2.*vm.*si2./(Km + si2);
    % Rate of metabolism into energy
    nuimp2= 2*et2.*vt.*s0./(Kt + s0);
    

% Plot of growth rate of strain 1 against time
figure(1)
plot(t,lam1,'LineWidth',2)  
xlabel('Time (minutes)')
ylabel('Growth Rate')
title('Growth Rate of Strain 1')
xlim([0 5e4]);

% Plot of growth rate of strain 2 against time
figure(2)
plot(t,lam2,'LineWidth',2)
xlabel('Time (minutes)')
ylabel('Growth Rate')
title('Growth Rate of Strain 2')
xlim([0 5e4]);

% Plots of both growth rates against time
figure(3)
plot(t,lam1,'LineWidth',2)
hold on
plot(t,lam2,'LineWidth',2)
legend('Growth rate of strain 1', 'Growth rate of strain 2')
xlabel('Time (minutes)')
ylabel('Growth Rate')
xlim([0 3000]);

% Plot of external and internalised nutrient
figure(4)
yyaxis left
plot(t,s0,'LineWidth',1.5)
xlabel('Time (minutes)')
ylabel('External Nutrient Concentration (molecules)')
yyaxis right
plot(t,si1,'LineWidth',1.5)
hold on
plot(t,si2,'LineWidth',1.5)
ylabel('Internal Nutrient Concentration (molecules)')
legend('External nutrient','Internalised nutrient of strain 1','Internalised nutrient of strain 2')
%xlim([0 200])
%xlim([0 0.002753728490607])


%%%%%%%%% Other plots
% figure(5)
% plot(t,si1,'LineWidth',1.5)
% hold on
% plot(t,si2,'LineWidth',1.5)
% xlabel('Time (minutes)')
% ylabel('Internal Nutrient Concentration (molecules)')
% legend('Internalised nutrient of strain 1','Internalised nutrient of strain 2')
% %xlim([0 20])

% figure(6)
% loglog(t,lam1,'LineWidth',1.5)
% xlim([10e-2 10e2]);
% 
% figure(7)
% loglog(t,lam2,'LineWidth',1.5)
% xlim([10e-2 10e2]);
% 
% figure(8)
% plot(t,n1,'LineWidth',1.5)
% hold on
% plot(t,n2,'LineWidth',1.5)
% xlabel('Time (minutes)')
% ylabel('Number of Cells')
% 
% figure(9)
% loglog(t,si1,'LineWidth',1.5)
% hold on
% loglog(t,si2,'LineWidth',1.5)
% xlabel('Time (minutes)')
% ylabel('Internal Nutrient Concentration (molecules)')
% legend('Internalised nutrient of strain 1','Internalised nutrient of strain 2')
% %xlim([0 20])