%
% This file initializes parameters and calls the solver routine.

% Adapted from "The impact of cellular trade-offs on synthetic genetic 
% circuits across cell populations", MSc Research Report, 2017

% Loops the solver over a range of external nutrient values to plot Monod's
% growth law, with two strains in a shared environment


%%
% parameters
thetar= 427;
k_cm= 0.00599;
s0= linspace(10,4e4,200);
%s0= 1.0e4;
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
ns= 0.5;
Kgamma= 3.0e8;
% ns=100;
% Kgamma=3.0e8;
parameters= [thetar k_cm s0 gmax cl thetax Kt M we Km vm nx Kq vt wr wq nq nr ns Kgamma];

% define rate constants
b= 0;
dm= 0.1;
kb= 1;
ku= 1.0;
f= cl*k_cm;
rates= [b dm kb ku f];

% define initial conditions
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

init= [rmr1_0 em1_0 rmq1_0 rmt1_0 et1_0 rmm1_0 mt1_0 mm1_0 q1_0 si1_0 mq1_0 mr1_0 r1_0 a1_0 rmr2_0 em2_0 rmq2_0 rmt2_0 et2_0 rmm2_0 mt2_0 mm2_0 q2_0 si2_0 mq2_0 mr2_0 r2_0 a2_0];

% call solver routine 
t0= 0;
tf= 1e12;


for i=1:length(s0)
    
    
    parameters(3) = s0(i);
    [t,y]= ode15s(@(t,y) cellmodel_odes_2_edit(t, y, rates, parameters), [t0 tf], init);
    
    
    gamma1= gmax*y(end,14)/(Kgamma + y(end,14));
    Gam1(i)= gamma1;
    % total translation rate
	ttrate1= (y(end,1) + y(end,3) + y(end,4) + y(end,6)).*gamma1;
    lam1(i)= ttrate1/M;
    si1(i)= y(end,10);
    et1(i)= y(end,5);
    nuimp1(i)= y(end,5)*vt*s0(i)/(Kt + s0(i));
    nucat1(i)= y(end,2)*vm*si1(i)/(Km + si1(i));
    % Protein concentration
    pro1(i) = y(end,2)+y(end,5)+y(end,9)+y(end,13);
    
    gamma2= gmax*y(end,28)/(Kgamma + y(end,28));
    Gam2(i)= gamma2;
    % total translation rate
	ttrate2= (y(end,15) + y(end,17) + y(end,18) + y(end,20)).*gamma2;
    lam2(i)= ttrate2/M;
    si2(i)= y(end,24);
    et2(i)= y(end,19);
    nuimp2(i)= y(end,16)*vt*s0(i)/(Kt + s0(i));
    nucat2(i)= y(end,19)*vm*si2(i)/(Km + si2(i));
    % Protein concentration 
    pro2(i) = y(end,16)+y(end,19)+y(end,23)+y(end,27);
    
end

figure(1)
plot(s0,lam1,'LineWidth',2)
hold on
plot(s0,lam2,'LineWidth',2)
hold off
xlabel('External Nutrient (molecules)')
ylabel('Growth Rate (1/min)')
legend('Strain 1','Strain 2')


figure(2)
plot(lam1, pro1,'LineWidth',2)
hold on
plot(lam2, pro2,'LineWidth',2)
xlabel('Growth rate (1/min)')
ylabel('Protein Concentration (molecules)')




% figure(2)
% plot(s0,nuimp2)
% hold on
% plot(s0,nuimp1)
% hold off
% 
% figure(3)
% plot(s0,nucat2)
% hold on
% plot(s0,nucat1)

