%% 
% Supplementary material 
%
% Adapted From Weisse et al., "A mechanistic link between cellular trade-offs, 
% gene expression and growth", PNAS, 2015
%
% This file initializes parameters and calls the solver routine.

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
ns= 0.5;
Kgamma= 7;
% Define vector of parameters
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
r_0= 10.0;
a_0= 1000.0;
% Define vector of inital conditions
init= [rmr_0 em_0 rmq_0 rmt_0 et_0 rmm_0 mt_0 mm_0 q_0 si_0 mq_0 mr_0 r_0 a_0];

% call solver routine 
t0= 0;
tf= 1e9;
[t,y]= ode15s(@(t,y) cellmodel_odes(t, y, rates, parameters), [t0 tf], init);
% Defines which column in y represents which variables corresponding to
% their ode in the solver 
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
