%%%
% Adapted from "The impact of cellular trade-offs on synthetic genetic 
% circuits across cell populations", MSc Research Report, 2017

% Loops the solver over a range of external nutrient values to plot Monod's
% growth law


clear all

%%
% parameters
thetar= 427;
k_cm= 0.00599;
s0= linspace(10,10e4,400);     % Creates array of external nutrient values
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

% Alternative parameters
ns=100;
Kgamma=3.0e08;

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
% Defines vector of intial conditions
init= [rmr_0 em_0 rmq_0 rmt_0 et_0 rmm_0 mt_0 mm_0 q_0 si_0 mq_0 mr_0 r_0 a_0];


% call solver routine 
t0= 0;
tf= 1e12;

for i=1:length(s0)      %loops over range of external nutrient values
    
    %clear y
    %clear t
    
    parameters(3) = s0(i);
    [t,y]= ode15s(@(t,y) cellmodel_odes(t, y, rates, parameters), [t0 tf], init);

    %Kgamma= gmax/Kp;
    gamma= gmax*y(end,14)/(Kgamma + y(end,14));
    Gam(i)= gamma;
    % total translation rate
	ttrate= (y(end,1) + y(end,3) + y(end,4) + y(end,6))*gamma;
    % Growth rate
    lam(i)= ttrate/M;
    si(i)= y(end,10);
    et(i)= y(end,5);
    % nutrient import rate
    nuimp(i)= y(end,5)*vt*s0(i)/(Kt + s0(i));
    % Nutrient metabolism rate
    nucat(i)= y(end,2)*vm*si(i)/(Km + si(i));
    
end    

figure(1)       % Plot of external nutrient against growth rate
plot(s0,lam)

figure(2)       % Plot of external nutrient against nutrient import rate
plot(s0,nuimp)

figure(3)       % Plot of external nutrient against metabolism of nutrient
plot(s0,nucat)

