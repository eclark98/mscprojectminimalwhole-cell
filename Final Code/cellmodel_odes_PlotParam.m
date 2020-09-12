%% 
% Supplementary material 
%
% From Weisse et al., "A mechanistic link between cellular trade-offs, 
% gene expression and growth", PNAS, 2015
%
% This file initializes parameters, calls the solver routine and plots
% concentrations and transcription rates

%%
% parameters
thetar= 427;
k_cm= 0.00599;
s0= 1.0e4;
% s0= 1.0e16;
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
ns=100;
% Kgamma= 7;
Kgamma=3.0e8;
% define vectors of parameters
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
% define vectors of intial conditions
init= [rmr_0 em_0 rmq_0 rmt_0 et_0 rmm_0 mt_0 mm_0 q_0 si_0 mq_0 mr_0 r_0 a_0];

% call solver routine 
t0= 0;
tf= 1e9;
[t,y]= ode15s(@(t,y) cellmodel_odes(t, y, rates, parameters), [t0 tf], init);
% defines which column in y array is which variable
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


% %%%Plotting concentration of variables over time 
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
 grid off
legend('Intracellular Nutrients','Energy', 'Ribosomes', 'q protein', 'Transporter enzymes', 'Metabolic enzymes')
xlabel('Time (minutes)')
ylabel('Concentration (molecules)')
xlim([1e0 1e10]);
hold off


    % defined from paper
    % translation elongation rate
    gamma= gmax*a/(Kgamma + y(:,14));
    % total translation rate
	ttrate= (y(:,3) + y(:,1) + y(:,4) + y(:,6)).*gamma;
    % growth/dilation rate
	lam= ttrate./M;
    %lam= lam(:,1);
	fr= nr*(y(:,13) + y(:,1) + y(:,4) + y(:,6) + y(:,3)) / ( nr*(y(:,13) + y(:,1) + y(:,4) + y(:,6) + y(:,3)) + nx * (y(:,9) + y(:,5) + y(:,2)));
	% rate of metabloism of nutrient
    nucat= y(:,2).*vm.*y(:,10)./(Km + y(:,10));
    % rate of import of nutrient
    nuimp= (y(:,5)*vt*s0/(Kt + s0));
    
% relative transcription rates    
r_rate = ((y(:,1).*gamma))./ttrate;
q_rate = (y(:,3).*gamma)./ttrate;
e_rate = ((y(:,4) + y(:,6)).*gamma)./ttrate;
  

%%%%%%% Plots rates
% figure(2)
% loglog(t,lam,'LineWidth',2);
% hold on
% loglog(t,si,'LineWidth',2);
% hold on
% loglog(t,nucat,'LineWidth',2);
% hold on 
% loglog(t,nuimp,'LineWidth',2);
% legend('Growth Rate','Internalised Nutrient','Nutrient Metabolism Rate', 'Nutrient Import Rate')
% xlabel('Time (minutes)')
% ylabel('Rate (min^{-1})')
% hold off

%%%%%%%%% Plot of rates
figure(3)
loglog(t,lam,'LineWidth',2);
hold on
loglog(t,si,'LineWidth',2);
hold on
loglog(t,nucat,'LineWidth',2);
hold on 
loglog(t,nuimp,'LineWidth',2);
legend('Growth Rate','Internalised Nutrient','Nutrient Metabolism Rate', 'Nutrient Import Rate')
xlabel('Time (minutes)')
ylabel('Rate (min^{-1})')
hold off


% figure(4)
% loglog(a,r_rate,'LineWidth',2);
% hold on
% loglog(a,q_rate,'LineWidth',2);
% hold on
% loglog(a,e_rate,'LineWidth',2);
% legend('Ribosomes','q-proteins','Enzymes')
% xlabel('Energy (molecules)')
% ylabel('Relative Transcription Rate')
% hold off




