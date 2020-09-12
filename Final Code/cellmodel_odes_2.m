%%
% Code adapted from Weisse et al., "A mechanistic link between 
% cellular trade-offs, gene expression and growth", PNAS, 2015
% and a previous masters student code from their project 

% This file implements the right hand side of the ODE with external 
% nutrient and population growth with two strains of cells in the same 
% environment.


%%
function dydt= cellmodel_odes_2(t, y, rates, parameters)
    
    % location of rate constants in vector 
	b= rates(1);
	dm= rates(2);
	kb= rates(3);
	ku= rates(4);
	f= rates(5);

    % location of parameters in vector 
	thetar= parameters(1);
	k_cm= parameters(2);
	s0= parameters(3);
	gmax= parameters(4);
	cl= parameters(5);
	thetax= parameters(6);
	Kt= parameters(7);
	M= parameters(8);
	we= parameters(9);
	Km= parameters(10);
	vm= parameters(11);
	nx= parameters(12);
	Kq= parameters(13);
	vt= parameters(14);
	wr= parameters(15);
	wq= parameters(16);
	nq= parameters(17);
	nr= parameters(18);
	ns= parameters(19);
    Kgamma= parameters(20);

    % rm? means conc. of complex of ribosome and mrna species m?
    % mr,mt,mm,mq means conc. of mrna of either r,t,m or q proteins respectively
    % r,et,em,q,p means conc. of r,et,em,q or p proteins respectively
    
    % Strain 1 
	rmr1= y(1);
	em1= y(2);
	rmq1= y(3);
	rmt1= y(4);
	et1= y(5);
	rmm1= y(6);
	mt1= y(7);
	mm1= y(8);
	q1= y(9);
	si1= y(10);
	mq1= y(11);
	mr1= y(12);
	r1= y(13);
	a1= y(14);
    
    % Strain 2
    rmr2= y(15);
	em2= y(16);
	rmq2= y(17);
	rmt2= y(18);
	et2= y(19);
	rmm2= y(20);
	mt2= y(21);
	mm2= y(22);
	q2= y(23);
	si2= y(24);
	mq2= y(25);
	mr2= y(26);
	r2= y(27);
	a2= y(28);

    
	%Kgamma= gmax/Kp;
    %Kgamma= 3e8;
    
    
    %%%%%% Strain 1
    % Translation elongation rate
	gamma1= gmax*a1/(Kgamma + a1);
    % Total translation rate = nx*(gamma/nx*rmx) = gamma*rmx
	ttrate1= (rmq1 + rmr1 + rmt1 + rmm1)*gamma1;
    % Dilution rate = growth rate
	lam1= ttrate1/M;
	fr1= nr*(r1 + rmr1 + rmt1 + rmm1 + rmq1) / ( nr*(r1 + rmr1 + rmt1 + rmm1 + rmq1) + nx * (q1 + et1 + em1));
	% Rate of import of nutrients
    nucat1= em1*vm*si1/(Km + si1);
    
    
    %%%%%% Strain 2
    % Translation elongation rate
    gamma2= gmax*a2/(Kgamma + a2);
    % Total translation rate = nx*(gamma/nx*rmx) = gamma*rmx
	ttrate2= (rmq2 + rmr2 + rmt2 + rmm2)*gamma2;
    % Dilution rate = growth rate
	lam2= ttrate2/M;
	fr2= nr*(r2 + rmr2 + rmt2 + rmm2 + rmq2) / ( nr*(r2 + rmr2 + rmt2 + rmm2 + rmq2) + nx * (q2 + et2 + em2));
	% Rate of import of nutrients
    nucat2= em2*vm*si2/(Km + si2);

    
	dydt(size(y,1),1)= 0;
    % rate of change of... strain 1
	dydt(1)= +kb*r1*mr1-ku*rmr1-gamma1/nr*rmr1-f*rmr1-lam1*rmr1;    % rmr1
	dydt(2)= +gamma1/nx*rmm1-lam1*em1;                              % em1
	dydt(3)= +kb*r1*mq1-ku*rmq1-gamma1/nx*rmq1-f*rmq1-lam1*rmq1;    % rmq1
	dydt(4)= +kb*r1*mt1-ku*rmt1-gamma1/nx*rmt1-f*rmt1-lam1*rmt1;    % rmt1
	dydt(5)= +gamma1/nx*rmt1-lam1*et1;                              % et1
	dydt(6)= +kb*r1*mm1-ku*rmm1-gamma1/nx*rmm1-f*rmm1-lam1*rmm1;    % rmm1
	dydt(7)= +(we*a1/(thetax + a1))+ku*rmt1+gamma1/nx*rmt1-kb*r1*mt1-dm*mt1-lam1*mt1; % mt1
	dydt(8)= +(we*a1/(thetax + a1))+ku*rmm1+gamma1/nx*rmm1-kb*r1*mm1-dm*mm1-lam1*mm1; % mm1
	dydt(9)= +gamma1/nx*rmq1-lam1*q1;                                                 % q1
	dydt(10)= +(et1*vt*s0/(Kt + s0))-nucat1-lam1*si1;                                 % si1
	dydt(11)= +(wq*a1/(thetax + a1)/(1 + (q1/Kq)^nq))+ku*rmq1+gamma1/nx*rmq1-kb*r1*mq1-dm*mq1-lam1*mq1; % mq1
	dydt(12)= +(wr*a1/(thetar + a1))+ku*rmr1+gamma1/nr*rmr1-kb*r1*mr1-dm*mr1-lam1*mr1;                  % mr1
	dydt(13)= +ku*rmr1+ku*rmt1+ku*rmm1+ku*rmq1+gamma1/nr*rmr1+gamma1/nr*rmr1+gamma1/nx*rmt1+gamma1/nx*rmm1+gamma1/nx*rmq1-kb*r1*mr1-kb*r1*mt1-kb*r1*mm1-kb*r1*mq1-lam1*r1; % r1
	dydt(14)= +ns*nucat1-ttrate1-lam1*a1; % a1
    
    % rate of change of... strain 2
	dydt(15)= +kb*r2*mr2-ku*rmr2-gamma2/nr*rmr2-f*rmr2-lam2*rmr2;   % rmr2
	dydt(16)= +gamma2/nx*rmm2-lam2*em2;                             % em2
	dydt(17)= +kb*r2*mq2-ku*rmq2-gamma2/nx*rmq2-f*rmq2-lam2*rmq2;   % rmq2
	dydt(18)= +kb*r2*mt2-ku*rmt2-gamma2/nx*rmt2-f*rmt2-lam2*rmt2;   % rmt2
	dydt(19)= +gamma2/nx*rmt2-lam2*et2;                             % et2
	dydt(20)= +kb*r2*mm2-ku*rmm2-gamma2/nx*rmm2-f*rmm2-lam2*rmm2;   % rmm2
	dydt(21)= +(we*a2/(thetax + a2))+ku*rmt2+gamma2/nx*rmt2-kb*r2*mt2-dm*mt2-lam2*mt2; % mt2
	dydt(22)= +(we*a2/(thetax + a2))+ku*rmm2+gamma2/nx*rmm2-kb*r2*mm2-dm*mm2-lam2*mm2; % mm2
	dydt(23)= +gamma2/nx*rmq2-lam2*q2;                                                 % q2
	dydt(24)= +(et2*vt*s0/(Kt + s0))-nucat2-lam2*si2;                                  % si2
	dydt(25)= +(wq*a2/(thetax + a2)/(1 + (q2/Kq)^nq))+ku*rmq2+gamma2/nx*rmq2-kb*r2*mq2-dm*mq2-lam2*mq2; % mq2
	dydt(26)= +(wr*a2/(thetar + a2))+ku*rmr2+gamma2/nr*rmr2-kb*r2*mr2-dm*mr2-lam2*mr2;                  % mr2
	dydt(27)= +ku*rmr2+ku*rmt2+ku*rmm2+ku*rmq2+gamma2/nr*rmr2+gamma2/nr*rmr2+gamma2/nx*rmt2+gamma2/nx*rmm2+gamma2/nx*rmq2-kb*r2*mr2-kb*r2*mt2-kb*r2*mm2-kb*r2*mq2-lam2*r2; % r2
	dydt(28)= +ns*nucat2-ttrate2-lam2*a2; % a2