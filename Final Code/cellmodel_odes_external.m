%% 
% Supplementary material 
%
% Adapted from Weisse et al., "A mechanistic link between cellular trade-offs, 
% gene expression and growth", PNAS, 2015
%
% This file implements the right hand side of the ODE with
% external nutrient and population growth.

%%
function dydt= cellmodel_odes_external(t, y, rates, parameters)

    % location of rate constants in vector 
	b= rates(1);
	dm= rates(2);
	kb= rates(3);
	ku= rates(4);
	f= rates(5);
    ds= rates(6);
    dn= rates(7);
    kin= rates(8);

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

    % define location of variables
    % rm? means conc. of complex of ribosome and mrna species m?
    % mr,mt,mm,mq means conc. of mrna of either r,t,m or q proteins respectively
    % r,et,em,q,p means conc. of r,et,em,q or p proteins respectively
    
	rmr= y(1);
	em= y(2);
	rmq= y(3);
	rmt= y(4);
	et= y(5);
	rmm= y(6);
	mt= y(7);
	mm= y(8);
	q= y(9);
	si= y(10);
	mq= y(11);
	mr= y(12);
	r= y(13);
	a= y(14);
    
    n=y(15);
    s0=y(16);

    % Translation elongation rate
	gamma= gmax*a/(Kgamma + a);
    % Total translation rate
	ttrate= (rmq + rmr + rmt + rmm)*gamma;
    % Growth rate
	lam= ttrate/M;
	fr= nr*(r + rmr + rmt + rmm + rmq) / ( nr*(r + rmr + rmt + rmm + rmq) + nx * (q + et + em));
	% Nutrient metabolism rate
    nucat= em*vm*si/(Km + si);
    % Nutrient import rate
    nuimp= et*vt*s0/(Kt + s0);
    
    % ODEs taken from paper
	dydt(size(y,1),1)= 0;
    % rate of change of....
	dydt(1)= +kb*r*mr-ku*rmr-gamma/nr*rmr-f*rmr-lam*rmr;    % rmr
	dydt(2)= +gamma/nx*rmm-lam*em;                          % em
	dydt(3)= +kb*r*mq-ku*rmq-gamma/nx*rmq-f*rmq-lam*rmq;    % rmq
	dydt(4)= +kb*r*mt-ku*rmt-gamma/nx*rmt-f*rmt-lam*rmt;    % rmt
	dydt(5)= +gamma/nx*rmt-lam*et;                          % et
	dydt(6)= +kb*r*mm-ku*rmm-gamma/nx*rmm-f*rmm-lam*rmm;    % rmm
	dydt(7)= +(we*a/(thetax + a))+ku*rmt+gamma/nx*rmt-kb*r*mt-dm*mt-lam*mt; % mt
	dydt(8)= +(we*a/(thetax + a))+ku*rmm+gamma/nx*rmm-kb*r*mm-dm*mm-lam*mm; % mm
	dydt(9)= +gamma/nx*rmq-lam*q;                                           % q
	dydt(10)= +(et*vt*s0/(Kt + s0))-nucat-lam*(si);                         % si
	dydt(11)= +(wq*a/(thetax + a)/(1 + (q/Kq)^nq))+ku*rmq+gamma/nx*rmq-kb*r*mq-dm*mq-lam*mq; % mq
	dydt(12)= +(wr*a/(thetar + a))+ku*rmr+gamma/nr*rmr-kb*r*mr-dm*mr-lam*mr;                 % mr
	dydt(13)= +ku*rmr+ku*rmt+ku*rmm+ku*rmq+gamma/nr*rmr+gamma/nr*rmr+gamma/nx*rmt+gamma/nx*rmm+gamma/nx*rmq-kb*r*mr-kb*r*mt-kb*r*mm-kb*r*mq-lam*r; % r
	dydt(14)= +ns*nucat-ttrate-lam*a; % a
    
    % Additions for dynamic external nutrient
    dydt(15)= +lam*n-dn*n;          % n
    dydt(16)= kin-nuimp*n-ds*s0;    % s0
    
    