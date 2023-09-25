function sigma = exc_7774_O2(E)
% EXC_7774_O2 - Xs = exc_7774_O2(E)
% 
% Electron-impact excitation of OI 3s 5S - 3p 5P 7774A transition
% in the dissociative ionization/excitation of O2
% Erdman and Zipf, J. Chem Phys. 87, 1987, p.4540
% 
% 	typed in by NI August 2006
% parent       O2
% products     0    (emission)
% threshold    15.9   (???) Seems OK, 10.74 + 5.15 (O2-bond-energy)
% units        eV
% valid from   0.0   Remove this eventually - or replace with other comments ...
% valid to     0.0
% points       56
EnX =  [16	0.044
	17	0.35
	18	0.618
	19	0.888
	20	1.10
	22	1.52
	24	1.75
	26	1.91
	28	2.05
	30	2.16
	35	2.40
	40	2.62
	50	3.19
	60	3.83
	70	3.99
	80	4.07
	90	4.31
	100	4.23
	125	3.87
	150	3.42
	175	3.06
	200	2.76
	225	2.56
	250	2.34
	275	2.20
	300	2.04
	325	1.94
	350	1.82
	375	1.72
	400	1.65
	425	1.58
	450	1.50
	475	1.44
	500	1.39
	600	1.20
	700	1.06
	800	0.954
	900	0.868
	1000	0.797
	1200	0.687
	1400	0.605
	1600	0.542
	1800	0.492
	2000	0.451
	2200	0.416
	2400	0.387
	2600	0.362
	2800	0.340
	3000	0.321
	4000	0.251
	5000	0.208
	6000	0.178
	7000	0.156
	8000	0.139
	9000	0.125
	10000	0.114];


sigma = exp(interp1(log(EnX(:,1)),log(EnX(:,2)),log(E),'pchip'));
sigma(E<16.14) = 0;
sigma = sigma*1e-22;
