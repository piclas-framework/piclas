%-------------------------------------------------------------------------%
% calculate equilibirum constant by means of statistical mechanic
% is only for Ionization of Ar
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
PI      = 3.14159265358979323846;  % -
kappaB  = 1.380650524e-23;         % J / K
NL      = 2.6867805e25;            % m^-3
PlanckC = 6.62606957e-34;          % J s
q     = 1.602176565E-19; % J / eV
%-------------------------------------------------------------------------%
% electron
mass1        = 9.1093826E-31;     % kg
% at 300 k
dref1        = 2.8179E-15;      % m
omega1       = 0.7;          % -
Tref1        = 300;          % Kelvin
% at 1000 K
%dref1        = 5.64E-15;      % m
%omega1       = 0.60;          % -
%Tref1        = 1000;          % Kelvin
nel(1)       = 1;          % -
gel(1,1)     = 2;         % -
elcharT(1,1) = 0;  % Kelvin
% Argon
mass2        = 6.631368E-26;     % kg
dref2        = 3.0E-10;      % m
omega2       = 0.7;          % -
Tref2        = 300;          % Kelvin
data=importdata('Ar.txt');
nel(2)       = size(data,1) ;    % -
for ii = 1:nel(2)
    gel(2,ii)     = data(ii,1);
    elcharT(2,ii) = data(ii,2);
end
% Argon I
mass3        = 6.6312769E-26; % kg
dref3        = 3.0E-10;      % m
omega3       = 0.7;          % -
Tref3        = 300;          % Kelvin
clear data;
data=importdata('ArIon.csv');
nel(3)       = size(data,1) ;    % -
for ii = 1:nel(3)
    gel(3,ii)     = data(ii,1);
    elcharT(3,ii) = data(ii,2);
end
%-------------------------------------------------------------------------%
epsilon = 1;
Trefm   = 0.5*( Tref1 + Tref3 );
omegam  = 0.7;
massr   = mass1 * mass3 / ( mass3 + mass1 );
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% analytical recombination rate
kr_ana=kf./Ksm;
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
for mm = 1:size(Temp,2)
    coeff  = ( sqrt(PI) /( 2 * epsilon ) )      * ...
             ( dref1 + dref3 )^2    * ...
             ( Temp(mm) / Trefm)^( 1-omegam )   * ...
             sqrt( 2 * kappaB * Trefm/ massr )  ;

     Vref = PI/6 * ( dref1 + dref1 + dref3 )^3;
     Prob=kr_ana/Vref/coeff;
     a=248959927653.1;
     b=-1.875975049158;
     %a=4.7005072818987e+12;
     %b=-2.181468813430;
      %a=1.2e10; %4.5e11
      %b=-2.2;
     %a= 1.1e+18;
     %b= -3.4;
     kr_fit(mm)=coeff* Vref*a*Temp(mm)^b;
 end
%-------------------------------------------------------------------------%
disp('Prob for reverse reaction calculated.');
%-------------------------------------------------------------------------%

% fit
% xin=reshape(Temp,size(Temp,2),1);
% yin=reshape(Prob,size(Prob,2),1);
% Templaw =@(a,b,x) a*x.^b;
% t1=fit(xin,yin,Templaw);

%-------------------------------------------------------------------------%
figure;
semilogy(Temp,kr_ana,Temp,kr_fit)
l1=legend('analytisch','fit')
set(l1,'Orientation','horizontal','Location','NorthOutside','Box','off');
figure;
semilogy(Temp,Ksm,Temp,kf./kr_fit)
l1=legend('analytisch','fit')
set(l1,'Orientation','horizontal','Location','NorthOutside','Box','off');
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
disp('kr calculated.');
%-------------------------------------------------------------------------%


