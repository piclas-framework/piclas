%-------------------------------------------------------------------------%
% calculate equilibirum constant by means of statistical mechanic
% is only for Ionization of Ar
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% Script written by P. Ortwein.
% Last modified:      2014-20-3 
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% temperatur range
Temp=[10000:100:60000];
Eactive = 182829.95097605;
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% physical constants
PI      = 3.14159265358979323846;  % -
kappaB  = 1.380650524e-23;         % J / K
NL      = 2.6867805e25;            % m^-3
PlanckC = 6.62606957e-34;          % J s
q     = 1.602176565E-19 % J / eV
%-------------------------------------------------------------------------%

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
Tref1        = 1000;          % Kelvin
nel(1)       = 1;          % -
gel(1,1)     = 2;         % -
elcharT(1,1) = 0;  % Kelvin
% Argon
mass2        = 6.631368E-26;     % kg
dref2        = 3.0E-10;      % m
omega2       = 0.7;          % -
Tref2        = 300;          % Kelvin
data=importdata('Ar.csv');
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

%-------------------------------------------------------------------------%
% initialize arrays
Ksm  = zeros(1,size(Temp,2));
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
for mm = 1:size(Temp,2)
%         mm = 1
% argon
    qtr_ar=( ( 2 * PI * mass2 / PlanckC) * ...
                 ( kappaB / PlanckC )    * ...
                 ( Temp(mm) ) )^1.5;
    qel(2)=0.; 
    for kk = 1:nel(2);
       qel(2) = qel(2) + ...
       gel(2,kk) * exp( - elcharT(2,kk) / Temp(mm) );
    end
    q_ar = qtr_ar * qel(2);
% argon ion
    qtr_arion=( ( 2 * PI * mass3 / PlanckC) * ...
                 ( kappaB / PlanckC )    * ...
                 ( Temp(mm) ) )^1.5;
    qel(3)=0.; 
    for kk = 1:nel(3);
       qel(3) = qel(3) + ...
       gel(3,kk) * exp( - elcharT(3,kk) / Temp(mm) );
    end
    q_arion = qtr_arion * qel(3);
% electron
    qtr_e=( ( 2 * PI * mass1 / PlanckC) * ...
                 ( kappaB / PlanckC )    * ...
                 ( Temp(mm) ) )^1.5;
    qel(1)=0.; 
    for kk = 1:nel(1);
       qel(1) = qel(1) + ...
       gel(1,kk) * exp( - elcharT(1,kk) / Temp(mm) );
    end 
    q_e = qtr_e * qel(1);         
    Ksm(mm) = q_arion*q_e*q_e / ( q_ar*q_e  )*exp( - Eactive / Temp(mm) ); 
end
%-------------------------------------------------------------------------%
figure;
semilogy(Temp,Ksm);
legend('statistical mechanics');
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
disp('Ksm calculated.');
%-------------------------------------------------------------------------%


