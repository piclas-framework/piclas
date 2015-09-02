%-------------------------------------------------------------------------%
% impact ionization kf
%-------------------------------------------------------------------------%
Temp=[10000:100:60000];
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
%-------------------------------------------------------------------------%
epsilon = 1;
Trefm   = 0.5*( Tref1 + Tref2 );
omegam  = 0.7;
massr   = mass1 * mass2 / ( mass1 + mass2 );
imax    = nel(2);
Eactive = 182829.95097605; %15.9596*q/kappaB; %182829.95097605; 
%-------------------------------------------------------------------------%
for mm = 1:size(Temp,2)
    coeff  = ( sqrt(PI) /( 2 * epsilon ) )      * ...
             ( dref1 + dref2 )^2    * ...
             ( Temp(mm) / Trefm)^( 1-omegam )   * ...
             sqrt( 2 * kappaB * Trefm/ massr )  ;
    summ = 0.0;
    summ2 = 0.0;
    for ii = 0:imax-1
        a1 = 2.5 - omegam;
        a2 = (Eactive - elcharT(2,ii+1) )/ Temp(mm);
        summ = summ + ...
               gammainc(a2,a1,'upper') *...
               gel(2,ii+1)*exp( - elcharT(2,ii+1) / Temp(mm) );
        summ2 = summ2 + ...
                gel(2,ii+1)*exp( - elcharT(2,ii+1)/Temp(mm));
    end
    kf(mm) = coeff * summ / summ2;
end
%-------------------------------------------------------------------------%
Annaloro = importdata('Annaloro_Arg.txt');
DSMC     = importdata('DSMC_Arg.txt');
figure;
p1=semilogy(Temp,kf,Annaloro(:,1),Annaloro(:,2),DSMC(:,1),DSMC(:,2))
set(p1(3),'Color',[0 0 1],'LineStyle','None','Marker','o');
l1=legend(p1,'Analytic','Annaloro','DSMC')
set(l1,'Orientation','horizontal','Location','NorthOutside',...
    'Box','off');
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
