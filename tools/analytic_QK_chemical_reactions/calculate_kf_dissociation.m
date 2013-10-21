function [ kf ] = calculate_kf_dissociation( Temp , id, Eactive )
%-------------------------------------------------------------------------%
% matlab function calculating the forward rate coefficient for 
% dissociation and recombination reactions
% looks like the following
% AB + M --> A + B + M
%
% Script written by P. Ortwein.
% Last modified:      2012-21-5 
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% include physical constants ( if not called)
physical_constants;
species_data;
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% initialize arrays
kf      = zeros(1,size(Temp,2));
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% precalculation and constants
epsilon = 1;
Trefm   = 0.5*( Tref(id(1)) + Tref(id(2)) );
omegam  = 0.5*( omega(id(1)) + omega(id(2)) );
massr   = mass(id(1)) * mass(id(2)) / ( mass(id(1)) + mass(id(2)) );
imax    = floor( Eactive / Theta_vib(id(1)) );
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
for mm = 1:size(Temp,2)
    if id(1) == id(2)
        epsilon = 2;
    end
    coeff  = ( sqrt(PI) /( 2 * epsilon ) )      * ...
             ( dref(id(1)) + dref(id(2)) )^2    * ...
             ( Temp(mm) / Trefm)^( 1-omegam )   * ...
             sqrt( 2 * kappaB * Trefm/ massr )  ;
    summ = 0.0;
    zvib = 1 / ( 1 - exp( - Theta_vib(id(1)) / Temp(mm) ) );
    for ii = 0:imax
        a1 = 2.5 - omegam;
        a2 = (Eactive - ii * Theta_vib(id(1)) )/ Temp(mm);
        summ = summ + ...
               gammainc(a2,a1,'upper') *...
               exp( -ii * Theta_vib(id(1)) / Temp(mm) );
    end
    kf (mm) = coeff * summ / zvib;
end
%-------------------------------------------------------------------------%
disp('kf for dissociation calculated.');
%-------------------------------------------------------------------------%
end

