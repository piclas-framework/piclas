function [ Prob ] = calculate_kr_prob(Temp, id, kr )
%-------------------------------------------------------------------------%
% matlab function calculation the probability
% required to obtain a fit for the recombination and exothermic 
% exchange reactions
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
Prob      = zeros(1,size(Temp,2));
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% precalculation and constants
epsilon = 1;
Trefm   = 0.5*( Tref(id(3)) + Tref(id(4)) );
omegam  = 0.5*( omega(id(3)) + omega(id(4)) );
massr   = mass(id(3)) * mass(id(4)) / ( mass(id(3)) + mass(id(4)) );
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
for mm = 1:size(Temp,2)
    if id(3) == id(4)
        epsilon = 2;
    end
    coeff  = ( sqrt(PI) /( 2 * epsilon ) )      * ...
             ( dref(id(3)) + dref(id(4)) )^2    * ...
             ( Temp(mm) / Trefm)^( 1-omegam )   * ...
             sqrt( 2 * kappaB * Trefm/ massr )  ;
     if id(5) == 0
         Prob(mm) = kr(mm) / ( coeff );
     else
         Vref = PI/6 * ( dref(id(3)) + dref(id(4)) + dref(id(5)) )^3;
         Prob(mm) = kr(mm) / ( coeff * Vref );
     end
 end
%-------------------------------------------------------------------------%
disp('Prob for reverse reaction calculated.');
%-------------------------------------------------------------------------%
end

