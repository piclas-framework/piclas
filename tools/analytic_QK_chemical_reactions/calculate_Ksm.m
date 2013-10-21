function [ Ksm ] = calculate_Ksm( Temp , id, Eactive )
%-------------------------------------------------------------------------%
% matlab function calculating the equilibirum constant based on 
% statistical mechanics
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
Ksm  = zeros(1,size(Temp,2));
qtr  = zeros(1,4);
qvib = zeros(1,4);
qrot = zeros(1,4);
qel  = zeros(1,4);
q    = zeros(1,4);
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
if id(5) == 0
    for mm = 1:size(Temp,2)
%          mm = 1;
         for LL = 1:4 
            qtr(LL) =( ( 2 * PI * mass(id(LL)) / PlanckC) * ...
                     ( kappaB / PlanckC )             * ...
                     ( Temp(mm) ) )^1.5;
            if ( nvib(id(LL)) > 0 )
               qvib(LL) = 1 / ( 1 - exp(-Theta_vib(id(LL)) / Temp(mm) )   );
               qrot(LL) = Temp(mm) / ( delta(id(LL)) * Theta_rot(id(LL)) );
            else
               qvib(LL) = 1;
               qrot(LL) = 1;
            end 
            qel(LL) = 1;
            if ( nel(id(LL)) > 0 )
               for kk = 1:nel(id(LL));
                   qel(LL) = qel(LL) + ...
                   gel(id(LL),kk) * exp( - elcharT(id(LL),kk) / Temp(mm) );
               end
            end
            q(LL) = qtr(LL) * qvib(LL) * qrot(LL) * qel(LL);
        end
        Ksm(mm) = q(3)*q(4) / ( q(1) * q(2) ) * exp( - Eactive / Temp(mm) ); 
     end
else
    for mm = 1:size(Temp,2)
%         mm = 1
        for LL = 1:5
            qtr(LL) =( ( 2 * PI * mass(id(LL)) / PlanckC) * ...
                     ( kappaB / PlanckC )             * ...
                     ( Temp(mm) ) )^1.5;
            if ( nvib(id(LL)) > 0 )
               qvib(LL) = 1 / ( 1 - exp(-Theta_vib(id(LL)) / Temp(mm) )   );
               qrot(LL) = Temp(mm) / ( delta(id(LL)) * Theta_rot(id(LL)) );
            else
               qvib(LL) = 1;
               qrot(LL) = 1;
            end 
            qel(LL) = 1;
            if ( nel(id(LL)) > 0 )
               for kk = 1:nel(id(LL));
                   qel(LL) = qel(LL) + ...
                   gel(id(LL),kk) * exp( - elcharT(id(LL),kk) / Temp(mm) );
               end
            end
            q(LL) = qtr(LL) * qvib(LL) * qrot(LL) * qel(LL);
        end
        Ksm(mm) = q(3)*q(4)*q(5) / ( q(1) * q(2) ) * exp( - Eactive / Temp(mm) ); 
    end
end
%-------------------------------------------------------------------------%
disp('Ksm calculated.');
%-------------------------------------------------------------------------%
end

