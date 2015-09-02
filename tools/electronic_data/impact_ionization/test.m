%------------------------------------------------------
% validation of QK-Theory
%------------------------------------------------------

disp('hello')

% Energy equation at the beginning:

% Energy at the beginning: Ekin_Ar+ + Ekin_el + Eel_Ar+ 
% Ekin_el= 3/2 *k*N*T = 1/2 m v^2

kB=1.3807e-23;
% or
% EkinArIon
N0=1e5
EkinArIon0=3./2.*kB*N0*1000;
EkinEl0=3./2.*kB*N0*50000;

EelAr0=kB*N0*