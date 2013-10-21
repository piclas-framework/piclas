%-------------------------------------------------------------------------%
% This file contains all the species data
% currently five species are listed:
% species order like in this file
% N2 , N, O2, O, NO
% reference values for the species are the NIST data base and the small
% table in Helrich (original form Mcxxx in 1972)
%
% Script written by P. Ortwein.
% Last modified:      2012-21-5 
%-------------------------------------------------------------------------%

number_of_species = 5;
nvib       = zeros(1,number_of_species);
delta      = zeros(1,number_of_species);
delta      = zeros(1,number_of_species);
Theta_rot  = zeros(1,number_of_species);
Theta_vib  = zeros(1,number_of_species);
nel        = zeros(1,number_of_species);
mass       = zeros(1,number_of_species);
gel        = zeros(number_of_species,5);
elcharT    = zeros(number_of_species,5);

%-------------------------------------------------------------------------%
% species entries 
%-------------------------------------------------------------------------%
% molecular nitrogen - N2
spec_name (1)    = {'N2'};
mass      (1)    = 4.652E-26;   % kg
dref      (1)    = 4.07E-10;    % m
omega     (1)    = 0.7;         % -
Tref      (1)    = 288;         % Kelvin
nvib      (1)    = 1;           % -
Theta_vib (1)    = 3395;        % Kelvin
delta     (1)    = 2;           % -
Theta_rot (1)    = 2.8757;      % Kelvin
nel       (1)    = 0;           % -
%-------------------------------------------------------------------------%
% atomic nitrogen - N
spec_name (2)    = {'N'};
mass      (2)    = 2.326E-26;   % kg
dref      (2)    = 3.00E-10;    % m
omega     (2)    = 0.7;         % -
Tref      (2)    = 300;         % Kelvin
nvib      (2)    = 0;           % -
nel       (2)    = 5;           % -
gel       (2,1)  = 4;           % -
elcharT   (2,1)  = 0;           % Kelvin
gel       (2,2)  = 6;           % -
elcharT   (2,2)  = 27658.8172;  % Kelvin
gel       (2,3)  = 4;           % -
elcharT   (2,3)  = 27671.3529;  % Kelvin
gel       (2,4)  = 2;           % -
elcharT   (2,4)  = 41491.4256;  % Kelvin
gel       (2,5)  = 4;           % - 
elcharT   (2,5)  = 41491.4102;  % Kelvin
%-------------------------------------------------------------------------%
% moleculer oxygen - O2
spec_name (3)    = {'O2'};
mass      (3)    = 5.312E-26;   % kg
dref      (3)    = 3.96E-10;    % m
omega     (3)    = 0.7;         % -
Tref      (3)    = 288;         % Kelvin
nvib      (3)    = 1;           % -
Theta_vib (3)    = 2239;        % Kelvin
delta     (3)    = 2;           % -
Theta_rot (3)    = 2.0801;      % Kelvin
nel       (3)    = 0;           % -
%-------------------------------------------------------------------------%
% atomic oxygen  - O
spec_name (4)    = {'O'};
mass      (4)    = 2.656E-26;   % kg
dref      (4)    = 3.00E-10;    % m
omega     (4)    = 0.7;         % -
Tref      (4)    = 288;         % Kelvin
nvib      (4)    = 0;           % -
nel       (4)    = 5;           % -
gel       (4,1)  = 5;           % -
elcharT   (4,1)  = 0;           % Kelvin
gel       (4,2)  = 3;           % -
elcharT   (4,2)  = 227.7006;    % Kelvin
gel       (4,3)  = 1;           % -
elcharT   (4,3)  = 326.5874;    % Kelvin
gel       (4,4)  = 5;           % -
elcharT   (4,4)  = 22829.5725;  % Kelvin
gel       (4,5)  = 1;           % - 
elcharT   (4,5)  = 48618.4102;  % Kelvin
%-------------------------------------------------------------------------%
% nitric oxide - NO
spec_name (5)    = {'NO'};
mass      (5)    = 4.983E-26;   % kg
dref      (5)    = 4.00E-10;    % m
omega     (5)    = 0.7;         % -
Tref      (5)    = 288;         % Kelvin
nvib      (5)    = 1;           % -
Theta_vib (5)    = 2817;        % Kelvin
delta     (5)    = 1;           % -
Theta_rot (5)    = 2.45291;     % Kelvin
nel       (5)    = 2;           % -
% following values are from birds program
gel       (5,1)  = 2;           % -  
elcharT   (5,1)  = 0;           % Kelvin
gel       (5,2)  = 2;           % -
elcharT   (5,2)  = 174.2;       % Kelvin
%-------------------------------------------------------------------------%

