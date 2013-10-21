%-------------------------------------------------------------------------%
% TMain matlab script 
% This small script controls the creation of the fitting data base
% 
% Feel free to use and modify. ( But only this small peace of code :P)
%
% Script written by P. Ortwein.
% Last modified:      2012-21-5 
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% Reaction Database
reaction_database;
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% Temperature range 
Temp =1000:100:20000;
%-------------------------------------------------------------------------%

for iReac = 1: number_of_reactions
    % forward reaction rate
    if iReac <= 15
        kf  = calculate_kf_dissociation( Temp, reaction(iReac).id, ...
                                            reaction(iReac).Eactiv );
    else
        kf = calculate_kf_exchange( Temp , reaction(iReac).id, ...
                                           reaction(iReac).Eactiv );
    end
    % equilibrium constat
    Ksm = calculate_Ksm( Temp, reaction(iReac).id, ...
                         reaction(iReac).Eactiv );
    % reverse rate
    kr = kf ./ Ksm;
%     semilogy( Temp, kr);
    % calculate reaction probability
    Prob = calculate_kr_prob( Temp,reaction(iReac).id, kr );
    filename = strcat('kf_',reaction(iReac).name,'.csv');
    write_datafile(filename, Temp, kf);
    filename = strcat('kr_',reaction(iReac).name,'.csv');
    write_datafile(filename, Temp, kr);    
    filename = strcat('Ksm_',reaction(iReac).name,'.csv');
    write_datafile(filename, Temp, Ksm);    
    % -------
    clear kf;
    clear kr;
    clear ksm;
    clear Prob;
end

%-------------------------------------------------------------------------%
disp('All files for fitting reactions written.')
disp('Q-K Bingo can start')
%-------------------------------------------------------------------------%