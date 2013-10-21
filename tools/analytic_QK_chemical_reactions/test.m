%-------------------------------------------------------------------------%
% Testing skript
%
% Script written by P. Ortwein.
% Last modified:      2012-21-5 
%-------------------------------------------------------------------------%

reaction_database;

Temp =1000:10:20000;

% kf3 = calculate_kf_dissociation( Temp, reaction(16).id, reaction(16).Eactiv );
kf3 = calculate_kf_exchange( Temp, reaction(16).id, reaction(16).Eactiv );

Ksm3 = calculate_Ksm( Temp, reaction(16).id, reaction(16).Eactiv );

kr3 = kf3 ./ Ksm3;

% semilogy( Temp, Ksm3);


Temp  = reshape(Temp,size(Temp,2),size(Temp,1));
% kf3  = reshape(kf3,size(kf3,2),size(kf3,1));
kr3  = reshape(kr3,size(kr3,2),size(kr3,1));
kr = [Temp kr3];
% semilogy(Temp,kr)
% kr = [Temp kr];
% Prob = calculate_kr_prob( Temp,reaction(1).id, kr3 );

% filename = strcat('Prob_',reaction(1).name,'.csv');

% write_datafile(filename, Temp, Prob);