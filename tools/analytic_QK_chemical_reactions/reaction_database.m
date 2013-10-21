%-------------------------------------------------------------------------%
% matlab database of chemical reactions
% This database contains only dissociation and endothermic exchange
% reactions. All data is used to calculate the reverse reaction rate
% coefficient. These rate coefficients are used to calculate the two
% coefficients for recombination a,b or the three coefficients for the
% exothermic exchange reactions a,b and eb
%
% use with care
%
% Reaction references:
% Reaction  1 - 10 : Bird 2001
% Reaction 11 - 15 : Bird 1985
% Reaction 16      : Bird 2001
% Reaction 17      : ?? Taken from DSMCSpecies.ini
%
% Script written by P. Ortwein.
% Last modified:      2012-21-5 
%-------------------------------------------------------------------------%

number_of_reactions = 17;

%-------------------------------------------------------------------------%
% Only dissociation
%-------------------------------------------------------------------------%
% Reaction 1
% N2 + N2 --> N + N + N2
reaction(1).name    = 'N2_plus_N2_to_N_plus_N_plus_N2';
reaction(1).type    = 'dissociation';
reaction(1).id      = [1 1 2 2 1];
reaction(1).Eactiv  = 113200;
%-------------------------------------------------------------------------%
% Reaction 2
% N2 + O2 --> N + N + O2
reaction(2).name    = 'N2_plus_O2_to_N_plus_N_plus_O2';
reaction(2).type    = 'dissociation';
reaction(2).id      = [1 3 2 2 3];
reaction(2).Eactiv  = 113200;
%-------------------------------------------------------------------------%
% Reaction 3
% N2 + NO --> N + N + NO
reaction(3).name    = 'N2_plus_NO_to_N_plus_N_plus_NO';
reaction(3).type    = 'dissociation';
reaction(3).id      = [1 5 2 2 5];
reaction(3).Eactiv  = 113200;
%-------------------------------------------------------------------------%
% Reaction 4
% N2 + N  --> N + N + N
reaction(4).name    = 'N2_plus_N_to_N_plus_N_plus_N';
reaction(4).type    = 'dissociation';
reaction(4).id      = [1 2 2 2 2];
reaction(4).Eactiv  = 113200;
%-------------------------------------------------------------------------%
% Reaction 5
% N2 + O  --> N + N + O
reaction(5).name    = 'N2_plus_O_to_N_plus_N_plus_O';
reaction(5).type    = 'dissociation';
reaction(5).id      = [1 4 2 2 4];
reaction(5).Eactiv  = 113200;
%-------------------------------------------------------------------------%
% Reaction 6
% O2 + N2 --> O + O + N2
reaction(6).name    = 'O2_plus_N2_to_O_plus_O_plus_N2';
reaction(6).type    = 'dissociation';
reaction(6).id      = [3 1 4 4 1];
reaction(6).Eactiv  = 59360;
%-------------------------------------------------------------------------%
% Reaction 7
% O2 + O2 --> O + O + O2
reaction(7).name    = 'O2_plus_O2_to_O_plus_O_plus_O2';
reaction(7).type    = 'dissociation';
reaction(7).id      = [3 3 4 4 3];
reaction(7).Eactiv  = 59360;
%-------------------------------------------------------------------------%
% Reaction 8
% O2 + NO --> O + O + NO
reaction(8).name    = 'O2_plus_NO_to_O_plus_O_plus_NO';
reaction(8).type    = 'dissociation';
reaction(8).id      = [3 5 4 4 5];
reaction(8).Eactiv  = 59360;
%-------------------------------------------------------------------------%
% Reaction 9
% O2 + N  --> O + O + N 
reaction(9).name    = 'O2_plus_N_to_O_plus_O_plus_N';
reaction(9).type    = 'dissociation';
reaction(9).id      = [3 2 4 4 2];
reaction(9).Eactiv  = 59360;
%-------------------------------------------------------------------------%
% Reaction 10
% O2 + O  --> O + O + O
reaction(10).name   = 'O2_plus_O_to_O_plus_O_plus_O';
reaction(10).type   = 'dissociation';
reaction(10).id     = [3 4 4 4 4];
reaction(10).Eactiv = 59360;
%-------------------------------------------------------------------------%
% Reaction 11
% NO + N2  --> N + O + N2
reaction(11).name   = 'NO_plus_N2_to_N_plus_O_plus_N2';
reaction(11).type   = 'dissociation';
reaction(11).id     = [5 1 2 4 1];
reaction(11).Eactiv = 75500;
%-------------------------------------------------------------------------%
% Reaction 12
% NO + O2  --> N + O + O2
reaction(12).name   = 'NO_plus_O2_to_N_plus_O_plus_O2';
reaction(12).type   = 'dissociation';
reaction(12).id     = [5 3 2 4 3];
reaction(12).Eactiv = 75500;
%-------------------------------------------------------------------------%
% Reaction 13
% NO + NO  --> N + O + NO
reaction(13).name   = 'NO_plus_NO_to_N_plus_O_plus_NO';
reaction(13).type   = 'dissociation';
reaction(13).id     = [5 5 2 4 5];
reaction(13).Eactiv = 75500;
%-------------------------------------------------------------------------%
% Reaction 14
% NO + N   --> N + O + N
reaction(14).name   = 'NO_plus_N_to_N_plus_O_plus_N';
reaction(14).type   = 'dissociation';
reaction(14).id     = [5 2 2 4 2];
reaction(14).Eactiv = 75500;
%-------------------------------------------------------------------------%
% Reaction 15
% NO + O   --> N + O + O
reaction(15).name   = 'NO_plus_O_to_N_plus_O_plus_O';
reaction(15).type   = 'dissociation';
reaction(15).id     = [5 4 2 4 4];
reaction(15).Eactiv = 75500;
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% Only exchange
%-------------------------------------------------------------------------%
% Reaction 16
% N2 + O  --> NO + N
reaction(16).name    = 'N2_plus_O_to_NO_plus_N';
reaction(16).type    = 'exchange';
reaction(16).id      = [1 4 5 2 0];
reaction(16).Eactiv  = 42938;         % 42938 - Park01 37500 - Park85         
%-------------------------------------------------------------------------%
% Reaction 17
% O2 + N  --> NO + O
reaction(17).name    = 'O2_plus_N_to_NO_plus_O';
reaction(17).type    = 'exchange';
reaction(17).id      = [3 2 5 4 0];
reaction(17).Eactiv  = 3600;
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% all forward reation rates are specified.
%-------------------------------------------------------------------------%