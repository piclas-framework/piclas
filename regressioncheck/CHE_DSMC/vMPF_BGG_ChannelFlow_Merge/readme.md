# vMPF: Channel flow with particle merge
* Channel flow with a sudden volume change (mortar interface) and particle merge to 1000 particles
* Number density should remain constant across interface
* New particle weights after merge cannot distributed evenly as density distribution due to flow across cell has to be reproduced
* Example: Many particles enter the large cell with a low MPF, as they travel through the cell their MPF increases and only a few particles with a high MPF leave the cell