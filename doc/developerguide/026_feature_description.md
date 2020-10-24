\hypertarget{feature_description}{}

# Feature Description \label{chap:feature_description}

This chapter describes the structure and logic behind some of the features of PICLas.

## Chemistry Module

* Three available chemistry models
  * Quantum Kinetic (QK)
  * Cross-sections (XSec)
  * Total Collision Energy (TCE)

The model of each reaction can be chosen separately. If a collision pair has multiple reaction paths (e.g. CH3 + H, two possible dissociation reactions and a recombination), the reaction paths with the same model are treated together. If a reaction path is selected, the reaction is performed and the following routines of the chemistry module are not performed. Thus, if its determined that a reaction path with a QK reaction is to be performed, the other reaction paths will not be considered.

1. Calculate the collision probability
   a. XSec: Calculate the reaction probability as it has to be added to the collision probability (if no collision cross-section data is used)
2. Loop over colliding pairs
  1. DSMC_perform_collision: Determined whether CollisMode = 3 (required for chemical reactions)
     1. ReactionDecision
        1. Get the number of reaction paths available for the specific collision pair
        2. Check whether any QK reactions are possible or calculate the reaction probability for TCE; XSec reaction probability have to be calculated during the collision probability calculation
     2. Treat QK reaction paths (if a reaction path is chosen, leave the routine)
     3. Treat XSec reaction paths (if a reaction path is chosen, leave the routine)
     4. Treat TCE reaction paths