# PIC-MCC - Electron beam acceleration and ionization of background gas with species-specific time step
* Inserting electrons with a zero velocity and at a low temperature
* Acceleration with a potential difference of 20 kV
* Ionization of background gas through MCC
* Species-specific time step to push electrons with dt = 1e-12s but not included in MCC
  * ManualTimeStep  = 1.0000E-11
  * Part-Species1-TimeStepFactor = 0.1
  * Part-VariableTimeStep-DisableForMCC = T
* Reference file was generated using a regular time step of dt = 1e-12s
* Comparing the number density of the electron and ions in the channel