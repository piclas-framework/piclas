# MCC modelling of back-scattering collisions
* Cold (no thermal velocity component) beam of ions traverses channel, where elastic and back-scattering collisions are possible
* Selected input velocity corresponds to a cross-section, where 0.88 of collisions are back-scattering and a total collision probability of 0.0131 -> p_back = 0.88 * 0.0131
* Number of potential interactions until the particles reach the end of the domain is around n=16, resulting in a probability of 0.171 that a back-scattering event occurs (1-(p-1)^n)
* Mass flow at the output should be reduced to 0.0171 kg/s (from 0.0206 kg/s) accordingly as the particles are decelerated and the domain fills up
* Comparing outgoing flux (1/s) with reference, which corresponds to roughly 0.017 kg/s