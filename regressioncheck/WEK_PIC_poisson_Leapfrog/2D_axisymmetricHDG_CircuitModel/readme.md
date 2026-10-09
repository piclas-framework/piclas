# PIC Circuit model: 2D L-shapred setup with dielectric
- Minimum working example for the circuit model given by Kentaro Hara et al., Effects of macroparticle
  weighting in axisymmetric particle-in-cell Monte Carlo collision simulations (https://doi.org/10.1088/1361-6595/acb28b)
- Anode electric potential boundary condition $V_{anode}$ with sinusiodal electric potential $V_{rf}=V_a\sin(\omega t + \phi)$
  that is connected to a capacitor, which leads to a bias voltage $V_{c}$ build-up on the capacitor: 

      V_{anode} = V_{rf} + V_{c}