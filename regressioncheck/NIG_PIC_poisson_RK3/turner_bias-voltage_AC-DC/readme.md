# 1D turner setup and bias voltage feature and power control
- He plasma between two electrodes one of which is usually AC, but for this bias-voltage test case it is set to AC or DC starting at 0V
- bias voltage for
  - DC: BCType=50 (starting from 0V)
  - AC: BCType=51 (requiring a refstate to define the sinusoidal function, starting from 10.)
  - AC+power control: BCType=52 which utilizes a bias voltage offset and cos function amplitude via power control (starting from 1V)
- ion excess on both electrodes generates a positive bias voltage and electron excess a negative potential
- comparison of PartAnalyze.csv with reference file and integration of the bias voltage over time from SurfaceAnalyze.csv
