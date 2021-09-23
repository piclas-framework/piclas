# Multi-species background gas with chemistry
* Reservoir simulation of high-temperature air to test the TCE implementation and a multi-species background gas
* N2 and O2 mixture as constant background gas
* N and O are inserted intially
* Reactions paths:
  *  N2 + N/O -> N + N + N/O (Reaction 1 + 4)
  *  O2 + N/O -> O + O + N/O (Reaction 2 + 5)
  *  N2 + O   -> NO + N (Reaction 3)
  *  Reaction 4 + 5 are added through non-reactives array and should have the same rates as reaction 1 and 3, respectively
*  Comparison of the output (Database_Ttrans_XXXXXK.csv) with the theoretical Arrhenius rates shows excellent agreement, theoretical values given below:

| Temperature [K] | Reaction 1 [m3/s] | Reaction 2 [m3/s] | Reaction 3 [m3/s] |
| :-------------: | :---------------: | :---------------: | :---------------: |
|      10000      |     2.405E-19     |     4.327E-17     |     2.284E-18     |
|      15000      |     5.472E-18     |     1.712E-16     |     5.477E-18     |
|      20000      |     2.278E-17     |     2.997E-16     |     7.790E-18     |
|      25000      |     4.945E-17     |     3.888E-16     |     9.150E-18     |
|      30000      |     7.857E-17     |     4.398E-16     |     9.849E-18     |