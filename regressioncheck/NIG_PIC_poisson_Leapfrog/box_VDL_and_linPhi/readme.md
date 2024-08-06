# parallel plate testcase
* Single electron is accelerated according to electric field in plate capacitor
* Coupled Power output is compared with following analytical solution
* Constant direct current => constant acceleration:
  * $a=\frac{Eq}{m}$
* Time derivation => velocity:
  * $v=\frac{Eq}{m}t$
* Time derivation => distance:
  * $s=\frac{1}{2}\frac{Eq}{m}t^2$
* Coupled Power:
  * $P=vF=vEq=\frac{E^2q^2}{m}t$

* parameter.ini for activating coupled Power output

          CalcCoupledPower = T
