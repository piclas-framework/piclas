# parallel plate AC testcase
* Single electron is accelerated according to alternating electric field in plate capacitor
* Coupled Power output is compared with following analytical solution
* Alternating current => alternating acceleration:
  * $a=\frac{Eq}{m}\sin(\omega t)$
* Time derivation => velocity:
  * $v=-\frac{1}{\omega}\frac{Eq}{m}\cos(\omega t)$
* Time derivation => distance:
  * $s=-\frac{1}{\omega^2}\frac{Eq}{m}\sin(\omega t)$
* Coupled Power:
  * $P=vF=vEq=-\frac{1}{\omega}\frac{E^2q^2}{m}\sin(\omega t)\cos(\omega t)$

* parameter.ini for activating coupled Power output

          CalcCoupledPower = T 
