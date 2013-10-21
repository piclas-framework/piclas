# Script to vary the temperature for the reservoir simulationO
# this is an equilibirium simulation, thus, temp-vib and temp-rot are changed

# setting start value, end and temperature step
$temperature       = 8000.;
$temperature_end   = 10000.;
$delta_temperature = 500.;


while ($temperature <= $temperature_end) {
    print "Changing temperature to $temperature\n";
    `cp ./DSMCQK.ini_old ./DSMCQK.ini`;
    `perl -pi -e 's/PLACEHOLDER/$temperature/g' ./DSMCTest.ini`;
    system("./flexi DSMCQK.ini DSMCSpecies_QK.ini");
    $temperature = $temperature + $delta_temperature;
}
print "*  ...finish!\n";
