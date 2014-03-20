# Script to vary the temperature for the reservoir simulationO
# this is an equilibirium simulation, thus, temp-vib and temp-rot are changed

# setting start value, end and temperature step
$temperature       = 35000.;
$temperature_end   = 60000.;
$delta_temperature = 5000.;


while ($temperature <= $temperature_end) {
    print "Changing temperature to $temperature\n";
    `cp ./DSMCTest.ini_old ./DSMCTest.ini`;
    `perl -pi -e 's/PLACEHOLDER/$temperature/g' ./DSMCTest.ini`;
    system("./flexi DSMCTest.ini DSMCSpecies.ini");
    $temperature = $temperature + $delta_temperature;
}
print "*  ...finish!\n";
