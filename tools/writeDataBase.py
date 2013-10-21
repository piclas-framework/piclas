#Set Up the visualization
OpenDatabase("FILENAME_PLACEHOLDER.dat")
AddPlot("Pseudocolor", "ElectricFieldX")
DrawPlots() 
#export the database
dbAtts = ExportDBAttributes()
dbAtts.db_type = "VTK"
dbAtts.filename = "FILENAME_PLACEHOLDER"
dbAtts.variables = ("ElectricFieldX","ElectricFieldY","ElectricFieldZ","MagneticFieldX","MagneticFieldY","MagneticFieldZ","Phi","Psi")
ExportDatabase(dbAtts)
exit()
