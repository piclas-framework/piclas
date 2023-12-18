# Cell-local insertion with variable MPF
* Initial insertion with cell_local in the whole domain on a non-conform mesh with different cell volumes
* Utilization of vMPF and a SplitThreshold to insert 500/1000 particles per cell at a constant number density
* Comparison of the number density after one iteration, reference file was generated with 10000 particles per cell
* Without vMPF, less particles are inserted and the number density fluctuates stronger, failing the comparison
* Additionally, the total particle number is compared which should correspond to a multiple between number of cells (320) and the split threshold