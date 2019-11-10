import numpy as np
import h5py

datatype = np.dtype(np.float64)

data_Ar = np.genfromtxt('Ar-e_effective.txt',
                    delimiter='\t',
                    skip_header=1,
                    #names=True,
                    dtype=datatype,
                    )

data_He = np.genfromtxt('He-e_effective.txt',
                    delimiter='\t',
                    skip_header=1,
                    #names=True,
                    dtype=datatype,
                    )

hdf = h5py.File('MCC_Database.h5', 'w')
hdf.attrs['Info'] = 'Collision cross-section database for MCC-based probability calculation with PICLas. First column is the collision energy in [eV], second column is the cross-section in [m^2]'

dataset1 = hdf.create_dataset('Ar-electron', data=data_Ar)
dataset1.attrs['Type'] = 'effective'
dataset1.attrs['Source'] = 'Phelps database, www.lxcat.net, retrieved on October 27, 2019. COMMENT: Yamabe, Buckman, and Phelps, Phys. Rev. 27, 1345 (1983). Revised Oct 1997.'

dataset2 = hdf.create_dataset('He-electron', data=data_He)
dataset2.attrs['Type'] = 'effective'
dataset2.attrs['Source'] = 'Phelps database, www.lxcat.net, retrieved on October 27, 2019. COMMENT: He Momentum Transfer - FROM Crompton, Elford, and Jory, Australian J. Phys. 20, 369 (1967) Crompton, Elford, and Robertson, ibid 23, 667 (1970) Miloy and Crompton,Phys. Rev. A 15, 1847 (1977) AT LOW ENERGY, and from Hayashi, Institute of Plasma Physics Report No. IPPJ-AM-19, (1981) AT HIGH ENERGIES.'

hdf.close()