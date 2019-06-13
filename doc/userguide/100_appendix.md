\hypertarget{appendix}{}

# Appendix \label{chap:appendix}

## Validated compiler combinations

| User |  System  |    Compiler     |  HDF5  |       MPI        |  CMake   |                       Notes                       |
| ---- | :------: | :-------------: | :----: | :--------------: | :------: | :-----------------------------------------------: |
| PO   |  Laptop  |    gnu4.9.2     | 1.8.16 |  openmpi-1.8.4   |  3.4.3   | gnu-sanitizer not working with DSMC, memory leak. |
|      | giganto  |    intel15.0    | 1.8.16 |  openmpi-1.8.2   | 2.8.12.2 |                    no autolist                    |
|      | hazelhen | intel15.0.4.223 | 1.8.14 | cray-mpich-7.3.1 |  3.4.2   |                   manual tecio                    |
|      |          | cray-libsci13.3 |  cray  |                  |          |                                                   |
|      |  Laptop  |    gnu5.2.0     | 1.8.16 |  openmpi-1.10.1  |  3.4.3   | gnu-sanitizer not working with DSMC, memory leak. |
|      |  Laptop  |    gnu7.3.+     | patch1 |  openmpi-3.0.0   |  3.10.+  |       Requires HDF_ROOT instead of HDF5_DIR       |
| SC   |  Laptop  |    gnu4.8.4     | 1.8.16 |  openmpi-1.6.5   |  3.2.2   |                                                   |
|      |  Laptop  |    gnu5.4.0     | 1.8.18 |  openmpi-1.8.8   |  3.5.1   |                                                   |
|      | hazelhen | intel15.0.4.223 | 1.8.14 | cray-mpich-7.3.1 |  3.4.2   |   set tecio path by hand (copy from old PICLas)   |
| WR   |  Laptop  |    gnu5.2.0     | 1.8.16 |  openmpi-1.10.0  |  3.4.3   |  linking only works with gnu5.2.0    --> solved   |
|      |  Laptop  |    gnu4.8.4     | 1.8.16 |  openmpi-1.10.0  |  3.4.3   |                                                   |
|      |  Laptop  |    gnu4.8.4     | 1.8.14 |  openmpi-1.6.5   |  3.4.3   |                                                   |
|      |  Laptop  |    gnu4.8.4     | 1.8.16 |  openmpi-1.6.5   |  3.4.3   |                                                   |
|      |  Laptop  |    gnu5.2.0     | 1.8.16 |  openmpi-1.6.5   |  3.4.3   |                                                   |
|      |  Laptop  |   intel15.0.4   | 1.8.16 |  openmpi-1.10.0  |  3.4.3   |                                                   |

## **pvconnect** script

```sh
#!/bin/bash
###############################################################################
# pvconnect
# --------------
# The script opens a ssh-tunnel to a given compute node via a mom node and 
# launches a paraview client with the correct -url option
#
# Developed by:
# --------------
# Ralf Schneider <schneider@hlrs.de>
#
# Last edited by Ralf Schneider - 22 April 2015
# 
# Edited by Nico Krais - 15.10.2018:
#   - Adopt to usage on a local machine, do not use virtualGL and do not
#     require nettest to search for an open port
#   - REMEMBER TO SET THE PATH TO YOUR PARAVIEW BINARY!!!
#
###############################################################################
#



usage()
{
	echo
	echo "USAGE: $0"
	echo "       -pvs pvserver[:port] -via host"
	echo
	echo "-pvs = ParaView server to connect to. Either the hostname alone or "
	echo "       hostname:port can be given. The hostname:port combination is"
	echo "       normally returned the pvserver by the comment:"
	echo "       Accepting connection(s): hostname:port"
	echo "       if no port is given the default port 11111 is used"
	echo "-via = Hostname via which to connect to the server. If the server was"
	echo "       launched by aprun this should be the name of the mom node"
	echo "       on which the aprun was executed."
	echo
	exit $1
}
#
if [ $# -eq 0 ]; then
	usage 0
fi
#
#module load tools/VirtualGL
#
port=9999
#
# -----------------------------------------------------------------------------
# Parse Arguments -------------------------------------------------------------
#
while [ $# -gt 0 ]
do
	case "$1" in
	-pvs*) acc_con=$2; shift ;;
	-via*) VIA=$2; shift ;;
	-help*) usage 0;;
   --help*) usage 0;;
	*) break ;;
	esac
	shift
done
#
if [ -z $acc_con ]; then
   echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
   echo "!! Argument for ParaView server is missing !!"
   echo "!! Please specify the -pvs option          !!"
   echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
   exit 1
fi
#
if [ -z $VIA ]; then
   echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
   echo "!! Found no host to connect through in arguments !!"
   echo "!! Please specify the -via option                !!"
   echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
   exit 1
fi
#
# -----------------------------------------------------------------------------
# Get target and port ---------------------------------------------------------
#
#
TARGET=`echo $acc_con | cut -d: -f1`
PVP=`echo $acc_con | cut -d: -f2`
#
if [ -z $PVP ]; then
   PVP=11111
   echo " "
   echo " Found no ParaView server port in arguments"
   echo " Using the default one: 11111"
fi
#
# -----------------------------------------------------------------------------
# Open tunnel -----------------------------------------------------------------
#
ssh -N -L $port:$TARGET:$PVP $VIA &
sshpid=$!
#
sleep 10
#
echo " "
echo " Opened ssh-tunnel with process id : ${sshpid}"
echo " Tunnel command was                : ssh -N -L $port:$TARGET:$PVP $VIA &"
#
# -----------------------------------------------------------------------------
# Launch paraview client ------------------------------------------------------
#
echo " "
echo " Launching paraview client : paraview -url=cs://localhost:$port"
/PATH/TO/PARAVIEW/bin/paraview --mpi -url=cs://localhost:$port
#
kill ${sshpid}
#
```

