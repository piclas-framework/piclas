\hypertarget{appendix}{}

# Appendix \label{chap:appendix}

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

