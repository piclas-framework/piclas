### Check-Out from Gitlab.com

For the Cray XC40 Hazelhen cluster, direct access to gitlab.com is blocked. Instead, 
a ssh-tunnel has to be used. The required steps are the following

1. Login to a Server/Computer with direct access to the Hazelhen cluster, e.g. 
    * iagpc151
    * IP-Address from the range of the IAG
    * Certain IP-Addresses within Edurom of University of Stuttgart
2. Use a SSH-Tunnel to connect to the Cluster
 
* **Version 1** 
      
    Login with SSH-tunnel and clone repository
    from git@github.com:piclas/piclas.git

```
        ssh -R 7777:gitlab.com:22 USER@hazelhen.hww.de
        git clone ssh://git@localhost:7777/piclas/piclas.git
```
* **Version 2** 
      
    Login with SSH-tunnel and clone repository
    from git://github.com/piclas/piclas.git

```
        ssh -R 7777:gitlab.com:9418 USER@hazelhen.hww.de
        git clone git://localhost:7777/piclas/piclas.git
```


### Building with ccmake on hazelhen

For building on the hazelhen-cluster, the following steps are needed:

* include in .bashrc or .profile:
        
		module unload PrgEnv-cray
        module load PrgEnv-intel
        module load cray-hdf5-parallel
        module load tools/cmake

* Profiling with Craypad

    * Compile PICLas with 
'''    
        module load perftools-base && module load perftools-lite && export CRAYPAT_LITE=event_profile
'''    
    * Run PICLas with normal submit script
    * Program has to finish normally! Enough time during execution. Note, that the profiled version is slower, hence, 
      the testqueue is maybe to short. 
    * Visualize the *.app2 files 
