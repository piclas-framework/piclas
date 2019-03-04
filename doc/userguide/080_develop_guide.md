\hypertarget{develop_guide}{}

# Development guidelines \label{chap:develop_guide}

This chapter contains information about the development process and other issues concerning Git (GitLab/GitHub).

## Development process

Naming convention for branches, workflow for development, milestones etc.

After the successful completion of all regression checks (check-in, nightly, weekly), the master.dev branch can be merged into the master.

## Release and deploy

### Collaborative Numerics Group

The master branch of development group can be merged after the succesful regression check with the master of the collaborative group. For this purpose, the collaborative repository can be added as a remote

       git remote add remote_name git@gitlab.com:collaborative-numerics-group/piclas/piclas.git

Now you can checkout the most recent version of the master branch of the collaborative-numerics-group and create a local branch with that version (a simple checkout will create a detached HEAD state)

       git fetch
       git checkout -b branch_name remote_name/master

The master branch of the development repository can now be merged into the newly created branch. Make sure to have the most recent version of the master branch (of the development repository) as well.

       git merge origin/master

Finally, the changes can be pushed from the *branch_name* to the master of collaborative-numerics-group

       git push remote_name master

If a tag has also been created, it should be pushed separately.

       git push remote_name tag_name

### GitHub

Upon completion of a milestone leading to tagged version, the tag should be deployed to GitHub.

## Cloning and compiling at the HLRS \label{sec:cloninghlrs}

Unfortunately, the GitHub server is not available on machines at the HLRS, such as the Hazelhen, due to restricted internet access. The workaround is to use ssh tunnels to access the GitHub repositories. Note that the reomte repositories hosted at teh GitLab at the Institute of Aerodynamics and Gasdynamics (IAG), no ssh tunnel is required and cloning works straight forwardly.

The following instructions to access the GitHub repositories on HLRS machines is taken from the HLRS wickie page, see [https://wickie.hlrs.de/platforms/index.php/Secure_Shell_ssh#Git](https://wickie.hlrs.de/platforms/index.php/Secure_Shell_ssh#Git).

### HTTPS

Unfortunately, just using a SSH tunnel as with the SSH and git protocols is not sufficient in this case. Instead, one has to connect via an additional SOCKS proxy on a machine that has unlimited access to the internet, e.g. your local machine.

In order to do so, establish a proxy by using a special feature of OpenSSH: 

       ssh -N -D 1080 localhost

This will establish some kind of a "loopback" SSH connection from your local machine to itself which will not execute any command (-N) but act as an SOCKS proxy on port 1080 (-D 1080).

On a second shell, now login to the desired HWW-system and forward a port on the remote machine (e.g. 7777) to the port on your local machine where the newly established SOCKS proxy is listening on (1080): 

       ssh -R 7777:localhost:1080 <system-name>.hww.de

By doing so, you have a SOCKS proxy listening on port 7777 of the HWW-system. Hence you can use this proxy for accessing remote git repositories. Unfortunately, the default versions of git installed on the HWW-systems are not capable of doing this. You hence have to load an appropriate version first: 

       module load tools/git

In order to use the proxy, you can now add "-c https.proxy='socks5://localhost:7777'" to your git commands, e.g.:

       git -c https.proxy='socks5://localhost:7777' clone https://github.com/piclas-framework/piclas.git

In order to avoid typing this in every git call, you can also set the respective port to be used whenever git talks to a remote repository via HTTPS by

       git config --global https.proxy 'socks5://localhost:7777'

Unfortunately, to connect with GitHub for pulling or pushing, the connection to Hazelhen has to be done via the ssh tunnel.
