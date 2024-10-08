# Merge Request To-Do

* [ ] Update of piclas version numbers using the [updatePiclasVersion.sh](https://github.com/piclas-framework/piclas/blob/master/tools/updatePiclasVersion.sh) tool by running "./tools/updatePiclasVersion.sh 3.2.1" in the top level directory
  * [ ] `./src/globals/globals_vars.f90` (`MajorVersion`, `MinorVersion` and `PatchVersion`)
  * [ ] `.github/workflows/cmake-ninja.yml` (`piclas-binaries-vX.X.X`)
* [ ] Update the [Prerequisites](https://piclas.readthedocs.io/en/latest/userguide/installation.html#prerequisites) (e.g. Ubuntu or other Linux distros)
* [ ] AppImage
  * [ ] Copy the repository to the testing GitHub repository and test the AppImage compilation process before doing the real GitHub release
  * [ ] Update the prerequisites list
    *  Release 1.0.0 - 3.3.0: glibc 2.17 + OpenMPI X.X.X
    *  Release 3.3.0 - X.X.X: glibc 2.18 [9b09c795](https://piclas.boltzplatz.eu/piclas/piclas/-/commit/9b09c7957800915cbdf5ecc4a0d8ba43993060da) + OpenMPI X.X.X
  * [ ] Update the MPI/HDF5/PETSc library versions in the [GitHub workflow](https://github.com/piclas-framework/piclas/blob/master/.github/workflows/cmake-ninja.yml) if required
    * [ ] Push to the feature branch with the changed cmake-ninja.yml file that features an updated version of MPI/HDF5/PETSc for building the AppImage to the testing repository via
        ```
        cd ~
        git clone -b feature.branch.name --single-branch git@piclas.boltzplatz.eu:piclas/piclas.git piclas-testing
        cd piclas-testing
        git remote add piclas-testing git@github.com:scopplestone/piclas-testing.git
        git config --list
        git push --force piclas-testing feature.branch.name
        ```
      and go to the [GitHub Workflows](https://github.com/scopplestone/piclas-testing/actions) page and check if the build is running
    * [ ] Test functionality of the AppImage by running the respective executable and test cases (reggies)
* [ ] Update the [Required Libraries](https://piclas.readthedocs.io/en/latest/userguide/installation.html#required-libraries) from the reggie server (see gitlab CI/CD) for mpich and OpenMPI
  *  Release 3.3.0: [openmpi](https://piclas.boltzplatz.eu/piclas/piclas/-/jobs/675270), [mpich](https://piclas.boltzplatz.eu/piclas/piclas/-/jobs/674955)
* [ ] Release notes

## Regression Testing

* [ ] Check-in
* [ ] Nightly
* [ ] Weekly

# Release Notes

## Release 3.X.X

### Breaking/Parameter changes

*

### Documentation/Tools/Regression testing

*

### Features

*

### Improvements

*

### Fixes

*
