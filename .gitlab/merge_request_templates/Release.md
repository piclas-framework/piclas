# Merge Request To-Do

* [ ] Update the piclas version numbers using the [updatePiclasVersion.sh](https://github.com/piclas-framework/piclas/blob/master/tools/updatePiclasVersion.sh) tool by running
    ```
    ./tools/updatePiclasVersion.sh 3.2.1
    ```
  in the top level directory. Adjust the version number argument to the current. Check if the following two files have been correctly changed/updated:
  * [ ] Check that `./src/globals/globals_vars.f90` (`MajorVersion`, `MinorVersion` and `PatchVersion`) is updated.
  * [ ] Check that `.github/workflows/cmake-ninja.yml` (`piclas-binaries-vX.X.X`) is updated.
* [ ] Update the prerequisite table for compiling piclas under [Required Libraries](https://piclas.readthedocs.io/en/latest/userguide/installation.html#required-libraries) by checking the versions that are currently used by the reggie server (see gitlab CI/CD) for mpich and OpenMPI. For example, see Release 3.3.0 pipelines for [openmpi](https://piclas.boltzplatz.eu/piclas/piclas/-/jobs/675270) and [mpich](https://piclas.boltzplatz.eu/piclas/piclas/-/jobs/674955). Check the most recent pipelines running for `master.dev` and use these versions in the table.
* [ ] AppImage: Update the library versions and push the most recent version of the `master.dev` repository to the testing GitHub repository and make sure that the AppImage compilation process succeeds before doing the real GitHub release.
  * [ ] Update the GCC/MPI/HDF5/PETSc library versions used in the GitHub workflow by adjusting the versions in [cmake-ninja.yml](https://piclas.boltzplatz.eu/piclas/piclas/-/blob/master.dev/.github/workflows/cmake-ninja.yml) with which the AppImage is built.
  * [ ] Update the information in the user guide regarding the versions of the glibc (currently the OS version is used) and OpenMPI dependencies in the [AppImage dependeny table](https://piclas.readthedocs.io/en/latest/userguide/installation.html#appimage-executable-download) that are required for the AppImage. Check the following examples and linked commits on how this is done:
    *  Release 1.0.0 - 3.3.0: glibc 2.17 + OpenMPI X.X.X
    *  Release 3.3.0 - X.X.X: glibc 2.18 [9b09c795](https://piclas.boltzplatz.eu/piclas/piclas/-/commit/9b09c7957800915cbdf5ecc4a0d8ba43993060da) + OpenMPI X.X.X
  * [ ] Push the feature branch `feature.branch.name`, which is usually the `master.dev` branch, with the changed `cmake-ninja.yml` file that features an updated version of the MPI/HDF5/PETSc libraries for building the AppImage to the testing repository via gitlab Pipelines. Select "New Pipeline" and set "Run for branch name or tag" to the required `feature.branch.name` and supply `DO_CREATE_APPIMAGE` as "Input variable key" and set it to `T` for "Input variable value". Then go to the [piclas-testing GitHub Workflows](https://github.com/scopplestone/piclas-testing/actions) page and check if the build is running. Fix any errors that might occur during the build process and repeat the process.
  * [ ] When the workflow action pipeline has successfully created the AppImages, test their functionality by running a series of tests with the executables (check-in reggies).
    * Download AppImage zip file "piclas-binaries-vX.X.X" from the newest workflow run under [piclas-testing GitHub Workflows](https://github.com/scopplestone/piclas-testing/actions) (under "Artifacts").
    * Check correct (updated) version output and commit hash.
    * Check AppImage integrity by running `md5sum -c md5sum.txt`.
    * Select "New Pipeline" and set "Run for branch name or tag" to the required `feature.branch.name` and supply `DO_TEST_APPIMAGE` as "Input variable key" and set it to the latest artifact (get the name from  [piclas-testing GitHub Workflows](https://github.com/scopplestone/piclas-testing/actions), e.g., `piclas-binaries-v3.5.0` for "Input variable value".
* [ ] Draft the GitHub Release notes for creating the actual release that will be displayed under [Releases](https://github.com/piclas-framework/piclas/releases).

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
