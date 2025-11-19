# Merge Request To-Do

## 1. Update version number
* [ ] Update the piclas version numbers using the [updatePiclasVersion.sh](https://github.com/piclas-framework/piclas/blob/master/tools/updatePiclasVersion.sh) tool by running
    ```
    ./tools/updatePiclasVersion.sh 3.2.1
    ```
  in the top level directory. Adjust the version number argument to the current. Check if the following two files have been correctly changed/updated:
  * [ ] Check that `./src/globals/globals_vars.f90` (`MajorVersion`, `MinorVersion` and `PatchVersion`) is updated.
  * [ ] Check that `.github/workflows/cmake-ninja.yml` (`piclas-binaries-vX.X.X`) is updated.

## 2. Prerequisites
* [ ] Update the prerequisite table for compiling piclas under [Required Libraries](https://piclas.readthedocs.io/en/latest/userguide/installation.html#required-libraries) by checking the versions that are currently used by the reggie server (see gitlab CI/CD) for mpich and OpenMPI. For example, see Release 3.3.0 pipelines for [openmpi](https://piclas.boltzplatz.eu/piclas/piclas/-/jobs/675270) and [mpich](https://piclas.boltzplatz.eu/piclas/piclas/-/jobs/674955). Check the most recent pipelines running for `master.dev` and use these versions in the table.

## 3. AppImage
Update the library versions and push the most recent version of the `master.dev` repository to the testing GitHub repository and make sure that the AppImage compilation process succeeds before doing the real GitHub release.
  * [ ] Update the GCC/MPI/HDF5/PETSc library versions used in the GitHub workflow by adjusting the versions in [cmake-ninja.yml](https://piclas.boltzplatz.eu/piclas/piclas/-/blob/master.dev/.github/workflows/cmake-ninja.yml) with which the AppImage is built.
  * [ ] Update the information in the user guide regarding the versions of the glibc (currently the OS version is used) and OpenMPI dependencies in the [AppImage dependeny table](https://piclas.readthedocs.io/en/latest/userguide/installation.html#appimage-executable-download) that are required for the AppImage. Check the following examples and linked commits on how this is done:
    *  Release 1.0.0 - 3.3.0: glibc 2.17 + OpenMPI X.X.X
    *  Release 3.3.0 - X.X.X: glibc 2.18 [9b09c795](https://piclas.boltzplatz.eu/piclas/piclas/-/commit/9b09c7957800915cbdf5ecc4a0d8ba43993060da) + OpenMPI X.X.X
  * [ ] Push the feature branch `feature.branch.name`, which is usually the `master.dev` branch, with the changed `cmake-ninja.yml` file that features an updated version of the MPI/HDF5/PETSc libraries for building the AppImage to the testing repository via gitlab Pipelines. Select [New Pipeline](https://piclas.boltzplatz.eu/piclas/piclas/-/pipelines/new) and set "Run for branch name or tag" to the required `feature.branch.name` and supply `DO_CREATE_APPIMAGE` as "Input variable key" and set it to `T` for "Input variable value".
    ```
    DO_CREATE_APPIMAGE = T
    ```
    * [ ] Check the status of the pipeline, which pushes the`feature.branch.name` branch to the GitHub repository [piclas-testing](https://github.com/scopplestone/piclas-testing) in the final gitlab stage. The AppImage is built on GitHub.
    * [ ] When the previous step is completed, go to the [piclas-testing GitHub Workflows](https://github.com/scopplestone/piclas-testing/actions) page and check if the build is running. Fix any errors that might occur during the build process and repeat the process. Continue with the next step only and as soon as the AppImage has been built using the GitHub actions.
  * When the workflow action pipeline has successfully created the AppImages, test their functionality by running a series of tests with the executables (check-in reggies). The previous step must be completed successfully before the following steps can be started.
    - [ ] Download AppImage zip file "piclas-binaries-vX.X.X" from the newest workflow run under
          [piclas-testing GitHub Workflows](https://github.com/scopplestone/piclas-testing/actions) (under "Artifacts").
    - [ ] Run one of the AppImage executables and check if the updated *version number* and *commit hash* are correct.
    - [ ] Check AppImage integrity by running `md5sum -c md5sum.txt` in the directory where the AppImage executables and the
          md5sum.txt have been extracted.
    - [ ] Select [New Pipeline](https://piclas.boltzplatz.eu/piclas/piclas/-/pipelines/new) and set "Run for branch name or tag" to the required `feature.branch.name` and supply `DO_TEST_APPIMAGE` as "Input variable key" and set it to the latest artifact (get the name "piclas-binaries-vX.X.X" from [piclas-testing GitHub Workflows](https://github.com/scopplestone/piclas-testing/actions), e.g., `piclas-binaries-v3.5.0` for "Input variable value".

    ```
    DO_TEST_APPIMAGE = piclas-binaries-vX.X.X
    ```

## 4. Draft GitHub Release Notes
Draft the GitHub Release notes for creating the official release
* [ ] Press the button **Draft a new Release** under [Releases](https://github.com/piclas-framework/piclas/releases) and leave the *Select tag* empty as the repository is pushed
      to GitHub in the following. Copy the release notes from below. Do not publish the release yet, instead, press the *Save draft* button to store the draft.

## 5. Regression Testing

* [ ] Check-in
* [ ] Nightly
* [ ] Weekly

## 6. Complete MR and create Gitlab Release
 - [ ] Merge this MR by pushing the `Merge` button
 - [ ] Create a new gilab release and tag vX.X.X (using the master branch) under
   [Releases](https://piclas.boltzplatz.eu/piclas/piclas/-/releases/new) (note the required permissions to be able to select the
   master branch) by copying the release notes from below
 - [ ] Check that the new tag is available under [Tags](https://piclas.boltzplatz.eu/piclas/piclas/-/tags) using the format vX.X.X

## 7. Deploy master branch to GitHub
- [ ] When the `Merge` button for this merge request is pushed, start the deployment of the *master* branch to GitHub via
      [New Pipeline](https://piclas.boltzplatz.eu/piclas/piclas/-/pipelines/new) and set the pipeline variable

```
DO_DEPLOY = T
```
- [ ] Check that the code base under [GitHub piclas](https://github.com/piclas-framework/piclas) has been updated

## 8. Create tag on release on GitHub
- [ ] Go to the release draft under [Releases](https://github.com/piclas-framework/piclas/releases) and create a new tag using the
  format vX.X.X and [use the same name as before in gitlab](https://piclas.boltzplatz.eu/piclas/piclas/-/tags)
- [ ] Publish the GitHub Release Notes (the commit hashes in the release notes should now also be hyperlinks in GitHub)

# Release Notes

## Release 3.X.X

### Breaking/Parameter Changes

*

### Documentation/Tools/Regression Testing

*

### Features

*

### Improvements

*

### Fixes

*