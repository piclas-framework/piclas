# Release & Deployment

## Release Procedure

A new release version of PICLas is created from the **master** repository, which requires a merge of the current **master.dev**
branch into the **master** branch. The corresponding merge request should be associated with a release milestone (e.g. *Release
1.X.X*). Within this specific merge request, the
template `Release` is chosen, which contains the to-do list as well as the template for the release notes.
After the successful completion of all to-do's and regression checks (check-in, nightly, weekly), the **master.dev** branch can be
merged into the **master** branch.

### Release Tag

A new release tag can be created through the web interface ([Repository -> Tags](https://piclas.boltzplatz.eu/piclas/piclas/tags) -> New tag)
and as the `Tag name`, the new version number is used, e.g.,

    v1.X.X

The tag is then created from the **master** branch repository and the `Message` is left empty. The release notes, which were created within the release merge request and are, shall be copied into the `Release notes` text box, not including the `# Release notes` heading.

### GitHub

Finally, the release tag can be deployed to GitHub. This can be achieved by running the `Deploy` script in the
[CI/CD -> Schedules](https://piclas.boltzplatz.eu/piclas/piclas/pipeline_schedules) web interface. At the moment, the respective tag and the
release have to be created manually on GitHub through the web interface with **piclas-framework** account. The releases are
accessed through [Releases](https://github.com/piclas-framework/piclas/releases) and a new release (including the tag) can be
created with `Draft a new release`. The tag version should be set as before (`v1.X.X`) and the release title accordingly
(`Release 1.X.X`). The release notes (not including the `# Release notes` heading) can be copied from the GitLab release.

## Automatic Deployment to GitHub

1. Add the required ssh key to the deploy keys on the respective platform (e.g. github)
1. Clone a code from the platform to update the list of known hosts. Do not forget to copy the
    information to the correct location for the runner to have access to the platform
    ```
    sudo cp~/.ssh/.ssh/known_hosts /var/lib/gitlab-runner/.ssh/known_hosts
    ```
    This might have to be performed via the gitlab-runner user, which can be accomplished by
    executing the following command
    ```
    sudo -u gitlab-runner git clone git@github.com:piclas-framework/piclas.git piclas_github
    ```
1. PICLas deployment is performed by the gitlab runner in the *deployment stage*
    ```
    github:
      stage: deploy
      tags:
        - withmodules-concurrent
      script:
        - if [ "$[[ inputs.DO_DEPLOY ]]" == "false" ]; then exit ; fi
        - rm -rf piclas_github || true ;
        - git clone -b master --single-branch git@piclas.boltzplatz.eu:piclas/piclas.git piclas_github ;
        - cd piclas_github ;
        - git remote add piclas-framework git@github.com:piclas-framework/piclas.git ;
        - git push --force --follow-tags piclas-framework master ;
    ```

## AppImage Executable

### piclas

Navigate to the piclas repository and create a build directory

    mkdir build && cd build

and compile piclas using the following cmake flags

    cmake .. -DCMAKE_INSTALL_PREFIX=/usr

and then

    make install DESTDIR=AppDir

or when using Ninja run

    DESTDIR=AppDir ninja install

Then create an AppImage (and subsequent paths) directory in the build folder

    mkdir -p AppDir/usr/share/icons/

and copy the piclas logo into the icons directory

    cp ../docs/logo.png AppDir/usr/share/icons/piclas.png

A desktop file should already exist in the top-level directory containing

    [Desktop Entry]
    Type=Application
    Name=piclas
    Exec=piclas
    Comment=PICLas is a flexible particle-based plasma simulation suite.
    Icon=piclas
    Categories=Development;
    Terminal=true

Next, download the AppImage executable

    curl -L -O https://github.com/linuxdeploy/linuxdeploy/releases/download/continuous/linuxdeploy-x86_64.AppImage

and make it executable

    chmod +x linuxdeploy-x86_64.AppImage

Then run

    ./linuxdeploy-x86_64.AppImage --appdir AppDir --output appimage --desktop-file=../.github/workflows/piclas.desktop

The executable should be created in the top-level directory, e.g.,

    piclas-ad6830c7a-x86_64.AppImage

### piclas2vtk and other tools

Other tools such as `piclas2vtk` and `superB` etc. are also included in the AppImage container and can be extracted via

    ./piclas-ad6830c7a-x86_64.AppImage --appimage-extract

The tools are located under `./squashfs-root/usr/bin/`.
To make one of those tools the main application of the AppImage, remove the AppDir folder

    rm -rf AppDir

and then

    make install DESTDIR=AppDir

or when using Ninja run

    DESTDIR=AppDir ninja install

and change the following settings, e.g., for `piclas2vtk`

    PROG='piclas2vtk'
    cp ../.github/workflows/piclas.desktop ${PROG}.desktop
    mkdir -p AppDir/usr/share/icons/
    cp ../docs/logo.png AppDir/usr/share/icons/${PROG}.png
    sed -i -e "s/Name=.*/Name=${PROG}/" ${PROG}.desktop
    sed -i -e "s/Exec=.*/Exec=${PROG}/" ${PROG}.desktop
    sed -i -e "s/Icon=.*/Icon=${PROG}/" ${PROG}.desktop
    ./linuxdeploy-x86_64.AppImage --appdir AppDir --output appimage --desktop-file=${PROG}.desktop

This should create

    piclas2vtk-ad6830c7a-x86_64.AppImage

### Troubleshooting

If problems occur when executing the AppImage, check the [troubleshooting]( https://docs.appimage.org/user-guide/troubleshooting/index.html)
section for possible fixes.
