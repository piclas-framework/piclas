# Building the AppImage Executable

## piclas

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

## piclas2vtk and other tools

Other tools such as `piclas2vtk` and `superB` etc. are also included in the AppImage container and can be extracted via

    ./piclas-ad6830c7a-x86_64.AppImage --appimage-extract

The tools are located under `./squashfs-root/usr/bin/`.
To make on those tools the main application of the AppImage, remove the AppDir folder

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

## Troubleshooting

If problems occur when executing the AppImage, check the [troubleshooting]( https://docs.appimage.org/user-guide/troubleshooting/index.html)
section for possible fixes.
