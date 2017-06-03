sudo mkdir -p /usr/local/include
sudo ln -s /usr/X11/include/freetype2/freetype /usr/local/include/freetype
sudo ln -s /usr/X11/include/ft2build.h /usr/local/include/ft2build.h
sudo ln -s /usr/X11/include/png.h /usr/local/include/png.h
sudo ln -s /usr/X11/include/pngconf.h /usr/local/include/pngconf.h
sudo ln -s /usr/X11/include/pnglibconf.h /usr/local/include/pnglibconf.h
sudo mkdir -p /usr/local/lib
sudo ln -s /usr/X11/lib/libfreetype.dylib /usr/local/lib/libfreetype.dylib
sudo ln -s /usr/X11/lib/libpng.dylib /usr/local/lib/libpng.dylib

