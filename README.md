
          ___________________   _____    _________
         /  _____/\_   _____/  /     \  /   _____/
        /   \  ___ |    __)_  /  \ /  \ \_____  \ 
        \    \_\  \|        \/    Y    \/        \
         \______  /_______  /\____|__  /_______  /
                \/        \/         \/        \/ 
             is an Extensible Molecular Simulator
 

# Build and install

GEMS use either meson (recommended) or autotools build systems. 
To build with meson run:

    meson setup build
    meson install -C build

Run `meson dist` to create a distribution. Installation prefix, debug flags and
other variables can be modified following meson way. For instance:

    meson setup build --prefix=$PWD/usr --reconfigure && meson install -C build/
 
Instead, autotools can be used by running:

    autoreconf -fi; # Only if no configure script is given
    ./configure
    make

Run `make dist` to create a distribution. Further compiling options are handle
following autotools way. For instance:

    export PATH=/share/apps/gcc/6.2.0/bin/:$PATH
    export LD_LIBRARY_PATH=/share/apps/gcc/6.2.0/lib64:$LD_LIBRARY_PATH
    export FC=gfortran
    export FCFLAGS=-fno-use-linker-plugin
	./configure --disable-openmp FCFLAGS='-Ofast'

See `./configure --help` for more information.

## Dependencies

Mandatory dependencies are: 

- Fortran Preprocesor Templates for Dynamic Data Structures (FPT-DDS) 
  url: (https://github.com/alexispaz/FortranTemplates)

- Fortran 90 function parser v1.1
  url: (https://github.com/alexispaz/fparser)

Meson build system will automatically download and install these dependencies
if are no already present in the environment.

# About

GEMS code is hosted in [github](https://github.com/alexispaz/GEMS).

Copyright notices and license information for the different files used in the
GEMS project can be found in the ABOUT file that follows a *similar* format
to the [Debian's COPYRIGHT file format](https://www.debian.org/doc/packaging-manuals/copyright-format/1.0/).

Files in [lib/](lib/) are third party software distributed with none or quite
small changes. Licenses and copyright notices can be found in there as
presented by their authors.
 
## License

GEMS is an Extensible Molecular Simulator. 
Copyright (C) 2020  Sergio Alexis Paz

GEMS is free software: you can redistribute it and/or modify it under the terms
of the GNU General Public License as published by the Free Software Foundation,
   either version 3 of the License, or (at your option) any later version.

GEMS is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
GEMS. If not, see <https://www.gnu.org/licenses/>.

