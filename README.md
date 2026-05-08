
          ___________________   _____    _________
         /  _____/\_   _____/  /     \  /   _____/
        /   \  ___ |    __)_  /  \ /  \ \_____  \ 
        \    \_\  \|        \/    Y    \/        \
         \______  /_______  /\____|__  /_______  /
                \/        \/         \/        \/ 
             is an Extensible Molecular Simulator
 

# Build and install

First downlaod GEMS from repo:
 
    git clone https://github.com/alexispaz/gems.git
    cd gems  

Then, use either meson (recommended) or autotools systems to build an install.
Optionally, these can be also handle with spack.

## Using spack

You can add the package from the local repo:

    cd spack/
    spack repo add .

Then install and activate the package in your preferred spack way.

    spack install gems
    spack load gems

## Using meson

Configure, build an install

    meson setup build
    meson install -C build

Then, run `meson dist` to create a distribution. Installation prefix, debug
flags and other variables can be modified during configuration. For instance:

    meson setup build --prefix=$PWD/usr --reconfigure && meson install -C build/

All together with debug flags would be something like:
  
    meson setup build --prefix=$(realpath ./usr) --buildtype=debug \
    -Df_args='-fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace -g3 -ggdb3 -Wextra' \
    && meson install -C build/  

## Using autotools

Configure, build an install


    autoreconf -fi; # Only if no configure script is given
    ./configure
    make
    make install

Run `make dist` to create a distribution. Further compiling options are handle
following autotools way. For instance:

    export PATH=/share/apps/gcc/6.2.0/bin/:$PATH
    export LD_LIBRARY_PATH=/share/apps/gcc/6.2.0/lib64:$LD_LIBRARY_PATH
    export FC=gfortran
    export FCFLAGS=-fno-use-linker-plugin
	./configure --disable-openmp FCFLAGS='-Ofast'

See `./configure --help` for more information.

_Warning: If you switch between Autotools and Meson, remember that Autotools
places compiled .mod files inside the src/ directory, while Meson expects the
source tree to be clean. Leftover .mod files in src/ can cause incorrect
dependency resolution and misleading compilation errors._

## Dependencies

Mandatory dependencies are: 

- Fortran Preprocesor Templates for Dynamic Data Structures (FPT-DDS) 
  url: (https://github.com/alexispaz/FortranTemplates)

- Fortran 90 function parser v1.1
  url: (https://github.com/alexispaz/fparser)

Meson build system will automatically download and install these dependencies
if are no already present in the environment.

## Documentation

- [Architecture](docs/ARCHITECTURE.md)
- [Memory model](docs/MEMORY_MODEL.md)
- [Input language](docs/INPUT_LANGUAGE.md)
- [CLI reference](docs/CLI_REFERENCE.md)
- [Developer guide](docs/DEVELOPER_GUIDE.md)

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
