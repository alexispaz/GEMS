
          ___________________   _____    _________
         /  _____/\_   _____/  /     \  /   _____/
        /   \  ___ |    __)_  /  \ /  \ \_____  \ 
        \    \_\  \|        \/    Y    \/        \
         \______  /_______  /\____|__  /_______  /
                \/        \/         \/        \/ 
             is an Extensible Molecular Simulator
 

# Compilation 

GEMS uses autotools to configure the compilation. If you have download a
tarball of GEMS you can compile by running:

    ./configure
    make

Instead, if you want to clone from GitHub then:

    git clone git@github.com:alexispaz/GEMS.git
    cd GEMS
    autoreconf -fi
    ./configure
    make

Further compiling options are handle
in a standar way. For instance:

    export PATH=/share/apps/gcc/6.2.0/bin/:$PATH
    export LD_LIBRARY_PATH=/share/apps/gcc/6.2.0/lib64:$LD_LIBRARY_PATH
    export FC=gfortran
    export FCFLAGS=-fno-use-linker-plugin
	./configure --disable-openmp FCFLAGS='-Ofast'

See `./configure --help` for more information

To create a tar ball to run in other computers without autotools use:
  
    make dist

then you only need to configure and make (i.e. `./configure; make`).

## Git subtrees

GEMS uses the following `git subtrees`: 

- Fortran Preprocesor Templates for Dynamic Data Structures (FPT-DDS) 
  url: (https://github.com/alexispaz/FortranTemplates)

- Fortran 90 function parser v1.1
  url: (https://github.com/alexispaz/fparser)

The subtrees were created by:

	git remote add fpt git@github.com:alexispaz/FortranTemplates.git
    git subtree add --prefix lib/fpt fpt master --squash

	git remote add fparser git@github.com:alexispaz/fparser.git
    git subtree add --prefix lib/fparser fparser master --squash

After clone GEMS, to keep the remote relation with the sub projects you can
run:

	git remote add fpt git@github.com:alexispaz/FortranTemplates.git
	git remote add fparser git@github.com:alexispaz/fparser.git
	git fetch fpt
	git fetch fparser

This allows to pull sub projects updates like:

	git pull -s subtree fpt master

# About

GEMS code is hosted in_ [github](https://github.com/alexispaz/GEMS).

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

