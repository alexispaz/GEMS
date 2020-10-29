
# We are not following the standard layout of a GNU proyect
AUTOMAKE_OPTIONS = foreign

ACLOCAL_AMFLAGS = -I m4

EXTRA_DIST = COPYRIGHT COPYING README.md examples 

if SPRNG
  sprng_dir =  lib/sprng5
endif

# SUBDIRS = src doc examples scripts
SUBDIRS = lib/fparser $(sprng_dir) src

# DIST_SUBDIRS = ${SUBDIRS} lib/fpt-dds-a25f

# If SUBDIRS is defined conditionally using Automake conditionals, Automake
# will define DIST_SUBDIRS automatically from the possibles values of SUBDIRS
# in all conditions.
# If SUBDIRS contains AC_SUBST variables, DIST_SUBDIRS will not be defined
# correctly because Automake does not know the possible values of these
# variables. In this case DIST_SUBDIRS needs to be defined manually.
# DIST_SUBDIRS = 

# In order to prevent recursion in some non-configured directory you must
# therefore ensure that this directory does not appear in DIST_SUBDIRS (and
# SUBDIRS).


# dist_doc_DATA = README.md

# vimdir = $(datadir)/vim
# vim_DATA = gems.vim

