SUBDIRS = lfg lcg lcg64 cmrg mlfg pmlcg .

bin_SCRIPTS = checksprng timesprng 
CLEANFILES = $(bin_SCRIPTS)
dist_data_DATA = insertlib gencheck gentime .insertlib
EXTRA_DIST = insertlib gencheck gentime .insertlib
CHKDIR = $(top_srcdir)/check

checksprng: gencheck
	rm -f checksprng
	chmod u+x ./gencheck
	@./gencheck
	chmod u+x ./checksprng

timesprng: gentime
	rm -f timesprng
	chmod u+x ./gentime
	@./gentime
	chmod u+x ./timesprng
