#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_REVISION(Revision: 24 SEP 2012)

AC_INIT(nonadmd, 0.1.41-alpha, daniel.montemayor@ucd.ie)
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_HEADERS([config.h])
AM_INIT_AUTOMAKE

# Checks for programs.
AC_CHECK_PROGS([FPP],[fpp],[:])
if test "$FPP" = :; then
   AC_MSG_ERROR([This package needs fpp.])
fi
AC_CHECK_PROGS([PWD],[pwd],[:])
if test "$PWD" = :; then
   AC_MSG_ERROR([This package needs pwd.])
fi
AC_DEFINE_UNQUOTED([EXEDIR],["$PWD"],[Build directory.])
AC_SUBST([EXEDIR],[$PWD])

AC_DEFINE_UNQUOTED([EXAMPLESDIR],["$srcdir/share/examples/"],[Examples directory.])
AC_SUBST([EXAMPLESDIR],[$PWD])

# Check for source files
AC_CONFIG_SRCDIR([src/system_driver.f90])
AC_CONFIG_SRCDIR([src/quantum.f90])
AC_CONFIG_SRCDIR([src/classical.f90])
AC_CONFIG_SRCDIR([src/coupling.f90])
AC_CONFIG_SRCDIR([src/hs.f90])
AC_CONFIG_SRCDIR([src/hb.f90])
AC_CONFIG_SRCDIR([src/hc.f90])
AC_CONFIG_SRCDIR([src/subsystem1.in])
AC_CONFIG_SRCDIR([src/subsystem2.in])
AC_CONFIG_SRCDIR([src/subsystem3.in])
AC_CONFIG_SRCDIR([src/subsystem4.in])
AC_CONFIG_SRCDIR([src/subsystem5.in])
AC_CONFIG_SRCDIR([src/subsystem6.in])
AC_CONFIG_SRCDIR([src/subsystem7.in])
AC_CONFIG_SRCDIR([src/subsystem8.in])
AC_CONFIG_SRCDIR([src/subsystem9.in])
AC_CONFIG_SRCDIR([src/subsystemcap.in])
AC_CONFIG_SRCDIR([/etc/qsmods])
AC_CONFIG_SRCDIR([/etc/qstypes])
AC_CONFIG_SRCDIR([/etc/csmods])
AC_CONFIG_SRCDIR([/etc/cstypes])
AC_CONFIG_SRCDIR([/etc/cpmods])
AC_CONFIG_SRCDIR([/etc/cptypes])
AC_CONFIG_SRCDIR([/libexec/sandwich])
AC_CONFIG_SRCDIR([/libexec/sandwich2])
AC_CONFIG_SRCDIR([/libexec/sandwich3])
AC_CONFIG_SRCDIR([/libexec/sandwich4])
# Checks for libraries.

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.

AC_PROG_FC([ifort,mpif90,gfortran])
AC_PROG_CC
AC_PROG_F77

# Checks for library functions.

AC_CONFIG_FILES([Makefile src/Makefile])

# append quantum subsystem classes and create quantum.f90
AC_CONFIG_COMMANDS([quantum_subsystem],[sed -e "s/subsys/quantum/g" $ac_top_srcdir/src/subsystem1.in | sed -e "s/hprim/hs/g" | sed -e "s/suffix/qs/g" > $ac_top_srcdir/src/quantum.in
$ac_top_srcdir/libexec/sandwich '  use ' $ac_top_srcdir/etc/qsmods '_class' >> $ac_top_srcdir/src/quantum.in
sed -e "s/subsys/quantum/g" $ac_top_srcdir/src/subsystem2.in | sed -e "s/hprim/hs/g" | sed -e "s/suffix/qs/g" >> $ac_top_srcdir/src/quantum.in
$ac_top_srcdir/libexec/sandwich2 '     type(' $ac_top_srcdir/etc/qstypes ')::' ' '>> $ac_top_srcdir/src/quantum.in
sed -e "s/subsys/quantum/g" $ac_top_srcdir/src/subsystem3.in | sed -e "s/hprim/hs/g" | sed -e "s/suffix/qs/g" >> $ac_top_srcdir/src/quantum.in
$ac_top_srcdir/libexec/sandwich4 "    case('" $ac_top_srcdir/etc/qstypes "')\n       if(trim(infile).NE.'unknown')then;call new(this%" ",trim(infile));else;call new(this%" ");end if\n       this%hs => this%" "%hs" >> $ac_top_srcdir/src/quantum.in
sed -e "s/subsys/quantum/g" $ac_top_srcdir/src/subsystem4.in | sed -e "s/hprim/hs/g" | sed -e "s/suffix/qs/g" >> $ac_top_srcdir/src/quantum.in
$ac_top_srcdir/libexec/sandwich '    call kill(this%' $ac_top_srcdir/etc/qstypes ')' >> $ac_top_srcdir/src/quantum.in
sed -e "s/subsys/quantum/g" $ac_top_srcdir/src/subsystem5.in | sed -e "s/hprim/hs/g" | sed -e "s/suffix/qs/g" >> $ac_top_srcdir/src/quantum.in
$ac_top_srcdir/libexec/sandwich2 "    case('" $ac_top_srcdir/etc/qstypes "');call update(this%" ")" >> $ac_top_srcdir/src/quantum.in
sed -e "s/subsys/quantum/g" $ac_top_srcdir/src/subsystem6.in | sed -e "s/hprim/hs/g" | sed -e "s/suffix/qs/g" >> $ac_top_srcdir/src/quantum.in
$ac_top_srcdir/libexec/sandwich2 "    case('" $ac_top_srcdir/etc/qstypes "');call resample(this%" ")" >> $ac_top_srcdir/src/quantum.in
sed -e "s/subsys/quantum/g" $ac_top_srcdir/src/subsystem7.in | sed -e "s/hprim/hs/g" | sed -e "s/suffix/qs/g" >> $ac_top_srcdir/src/quantum.in

$ac_top_srcdir/libexec/sandwich "    case('" $ac_top_srcdir/etc/qstypes "')" >> $ac_top_srcdir/src/quantum.in
$ac_top_srcdir/libexec/sandwich "       if(withmsg)call display(this%" $ac_top_srcdir/etc/qstypes ",msg)" >> $ac_top_srcdir/src/quantum.in
$ac_top_srcdir/libexec/sandwich "       if(.not.withmsg)call display(this%" $ac_top_srcdir/etc/qstypes ")" >> $ac_top_srcdir/src/quantum.in

sed -e "s/subsys/quantum/g" $ac_top_srcdir/src/subsystem8.in | sed -e "s/hprim/hs/g" | sed -e "s/suffix/qs/g" >> $ac_top_srcdir/src/quantum.in
$ac_top_srcdir/libexec/sandwich2 "    case('" $ac_top_srcdir/etc/qstypes "');call save(this%" ",file//'.qs')" >> $ac_top_srcdir/src/quantum.in
sed -e "s/subsys/quantum/g" $ac_top_srcdir/src/subsystem9.in | sed -e "s/hprim/hs/g" | sed -e "s/suffix/qs/g" >> $ac_top_srcdir/src/quantum.in
$ac_top_srcdir/libexec/sandwich "    case('" $ac_top_srcdir/etc/qstypes "')" >> $ac_top_srcdir/src/quantum.in
sed -e "s/subsys/quantum/g" $ac_top_srcdir/src/subsystemcap.in | sed -e "s/hprim/hs/g" | sed -e "s/suffix/qs/g" >> $ac_top_srcdir/src/quantum.in
mv $ac_top_srcdir/src/quantum.in ./src/quantum.f90])

# append classical subsystem classes and create classical.f90
AC_CONFIG_COMMANDS([classical_subsystem],[sed -e "s/subsys/classical/g" $ac_top_srcdir/src/subsystem1.in | sed -e "s/hprim/hb/g" | sed -e "s/suffix/cs/g" > $ac_top_srcdir/src/classical.in
$ac_top_srcdir/libexec/sandwich '  use ' $ac_top_srcdir/etc/csmods '_class' >> $ac_top_srcdir/src/classical.in
sed -e "s/subsys/classical/g" $ac_top_srcdir/src/subsystem2.in | sed -e "s/hprim/hb/g" | sed -e "s/suffix/cs/g" >> $ac_top_srcdir/src/classical.in
$ac_top_srcdir/libexec/sandwich2 '     type(' $ac_top_srcdir/etc/cstypes ')::' ' '>> $ac_top_srcdir/src/classical.in
sed -e "s/subsys/classical/g" $ac_top_srcdir/src/subsystem3.in | sed -e "s/hprim/hb/g" | sed -e "s/suffix/cs/g" >> $ac_top_srcdir/src/classical.in
$ac_top_srcdir/libexec/sandwich4 "    case('" $ac_top_srcdir/etc/cstypes "')\n       if(trim(infile).NE.'unknown')then;call new(this%" ",trim(infile));else;call new(this%" ");end if\n       this%hb => this%" "%hb" >> $ac_top_srcdir/src/classical.in
sed -e "s/subsys/classical/g" $ac_top_srcdir/src/subsystem4.in | sed -e "s/hprim/hb/g" | sed -e "s/suffix/cs/g" >> $ac_top_srcdir/src/classical.in
$ac_top_srcdir/libexec/sandwich '    call kill(this%' $ac_top_srcdir/etc/cstypes ')' >> $ac_top_srcdir/src/classical.in
sed -e "s/subsys/classical/g" $ac_top_srcdir/src/subsystem5.in | sed -e "s/hprim/hb/g" | sed -e "s/suffix/cs/g" >> $ac_top_srcdir/src/classical.in
$ac_top_srcdir/libexec/sandwich2 "    case('" $ac_top_srcdir/etc/cstypes "');call update(this%" ")" >> $ac_top_srcdir/src/classical.in
sed -e "s/subsys/classical/g" $ac_top_srcdir/src/subsystem6.in | sed -e "s/hprim/hb/g" | sed -e "s/suffix/cs/g" >> $ac_top_srcdir/src/classical.in
$ac_top_srcdir/libexec/sandwich2 "    case('" $ac_top_srcdir/etc/cstypes "');call resample(this%" ")" >> $ac_top_srcdir/src/classical.in
sed -e "s/subsys/classical/g" $ac_top_srcdir/src/subsystem7.in | sed -e "s/hprim/hb/g" | sed -e "s/suffix/cs/g" >> $ac_top_srcdir/src/classical.in

$ac_top_srcdir/libexec/sandwich "    case('" $ac_top_srcdir/etc/cstypes "')" >> $ac_top_srcdir/src/classical.in
$ac_top_srcdir/libexec/sandwich "       if(withmsg)call display(this%" $ac_top_srcdir/etc/cstypes ",msg)" >> $ac_top_srcdir/src/classical.in
$ac_top_srcdir/libexec/sandwich "       if(.not.withmsg)call display(this%" $ac_top_srcdir/etc/cstypes ")" >> $ac_top_srcdir/src/classical.in

sed -e "s/subsys/classical/g" $ac_top_srcdir/src/subsystem8.in | sed -e "s/hprim/hb/g" | sed -e "s/suffix/cs/g" >> $ac_top_srcdir/src/classical.in
$ac_top_srcdir/libexec/sandwich2 "    case('" $ac_top_srcdir/etc/cstypes "');call save(this%" ",file//'.cs')" >> $ac_top_srcdir/src/classical.in
sed -e "s/subsys/classical/g" $ac_top_srcdir/src/subsystem9.in | sed -e "s/hprim/hb/g" | sed -e "s/suffix/cs/g" >> $ac_top_srcdir/src/classical.in
$ac_top_srcdir/libexec/sandwich "    case('" $ac_top_srcdir/etc/cstypes "')" >> $ac_top_srcdir/src/classical.in
sed -e "s/subsys/classical/g" $ac_top_srcdir/src/subsystemcap.in | sed -e "s/hprim/hb/g" | sed -e "s/suffix/cs/g" >> $ac_top_srcdir/src/classical.in
mv $ac_top_srcdir/src/classical.in ./src/classical.f90])

# append coupling subsystem classes and create coupling.f90
AC_CONFIG_COMMANDS([coupling_subsystem],[sed -e "s/subsys/coupling/g" $ac_top_srcdir/src/subsystem1.in | sed -e "s/hprim/hc/g" | sed -e "s/suffix/cp/g" > $ac_top_srcdir/src/coupling.in
$ac_top_srcdir/libexec/sandwich '  use ' $ac_top_srcdir/etc/cpmods '_class' >> $ac_top_srcdir/src/coupling.in
echo '  use quantum_class' >> $ac_top_srcdir/src/coupling.in 
echo '  use classical_class' >> $ac_top_srcdir/src/coupling.in 
sed -e "s/subsys/coupling/g" $ac_top_srcdir/src/subsystem2.in | sed -e "s/hprim/hc/g" | sed -e "s/suffix/cp/g" >> $ac_top_srcdir/src/coupling.in
$ac_top_srcdir/libexec/sandwich2 '     type(' $ac_top_srcdir/etc/cptypes ')::' ' '>> $ac_top_srcdir/src/coupling.in
sed -e "s/subsys/coupling/g" $ac_top_srcdir/src/coupling3.in | sed -e "s/hprim/hc/g" | sed -e "s/suffix/cp/g" >> $ac_top_srcdir/src/coupling.in
$ac_top_srcdir/libexec/sandwich4 "    case('" $ac_top_srcdir/etc/cptypes "')\n       if(trim(infile).NE.'unknown')then;call new(this%" ",qs,cs,trim(infile));else;call new(this%" ",qs,cs);end if\n       this%hc => this%" "%hc" >> $ac_top_srcdir/src/coupling.in
sed -e "s/subsys/coupling/g" $ac_top_srcdir/src/subsystem4.in | sed -e "s/hprim/hc/g" | sed -e "s/suffix/cp/g" >> $ac_top_srcdir/src/coupling.in
$ac_top_srcdir/libexec/sandwich '    call kill(this%' $ac_top_srcdir/etc/cptypes ')' >> $ac_top_srcdir/src/coupling.in
sed -e "s/subsys/coupling/g" $ac_top_srcdir/src/subsystem5.in | sed -e "s/hprim/hc/g" | sed -e "s/suffix/cp/g" >> $ac_top_srcdir/src/coupling.in
$ac_top_srcdir/libexec/sandwich2 "    case('" $ac_top_srcdir/etc/cptypes "');call update(this%" ")" >> $ac_top_srcdir/src/coupling.in
sed -e "s/subsys/coupling/g" $ac_top_srcdir/src/subsystem6.in | sed -e "s/hprim/hc/g" | sed -e "s/suffix/cp/g" >> $ac_top_srcdir/src/coupling.in
$ac_top_srcdir/libexec/sandwich2 "    case('" $ac_top_srcdir/etc/cptypes "');call resample(this%" ")" >> $ac_top_srcdir/src/coupling.in
sed -e "s/subsys/coupling/g" $ac_top_srcdir/src/subsystem7.in | sed -e "s/hprim/hc/g" | sed -e "s/suffix/cp/g" >> $ac_top_srcdir/src/coupling.in

$ac_top_srcdir/libexec/sandwich "    case('" $ac_top_srcdir/etc/cptypes "')" >> $ac_top_srcdir/src/coupling.in
$ac_top_srcdir/libexec/sandwich "       if(withmsg)call display(this%" $ac_top_srcdir/etc/cptypes ",msg)" >> $ac_top_srcdir/src/coupling.in
$ac_top_srcdir/libexec/sandwich "       if(.not.withmsg)call display(this%" $ac_top_srcdir/etc/cptypes ")" >> $ac_top_srcdir/src/coupling.in

sed -e "s/subsys/coupling/g" $ac_top_srcdir/src/subsystem8.in | sed -e "s/hprim/hc/g" | sed -e "s/suffix/cp/g" >> $ac_top_srcdir/src/coupling.in
$ac_top_srcdir/libexec/sandwich2 "    case('" $ac_top_srcdir/etc/cptypes "');call save(this%" ",file//'.cp')" >> $ac_top_srcdir/src/coupling.in
sed -e "s/subsys/coupling/g" $ac_top_srcdir/src/subsystem9.in | sed -e "s/hprim/hc/g" | sed -e "s/suffix/cp/g" >> $ac_top_srcdir/src/coupling.in
$ac_top_srcdir/libexec/sandwich "    case('" $ac_top_srcdir/etc/cptypes "')" >> $ac_top_srcdir/src/coupling.in
sed -e "s/subsys/coupling/g" $ac_top_srcdir/src/subsystemcap.in | sed -e "s/hprim/hc/g" | sed -e "s/suffix/cp/g" >> $ac_top_srcdir/src/coupling.in
mv $ac_top_srcdir/src/coupling.in ./src/coupling.f90])

# create system driver in build directory
#AC_CONFIG_COMMANDS([create_system_driver],[cp $ac_top_srcdir/src/system_driver_template.f90 src/system_driver.f90])
AC_CONFIG_COMMANDS([create_system_driver],[cp $ac_top_srcdir/src/system_driver.f90 src/system_driver.f90])

AC_OUTPUT
