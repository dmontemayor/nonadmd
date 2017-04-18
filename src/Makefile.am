bin_PROGRAMS = App

App_FCFLAGS = -O3 -$(FPP) -I $(EXEDIR)
#App_FCFLAGS = -O3 -I $(EXEDIR)

corefiles = utils.f90 subroutines.f90

qsfiles = dipoles.f90 #MorseOscillator.f90

csfiles = harmonicbath.f90

cpfiles = bilinear.f90 #DiagonalCaldeiraLeggett.f90 #DMBLcoupling.f90 #CaldeiraLeggett.f90

thefiles = PLDM.f90 #islandmapserial.f90 #verlet.f90 #LDM.f90

tools = molreader.f90 spectrometer.f90 praxis.f

#custom = custom/JSB850.f90 custom/JSB800.f90 custom/JSLH2.f90 custom/LH2.f90

testfiles = #tests.f90

#App_SOURCES = $(corefiles) hs.f90 $(qsfiles) quantum.f90 hb.f90 $(csfiles) classical.f90 hc.f90 $(cpfiles) coupling.f90 $(thefiles) $(tools) $(custom) $(testfiles) system_driver.f90

App_SOURCES = $(corefiles) hs.f90 $(qsfiles) quantum.f90 hb.f90 $(csfiles) classical.f90 hc.f90 $(cpfiles) coupling.f90 $(thefiles) $(tools) $(testfiles) system_driver.f90

EXTRA_DIST = subsystem1.in subsystem2.in subsystem3.in subsystem4.in subsystem5.in subsystem6.in subsystem7.in subsystem8.in subsystem9.in subsystemcap.in coupling3.in TEMPLATE.f90

