
sinclude ../Makefile.tests_common

CLEANCOMMAND+= wstar_units

#Define this tests level
MY_TEST_LEVEL=1
VAL:= $(shell expr ${MY_TEST_LEVEL} \<= ${TEST_LEVEL})
ifeq "$(VAL)" "0"
all:
	@echo "==> Skipping test in ${PWD}/$(TEST_FOLD) (Test_level ${TEST_LEVEL}, This test ${MY_TEST_LEVEL})."
else
all: wstar_units 
	$(call RUNCOMMAND,wstar_units,wstar_units.list)
endif

wstar_units.o:  ${TEST_MODS} 

wstar_units: wstar_units.o $(TEST_MODS)
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

