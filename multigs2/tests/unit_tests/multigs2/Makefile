
sinclude ../Makefile.tests_common

CLEANCOMMAND+= test_multigs2

#Define this tests level
MY_TEST_LEVEL=1
VAL:= $(shell expr ${MY_TEST_LEVEL} \<= ${TEST_LEVEL})
ifeq "$(VAL)" "0"
all:
	@echo "==> Skipping test in ${PWD}/$(TEST_FOLD) (Test_level ${TEST_LEVEL}, This test ${MY_TEST_LEVEL})."
else
all: test_multigs2 
	$(call RUNCOMMAND,test_multigs2,test_multigs2.in)

endif

#$(info "TEST_MODS=$(TEST_MODS)")

test_multigs2.o:  ${TEST_MODS} 

test_multigs2: test_multigs2.o $(TEST_MODS)
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

