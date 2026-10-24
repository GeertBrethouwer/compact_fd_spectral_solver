# Accuracy tests of hichan's compact finite differences (src/compact.f90).
# HICHAN points to a checkout of https://github.com/GeertBrethouwer/hichan
HICHAN ?= ../hichan
FC      = gfortran
FFLAGS ?= -O2
SRC     = $(HICHAN)/src
B       = build

TESTS = test_global test_interior test_wall

all: $(addprefix $(B)/,$(TESTS)) $(B)/test_wall_quad

$(B)/%: tests/%.f90 tests/stub_param.f90 $(SRC)/compact.f90 $(SRC)/types.f90
	mkdir -p $(B)/mod_$*
	$(FC) $(FFLAGS) -J$(B)/mod_$* -o $@ $(SRC)/types.f90 tests/stub_param.f90 $(SRC)/compact.f90 $<

# same wall test in quad precision, to separate round-off from truncation error
$(B)/test_wall_quad: tests/test_wall.f90 tests/stub_param.f90 $(SRC)/compact.f90 $(SRC)/types.f90
	mkdir -p $(B)/mod_quad
	sed 's/rp = dp/rp = selected_real_kind(30)/' $(SRC)/types.f90 > $(B)/types_quad.f90
	sed 's/1.e-17_rp/1.e-34_rp/' $(SRC)/compact.f90 > $(B)/compact_quad.f90
	$(FC) $(FFLAGS) -J$(B)/mod_quad -o $@ $(B)/types_quad.f90 tests/stub_param.f90 $(B)/compact_quad.f90 $<

run: all
	@for t in $(TESTS) test_wall_quad; do echo "=== $$t"; ./$(B)/$$t; done

clean:
	rm -rf $(B)

.PHONY: all run clean
