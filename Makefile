# Accuracy tests of hichan's compact finite differences (src/compact.f90).
# HICHAN: the copy of https://github.com/GeertBrethouwer/hichan in hichan/ (or another checkout)
HICHAN ?= hichan
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

# D1, D1m vs summation-by-parts operators (needs python3 with numpy and scipy)
PYTHON ?= python3
$(B)/ops/done: $(B)/dump_ops
	mkdir -p $(B)/ops && cd $(B)/ops && ../dump_ops && touch done

sbp: $(B)/ops/done
	$(PYTHON) tests/compare_sbp.py $(B)/ops

# Corbino-Castillo mimetic operators; reads the coefficient tables from MOLE (GPL-3.0) at run time
mimetic: $(B)/ops/done
	$(PYTHON) tests/compare_mimetic.py $(B)/ops $(B)/mole

clean:
	rm -rf $(B)

.PHONY: all run sbp mimetic clean
