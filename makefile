test_uel: test_uel.f90 uel_phase_field.f
	gfortran test_uel.f90 uel_phase_field.f -ffpe-trap=invalid,zero,overflow -cpp -g -Wall -Wextra -Wno-tabs -o test_uel
