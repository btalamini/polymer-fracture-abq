test_uel: test_uel.f90 uel_phase_field.f
	gfortran test_uel.f90 uel_phase_field.f -ffpe-trap=invalid,zero,overflow -cpp -g -Wall -Wextra -Wno-tabs -o test_uel

test_pl_stress: test_uel_polymer_plane_stress.f90 uel_phase_field_polymer_plane_stress.f
	gfortran test_uel_polymer_plane_stress.f90 uel_phase_field_polymer_plane_stress.f -ffpe-trap=invalid,zero,overflow -cpp -g -Wall -Wextra -Wno-tabs -o test_pl_stress
