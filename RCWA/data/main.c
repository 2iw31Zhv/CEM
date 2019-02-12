#include <S4.h>
#include <math.h>


int main(int argc, char *argv[]){
	double freq = 1;
	double hole_radius = 0.2;

	if(argc > 2){
		freq = atof(argv[1]);
		hole_radius = atof(argv[2]);
	}
	else if(argc > 1){
		freq = atof(argv[1]);
	}

	double eps[2];
	Simulation *S = (Simulation*)malloc(sizeof(Simulation));
	Simulation_Init(S);
	S->Lr[0] = 1;
	S->Lr[1] = 0;
	S->Lr[2] = 0;
	S->Lr[3] = 1;
	//S->options.use_polarization_basis = 1;
	Simulation_MakeReciprocalLattice(S);
	Simulation_SetNumG(S, 100);

	Material *M_air = Simulation_AddMaterial(S);
	eps[0] = 1; eps[1] = 0;
	Material_Init(M_air, "air", eps);
	Material *M_Si = Simulation_AddMaterial(S);
	eps[0] = 12.25; eps[1] = 0;
	Material_Init(M_Si, "Si", eps);

	Layer *L1 = Simulation_AddLayer(S);
	Layer_Init(L1, "above", 0, "air", NULL);

	Layer *L2 = Simulation_AddLayer(S);
	Layer_Init(L2, "slab", 0.5, "Si", NULL);

	double center[2] = {0,0};
	Simulation_AddLayerPatternCircle(S, L2, 0, center, hole_radius);

	Layer *L3 = Simulation_AddLayer(S);
	Layer_Init(L3, "below", 0, NULL, "above");

	double angle[2], pol_s[2], pol_p[2];
	angle[0] = 0;
	angle[1] = 0;
	pol_s[0] = 1; pol_s[1] = 0;
	pol_p[0] = 0; pol_p[1] = 0;
	Simulation_MakeExcitationPlanewave(S, angle, pol_s, pol_p, 0);
	S->omega[0] = 2*M_PI*freq;
	S->omega[1] = 0;

	double below_powers[4];
	//powers = {forward_real, backward_real, forward_imag, backward_imag};
	
	Simulation_GetPoyntingFlux(S, L3, 0, below_powers);
	printf("%g\t%g\n", freq, below_powers[0]);

	// int n4 = 4 * S->n_G;

	// std::complex<double> * S01 = (std::complex<double>*)malloc(sizeof(std::complex<double>) * n4 * n4);
	// std::complex<double> * S12 = (std::complex<double>*)malloc(sizeof(std::complex<double>) * n4 * n4);
	// std::complex<double> * S23 = (std::complex<double>*)malloc(sizeof(std::complex<double>) * n4 * n4);
	// std::complex<double> * S30 = (std::complex<double>*)malloc(sizeof(std::complex<double>) * n4 * n4);

	// Simulation_GetSMatrix(S, 0, 1, S01);
	// Simulation_GetSMatrix(S, 1, 2, S12);
	// Simulation_GetSMatrix(S, 2, 3, S23);
	// //Simulation_GetSMatrix(S, 3, -1, S23);

	// printf("S01:\n");
	// for (int i = 0; i < n4; ++i)
	// {
	// 	for (int j = 0; j < n4; ++j)
	// 	{
	// 		printf("%g+%gj\t", S01[i*n4 + j].real(), S01[i*n4 + j].imag());
	// 	}
	// 	printf(";\n");
	// }

	// printf("S12:\n");
	// for (int i = 0; i < n4; ++i)
	// {
	// 	for (int j = 0; j < n4; ++j)
	// 	{
	// 		printf("%g+%gj\t", S12[i*n4 + j].real(), S12[i*n4 + j].imag());
	// 	}
	// 	printf(";\n");
	// }


	// printf("S23:\n");
	// for (int i = 0; i < n4; ++i)
	// {
	// 	for (int j = 0; j < n4; ++j)
	// 	{
	// 		printf("%g+%gj\t", S23[i*n4 + j].real(), S23[i*n4 + j].imag());
	// 	}
	// 	printf(";\n");
	// }

	// printf("S30:\n");
	// for (int i = 0; i < n4; ++i)
	// {
	// 	for (int j = 0; j < n4; ++j)
	// 	{
	// 		printf("%g+%gj\t", S30[i*n4 + j].real(), S30[i*n4 + j].imag());
	// 	}
	// 	printf(";\n");
	// }

	// free(S01);
	// free(S12);
	// free(S23);
	// free(S30);

	Simulation_Destroy(S);
	free(S);
	return 0;
}
