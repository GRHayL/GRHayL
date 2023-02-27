#include "./NRPy_basic_defines.h"
#include "./NRPy_function_prototypes.h"
/*
 * Main initial data function.
 */
void initial_data(const char *initial_data_option, const paramstruct *restrict params, REAL *restrict xx[3], REAL *restrict evol_gfs,REAL *restrict auxevol_gfs) {
#include "./set_Cparameters.h"



    const char *option1 = "b1";
    const char *option2 = "b2";
    const char *option3 = "b3";
    const char *option4 = "b4";
    const char *option5 = "b5";
    const char *option6 = "ce";
    const char *option7 = "mr";
    const char *option8 = "la";

    if (strcmp(initial_data_option, option1) == 0) {
        Balsara1(params, xx, evol_gfs, auxevol_gfs);
    }

    else if (strcmp(initial_data_option, option2) == 0) {
        Balsara2(params, xx, evol_gfs, auxevol_gfs);
    }
    
    else if (strcmp(initial_data_option, option3) == 0) {
        Balsara3(params, xx, evol_gfs, auxevol_gfs);
    }
    
    else if (strcmp(initial_data_option, option4) == 0) {
        Balsara4(params, xx, evol_gfs, auxevol_gfs);
    }
    
    else if (strcmp(initial_data_option, option5) == 0) {
        Balsara5(params, xx, evol_gfs, auxevol_gfs);
    }
    
    else if (strcmp(initial_data_option, option6) == 0) {
        Cylindrical_Explosion(params, xx, evol_gfs, auxevol_gfs);
    }
    
    else if (strcmp(initial_data_option, option7) == 0) {
        Magnetic_Rotor(params, xx, evol_gfs, auxevol_gfs);
    }
    
    else if (strcmp(initial_data_option, option8) == 0) {
        Loop_Advection(params, xx, evol_gfs, auxevol_gfs);
    }    

    else {
        printf("ERROR: Invalid choice of initial data.\n");
        exit(1);
    }

    }
