# PLEASE DO NOT USE THESE FUNCTIONS; THEY HAVE BEEN DEPRECATED
import indexedexp as ixp
import reference_metric as rfm
import outputC as outC
import NRPy_param_funcs as par
import os

# Find the appropriate timestep for the CFL condition.
def add_to_Cfunction_dict_find_timestep():
    # Compute proper distance in all 3 directions.
    delxx = ixp.declarerank1("dxx", DIM=3)
    ds_drn = rfm.ds_dirn(delxx)

    ds_dirn_h = outC.outputC([ds_drn[0], ds_drn[1], ds_drn[2]], ["ds_dirn0", "ds_dirn1", "ds_dirn2"],"returnstring")

    desc="Find the CFL-constrained timestep"
    outC.add_to_Cfunction_dict(
        desc     =desc,
        c_type   ="REAL",
        name     ="find_timestep",
        params   ="const paramstruct *restrict params, REAL *restrict xx[3]",
        preloop  ="REAL dsmin = 1e38; // Start with a crazy high value... close to the largest number in single precision.",
        body     ="REAL ds_dirn0, ds_dirn1, ds_dirn2;\n"+ds_dirn_h+"""
#ifndef MIN
#define MIN(A, B) ( ((A) < (B)) ? (A) : (B) )
#endif
        // Set dsmin = MIN(dsmin, ds_dirn0, ds_dirn1, ds_dirn2);
        dsmin = MIN(dsmin,MIN(ds_dirn0,MIN(ds_dirn1,ds_dirn2)));
""",
        loopopts ="InteriorPoints,Read_xxs,DisableOpenMP",
        postloop ="return dsmin*CFL_FACTOR/wavespeed;\n")


def out_timestep_func_to_file(outfile):
    add_to_Cfunction_dict_find_timestep()
    with open(outfile, "w") as file:
        file.write(outC.outC_function_dict["find_timestep"])


############################
## TO BE DEPRECATED:
def set_Nxx_dxx_invdx_params__and__xx_h(outdir=".",grid_centering="cell"):
    if grid_centering not in ('cell', 'vertex'):
        print("rfm.set_Nxx_dxx_invdx_params__and__xx_h(): grid_centering = \""+grid_centering+"\" not supported!")
        sys.exit(1)

    with open(os.path.join(outdir,"set_Nxx_dxx_invdx_params__and__xx.h"),"w") as file:
        file.write(r"""
void set_Nxx_dxx_invdx_params__and__xx(const int EigenCoord, const int Nxx[3],
                                       paramstruct *restrict params, REAL *restrict xx[3]) {
    // Override parameter defaults with values based on command line arguments and NGHOSTS.
    params->Nxx0 = Nxx[0];
    params->Nxx1 = Nxx[1];
    params->Nxx2 = Nxx[2];
    params->Nxx_plus_2NGHOSTS0 = Nxx[0] + 2*NGHOSTS;
    params->Nxx_plus_2NGHOSTS1 = Nxx[1] + 2*NGHOSTS;
    params->Nxx_plus_2NGHOSTS2 = Nxx[2] + 2*NGHOSTS;
    // Step 0d: Set up space and time coordinates
    // Step 0d.i: Declare \Delta x^i=dxx{0,1,2} and invdxx{0,1,2}, as well as xxmin[3] and xxmax[3]:
#include "set_Cparameters.h"
    REAL xxmin[3],xxmax[3];
    if(EigenCoord == 0) {
""")
        for i in range(3):
            file.write("        xxmin["+str(i)+"] = "+str(rfm.xxmin[i])+";\n")
            file.write("        xxmax["+str(i)+"] = "+str(rfm.xxmax[i])+";\n")
        file.write("""
    } else { // if (EigenCoord == 1)
""")
        CoordSystem_orig = par.parval_from_str("reference_metric::CoordSystem")
        par.set_parval_from_str("reference_metric::CoordSystem",rfm.get_EigenCoord())
        rfm.reference_metric()
        for i in range(3):
            file.write("        xxmin["+str(i)+"] = "+str(rfm.xxmin[i])+";\n")
            file.write("        xxmax["+str(i)+"] = "+str(rfm.xxmax[i])+";\n")
        par.set_parval_from_str("reference_metric::CoordSystem",CoordSystem_orig)
        rfm.reference_metric()
        file.write("""
    }

    params->dxx0 = (xxmax[0] - xxmin[0]) / ((REAL)Nxx[0]);
    params->dxx1 = (xxmax[1] - xxmin[1]) / ((REAL)Nxx[1]);
    params->dxx2 = (xxmax[2] - xxmin[2]) / ((REAL)Nxx[2]);
    params->invdx0 = 1.0/params->dxx0;
    params->invdx1 = 1.0/params->dxx1;
    params->invdx2 = 1.0/params->dxx2;\n""")
        # The following capability was suggested by Terrence Pierre Jacques (Thanks!)
        cell_offset = "(1.0/2.0)" # default cell-centered offset
        cell_comment = "Cell-centered grid."
        if grid_centering == "vertex":
            cell_offset = "0.0"
            cell_comment = "Vertex-centered grid."
        file.write("""
    // Now that params.dxx{0,1,2} and params.invdxx{0,1,2} have been set,
    // Step 0d.iii: Set up uniform coordinate grids
    xx[0] = (REAL *)malloc(sizeof(REAL)*Nxx_plus_2NGHOSTS0);
    for(int j=0;j<Nxx_plus_2NGHOSTS0;j++)
        xx[0][j] = xxmin[0] + ((REAL)(j-NGHOSTS) + """+cell_offset+""")*params->dxx0; // """+cell_comment+"""
    xx[1] = (REAL *)malloc(sizeof(REAL)*Nxx_plus_2NGHOSTS1);
    for(int j=0;j<Nxx_plus_2NGHOSTS1;j++)
        xx[1][j] = xxmin[1] + ((REAL)(j-NGHOSTS) + """+cell_offset+""")*params->dxx1; // """+cell_comment+"""
    xx[2] = (REAL *)malloc(sizeof(REAL)*Nxx_plus_2NGHOSTS2);
    for(int j=0;j<Nxx_plus_2NGHOSTS2;j++)
        xx[2][j] = xxmin[2] + ((REAL)(j-NGHOSTS) + """+cell_offset+""")*params->dxx2; // """+cell_comment+"""
    //fprintf(stderr,"hey inside setxx: %e %e %e | %e %e\\n",xxmin[0],xxmin[1],xxmin[2],xx[0][0],params->dxx0);
}
""")

def xx_to_Cart_h(funcname,cparamsloc,outfile):
    # Arbitrary-coordinate NRPy+ file output, Part 1: output the conversion from (x0,x1,x2) to Cartesian (x,y,z)
    Cout = outC.outputC([rfm.xx_to_Cart[0],rfm.xx_to_Cart[1],rfm.xx_to_Cart[2]],
                        ["xCart[0]","xCart[1]","xCart[2]"],
                        "returnstring",params="preindent=1")

    with open(outfile, "w") as file:
        file.write("""
static inline void """+funcname+"""(const paramstruct *restrict params, REAL *restrict xx[3],const int i0,const int i1,const int i2, REAL xCart[3]) {
#include """+"\""+cparamsloc+"\""+"""
    REAL xx0 = xx[0][i0];
    REAL xx1 = xx[1][i1];
    REAL xx2 = xx[2][i2];\n"""+Cout+"}\n")
