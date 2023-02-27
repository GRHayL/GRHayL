# As documented in the NRPy+ tutorial module
#   Tutorial-RK_Butcher_Table_Generating_C_Code_new_way.ipynb,
#   this module will produce the required C codes for
#   allocating required memory Method of Lines (MoL) timestepping,
#   implementing MoL timestepping, and deallocating memory

# Authors: Brandon Clark
#          Zachariah B. Etienne (maintainer)
#          zachetie **at** gmail **dot* com

import sympy as sp  # Import SymPy, a computer algebra system written entirely in Python
import os, sys      # Standard Python modules for multiplatform OS-level functions
from MoLtimestepping.RK_Butcher_Table_Dictionary import Butcher_dict
from outputC import add_to_Cfunction_dict, indent_Ccode, outC_NRPy_basic_defines_h_dict, outputC, superfast_uniq  # NRPy+: Basic C code output functionality


# Check if Butcher Table is diagonal
def diagonal(key):
    Butcher = Butcher_dict[key][0]
    L = len(Butcher)-1  # Establish the number of rows to check for diagonal trait, all bust last row
    row_idx = 0  # Initialize the Butcher table row index
    for i in range(L):  # Check all the desired rows
        for j in range(1, row_idx):  # Check each element before the diagonal element in a row
            if Butcher[i][j] != sp.sympify(0):  # If any non-diagonal coeffcient is non-zero,
                                                # then the table is not diagonal
                return False
        row_idx += 1  # Update to check the next row
    return True


# Each MoL method has its own set of names for groups of gridfunctions,
#   aiming to be sufficiently descriptive. So for example a set of
#   gridfunctions that store "k_1" in an RK-like method could be called
#   "k1_gfs".
def generate_gridfunction_names(MoL_method = "RK4"):
    # Step 3.a: MoL gridfunctions fall into 3 overlapping categories:
    #           1) y_n=y_i(t_n) gridfunctions y_n_gfs, which stores data for the vector of gridfunctions y_i at t_n,
    #              the start of each MoL timestep.
    #           2) non-y_n gridfunctions, needed to compute the data at t_{n+1}. Often labeled with k_i in the name,
    #              these gridfunctions are *not* needed at the start of each timestep, so are available for temporary
    #              storage when gridfunctions needed for diagnostics are computed at the start of each timestep.
    #              These gridfunctions can also be freed during a regrid, to enable storage for the post-regrid
    #              destination y_n_gfs.
    #           3) Diagnostic output gridfunctions diagnostic_output_gfs & diagnostic_output_gfs2, which simply use
    #              the memory from auxiliary gridfunctions at one auxiliary time to compute diagnostics at t_n.

    # Here we specify which gridfunctions fall into each category, starting with the obvious: y_n_gridfunctions
    y_n_gridfunctions = "y_n_gfs"

    # Next the less-obvious, which depend on non-y_n_gfs
    non_y_n_gridfunctions_list = []

    # No matter the method we define gridfunctions "y_n_gfs" to store the initial data
    diagnostic_gridfunctions2_point_to = ""
    if diagonal(MoL_method) and "RK3" in MoL_method:
        non_y_n_gridfunctions_list.append("k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs")
        non_y_n_gridfunctions_list.append("k2_or_y_nplus_a32_k2_gfs")
        diagnostic_gridfunctions_point_to = "k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs"
        diagnostic_gridfunctions2_point_to = "k2_or_y_nplus_a32_k2_gfs"
    else:
        if not diagonal(MoL_method):  # Allocate memory for non-diagonal Butcher tables
            # Determine the number of k_i steps based on length of Butcher Table
            num_k = len(Butcher_dict[MoL_method][0])-1
            # For non-diagonal tables an intermediate gridfunction "next_y_input" is used for rhs evaluations
            non_y_n_gridfunctions_list.append("next_y_input_gfs")
            for i in range(num_k): # Need to allocate all k_i steps for a given method
                non_y_n_gridfunctions_list.append("k" + str(i + 1) + "_gfs")
            diagnostic_gridfunctions_point_to = "k1_gfs"
            if "k2_gfs" in non_y_n_gridfunctions_list:
                diagnostic_gridfunctions2_point_to = "k2_gfs"
            else:
                print("MoL WARNING: No gridfunction group available for diagnostic_output_gfs2")
        else:  # Allocate memory for diagonal Butcher tables, which use a "y_nplus1_running_total gridfunction"
            non_y_n_gridfunctions_list.append("y_nplus1_running_total_gfs")
            if MoL_method != 'Euler':  # Allocate memory for diagonal Butcher tables that aren't Euler
                # Need k_odd for k_1,3,5... and k_even for k_2,4,6...
                non_y_n_gridfunctions_list.append("k_odd_gfs")
                non_y_n_gridfunctions_list.append("k_even_gfs")
            diagnostic_gridfunctions_point_to  = "y_nplus1_running_total_gfs"
            diagnostic_gridfunctions2_point_to = "k_odd_gfs"
    non_y_n_gridfunctions_list.append("auxevol_gfs")

    return y_n_gridfunctions, non_y_n_gridfunctions_list,\
           diagnostic_gridfunctions_point_to, diagnostic_gridfunctions2_point_to


# add_to_Cfunction_dict_MoL_malloc() registers
#           MoL_malloc_y_n_gfs() and
#           MoL_malloc_non_y_n_gfs(), which allocate memory for
#           the indicated sets of gridfunctions
def add_to_Cfunction_dict_MoL_malloc(MoL_method, which_gfs):
    includes = ["NRPy_basic_defines.h", "NRPy_function_prototypes.h"]
    desc  = "Method of Lines (MoL) for \"" + MoL_method + "\" method: Allocate memory for \""+which_gfs+"\" gridfunctions\n"
    desc += "   * y_n_gfs are used to store data for the vector of gridfunctions y_i at t_n, at the start of each MoL timestep\n"
    desc += "   * non_y_n_gfs are needed for intermediate (e.g., k_i) storage in chosen MoL method\n"
    c_type = "void"

    y_n_gridfunctions, non_y_n_gridfunctions_list, diagnostic_gridfunctions_point_to, diagnostic_gridfunctions2_point_to = \
        generate_gridfunction_names(MoL_method = MoL_method)

    if which_gfs == "y_n_gfs":
        gridfunctions_list = [y_n_gridfunctions]
    elif which_gfs == "non_y_n_gfs":
        gridfunctions_list = non_y_n_gridfunctions_list
    else:
        print("ERROR: which_gfs = \"" + which_gfs + "\" unrecognized.")
        sys.exit(1)
    name = "MoL_malloc_" + which_gfs
    params = "const paramstruct *restrict params, MoL_gridfunctions_struct *restrict gridfuncs"
    body = "const int Nxx_plus_2NGHOSTS_tot = Nxx_plus_2NGHOSTS0*Nxx_plus_2NGHOSTS1*Nxx_plus_2NGHOSTS2;\n"
    for gridfunctions in gridfunctions_list:
        num_gfs = "NUM_EVOL_GFS"
        if gridfunctions == "auxevol_gfs":
            num_gfs = "NUM_AUXEVOL_GFS"
        body += "gridfuncs->" + gridfunctions + " = (REAL *restrict)malloc(sizeof(REAL) * " + num_gfs + " * Nxx_plus_2NGHOSTS_tot);\n"
    body += "\ngridfuncs->diagnostic_output_gfs  = gridfuncs->" + diagnostic_gridfunctions_point_to + ";\n"
    body += "\ngridfuncs->diagnostic_output_gfs2 = gridfuncs->" + diagnostic_gridfunctions2_point_to + ";\n"
    add_to_Cfunction_dict(
        includes=includes,
        desc=desc,
        c_type=c_type, name=name, params=params,
        body=indent_Ccode(body, "  "),
        rel_path_to_Cparams=os.path.join("."))

# single_RK_substep_input_symbolic() performs necessary replacements to
#   define C code for a single RK substep
#   (e.g., computing k_1 and then updating the outer boundaries)
def single_RK_substep_input_symbolic(commentblock, RHS_str, RHS_input_str, RHS_output_str, RK_lhss_list, RK_rhss_list,
                                     post_RHS_list, post_RHS_output_list, enable_SIMD=False,
                                     enable_griddata=False, gf_aliases="", post_post_RHS_string=""):
    return_str  = commentblock + "\n"
    if not isinstance(RK_lhss_list, list):
        RK_lhss_list = [RK_lhss_list]
    if not isinstance(RK_rhss_list, list):
        RK_rhss_list = [RK_rhss_list]

    if not isinstance(post_RHS_list, list):
        post_RHS_list = [post_RHS_list]
    if not isinstance(post_RHS_output_list, list):
        post_RHS_output_list = [post_RHS_output_list]

    indent = ""
    if enable_griddata:
        return_str += "{\n" + indent_Ccode(gf_aliases, "  ")
        indent = "  "
    # Part 1: RHS evaluation:
    return_str += indent_Ccode(str(RHS_str).replace("RK_INPUT_GFS",  str(RHS_input_str).replace("gfsL", "gfs")).
                               replace("RK_OUTPUT_GFS", str(RHS_output_str).replace("gfsL", "gfs"))+"\n", indent=indent)

    # Part 2: RK update
    if enable_SIMD:
        return_str += "#pragma omp parallel for\n"
        return_str += indent + "for(int i=0;i<Nxx_plus_2NGHOSTS0*Nxx_plus_2NGHOSTS1*Nxx_plus_2NGHOSTS2*NUM_EVOL_GFS;i+=SIMD_width) {\n"
    else:
        return_str += indent + "LOOP_ALL_GFS_GPS(i) {\n"
    vartype = "REAL"
    if enable_SIMD:
        vartype = "REAL_SIMD_ARRAY"
    RK_lhss_str_list = []
    for i, el in enumerate(RK_lhss_list):
        if enable_SIMD:
            RK_lhss_str_list.append(indent + "const REAL_SIMD_ARRAY __RHS_exp_" + str(i))
        else:
            RK_lhss_str_list.append(indent + str(el).replace("gfsL", "gfs[i]"))
    read_list = []
    for el in RK_rhss_list:
        for read in list(sp.ordered(el.free_symbols)):
            read_list.append(read)
    read_list_uniq = superfast_uniq(read_list)
    for el in read_list_uniq:
        if str(el) != "dt":
            if enable_SIMD:
                return_str += indent + "  const " + vartype + " " + str(el) + " = ReadSIMD(&" + str(el).replace("gfsL", "gfs[i]") + ");\n"
            else:
                return_str += indent + "  const " + vartype + " " + str(el) + " = " + str(el).replace("gfsL", "gfs[i]") + ";\n"
    if enable_SIMD:
        return_str += indent + "  const REAL_SIMD_ARRAY DT = ConstSIMD(dt);\n"
    preindent = "1"
    if enable_griddata:
        preindent = "2"
    kernel = outputC(RK_rhss_list, RK_lhss_str_list, filename="returnstring",
                     params="includebraces=False,preindent="+preindent+",outCverbose=False,enable_SIMD="+str(enable_SIMD))
    if enable_SIMD:
        return_str += kernel.replace("dt", "DT")
        for i, el in enumerate(RK_lhss_list):
            return_str += "  WriteSIMD(&" + str(el).replace("gfsL", "gfs[i]") + ", __RHS_exp_" + str(i) + ");\n"
    else:
        return_str += kernel

    return_str += indent + "}\n"

    # Part 3: Call post-RHS functions
    for post_RHS, post_RHS_output in zip(post_RHS_list, post_RHS_output_list):
        return_str += indent_Ccode(post_RHS.replace("RK_OUTPUT_GFS", str(post_RHS_output).replace("gfsL", "gfs")))

    if enable_griddata:
        return_str += "}\n"

    for post_RHS, post_RHS_output in zip(post_RHS_list, post_RHS_output_list):
        return_str += indent_Ccode(post_post_RHS_string.replace("RK_OUTPUT_GFS", str(post_RHS_output).replace("gfsL", "gfs")), "")

    return return_str


########################################################################################################################
# EXAMPLE
# ODE: y' = f(t,y), y(t_0) = y_0
# Starting at time t_n with solution having value y_n and trying to update to y_nplus1 with timestep dt

# Example of scheme for RK4 with k_1, k_2, k_3, k_4 (Using non-diagonal algorithm) Notice this requires storage of
# y_n, y_nplus1, k_1 through k_4

# k_1      = dt*f(t_n, y_n)
# k_2      = dt*f(t_n + 1/2*dt, y_n + 1/2*k_1)
# k_3      = dt*f(t_n + 1/2*dt, y_n + 1/2*k_2)
# k_4      = dt*f(t_n + dt, y_n + k_3)
# y_nplus1 = y_n + 1/3k_1 + 1/6k_2 + 1/6k_3 + 1/3k_4

# Example of scheme RK4 using only k_odd and k_even (Diagonal algroithm) Notice that this only requires storage

# k_odd     = dt*f(t_n, y_n)
# y_nplus1  = 1/3*k_odd
# k_even    = dt*f(t_n + 1/2*dt, y_n + 1/2*k_odd)
# y_nplus1 += 1/6*k_even
# k_odd     = dt*f(t_n + 1/2*dt, y_n + 1/2*k_even)
# y_nplus1 += 1/6*k_odd
# k_even    = dt*f(t_n + dt, y_n + k_odd)
# y_nplus1 += 1/3*k_even
########################################################################################################################
def add_to_Cfunction_dict_MoL_step_forward_in_time(MoL_method,
                                                   RHS_string = "", post_RHS_string = "", post_post_RHS_string="",
                                                   enable_rfm=False, enable_curviBCs=False, enable_SIMD=False,
                                                   enable_griddata=False):
    includes = ["NRPy_basic_defines.h", "NRPy_function_prototypes.h"]
    if enable_SIMD:
        includes += [os.path.join("SIMD", "SIMD_intrinsics.h")]
    desc  = "Method of Lines (MoL) for \"" + MoL_method + "\" method: Step forward one full timestep.\n"
    c_type = "void"
    name = "MoL_step_forward_in_time"
    if enable_griddata:
        params = "griddata_struct *restrict griddata, const REAL dt"
    else:
        params  = "const paramstruct *restrict params, "
        if enable_rfm:
            params += "const rfm_struct *restrict rfmstruct, "
        else:
            params += "REAL *restrict xx[3], "
        if enable_curviBCs:
            params += "const bc_struct *restrict bcstruct, "
        params += "MoL_gridfunctions_struct *restrict gridfuncs, const REAL dt"

    indent = ""  # We don't bother with an indent here.

    body = indent + "// C code implementation of -={ " + MoL_method + " }=- Method of Lines timestepping.\n\n"

    y_n_gridfunctions, non_y_n_gridfunctions_list, _throwaway, _throwaway2 = generate_gridfunction_names(MoL_method)
    if enable_griddata:
        gf_prefix = "griddata->gridfuncs."
    else:
        gf_prefix = "gridfuncs->"
    gf_aliases = """// Set gridfunction aliases from gridfuncs struct
REAL *restrict """ + y_n_gridfunctions + " = "+gf_prefix + y_n_gridfunctions + """;  // y_n gridfunctions
// Temporary timelevel & AUXEVOL gridfunctions:\n"""
    for gf in non_y_n_gridfunctions_list:
        gf_aliases += "REAL *restrict " + gf + " = "+gf_prefix + gf + ";\n"
    if enable_griddata:
        gf_aliases += "paramstruct *restrict params = &griddata->params;\n"
        if enable_rfm:
            gf_aliases += "const rfm_struct *restrict rfmstruct = &griddata->rfmstruct;\n"
        else:
            gf_aliases += "REAL * xx[3]; for(int ww=0;ww<3;ww++) xx[ww] = griddata->xx[ww];\n"
        if enable_curviBCs:
            gf_aliases += "const bc_struct *restrict bcstruct = &griddata->bcstruct;\n"
        for i in ["0", "1", "2"]:
            gf_aliases += "const int Nxx_plus_2NGHOSTS" + i + " = griddata->params.Nxx_plus_2NGHOSTS" + i + ";\n"

    if not enable_griddata:
        body += gf_aliases

    # Implement Method of Lines (MoL) Timestepping
    Butcher = Butcher_dict[MoL_method][0]  # Get the desired Butcher table from the dictionary
    num_steps = len(Butcher)-1  # Specify the number of required steps to update solution

    # Diagonal RK3 only!!!
    dt = sp.Symbol("dt", real=True)
    if diagonal(MoL_method) and "RK3" in MoL_method:
        #  In a diagonal RK3 method, only 3 gridfunctions need be defined. Below implements this approach.
        y_n_gfs = sp.Symbol("y_n_gfsL", real=True)
        k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs = sp.Symbol("k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfsL", real=True)
        k2_or_y_nplus_a32_k2_gfs = sp.Symbol("k2_or_y_nplus_a32_k2_gfsL", real=True)
        # k_1
        body += """
    // In a diagonal RK3 method like this one, only 3 gridfunctions need be defined. Below implements this approach.
    // Using y_n_gfs as input, k1 and apply boundary conditions\n"""
        body += single_RK_substep_input_symbolic(
            commentblock="""// -={ START k1 substep }=-
    // RHS evaluation:
    //  1. We will store k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs now as
    //     ...  the update for the next rhs evaluation y_n + a21*k1*dt
    // Post-RHS evaluation:
    //  1. Apply post-RHS to y_n + a21*k1*dt""",
            RHS_str=RHS_string,
            RHS_input_str=y_n_gfs,
            RHS_output_str=k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs,
            RK_lhss_list=[k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs],
            RK_rhss_list=[Butcher[1][1]*k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs*dt + y_n_gfs],
            post_RHS_list=[post_RHS_string], post_RHS_output_list=[k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs],
            enable_SIMD=enable_SIMD, enable_griddata=enable_griddata, gf_aliases=gf_aliases,
            post_post_RHS_string=post_post_RHS_string) + "// -={ END k1 substep }=-\n\n"

        # k_2
        body += single_RK_substep_input_symbolic(
            commentblock="""// -={ START k2 substep }=-
    // RHS evaluation:
    //    1. Reassign k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs to be the running total y_{n+1}; a32*k2*dt to the running total
    //    2. Store k2_or_y_nplus_a32_k2_gfs now as y_n + a32*k2*dt
    // Post-RHS evaluation:
    //    1. Apply post-RHS to both y_n + a32*k2 (stored in k2_or_y_nplus_a32_k2_gfs)
    //       ... and the y_{n+1} running total, as they have not been applied yet to k2-related gridfunctions""",
            RHS_str=RHS_string,
            RHS_input_str=k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs,
            RHS_output_str=k2_or_y_nplus_a32_k2_gfs,
            RK_lhss_list=[k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs, k2_or_y_nplus_a32_k2_gfs],
            RK_rhss_list=[Butcher[3][1]*(k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs - y_n_gfs)/Butcher[1][1] + y_n_gfs + Butcher[3][2]*k2_or_y_nplus_a32_k2_gfs*dt,
                          Butcher[2][2]*k2_or_y_nplus_a32_k2_gfs*dt + y_n_gfs],
            post_RHS_list=[post_RHS_string, post_RHS_string],
            post_RHS_output_list=[k2_or_y_nplus_a32_k2_gfs,
                                  k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs],
            enable_SIMD=enable_SIMD, enable_griddata=enable_griddata, gf_aliases=gf_aliases,
            post_post_RHS_string=post_post_RHS_string) + "// -={ END k2 substep }=-\n\n"

        # k_3
        body += single_RK_substep_input_symbolic(
            commentblock="""// -={ START k3 substep }=-
        // RHS evaluation:
        //    1. Add k3 to the running total and save to y_n
        // Post-RHS evaluation:
        //    1. Apply post-RHS to y_n""",
            RHS_str=RHS_string,
            RHS_input_str=k2_or_y_nplus_a32_k2_gfs, RHS_output_str=y_n_gfs,
            RK_lhss_list=[y_n_gfs],
            RK_rhss_list=[k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs + Butcher[3][3]*y_n_gfs*dt],
            post_RHS_list=[post_RHS_string],
            post_RHS_output_list=[y_n_gfs],
            enable_SIMD=enable_SIMD, enable_griddata=enable_griddata, gf_aliases=gf_aliases,
            post_post_RHS_string=post_post_RHS_string) + "// -={ END k3 substep }=-\n\n"
    else:
        y_n = sp.Symbol("y_n_gfsL", real=True)
        if not diagonal(MoL_method):
            for s in range(num_steps):
                next_y_input = sp.Symbol("next_y_input_gfsL", real=True)

                # If we're on the first step (s=0), we use y_n gridfunction as input.
                #      Otherwise next_y_input is input. Output is just the reverse.
                if s == 0:  # If on first step:
                    RHS_input = y_n
                else:  # If on second step or later:
                    RHS_input = next_y_input
                RHS_output = sp.Symbol("k" + str(s + 1) + "_gfs", real=True)
                if s == num_steps - 1:  # If on final step:
                    RK_lhs = y_n
                else:  # If on anything but the final step:
                    RK_lhs = next_y_input
                RK_rhs = y_n
                for m in range(s + 1):
                    k_mp1_gfs = sp.Symbol("k" + str(m + 1) + "_gfsL")
                    if Butcher[s + 1][m + 1] != 0:
                        if Butcher[s + 1][m + 1] != 1:
                            RK_rhs += dt * k_mp1_gfs*Butcher[s + 1][m + 1]
                        else:
                            RK_rhs += dt * k_mp1_gfs

                post_RHS = post_RHS_string
                if s == num_steps - 1:  # If on final step:
                    post_RHS_output = y_n
                else:  # If on anything but the final step:
                    post_RHS_output = next_y_input

                body += single_RK_substep_input_symbolic(
                    commentblock="// -={ START k" + str(s + 1) + " substep }=-",
                    RHS_str=RHS_string,
                    RHS_input_str=RHS_input, RHS_output_str=RHS_output,
                    RK_lhss_list=[RK_lhs], RK_rhss_list=[RK_rhs],
                    post_RHS_list=[post_RHS],
                    post_RHS_output_list=[post_RHS_output],
                    enable_SIMD=enable_SIMD, enable_griddata=enable_griddata, gf_aliases=gf_aliases,
                    post_post_RHS_string=post_post_RHS_string) + "// -={ END k" + str(s + 1) + " substep }=-\n\n"
        else:
            y_n = sp.Symbol("y_n_gfsL", real=True)
            y_nplus1_running_total = sp.Symbol("y_nplus1_running_total_gfsL", real=True)
            if MoL_method == 'Euler':  # Euler's method doesn't require any k_i, and gets its own unique algorithm
                body += single_RK_substep_input_symbolic(
                    commentblock=indent + "// ***Euler timestepping only requires one RHS evaluation***",
                    RHS_str=RHS_string,
                    RHS_input_str=y_n, RHS_output_str=y_nplus1_running_total,
                    RK_lhss_list=[y_n], RK_rhss_list=[y_n + y_nplus1_running_total*dt],
                    post_RHS_list=[post_RHS_string],
                    post_RHS_output_list=[y_n],
                    enable_SIMD=enable_SIMD, enable_griddata=enable_griddata, gf_aliases=gf_aliases,
                    post_post_RHS_string=post_post_RHS_string)
            else:
                for s in range(num_steps):
                    # If we're on the first step (s=0), we use y_n gridfunction as input.
                    # and k_odd as output.
                    if s == 0:
                        RHS_input  = sp.Symbol("y_n_gfsL", real=True)
                        RHS_output = sp.Symbol("k_odd_gfsL", real=True)
                    # For the remaining steps the inputs and ouputs alternate between k_odd and k_even
                    elif s % 2 == 0:
                        RHS_input = sp.Symbol("k_even_gfsL", real=True)
                        RHS_output = sp.Symbol("k_odd_gfsL", real=True)
                    else:
                        RHS_input = sp.Symbol("k_odd_gfsL", real=True)
                        RHS_output = sp.Symbol("k_even_gfsL", real=True)
                    RK_lhs_list = []
                    RK_rhs_list = []
                    if s != num_steps-1:  # For anything besides the final step
                        if s == 0:  # The first RK step
                            RK_lhs_list.append(y_nplus1_running_total)
                            RK_rhs_list.append(RHS_output*dt*Butcher[num_steps][s+1])

                            RK_lhs_list.append(RHS_output)
                            RK_rhs_list.append(y_n + RHS_output*dt*Butcher[s+1][s+1])
                        else:
                            if Butcher[num_steps][s+1] != 0:
                                RK_lhs_list.append(y_nplus1_running_total)
                                if Butcher[num_steps][s+1] != 1:
                                    RK_rhs_list.append(y_nplus1_running_total + RHS_output*dt*Butcher[num_steps][s+1])
                                else:
                                    RK_rhs_list.append(y_nplus1_running_total + RHS_output*dt)
                            if Butcher[s+1][s+1] != 0:
                                RK_lhs_list.append(RHS_output)
                                if Butcher[s+1][s+1] != 1:
                                    RK_rhs_list.append(y_n + RHS_output*dt*Butcher[s+1][s+1])
                                else:
                                    RK_rhs_list.append(y_n + RHS_output*dt)
                        post_RHS_output = RHS_output
                    if s == num_steps-1:  # If on the final step
                        if Butcher[num_steps][s+1] != 0:
                            RK_lhs_list.append(y_n)
                            if Butcher[num_steps][s+1] != 1:
                                RK_rhs_list.append(y_n + y_nplus1_running_total + RHS_output*dt*Butcher[num_steps][s+1])
                            else:
                                RK_rhs_list.append(y_n + y_nplus1_running_total + RHS_output*dt)
                        post_RHS_output = y_n
                    body += single_RK_substep_input_symbolic(
                        commentblock=indent + "// -={ START k" + str(s + 1) + " substep }=-",
                        RHS_str=RHS_string,
                        RHS_input_str=RHS_input, RHS_output_str=RHS_output,
                        RK_lhss_list=RK_lhs_list, RK_rhss_list=RK_rhs_list,
                        post_RHS_list=[post_RHS_string],
                        post_RHS_output_list=[post_RHS_output],
                        enable_SIMD=enable_SIMD, enable_griddata=enable_griddata, gf_aliases=gf_aliases,
                        post_post_RHS_string=post_post_RHS_string) + "// -={ END k" + str(s + 1) + " substep }=-\n\n"

    enableCparameters=True
    if enable_griddata:
        enableCparameters=False
    add_to_Cfunction_dict(
        includes=includes,
        desc=desc,
        c_type=c_type, name=name, params=params,
        body=indent_Ccode(body, "  "),
        enableCparameters=enableCparameters, rel_path_to_Cparams=os.path.join("."))


# add_to_Cfunction_dict_MoL_free_memory() registers
#           MoL_free_memory_y_n_gfs() and
#           MoL_free_memory_non_y_n_gfs(), which free memory for
#           the indicated sets of gridfunctions
def add_to_Cfunction_dict_MoL_free_memory(MoL_method, which_gfs):
    includes = ["NRPy_basic_defines.h", "NRPy_function_prototypes.h"]
    desc = "Method of Lines (MoL) for \"" + MoL_method + "\" method: Free memory for \"" + which_gfs + "\" gridfunctions\n"
    desc += "   - y_n_gfs are used to store data for the vector of gridfunctions y_i at t_n, at the start of each MoL timestep\n"
    desc += "   - non_y_n_gfs are needed for intermediate (e.g., k_i) storage in chosen MoL method\n"
    c_type = "void"

    y_n_gridfunctions, non_y_n_gridfunctions_list, _diagnostic_gridfunctions_point_to, \
        _diagnostic_gridfunctions2_point_to = generate_gridfunction_names(MoL_method=MoL_method)

    if which_gfs == "y_n_gfs":
        gridfunctions_list = [y_n_gridfunctions]
    elif which_gfs == "non_y_n_gfs":
        gridfunctions_list = non_y_n_gridfunctions_list
    else:
        print("ERROR: which_gfs = \"" + which_gfs + "\" unrecognized.")
        sys.exit(1)
    name = "MoL_free_memory_" + which_gfs
    params = "const paramstruct *restrict params, MoL_gridfunctions_struct *restrict gridfuncs"
    body = ""
    for gridfunctions in gridfunctions_list:
        body += "    free(gridfuncs->" + gridfunctions + ");\n"
    add_to_Cfunction_dict(
        includes=includes,
        desc=desc,
        c_type=c_type, name=name, params=params,
        body=indent_Ccode(body, "  "),
        rel_path_to_Cparams=os.path.join("."))


# Register MoL_gridfunctions_struct in NRPy_basic_defines
def NRPy_basic_defines_MoL_timestepping_struct(MoL_method="RK4"):
    y_n_gridfunctions, non_y_n_gridfunctions_list, _diagnostic_gridfunctions_point_to, \
        _diagnostic_gridfunctions2_point_to = generate_gridfunction_names(MoL_method=MoL_method)
    # Step 3.b: Create MoL_timestepping struct:
    indent = "  "
    Nbd = "typedef struct __MoL_gridfunctions_struct__ {\n"
    Nbd += indent + "REAL *restrict " + y_n_gridfunctions + ";\n"
    for gfs in non_y_n_gridfunctions_list:
        Nbd += indent + "REAL *restrict " + gfs + ";\n"
    Nbd += indent + "REAL *restrict diagnostic_output_gfs;\n"
    Nbd += indent + "REAL *restrict diagnostic_output_gfs2;\n"
    Nbd += "} MoL_gridfunctions_struct;\n"
    Nbd += """#define LOOP_ALL_GFS_GPS(ii) _Pragma("omp parallel for") \\
  for(int (ii)=0;(ii)<Nxx_plus_2NGHOSTS0*Nxx_plus_2NGHOSTS1*Nxx_plus_2NGHOSTS2*NUM_EVOL_GFS;(ii)++)\n"""

    outC_NRPy_basic_defines_h_dict["MoL"] = Nbd

    # Finally, register per-grid (griddata) and all-grids
    #   (commondata) data for this module. griddata and
    #   and commondata are declared inside NRPy_basic_defines.h.
    import grid as gri
    gri.glb_griddata_struct_list += [gri.glb_griddata(__name__, "MoL_gridfunctions_struct gridfuncs;")]


# Finally declare the master registration function
def register_C_functions_and_NRPy_basic_defines(MoL_method = "RK4",
            RHS_string =  "rhs_eval(Nxx,Nxx_plus_2NGHOSTS,dxx, RK_INPUT_GFS, RK_OUTPUT_GFS);",
            post_RHS_string = "apply_bcs(Nxx,Nxx_plus_2NGHOSTS, RK_OUTPUT_GFS);", post_post_RHS_string = "",
            enable_rfm=False, enable_curviBCs=False, enable_SIMD=False, enable_griddata=False):
    for which_gfs in ["y_n_gfs", "non_y_n_gfs"]:
        add_to_Cfunction_dict_MoL_malloc(MoL_method, which_gfs)
        add_to_Cfunction_dict_MoL_free_memory(MoL_method, which_gfs)
    add_to_Cfunction_dict_MoL_step_forward_in_time(MoL_method, RHS_string, post_RHS_string, post_post_RHS_string,
                                                   enable_rfm=enable_rfm, enable_curviBCs=enable_curviBCs,
                                                   enable_SIMD=enable_SIMD, enable_griddata=enable_griddata)
    NRPy_basic_defines_MoL_timestepping_struct(MoL_method=MoL_method)
