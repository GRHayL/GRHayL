# PLEASE DO NOT USE THESE FUNCTIONS; THEY HAVE BEEN DEPRECATED

# As documented in the NRPy+ tutorial module
#   Tutorial-Coutput__Parameter_Interface.ipynb
#   this core NRPy+ module is used for
#   initializing, storing, and recalling
#   parameters.

# Author: Zachariah B. Etienne
#         zachetie **at** gmail **dot* com

import sympy as sp                   # Import SymPy
import os, sys                       # Standard Python: OS-independent system functions
from collections import namedtuple   # Standard Python: Enable namedtuple data type
import re
import textwrap
import NRPy_param_funcs as par

## TO BE DEPRECATED BY NRPy_param_funcs_register_C_functions_and_NRPy_basic_defines() WITHIN OUTPUTC
def generate_Cparameters_Ccodes(directory=os.path.join(".")):
    # Step 1: Check that Cparams types are supported.
    for i in range(len(par.glb_Cparams_list)):
        partype = par.glb_Cparams_list[i].type
        if partype not in ('bool', '#define', 'char', 'int', 'REAL'):
            print("Error: parameter "+par.glb_Cparams_list[i].module+"::"+par.glb_Cparams_list[i].parname+" has unsupported type: \""
                  + par.glb_Cparams_list[i].type + "\"")
            sys.exit(1)

    # Step 2: Generate C code to declare C paramstruct;
    #         output to "declare_Cparameters_struct.h"
    #         We want the elements of this struct to be *sorted*,
    #         to ensure that the struct is consistently ordered
    #         for checkpointing purposes.
    Cparameters_struct = "typedef struct __paramstruct__ {\n"
    CCodelines = []
    for i in range(len(par.glb_Cparams_list)):
        if par.glb_Cparams_list[i].type != "#define":
            if par.glb_Cparams_list[i].type == "char":
                c_type = "char"
            else:
                c_type = par.glb_Cparams_list[i].type
            comment = "  // " + par.glb_Cparams_list[i].module + "::" + par.glb_Cparams_list[i].parname
            CCodelines.append("    " + c_type + " " + par.glb_Cparams_list[i].parname + ";" + comment + "\n")
    for line in sorted(CCodelines):
        Cparameters_struct += line
    Cparameters_struct += "} paramstruct;\n"
    with open(os.path.join(directory, "declare_Cparameters_struct.h"), "w") as file:
        file.write(Cparameters_struct)

    # Step 3: Generate C code to set all elements in
    #         C paramstruct to default values; output to
    #         "set_Cparameters_default.h"
    with open(os.path.join(directory, "set_Cparameters_default.h"), "w") as file:
        for i in range(len(par.glb_Cparams_list)):
            if par.glb_Cparams_list[i].type != "#define":
                c_output = "params." + par.glb_Cparams_list[i].parname
                comment = "  // " + par.glb_Cparams_list[i].module + "::" + par.glb_Cparams_list[i].parname
                if isinstance(par.glb_Cparams_list[i].defaultval, (bool, int, float)):
                    c_output += " = " + str(par.glb_Cparams_list[i].defaultval).lower() + ";" + comment + "\n"
                elif par.glb_Cparams_list[i].type == "char" and isinstance(par.glb_Cparams_list[i].defaultval, (str)):
                    c_output += " = \"" + str(par.glb_Cparams_list[i].defaultval).lower() + "\";" + comment + "\n"
                else:
                    c_output += " = " + str(par.glb_Cparams_list[i].defaultval) + ";" + comment + "\n"
                file.write(c_output)

    # Step 4: Generate C code to set C parameter constants
    #         (i.e., all ints != -12345678 and REALs != 1e300);
    #         output to filename "set_Cparameters.h" if enable_SIMD==False
    #         or "set_Cparameters-SIMD.h" if enable_SIMD==True
    # Step 4.a: Output non-SIMD version, set_Cparameters.h
    def gen_set_Cparameters(pointerEnable=True):
        returnstring = ""
        for i in range(len(par.glb_Cparams_list)):
            if par.glb_Cparams_list[i].type == "char":
                c_type = "char *"
            else:
                c_type = par.glb_Cparams_list[i].type

            pointer = "->"
            if pointerEnable==False:
                pointer = "."

            if not ((c_type == "REAL" and par.glb_Cparams_list[i].defaultval == 1e300) or c_type == "#define"):
                comment = "  // " + par.glb_Cparams_list[i].module + "::" + par.glb_Cparams_list[i].parname
                Coutput = "const "+c_type+" "+par.glb_Cparams_list[i].parname+" = "+"params"+pointer+par.glb_Cparams_list[i].parname + ";" + comment + "\n"
                returnstring += Coutput
        return returnstring

    with open(os.path.join(directory, "set_Cparameters.h"), "w") as file:
        file.write(gen_set_Cparameters(pointerEnable=True))
    with open(os.path.join(directory, "set_Cparameters-nopointer.h"), "w") as file:
        file.write(gen_set_Cparameters(pointerEnable=False))

    # Step 4.b: Output SIMD version, set_Cparameters-SIMD.h
    with open(os.path.join(directory, "set_Cparameters-SIMD.h"), "w") as file:
        for i in range(len(par.glb_Cparams_list)):
            if par.glb_Cparams_list[i].type == "char":
                c_type = "char *"
            else:
                c_type = par.glb_Cparams_list[i].type

            comment = "  // " + par.glb_Cparams_list[i].module + "::" + par.glb_Cparams_list[i].parname
            parname = par.glb_Cparams_list[i].parname
            if c_type == "REAL" and par.glb_Cparams_list[i].defaultval != 1e300:
                c_output =  "const REAL            NOSIMD" + parname + " = " + "params->" + par.glb_Cparams_list[i].parname + ";"+comment+"\n"
                c_output += "const REAL_SIMD_ARRAY " + parname + " = ConstSIMD(NOSIMD" + parname + ");"+comment+"\n"
                file.write(c_output)
            elif par.glb_Cparams_list[i].defaultval != 1e300 and c_type != "#define":
                c_output = "const "+c_type+" "+parname + " = " + "params->" + par.glb_Cparams_list[i].parname + ";"+comment+"\n"
                file.write(c_output)
