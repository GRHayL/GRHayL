# As documented in the NRPy+ tutorial module
#   Tutorial-Coutput__Parameter_Interface.ipynb
#   this core NRPy+ module is used for
#   initializing, storing, and recalling
#   parameters.

# Author: Zachariah B. Etienne
#         zachetie **at** gmail **dot* com

import sympy as sp                   # Import SymPy
import sys                           # Standard Python: OS-independent system functions
from collections import namedtuple   # Standard Python: Enable namedtuple data type
import re
import textwrap

wrapper = textwrap.TextWrapper(initial_indent="", subsequent_indent="  ", width=100)

glb_params_list = []  # = where we store NRPy+ parameters and default values of parameters. A list of named tuples
glb_paramsvals_list = []  # = where we store NRPy+ parameter values.
glb_param = namedtuple('glb_param', 'type module parname defaultval')

glb_Cparams_list = []  # = where we store C runtime parameters and default values of parameters. A list of named tuples
glb_Cparam = namedtuple('glb_Cparam', 'type module parname defaultval')

veryverbose = False


def initialize_param(param):
    if get_params_idx(param) == -1:
        glb_params_list.append(param)
        glb_paramsvals_list.append(param.defaultval)
    else:
        if veryverbose == True:
            print("initialize_param() minor warning: Did nothing; already initialized parameter " + param.module + "::" + param.parname)


def initialize_Cparam(param):
    if get_params_idx(param, Cparam=True) == -1:
        glb_Cparams_list.append(param)
    else:
        if veryverbose == True:
            print("initialize_Cparam() minor warning: Did nothing; already initialized parameter " + param.module + "::" + param.parname)


# Given the named tuple `input` and list of named tuples `params`,
#    defined according to namedtuple('param', 'type module name defaultval'),
#    where in the case of `input`, defaultval need not be set,
#    return the list index of `params` that matches `input`.
# On error returns -1
def get_params_idx(param, Cparam=False):
    # inspired by: https://stackoverflow.com/questions/2917372/how-to-search-a-list-of-tuples-in-python:
    if Cparam==False:
        lst = [i for i, v in enumerate(glb_params_list)
               if (param.type == "ignoretype" or param.type == v[0]) and param.module == v[1] and param.parname == v[2]]
    else:
        lst = [i for i, v in enumerate(glb_Cparams_list) if param.parname == v[2]]
    if lst == []:
        return -1  # No match found => error out!
    if len(lst) > 1:
        print("Error: Found multiple parameters matching " + str(param))
        sys.exit(1)
    return lst.pop()  # pop() returns the index


def get_params_value(param):
    idx = get_params_idx(param)
    if idx < 0:
        print("Error: could not find a parameter matching:", param)
        print("Full list of modules:\n", wrapper.fill(str(glb_params_list)))
        sys.exit(1)
    else:
        return glb_paramsvals_list[idx]


#
def idx_from_str(varname,modname=""):
    if "::" in varname:
        splitstring = re.split('::', varname)
        modname=splitstring[0]
        varname=splitstring[1]

    # inspired by: https://stackoverflow.com/questions/2917372/how-to-search-a-list-of-tuples-in-python:
    if modname == "":
        lst = [i for i, v in enumerate(glb_params_list) if v[2] == varname]
    else:
        lst = [i for i, v in enumerate(glb_params_list) if (v[1] == modname and v[2] == varname)]
    if lst == []:
        print("Error: Could not find a parameter matching \""+varname+"\" in \n",wrapper.fill(str(glb_params_list)))
        sys.exit(1)
    if len(lst) > 1:
        print("Error: Found more than one parameter named \""+varname+"\". Use get_params_value() instead.")
        sys.exit(1)
    return lst.pop()


def parval_from_str(string):
    return glb_paramsvals_list[idx_from_str(string)]


def set_parval_from_str(string, value):
    glb_paramsvals_list[idx_from_str(string)] = value

# parse_param_string__set__params_and_paramsvars:
# Summary: This function parses a string like this
# module::variablename=value
# into its component parts: module,variablename,value
# then hunts for module,variablename in the
# params list. When it finds the matching index in the
# list "idx", it sets paramsvals[idx] = value.
# Usage comment: This function parses both parameter file
#   inputs as well as parameter inputs from the command
#   line. You should set filename = "" if reading from the
#   command line. This ensures the error messages are
#   appropriate for the context.
def set_paramsvals_value(line, filename="", FindMainModuleMode=False):
    MainModuleFound = True
    if FindMainModuleMode == True:
        MainModuleFound = False
    # First remove carriage return and leading whitespace:
    stripped_line_of_text = line.strip()
    # Only process lines that do NOT start with a hash
    if not stripped_line_of_text.startswith("#"):
        # Valid lines take the form:
        #    module::variable = value # Comment here
        # Thus, using the delimiters "::", "=", and "#",
        # with the regex split command, the first three
        # items in single_param_def will be [module, variablename, value]
        single_param_def = re.split('::|=|#', stripped_line_of_text)

        # First verify that single_param_def has at least 3 parts:
        if len(single_param_def) < 3:
            if filename != "":
                print("Error: the line " + line + " in parameter file " + filename + " is not in the form")
                print("\"module::variable = value\"")
            else:
                print("Error: the command-line argument " + stripped_line_of_text + " is not in the form")
                print("\"module::variable=value\"   <-- NOTICE NO SPACES ALLOWED!")
            sys.exit(1)

        # Next remove all leading/trailing whitespace from single_param_def
        for i in range(len(single_param_def)):
            single_param_def[i] = single_param_def[i].strip()

        if FindMainModuleMode == False:
            # Next find the parameter in the params list and set
            # the corresponding element in paramsvals to the value.
            idx = get_params_idx(glb_param("ignoretype", single_param_def[0], single_param_def[1], "ignoredefval"))
            # If parameter is not found, print useful error message, then exit:
            if idx == -1:
                if filename != "":
                    print("Error: when reading line \"" + line + "\" in parameter file \"" + filename + "\":")
                else:
                    print("Error: when parsing command-line argument \"" + stripped_line_of_text + "\":")
                print("\t\tcould not find parameter \""+ single_param_def[1] + "\" in \""+single_param_def[0]+"\" module.")
                sys.exit(1)
            # If parameter is found at index idx, set paramsval[idx] to the value specified in the file.
            partype = glb_params_list[idx].type
            if partype == "bool":
                if single_param_def[2] == "True":
                    glb_paramsvals_list[idx] = True
                elif single_param_def[2] == "False":
                    glb_paramsvals_list[idx] = False
                else:
                    print("Error: \"bool\" type can only take values of \"True\" or \"False\"")
                    sys.exit(1)
            elif partype == "int":
                glb_paramsvals_list[idx] = int(single_param_def[2])
            elif partype in ('REAL', 'char', 'char *'):
                glb_paramsvals_list[idx] = single_param_def[2]
            else:
                print("Error: type \""+partype+"\" on variable \""+ glb_params_list[idx].parname +"\" is unsupported.")
                print("Supported types include: bool, int, REAL, REALARRAY, char, and char *")
                sys.exit(1)
        elif FindMainModuleMode == True and MainModuleFound == False:
            if single_param_def[0] == "NRPy" and single_param_def[1] == "MainModule":
                idx = get_params_idx(glb_param("ignoretype", single_param_def[0], single_param_def[1], "ignoredefval"))
                if idx == -1:
                    print("Critical error: NRPy::MainModule is uninitialized!")
                    sys.exit(1)
                glb_paramsvals_list[idx] = single_param_def[2]


def Cparameters(c_type, module, names, default_vals, assumption="Real"):
    output = []
    # if names is not a list, make it a list, to
    #      simplify the remainder of this routine.
    if not isinstance(names, list):
        names = [names]
    default_val_list = []
    if not isinstance(default_vals, list):
        for i in range(len(names)):
            default_val_list.append(default_vals)
    else:
        # If default_vals *is* a list, then make sure it has the same number of elements as "names".
        if len(default_vals) != len(names):
            print("Error in Cparameters(): Was provided a list of variables:\n"+str(names)+"\n")
            print("and a list of their default values:\n" + str(default_vals) + "\n")
            print("but the lists have different lengths (" + str(len(names)) +" != " + str(len(default_vals)) + ")\n")
            sys.exit(1)
        default_val_list = default_vals

    names = [name for name in names]
    for idx, name in enumerate(names):
        initialize_Cparam(glb_Cparam(c_type, module, name, default_val_list[idx]))
        if assumption == "Real":
            tmp = sp.Symbol(name, real=True)  # Assumes all Cparameters are real.
        elif assumption == "RealPositive":
            tmp = sp.Symbol(name, real=True, positive=True)  # Assumes all Cparameters are real and positive.
        else:
            print("Error: assumption "+str(assumption)+" not supported.")
            sys.exit(1)
        output.append(tmp)
    if len(names) == 1:
        return output[0]
    return output


def register_NRPy_basic_defines():
    import grid as gri
    # Register griddata_struct variables for this module,
    #   where griddata_struct is declared in NRPy_basic_defines.h
    gri.glb_griddata_struct_list += [gri.glb_griddata(__name__, "paramstruct params;")]
