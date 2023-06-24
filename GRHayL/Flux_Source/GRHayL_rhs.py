# Step 0: Add NRPy's directory to the path
# https://stackoverflow.com/questions/16780014/import-file-from-parent-directory
import os,sys,shutil

nrpy_dir_path = os.path.join("nrpy/")
if nrpy_dir_path not in sys.path:
    sys.path.append(nrpy_dir_path)

import IGM_All_fluxes as fl
import IGM_Characteristic_Speeds as chsp
import IGM_All_Source_Terms as st    # NRPy+: Generate general relativistic magnetohydrodynamics equations


includes = ["flux_source.h"]

def write_ghl_Cfunctions_to_dir(Ccodesrootdir, includes, formalism="ADM", tabulated=False, entropy=False):
    outCparams = "outCverbose=False,GoldenKernelsEnable=True"


#     st.Cfunction__GRMHD_SourceTerms(Ccodesrootdir,
#                                     includes=includes,
#                                     formalism=formalism,
#                                     outCparams=outCparams)

#     chsp.Cfunction__GRMHD_characteristic_speeds(Ccodesrootdir,
#                                                 includes=includes,
#                                                 formalism=formalism,
#                                                 outCparams=outCparams)

    fl.Cfunction__GRMHD_fluxes(Ccodesrootdir,
                               includes=includes,
                               formalism=formalism,
                               outCparams=outCparams,
                               tabulated=tabulated,
                               entropy=entropy)


def create_fresh_directory(name):
    shutil.rmtree(name, ignore_errors=True)
    os.mkdir(name)

if __name__ == '__main__':
    for dirname in ["hybrid", "hybrid_entropy", "tabulated", "tabulated_entropy"]:
        create_fresh_directory(dirname)
        entropy   = True if "entropy"   in dirname else False
        tabulated = True if "tabulated" in dirname else False
        write_ghl_Cfunctions_to_dir(dirname, includes, tabulated=tabulated, entropy=entropy)
