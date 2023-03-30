# Step 0: Add NRPy's directory to the path
# https://stackoverflow.com/questions/16780014/import-file-from-parent-directory
import os,sys

nrpy_dir_path = os.path.join("nrpy/")
if nrpy_dir_path not in sys.path:
    sys.path.append(nrpy_dir_path)

import IGM_All_fluxes as fl
import IGM_Characteristic_Speeds as chsp
import IGM_All_Source_Terms as st    # NRPy+: Generate general relativistic magnetohydrodynamics equations


includes = ["flux_source.h"]

def write_grhayl_Cfunctions_to_dir(Ccodesrootdir, includes, formalism="ADM", tabulated=False):
    outCparams = "outCverbose=False,GoldenKernelsEnable=True"


    st.Cfunction__GRMHD_SourceTerms(Ccodesrootdir,
                                    includes=includes,
                                    formalism=formalism,
                                    outCparams=outCparams)

    chsp.Cfunction__GRMHD_characteristic_speeds(Ccodesrootdir,
                                                includes=includes,
                                                formalism=formalism,
                                                outCparams=outCparams)

    fl.Cfunction__GRMHD_fluxes(Ccodesrootdir,
                               includes=includes,
                               formalism=formalism,
                               outCparams=outCparams,
                               tabulated=tabulated)



if __name__ == '__main__':

    tabulated=False
    if len(sys.argv) > 1:
        if sys.argv[-1].lower() == "tabulated":
            tabulated=True

    write_grhayl_Cfunctions_to_dir("./tabulated", includes, tabulated=tabulated)

    print("Finished printing C files!")
