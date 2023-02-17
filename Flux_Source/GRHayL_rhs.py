# Step 0: Add NRPy's directory to the path
# https://stackoverflow.com/questions/16780014/import-file-from-parent-directory
import os,sys

nrpy_dir_path = os.path.join("nrpy/")
if nrpy_dir_path not in sys.path:
    sys.path.append(nrpy_dir_path)

import IGM_All_fluxes as fl    
import IGM_Characteristic_Speeds as chsp 
import IGM_All_Source_Terms as st    # NRPy+: Generate general relativistic magnetohydrodynamics equations


includes = ["flux_source_terms.h"]

def write_grhayl_Cfunctions_to_dir(Ccodesrootdir, includes, formalism="ADM"):
    outCparams = "outCverbose=False,CSE_sorting=True"


    st.Cfunction__GRMHD_SourceTerms(Ccodesrootdir, 
                                    includes=includes,
                                    formalism=formalism)

    chsp.Cfunction__GRMHD_characteristic_speeds(Ccodesrootdir, 
                                                includes=includes,
                                                formalism=formalism,
                                                outCparams=outCparams)

    fl.Cfunction__GRMHD_fluxes(Ccodesrootdir, 
                               includes=includes,
                               formalism=formalism,
                               outCparams=outCparams)



if __name__ == '__main__':
    write_grhayl_Cfunctions_to_dir("./", includes)

    print("Finished printing C files!")
    