from numpy import linspace, meshgrid, log
from h5py import File as H5File

CGS_TO_CODE_LENGTH   = 6.77269222552442e-06
CGS_TO_CODE_TIME     = 2.03040204956746e+05
CGS_TO_CODE_DENSITY  = 1.61887093132742e-18
CGS_TO_CODE_PRESSURE = 1.80123683248503e-39
CGS_TO_CODE_ENERGY   = 1.11265005605362e-21

def set_table_quantity(i, rho, T, Y_e):
    if i==0:
        return (rho + Y_e + T - log(CGS_TO_CODE_PRESSURE)) / log(10.0) # Pressure
    elif i==1:
        return (rho + Y_e - T - log(CGS_TO_CODE_ENERGY)) / log(10.0) # eps
    elif i==2:
        return rho - Y_e + T # entropy
    elif i==3:
        return rho - Y_e - T # munu
    elif i==4:
        return (-rho + Y_e + T) / (CGS_TO_CODE_LENGTH*CGS_TO_CODE_LENGTH/CGS_TO_CODE_TIME/CGS_TO_CODE_TIME) # cs2
    elif i==5:
        return (-rho + Y_e - T) / (CGS_TO_CODE_ENERGY)  # deps/dT
    elif i==6:
        return (-rho - Y_e + T) / (CGS_TO_CODE_PRESSURE/CGS_TO_CODE_DENSITY) # dP/drho
    elif i==7:
        return (-rho - Y_e - T) / (CGS_TO_CODE_PRESSURE/CGS_TO_CODE_ENERGY) # dP/deps
    elif i==8:
        return 2*rho + Y_e + T # muhat
    elif i==9:
        return 2*rho + Y_e - T # mu_e
    elif i==10:
        return 2*rho - Y_e + T # mu_p
    elif i==11:
        return 2*rho - Y_e - T # mu_n
    elif i==12:
        return -2*rho + Y_e + T # X_a
    elif i==13:
        return -2*rho + Y_e - T # X_h
    elif i==14:
        return -2*rho - Y_e + T # X_n
    elif i==15:
        return -2*rho - Y_e - T # X_p
    elif i==16:
        return rho + 3*Y_e - T # Abar
    elif i==17:
        return rho - 3*Y_e + T # Zbar
    elif i==18:
        return rho - 3*Y_e - T # Gamma

def create_simple_table():

    # Set N_rho, N_T, N_Y_e
    N_rho   = 7
    N_T     = 5
    N_Y_e   = 3
    npoints = N_rho * N_T * N_Y_e
    dims    = (N_Y_e, N_T, N_rho)

    # Set rho, T, and Y_e
    table_log_rho = linspace(1, N_rho, N_rho)
    table_log_T   = linspace(1, N_T  , N_T  )
    table_Y_e     = linspace(1, N_Y_e, N_Y_e)

    Y_e, log_T, log_rho = meshgrid(table_Y_e, table_log_T, table_log_rho, indexing='ij')

    # Set energy shift
    energy_shift = 123e-45 / CGS_TO_CODE_ENERGY

    # Set table quantities
    i=0
    table_logpress  = set_table_quantity(i, log_rho, log_T, Y_e); i+=1
    table_logenergy = set_table_quantity(i, log_rho, log_T, Y_e); i+=1
    table_entropy   = set_table_quantity(i, log_rho, log_T, Y_e); i+=1
    table_munu      = set_table_quantity(i, log_rho, log_T, Y_e); i+=1
    table_cs2       = set_table_quantity(i, log_rho, log_T, Y_e); i+=1
    table_dedt      = set_table_quantity(i, log_rho, log_T, Y_e); i+=1
    table_dpdrhoe   = set_table_quantity(i, log_rho, log_T, Y_e); i+=1
    table_dpderho   = set_table_quantity(i, log_rho, log_T, Y_e); i+=1
    table_muhat     = set_table_quantity(i, log_rho, log_T, Y_e); i+=1
    table_mu_e      = set_table_quantity(i, log_rho, log_T, Y_e); i+=1
    table_mu_p      = set_table_quantity(i, log_rho, log_T, Y_e); i+=1
    table_mu_n      = set_table_quantity(i, log_rho, log_T, Y_e); i+=1
    table_Xa        = set_table_quantity(i, log_rho, log_T, Y_e); i+=1
    table_Xh        = set_table_quantity(i, log_rho, log_T, Y_e); i+=1
    table_Xn        = set_table_quantity(i, log_rho, log_T, Y_e); i+=1
    table_Xp        = set_table_quantity(i, log_rho, log_T, Y_e); i+=1
    table_Abar      = set_table_quantity(i, log_rho, log_T, Y_e); i+=1
    table_Zbar      = set_table_quantity(i, log_rho, log_T, Y_e); i+=1
    table_gamma     = set_table_quantity(i, log_rho, log_T, Y_e)

    table_log_rho = ( table_log_rho - log(CGS_TO_CODE_DENSITY) ) / log(10.0)
    table_log_T  /= log(10.0)

    # Now create the output table
    with H5File("simple_table.h5", "w") as f:
        # Create basic datasets
        f.create_dataset("pointsrho"   , shape=1    , dtype="i4", data=N_rho          )
        f.create_dataset("pointstemp"  , shape=1    , dtype="i4", data=N_T            )
        f.create_dataset("pointsye"    , shape=1    , dtype="i4", data=N_Y_e          )
        f.create_dataset("energy_shift", shape=1    , dtype="f8", data=energy_shift   )
        f.create_dataset("logrho"      , shape=N_rho, dtype="f8", data=table_log_rho  )
        f.create_dataset("logtemp"     , shape=N_T  , dtype="f8", data=table_log_T    )
        f.create_dataset("ye"          , shape=N_Y_e, dtype="f8", data=table_Y_e      )

        # Create derived datasets
        f.create_dataset("Abar"        , shape=dims , dtype="f8", data=table_Abar     )
        f.create_dataset("Xa"          , shape=dims , dtype="f8", data=table_Xa       )
        f.create_dataset("Xh"          , shape=dims , dtype="f8", data=table_Xh       )
        f.create_dataset("Xn"          , shape=dims , dtype="f8", data=table_Xn       )
        f.create_dataset("Xp"          , shape=dims , dtype="f8", data=table_Xp       )
        f.create_dataset("Zbar"        , shape=dims , dtype="f8", data=table_Zbar     )
        f.create_dataset("cs2"         , shape=dims , dtype="f8", data=table_cs2      )
        f.create_dataset("dedt"        , shape=dims , dtype="f8", data=table_dedt     )
        f.create_dataset("dpderho"     , shape=dims , dtype="f8", data=table_dpderho  )
        f.create_dataset("dpdrhoe"     , shape=dims , dtype="f8", data=table_dpdrhoe  )
        f.create_dataset("entropy"     , shape=dims , dtype="f8", data=table_entropy  )
        f.create_dataset("gamma"       , shape=dims , dtype="f8", data=table_gamma    )
        f.create_dataset("logenergy"   , shape=dims , dtype="f8", data=table_logenergy)
        f.create_dataset("logpress"    , shape=dims , dtype="f8", data=table_logpress )
        f.create_dataset("mu_e"        , shape=dims , dtype="f8", data=table_mu_e     )
        f.create_dataset("mu_n"        , shape=dims , dtype="f8", data=table_mu_n     )
        f.create_dataset("mu_p"        , shape=dims , dtype="f8", data=table_mu_p     )
        f.create_dataset("muhat"       , shape=dims , dtype="f8", data=table_muhat    )
        f.create_dataset("munu"        , shape=dims , dtype="f8", data=table_munu     )

if __name__ == '__main__':
    create_simple_table()

# string = """
# double set_table_quantity(
#       const int which_var,
#       const double rho,
#       const double Y_e,
#       const double T ) {
#   switch (which_var) {
#     if i==0:
#       return rho + Y_e + T # Pressure
#     elif i==1:
#       return rho + Y_e - T # eps
#     elif i==2:
#       return rho - Y_e + T # entropy
#     elif i==3:
#       return rho - Y_e - T # munu
#     elif i==4:
#       return -rho + Y_e + T # cs2
#     elif i==5:
#       return -rho + Y_e - T  # deps/dT
#     elif i==6:
#       return -rho - Y_e + T # dP/drho
#     elif i==7:
#       return -rho - Y_e - T # dP/deps
#     elif i==8:
#       return 2*rho + Y_e + T # muhat
#     elif i==9:
#       return 2*rho + Y_e - T # mu_e
#     elif i==10:
#       return 2*rho - Y_e + T # mu_p
#     elif i==11:
#       return 2*rho - Y_e - T # mu_n
#     elif i==12:
#       return -2*rho + Y_e + T # X_a
#     elif i==13:
#       return -2*rho + Y_e - T # X_h
#     elif i==14:
#       return -2*rho - Y_e + T # X_n
#     elif i==15:
#       return -2*rho - Y_e - T # X_p
#     elif i==16:
#       return rho + 3*Y_e - T # Abar
#     elif i==17:
#       return rho - 3*Y_e + T # Zbar
#     elif i==18:
#       return rho - 3*Y_e - T # Gamma
#     default:
#       grhayl_error("Invalid variable %d\\n", which_var);
#   }
# }\n"""

# string = string.replace("elif ", "case ")
# string = string.replace("if ", "case ")
# string = string.replace("i==", "")
# string = string.replace(" #", "; //")
# print(string)
