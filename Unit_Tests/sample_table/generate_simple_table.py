from numpy import linspace
from h5py import File as H5File

def set_table_quantity(i, npoints):
    return linspace(i*10, i*10 + 9, npoints)

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

    # Set energy shift
    energy_shift = 123e45

    # Set table quantities
    i=0
    table_logpress  = set_table_quantity(i, npoints); i+=1
    table_logenergy = set_table_quantity(i, npoints); i+=1
    table_entropy   = set_table_quantity(i, npoints); i+=1
    table_munu      = set_table_quantity(i, npoints); i+=1
    table_cs2       = set_table_quantity(i, npoints); i+=1
    table_dedt      = set_table_quantity(i, npoints); i+=1
    table_dpdrhoe   = set_table_quantity(i, npoints); i+=1
    table_dpderho   = set_table_quantity(i, npoints); i+=1
    table_muhat     = set_table_quantity(i, npoints); i+=1
    table_mu_e      = set_table_quantity(i, npoints); i+=1
    table_mu_p      = set_table_quantity(i, npoints); i+=1
    table_mu_n      = set_table_quantity(i, npoints); i+=1
    table_Xa        = set_table_quantity(i, npoints); i+=1
    table_Xh        = set_table_quantity(i, npoints); i+=1
    table_Xn        = set_table_quantity(i, npoints); i+=1
    table_Xp        = set_table_quantity(i, npoints); i+=1
    table_Abar      = set_table_quantity(i, npoints); i+=1
    table_Zbar      = set_table_quantity(i, npoints); i+=1
    table_gamma     = set_table_quantity(i, npoints)

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
