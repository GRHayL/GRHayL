from os import getcwd
from os.path import join
from sys import exit, argv
from h5py import File as H5File
from numpy import log10, linspace, meshgrid
from scipy.interpolate import RegularGridInterpolator

def check_limits(intable, rho_min, rho_max, T_min, T_max, Y_e_min, Y_e_max):
    # Set mins and maxs
    table_rho_min = 10**intable["logrho"][0]
    table_rho_max = 10**intable["logrho"][-1]
    table_T_min   = 10**intable["logtemp"][0]
    table_T_max   = 10**intable["logtemp"][-1]
    table_Y_e_min = intable["ye"][0]
    table_Y_e_max = intable["ye"][-1]

    if rho_min is None:
        rho_min = 10**intable["logrho"][0]
    else:
        if rho_min < table_rho_min:
            print(f'Error: rho_min exceeds table bounds ({rho_min:e} < {table_rho_min:e}).')
            exit(1)
    if rho_max is None:
        rho_max = 10**intable["logrho"][-1]
    else:
        if rho_max > table_rho_max:
            print(f'Error: rho_max exceeds table bounds ({rho_max:e} > {table_rho_max:e}).')
            exit(1)
    if T_min is None:
        T_min = 10**intable["logtemp"][0]
    else:
        if T_min < table_T_min:
            print(f'Error: T_max exceeds table bounds ({T_min:e} < {table_T_min:e}).')
            exit(1)
    if T_max is None:
        T_max = 10**intable["logtemp"][-1]
    else:
        if T_max > table_T_max:
            print(f'Error: T_max exceeds table bounds ({T_max:e} > {table_T_max:e}).')
            exit(1)
    if Y_e_min is None:
        Y_e_min = intable["ye"][0]
    else:
        if Y_e_min < table_Y_e_min:
            print(f'Error: Y_e_min exceeds table bounds ({Y_e_min:e} < {table_Y_e_min:e}).')
            exit(1)
    if Y_e_max is None:
        Y_e_max = intable["ye"][-1]
    else:
        if Y_e_max > table_Y_e_max:
            print(f'Error: Y_e_max exceeds table bounds ({Y_e_max:e} > {table_Y_e_max:e}).')
            exit(1)

    return rho_min, rho_max, T_min, T_max, Y_e_min, Y_e_max

def create_simple_table(intablename, outtablename,
                        N_rho, N_T, N_Y_e,
                        rho_min=None, rho_max=None,
                        T_min=None, T_max=None,
                        Y_e_min=None, Y_e_max=None):
    
    # First we read in the input table
    with H5File(intablename, "r") as intable:

        print(intable.keys())
        
        rho_min, rho_max, \
        T_min  , T_max  , \
        Y_e_min, Y_e_max = check_limits(intable, rho_min, rho_max, T_min, T_max, Y_e_min, Y_e_max)

        # Set basic quantities
        dims          = (N_Y_e, N_T, N_rho)
        table_log_rho = linspace(log10(rho_min), log10(rho_max), N_rho)
        table_log_T   = linspace(log10(T_min)  , log10(T_max)  , N_T)
        table_Y_e     = linspace(Y_e_min       , Y_e_max       , N_Y_e)

        # Now prepare for the interpolation
        src_pts   = (intable["ye"][:], intable["logtemp"][:], intable["logrho"][:])
        dest_pts  = meshgrid(table_Y_e, table_log_T, table_log_rho, indexing='ij')
        def interpolate(src_pts, dest_pts, table, key):
            interp = RegularGridInterpolator(src_pts, table[key][:])
            return interp((dest_pts[0], dest_pts[1], dest_pts[2]))

        # Compute table quantities
        table_Abar      = interpolate(src_pts, dest_pts, intable, "Abar"     )
        table_Xa        = interpolate(src_pts, dest_pts, intable, "Xa"       )
        table_Xh        = interpolate(src_pts, dest_pts, intable, "Xh"       )
        table_Xn        = interpolate(src_pts, dest_pts, intable, "Xn"       )
        table_Xp        = interpolate(src_pts, dest_pts, intable, "Xp"       )
        table_Zbar      = interpolate(src_pts, dest_pts, intable, "Zbar"     )
        table_cs2       = interpolate(src_pts, dest_pts, intable, "cs2"      )
        table_dedt      = interpolate(src_pts, dest_pts, intable, "dedt"     )
        table_dpderho   = interpolate(src_pts, dest_pts, intable, "dpderho"  )
        table_dpdrhoe   = interpolate(src_pts, dest_pts, intable, "dpdrhoe"  )
        table_entropy   = interpolate(src_pts, dest_pts, intable, "entropy"  )
        table_gamma     = interpolate(src_pts, dest_pts, intable, "gamma"    )
        table_logenergy = interpolate(src_pts, dest_pts, intable, "logenergy")
        table_logpress  = interpolate(src_pts, dest_pts, intable, "logpress" )
        table_mu_e      = interpolate(src_pts, dest_pts, intable, "mu_e"     )
        table_mu_n      = interpolate(src_pts, dest_pts, intable, "mu_n"     )
        table_mu_p      = interpolate(src_pts, dest_pts, intable, "mu_p"     )
        table_muhat     = interpolate(src_pts, dest_pts, intable, "muhat"    )
        table_munu      = interpolate(src_pts, dest_pts, intable, "munu"     )

        # Print basic information
        print("Input table information:")
        print(f'  Filename = {intablename}')
        print(f'  N_rho    = {intable["pointsrho"][0]}')
        print(f'  N_T      = {intable["pointstemp"][0]}')
        print(f'  N_Y_e    = {intable["pointsye"][0]}')
        print(f'  rho_min  = {10**intable["logrho"][0]:.15e} g/cm3')
        print(f'  rho_max  = {10**intable["logrho"][-1]:.15e} g/cm3')
        print(f'  T_min    = {10**intable["logtemp"][0]:.15e} MeV')
        print(f'  T_max    = {10**intable["logtemp"][-1]:.15e} MeV')
        print(f'  Y_e_min  = {intable["ye"][0]:.15e}')
        print(f'  Y_e_max  = {intable["ye"][-1]:.15e}')
    
        print("\nOutput table information:")
        print(f'  Filename = {outtablename}')
        print(f'  N_rho    = {N_rho}')
        print(f'  N_T      = {N_T}')
        print(f'  N_Y_e    = {N_Y_e}')
        print(f'  rho_min  = {rho_min:.15e} g/cm3')
        print(f'  rho_max  = {rho_max:.15e} g/cm3')
        print(f'  T_min    = {T_min:.15e} MeV')
        print(f'  T_max    = {T_max:.15e} MeV')
        print(f'  Y_e_min  = {Y_e_min:.15e}')
        print(f'  Y_e_max  = {Y_e_max:.15e}')

        # Write information to file
        with open("table_info.txt", "w") as f:
            f.write("Input table information:\n")
            f.write(f'  Filename = {intablename}'+"\n")
            f.write(f'  N_rho    = {intable["pointsrho"][0]}'+"\n")
            f.write(f'  N_T      = {intable["pointstemp"][0]}'+"\n")
            f.write(f'  N_Y_e    = {intable["pointsye"][0]}'+"\n")
            f.write(f'  rho_min  = {10**intable["logrho"][0]:.15e} g/cm3'+"\n")
            f.write(f'  rho_max  = {10**intable["logrho"][-1]:.15e} g/cm3'+"\n")
            f.write(f'  T_min    = {10**intable["logtemp"][0]:.15e} MeV'+"\n")
            f.write(f'  T_max    = {10**intable["logtemp"][-1]:.15e} MeV'+"\n")
            f.write(f'  Y_e_min  = {intable["ye"][0]:.15e}'+"\n")
            f.write(f'  Y_e_max  = {intable["ye"][-1]:.15e}'+"\n")

            f.write("\nOutput table information:\n")
            f.write(f'  Filename = {outtablename}'+"\n")
            f.write(f'  N_rho    = {N_rho}'+"\n")
            f.write(f'  N_T      = {N_T}'+"\n")
            f.write(f'  N_Y_e    = {N_Y_e}'+"\n")
            f.write(f'  rho_min  = {rho_min:.15e} g/cm3'+"\n")
            f.write(f'  rho_max  = {rho_max:.15e} g/cm3'+"\n")
            f.write(f'  T_min    = {T_min:.15e} MeV'+"\n")
            f.write(f'  T_max    = {T_max:.15e} MeV'+"\n")
            f.write(f'  Y_e_min  = {Y_e_min:.15e}'+"\n")
            f.write(f'  Y_e_max  = {Y_e_max:.15e}'+"\n")


    # Now create the output table
    with H5File(outtablename, "w") as f:
        # Create basic datasets
        f.create_dataset("pointsrho"   , shape=1    , dtype="i4", data=N_rho          )
        f.create_dataset("pointstemp"  , shape=1    , dtype="i4", data=N_T            )
        f.create_dataset("pointsye"    , shape=1    , dtype="i4", data=N_Y_e          )
        f.create_dataset("energy_shift", shape=1    , dtype="i4", data=N_Y_e          )
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
    argc = len(argv)
    
    if argc < 2:
        print("ERROR: Need at least an input table")
        exit(1)

    intablename  = argv[1]
    outtablename = intablename.split(".h5")[0]+"_simple.h5"
    outtablename = join(getcwd(), outtablename.split("/")[-1])
    if argc > 2:
        outtablename = argv[2]
        if argc > 3:
            print("Warning: ignoring arguments ",argv[2:])

    N_rho        = 13
    N_T          = 11
    N_Y_e        = 7
    # rho_min      = 1e-12
    # rho_max      = 1e-3
    # T_min        = 1e-2
    # T_max        = 1e+2
    # Y_e_min      = 0.05
    # Y_e_max      = 0.55
    create_simple_table(intablename, outtablename, N_rho, N_T, N_Y_e)#,
    #                    rho_min=rho_min, rho_max=rho_max,
    #                    T_min=T_min, T_max=T_max,
    #                    Y_e_min=Y_e_min, Y_e_max=Y_e_max)
