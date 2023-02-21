""" Assert SymPy Expression Equality """
# Author: Ken Sible
# Email:  ksible *at* outlook *dot* com

from UnitTesting.standard_constants import precision
from mpmath import mp, mpf, mpc, fabs, log10
import sympy as sp, sys, random

def expand_vardict(vardict):
    if all(not isinstance(vardict[var], list) for var in vardict):
        return vardict
    for var in vardict:
        if not isinstance(vardict[var], list):
            vardict[var] = [vardict[var]]
    return expand_vardict({(var + '[' + str(i) + ']') if len(vardict[var]) > 1 else var: vardict[var][i]
            for var in vardict for i in range(len(vardict[var]))})

def compute_value(symdict, replaced, reduced, factor):
    # determine the precision for floating point arithmetic
    mp.dps = factor * precision
    # substitute unique random number for each free symbol in symdict
    # into every common subexpression from CSE and into main expression
    for var, subexpr in replaced:
        for symbol in subexpr.free_symbols:
            subexpr = subexpr.subs(symbol, symdict[symbol])
        if subexpr.atoms(sp.pi):
            subexpr = subexpr.subs(sp.pi, mp.pi)
        symdict[var] = subexpr
    reduced = reduced[0]
    for symbol in reduced.free_symbols:
        reduced = reduced.subs(symbol, symdict[symbol])
        if reduced.atoms(sp.pi):
            reduced = reduced.subs(sp.pi, mp.pi)
    reduced = reduced.subs(sp.Function('NRPyAbs'), sp.Abs)
    value = mpc(sp.N(reduced)) if isinstance(reduced, complex) else mpf(reduced)
    mp.dps = precision
    return value

def update_vardict(vardict):
    # expand vardict into component mapping
    vardict = expand_vardict(vardict)
    # extract every free symbol present in vardict
    symdict = {symbol: None for var in vardict if isinstance(vardict[var], sp.Basic)
                   for symbol in vardict[var].free_symbols}
    for free_symbol in symdict:
        # seed random number generator with free_symbol hash value
        random.seed(hash(free_symbol))
        # update symdict with mapping: free_symbol -> unique random number
        symdict[free_symbol] = mpf(random.random())
    for var in vardict:
        # apply CSE to every expression in vardict
        replaced, reduced = sp.cse(vardict[var], order='none')
        # calculate value after substituting the unique random number
        # from each free symbol in symdict into every expression in vardict
        value = compute_value(symdict, replaced, reduced, factor=1)
        # double the precision (factor = 2) whenever value within range of zero
        if fabs(value) != mpf('0.0') and fabs(value) < 10 ** ((-2.0/3) * precision):
            _value = compute_value(symdict, replaced, reduced, factor=2)
            if fabs(_value) < 10 ** (-(4.0/3) * precision):
                value = mpf('0.0')
        # update vardict with mapping: variable -> (pseudo) unique number
        vardict[var] = value
    return vardict

def assert_equal(vardict_1, vardict_2, suppress_message=False):
    """ Assert SymPy Expression Equality

        >>> from sympy import sin, cos
        >>> from sympy.abc import x

        >>> assert_equal(sin(2*x), 2*sin(x)*cos(x))
        Assertion Passed!
        >>> assert_equal(cos(2*x), cos(x)**2 - sin(x)**2)
        Assertion Passed!

        >>> assert_equal(cos(2*x), 1 - 2*sin(x)**2)
        Assertion Passed!
        >>> assert_equal(cos(2*x), 1 + 2*sin(x)**2)
        Traceback (most recent call last):
        ...
        AssertionError

        >>> vardict_1 = {'A': sin(2*x), 'B': cos(2*x)}
        >>> vardict_2 = {'A': 2*sin(x)*cos(x), 'B': cos(x)**2 - sin(x)**2}
        >>> assert_equal(vardict_1, vardict_2)
        Assertion Passed!

        >>> assert_equal('(a^2 - b^2) - (a + b)*(a - b)', 0)
        Assertion Passed!
    """
    if not isinstance(vardict_1, dict):
        vardict_1 = {'': vardict_1}
    if not isinstance(vardict_2, dict):
        vardict_2 = {'': vardict_2}
    for var_1, var_2 in zip(vardict_1, vardict_2):
        if not isinstance(vardict_1[var_1], sp.Basic):
            vardict_1[var_1] = sp.sympify(vardict_1[var_1])
        if not isinstance(vardict_2[var_2], sp.Basic):
            vardict_2[var_2] = sp.sympify(vardict_2[var_2])
    # update each vardict with mapping: variable -> (pseudo) unique number
    vardict_1, vardict_2 = update_vardict(vardict_1), update_vardict(vardict_2)
    # assert whether SDA >= (2/3) * precision, implying expression equality
    for var_1, var_2 in zip(vardict_1, vardict_2):
        n_1, n_2 = vardict_1[var_1], vardict_2[var_2]
        if n_1 == n_2: continue
        E_rel = 2 * fabs(n_1 - n_2)/(fabs(n_1) + fabs(n_2))
        assert -log10(E_rel) + 1 >= (2.0/3) * precision
    if not suppress_message:
        print('Assertion Passed!')


# This function checks whether an expression evaluates to zero.
import hashlib
from outputC import cse_postprocess
from UnitTesting.cse_simplify_and_evaluate_sympy_expressions import calculate_value
def check_zero(expression, verbose=False):
    # Setting precision
    precision = 30 # 30 sig digits
    mp.dps = precision

    free_symbols_dict = dict()

    # Setting each variable in free_symbols_set to a random number in [0, 1) according to the hashed string
    # representation of each variable.
    for var in expression.free_symbols:
        # Make sure M_PI is set to its correct value, pi
        if str(var) == "M_PI":
            free_symbols_dict[var] = mp.mpf(pi)
        # Then make sure M_SQRT1_2 is set to its correct value, 1/sqrt(2)
        elif str(var) == "M_SQRT1_2":
            free_symbols_dict[var] = mp.mpf(1.0/sqrt(2.0))
        # All other free variables are set to random numbers
        else:
            # Take the variable [var], turn it into a string, encode the string, hash the string using the md5
            # algorithm, turn the hash into a hex number, turn the hex number into an int, set the random seed to
            # that int. This ensures each variable gets a unique but consistent value.
            random.seed(int(hashlib.md5(str(var).encode()).hexdigest(), 16))
            # Store the random value in free_symbols_dict as a mpf
            free_symbols_dict[var] = mp.mpf(random.random())
            # Warning: might slow Travis CI too much: logging.debug(' ...Setting '+str(var)+' to the random value: '+str(free_symbols_dict[var]))

    # Using SymPy's cse algorithm to optimize our value substitution
    replaced, reduced = cse_postprocess(sp.cse(expression, order='none'))

    # Warning: might slow Travis CI too much: logging.debug(' var = '+str(var)+' |||| replaced = '+str(replaced))

    # Calculate our result_value
    result_value = calculate_value(free_symbols_dict, replaced, reduced)

    # Check if result_value is near-zero, and double checking if it should be zero
    if mp.fabs(result_value) != mp.mpf('0.0') and mp.fabs(result_value) < 10 ** ((-2.0/3)*precision):
        if verbose:
            print("Found |result| (" + str(mp.fabs(result_value)) + ") close to zero. "
                         "Checking if indeed it should be zero.")
        new_result_value = calculate_value(free_symbols_dict, replaced, reduced, precision_factor=2)
        if mp.fabs(new_result_value) < 10 ** (-(4.0/3) * precision):
            if verbose:
                print("After re-evaluating with twice the digits of precision, |result| dropped to " +
                             str(new_result_value) + ". Setting value to zero")
            result_value = mp.mpf('0.0')

    # Store result_value in calculated_dict
    if result_value == mp.mpf('0.0'):
        return True
    return False

if __name__ == "__main__":
    import doctest
    sys.exit(doctest.testmod()[0])
