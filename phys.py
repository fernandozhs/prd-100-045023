# Copy
import copy

# SciPy ODE Solver and Integrator:
from scipy.integrate import odeint
from scipy.integrate import quad

# SciPy B-spline Interpolation:
from scipy.interpolate import interp1d
from scipy.interpolate import splrep

# SciPy Special Functions:
from scipy.special import erfc

# NumPy
import numpy as np

# Tools
from tools import *


def FrequencySample(k, m, Xi, xSample, sSample, s1Sample, K):
    """ Generates a sample of the frequency function `f(x)` along the Real axis for a given field mode.

    Args:
        k: a float indicating the wavenumber of the field mode under consideration.
        m: a float indicating the field mass.
        Xi: a float indicating the field coupling.
        xSample: a (n,)-dimensional NumPy array sampling the points `x` over the Real axis at which at which the mode frequency is to be evaluated.
        sSample: a (n,)-dimensional NumPy array sampling the scale factor `s(x)` along the Real axis.
        s1Sample: a (n,)-dimensional NumPy array sampling the first order derivative of the scale factor `s1(x)` along the Real axis.
        K: an integer equal to `-1`, `0`, or `1` indicating the background curvature.

    Returns:
        An (n,)-dimensional NumPy array sampling the frequency function `f(x)` at the points `xSample` along the Real axis for a the field mode of wavenumber `k`.
    """

    # Computes the a sample of the Hubble parameter `h(xSample)` from the input samples for the scale factor `xSample` and its derivative `s1Sample`.
    hSample = s1Sample/sSample

    # Obtains a sample the first order derivative of the Hubble parameter, `h1(xSample)`, from `hSample` and `xSample`.
    h1 = splrep(xSample[:-1], hSample[:-1], w=np.exp(-100*(xSample[:-1] - xSample[:-1].mean())**2), s=10)
    h1Sample = splev(xSample, h1, der=1)

    # Computes a sample of the Ricci scalar.
    ricciSample = 6.*(h1Sample + 2.*hSample**2 + K/sSample**2)

    # Defines the conformal factor.
    cFactor = Xi - 1.0/6.0

    # Generates the frequency function sample.
    fSample = np.sqrt((k/sSample)**2 + m**2 + cFactor*ricciSample - 0.5*(h1Sample + hSample**2/2.0))

    # Reuturns a sample of the frequency function for the field mode of wavenumber `k`.
    return fSample


def MinkowskiFrequencySample(k, m, sSample):
    """ Generates a sample of the Minkowski frequency function `f(x)` along the Real axis for a given field mode.

    Args:
        k: a float indicating the wavenumber of the field mode under consideration.
        m: a float indicating the field mass.
        sSample: a (n,)-dimensional NumPy array sampling the scale factor `s(x)` along the Real axis.

    Returns:
        An (n,)-dimensional NumPy array sampling the Minkowski frequency function `f(x)` along the Real axis for a the field mode of wavenumber `k`.
    """

    # Generates the intrinsic frequency function sample.
    mfSample = np.sqrt((k/sSample)**2 + m**2)

    # Reuturns a sample of the intrinsic frequency function for the field mode of wavenumber `k`.
    return mfSample


def EventTime(Root, f):
    """ Finds at which point a Stokes line associated with `f` crosses the Real axis.

    Args:
        Root: a NumPy array of the form `np.array([x, y])` specifying the location ` z = x + 1j*y` of a root of `f` on the complex plane. This root and its complex conjugate are connected by a Stokes line.
        f: a callable function with two arguments `f(x, y)` represeting the frequency function over the complex plane for the field mode of interest.

    Returns:
        A float characterizing the point at which the Stokes line crosses the Real axis.
    """

    # Defines the parameter values for which the Stokes line will be constructed. Notice that a small buffer of `1e-3` is added to `yRange` in order to ensure the correct Stokes line is being tracked.
    yRange = np.linspace(-Root[1] + 1e-3, 0., 2)

    # Solves for the Stokes line `x = x(y)`.
    Line = odeint(StokesLine, Root[0], yRange, args=(f,))

    # Returns the value of the `x` coordinate at which the Stokes line crosses the Real axis: `x(0)`.
    return Line[-1]


def EventAmplitude(Root, f):
    """ Determines the amplitude of a partile creation event in given field mode sourced by the crossing of a Stokes line.

    Args:
        Root: a NumPy array of the form `np.array([x, y])` specifying the location ` z = x + 1j*y` of a root of `f` on the complex plane.
        f: a callable function with two arguments `f(x, y)` represeting the frequency function over the complex plane for the field mode of interest.

    Returns:
        A float which characterizes the number of particles produced per comoving volume for a given field mode.
    """

    # Swaps the argument order of `f` to facilitate integration along the complex `y` axis.
    fSwap = lambda y, x: f(x, y)

    # Integrates `f` along a straight line connecting `Root` and the complex conjugate of that point.
    fIntegral = quad(fSwap, -Root[1], Root[1], args=(Root[0]))[0]

    # The event amplitude is finally given by:
    Amplitude = np.exp(-fIntegral)

    # Returns the number of particles produced per comoving volume for a given field mode.
    return Amplitude


def EventWidth(Root, x, f):
    """ Computes the width of a partile creation event in given field mode sourced by the crossing of a Stokes line.

    Args:
        Root: a NumPy array of the form `np.array([x, y])` specifying the location ` z = x + 1j*y` of a root of `f` on the complex plane. This root and its complex conjugate are connected by a Stokes line.
        x: a float which characterizes the `x` coordinate of the point at which the Stokes line crosses the Real axis.
        f: a callable function with two arguments `f(x, y)` represeting the frequency function over the complex plane for the field mode of interest.

    Returns:
        A float which locates the point at which the Stokes line crosses the Real axis.
    """

    # Swaps the argument order of `f` to facilitate integration along the complex `y` axis.
    fSwap = lambda y, x: f(x, y)

    # Integrates `f` along a path connecting `Root` on the upper complex conjugate to `x` on the Real axis.
    fIntegral = quad(fSwap, -Root[1], 0., args=(Root[0]))[0] + quad(f, Root[0], x, args=(0.))[0]

    # Computes and returns the width of a partile creation event.
    Width = 2.*f(x,0)/np.sqrt(2.*fIntegral)

    return Width


def EventPhase(x, x0, f):
    """ Computes the relative phase between partile creation events in given field mode.

    Args:
        x: a float which characterizes the `x` coordinate of the point at which the Stokes line crosses the Real axis.
        x0: a float which characterizes the `x0` coordinate of the point at which the reference Stokes line crosses the Real axis.
        f: a callable function with two arguments `f(x, y)` represeting the frequency function over the complex plane for the field mode of interest.

    Returns:
        A float which is the phase between a particle creation event and a reference particle creation event in a given field mode.
    """

    # Integrates the mode frequency function `f` along the Real axis interval `[x0, x]`.
    fIntegral = quad(f, x0, x, args=(0.))[0]

    # Computes and returns the phase between the partile creation events under consideration.
    Phase = np.exp(-2.j*fIntegral)

    return Phase


def BogolyubovBeta(x, modeRecord):
    """ Computes the Bogolyubov coefficient `Beta` for a given field mode from the information stored in `modeRecord`.

    Args:
        x: a float which characterizes the `x` time coordinate on the Real axis.
        modeRecord: the dictionary entry for a particular field mode in `fieldRecord`.

    Returns:
        The Bogolyubov coefficient `Beta(x)` per comoving volume at the time `x`.
    """

    # Initializes the Bogolyubov coefficient.
    Beta = 0.

    # Loops over each Stokes line crossing stored in `modeRecord`.
    for i in range(0, len(modeRecord['root'])):

        # Extracts the quantities characterizing the particle production event.
        eventAmplitude = modeRecord['eventAmplitude'][i]
        eventPhase = modeRecord['eventPhase'][i]
        eventTime = modeRecord['eventTime'][i]
        eventWidth = modeRecord['eventWidth'][i]

        # Incorporates the `i-th` Stokes line crossing to the Bogolyubov coefficient `Beta` of the field mode under consideration.
        Beta += -1.j*(eventAmplitude/2.)*erfc(eventWidth*(eventTime - x))*eventPhase

    # Returns the Bogolyubov coefficient `Beta` at `x`.
    return Beta


def BogolyubovAlpha(x, modeRecord):
    """ Computes the Bogolyubov coefficient `Alpha` for a given field mode from the information stored in `modeRecord`.

    Args:
        x: a float which characterizes the `x` time coordinate on the Real axis.
        modeRecord: the dictionary entry for a particular field mode in `fieldRecord`.

    Returns:
        The Bogolyubov coefficient `Alpha(x)` per comoving volume at the time `x`.
    """

    # Computes the Bogolyubov coefficient `Beta` at `x`.
    Beta= BogolyubovBeta(x, modeRecord)

    # Computes the Bogolyubov coefficient `Alpha` at `x`, up to a small undetermined phase.
    Alpha = np.sqrt(1 + np.abs(Beta)**2)

    # Returns the Bogolyubov coefficient `Alpha` at `x`.
    return Alpha


def AdiabaticR(x, fieldRecord, modeRecord):
    """ Computes the adiabatic quantity `R` for a given field mode from the information stored in `fieldRecord` and `modeRecord`.

    Args:
        x: a float which characterizes the `x` time coordinate on the Real axis.
        fieldRecord: a dictionary containing a dictionary entry for each field mode.
        modeRecord: the dictionary entry for a particular field mode in `fieldRecord`.

    Returns:
        The adiabatic quantity `R(x)` per comoving volume at the time `x`.
    """

    # Extracts the initial comoving number of adiabatic particles in the field mode under consideration.
    initialN = modeRecord['initialN']

    # Computes the Bogolyubov coefficients `Alpha` and `Beta` at `x`.
    Beta = BogolyubovBeta(x, modeRecord)
    Alpha = BogolyubovAlpha(x, modeRecord)

    # Constructs the function `R`.
    R = (1. + 2.*initialN)*np.real(Alpha*np.conjugate(Beta))

    # Returns the function `R(x)` computed at `x`.
    return R


def AdiabaticI(x, fieldRecord, modeRecord):
    """ Computes the adiabatic quantity `I` for a given field mode from the information stored in `fieldRecord` and `modeRecord`.

    Args:
        x: a float which characterizes the `x` time coordinate on the Real axis.
        fieldRecord: a dictionary containing a dictionary entry for each field mode.
        modeRecord: the dictionary entry for a particular field mode in `fieldRecord`.

    Returns:
        The adiabatic quantity `I(x)` per comoving volume at the time `x`.
    """

    # Extracts the initial comoving number of adiabatic particles in the field mode under consideration.
    initialN = modeRecord['initialN']

    # Computes the Bogolyubov coefficients `Alpha` and `Beta` at `x`.
    Beta = BogolyubovBeta(x, modeRecord)
    Alpha = BogolyubovAlpha(x, modeRecord)

    # Constructs the function `I`.
    I = (1. + 2.*initialN)*np.imag(Alpha*np.conjugate(Beta))

    # Returns the function `I(x)` computed at `x`.
    return I


def AdiabaticN(x, modeRecord):
    """ Assembles all particle creation events stored in `modeRecord` to create a function which tracks the time dependence of the number of adiabatic particles populating a given field mode.

    Args:
        x: a float which characterizes the `x` time coordinate on the Real axis.
        modeRecord: the dictionary entry for a particular field mode in `fieldRecord`.

    Returns:
        The total number of adiabatic particles `N(x)` per comoving volume at the time `x`.
    """

    # Extracts the initial comoving number of adiabatic particles in the field mode under consideration.
    initialN = modeRecord['initialN']

    # Computes the Bogolyubov coefficient `Beta` at `x`.
    Beta = BogolyubovBeta(x, modeRecord)

    # Computes the number of adiabatic particles populating the field mode of interest at time `x`.
    N = initialN + (1. + 2.*initialN)*np.abs(Beta)**2

    # Returns the total number of adiabatic particles at `x`.
    return N


def EnergyDensityN(x, s, fieldRecord, modeRecord, metricRecord):
    """ Computes the contribution to the the energy density due to the adiabatic particles populating a given field mode.

    Args:
        x: a float which characterizes the `x` time coordinate on the Real axis.
        s: a float specifying the scale factor `s(x)` at the time `x`.
        fieldRecord: a dictionary containing a dictionary entry for each field mode.
        modeRecord: the dictionary entry for a particular field mode in `fieldRecord`.
        metricRecord: a metric dictionary containing NumPy arrays which detail the evolution of the spacetime metric.

    Returns:
        The portion of the field energy density due to the adiabatic particles populating a given field mode.
    """

    # Extracts `Xi` and `K` from `fieldRecord` and `metricRecord`.
    Xi = fieldRecord['Xi']
    K = metricRecord['K']

    # Extracts `WSample`, `VSample`, and `mfSample` from `modeRecord`.
    WSample = modeRecord['phaseIntegralW']
    VSample = modeRecord['phaseIntegralV']
    mfSample = modeRecord['minkowskiFrequency']

    # Extracts `hSample` from `metricRecord`.
    hSample = metricRecord['s1']/metricRecord['s']

    # Creates a linear interpolation of the functions `W(x)`, `V(x)`, `mf(x)`, and `h(x)` over the Real axis.
    W = interp1d(fieldRecord['x'], WSample, fill_value="extrapolate")
    V = interp1d(fieldRecord['x'], VSample, fill_value="extrapolate")
    mf = interp1d(fieldRecord['x'], mfSample, fill_value="extrapolate")
    h = interp1d(metricRecord['x'], hSample, fill_value="extrapolate")

    # Constructs the function `rhoN(x)`.
    rhoN = ( W(x)**2 + mf(x)**2 + 0.25*(V(x) - h(x))**2 + (6.*Xi - 1.)*(h(x)*V(x) + K/s**2 - 2.*h(x)**2) )/W(x)

    # Returns `rhoN(x)` computed at the time `x`.
    return rhoN


def EnergyDensityR(x, s, fieldRecord, modeRecord, metricRecord):
    """ Computes the contribution to the the energy density due to the interference term `R` of a given field mode.

    Args:
        x: a float which characterizes the `x` time coordinate on the Real axis.
        s: a float specifying the scale factor `s(x)` at the time `x`.
        fieldRecord: a dictionary containing a dictionary entry for each field mode.
        modeRecord: the dictionary entry for a particular field mode in `fieldRecord`.
        metricRecord: a metric dictionary containing NumPy arrays which detail the evolution of the spacetime metric.

    Returns:
        The portion of the field energy density due to the interference term `R` of a given field mode.
    """

    # Extracts `Xi` and `K` from `fieldRecord` and `metricRecord`.
    Xi = fieldRecord['Xi']
    K = metricRecord['K']

    # Extracts `WSample`, `VSample`, and `mfSample` from `modeRecord`.
    WSample = modeRecord['phaseIntegralW']
    VSample = modeRecord['phaseIntegralV']
    mfSample = modeRecord['minkowskiFrequency']

    # Extracts `hSample` from `metricRecord`.
    hSample = metricRecord['s1']/metricRecord['s']

    # Creates a linear interpolation of the functions `W(x)`, `V(x)`, `mf(x)`, and `h(x)` over the Real axis.
    W = interp1d(fieldRecord['x'], WSample, fill_value="extrapolate")
    V = interp1d(fieldRecord['x'], VSample, fill_value="extrapolate")
    mf = interp1d(fieldRecord['x'], mfSample, fill_value="extrapolate")
    h = interp1d(metricRecord['x'], hSample, fill_value="extrapolate")

    # Constructs the function `rhoR(x)`.
    rhoR = ( -W(x)**2 + mf(x)**2 + 0.25*(V(x) - h(x))**2 + (6.*Xi - 1.)*(h(x)*V(x) + K/s**2 - 2.*h(x)**2) )/W(x)

    # Returns `rhoR(x)` computed at the time `x`.
    return rhoR


def EnergyDensityI(x, s, fieldRecord, modeRecord, metricRecord):
    """ Computes the contribution to the the energy density due to the interference term `I` of a given field mode.

    Args:
        x: a float which characterizes the `x` time coordinate on the Real axis.
        s: a float specifying the scale factor `s(x)` at the time `x`.
        fieldRecord: a dictionary containing a dictionary entry for each field mode.
        modeRecord: the dictionary entry for a particular field mode in `fieldRecord`.
        metricRecord: a metric dictionary containing NumPy arrays which detail the evolution of the spacetime metric.

    Returns:
        The portion of the field energy density due to the interference term `I` of a given field mode.
    """

    # Extracts `Xi` and `K` from `fieldRecord` and `metricRecord`.
    Xi = fieldRecord['Xi']
    K = metricRecord['K']

    # Extracts `VSample` from `modeRecord`.
    VSample = modeRecord['phaseIntegralV']

    # Extracts `hSample` from `metricRecord`.
    hSample = metricRecord['s1']/metricRecord['s']

    # Creates a linear interpolation of the functions `V(x)` and `h(x)` over the Real axis.
    V = interp1d(fieldRecord['x'], VSample, fill_value="extrapolate")
    h = interp1d(metricRecord['x'], hSample, fill_value="extrapolate")

    # Constructs the function `rhoI(x)`.
    rhoI = V(x) - h(x) + 2.*(6.*Xi - 1.)*h(x)

    # Returns `rhoI(x)` computed at the time `x`.
    return rhoI


def EnergyDensity(x, s, fieldRecord, metricRecord):
    """ Assembles all contributions to the field energy density to create a function which tracks the time dependence of the field total energy density.

    Args:
        x: a float which characterizes the `x` time coordinate on the Real axis.
        s: a float specifying the scale factor `s(x)` at the time `x`.
        fieldRecord: a dictionary containing a dictionary entry for each field mode.

    Returns:
        The total field energy density at the time `x` due to all particles populating field modes.
    """

    # Initializes the field energy density due to adiabatic particles.
    rho = 0.

    # Loops over all field modes being tracked by `fieldRecord`.
    for k in fieldRecord.keys()[:-5]:

        # Creates a copy of the `modeRecord` associated with the mode of wavenumber `k`.
        modeRecord = copy.deepcopy(fieldRecord[k])

        # Extracts the mode wavenumber from `modeRecord`.
        q = modeRecord['k']

        # Computes the adiabatic quantities `N`, `R`, and `I` for the mode `k` at time `x`.
        N = AdiabaticN(x, modeRecord)
        R = AdiabaticR(x, fieldRecord, modeRecord)
        I = AdiabaticI(x, fieldRecord, modeRecord)

        # Computes the adiabatic quantities `rhoN`, `rhoR`, and `rhoI` for the mode `k` at time `x`.
        rhoN = EnergyDensityN(x, s, fieldRecord, modeRecord, metricRecord)
        rhoR = EnergyDensityR(x, s, fieldRecord, modeRecord, metricRecord)
        rhoI = EnergyDensityI(x, s, fieldRecord, modeRecord, metricRecord)

        # Incorporates all contributions from the mode of wavenumber `k` to the the field energy density.
        rho += (q**2 * ( rhoN*N + rhoR*R + rhoI*I ))/(4. * s**3 * np.pi**2)

    # Returns the total field energy density at `x`.
    return rho


def EnergyDensityRI(x, s, fieldRecord, metricRecord):
    """ Computes all contributions to the field energy density due to the interference terms `R` and `I`.

    Args:
        x: a float which characterizes the `x` time coordinate on the Real axis.
        s: a float specifying the scale factor `s(x)` at the time `x`.
        fieldRecord: a dictionary containing a dictionary entry for each field mode.

    Returns:
        The total contribution to field energy density at the time `x` due to the interference terms `R` and `I`.
    """
    
    # Initializes the field energy density due to the interference terms `R` and `I`.
    rhoRI = 0.

    # Loops over all field modes being tracked by `fieldRecord`.
    for k in fieldRecord.keys()[:-5]:

        # Creates a copy of the `modeRecord` associated with the mode of wavenumber `k`.
        modeRecord = copy.deepcopy(fieldRecord[k])

        # Extracts the mode wavenumber from `modeRecord`.
        q = modeRecord['k']

        # Computes the adiabatic quantities `R`, and `I` for the mode `k` at time `x`.
        R = AdiabaticR(x, fieldRecord, modeRecord)
        I = AdiabaticI(x, fieldRecord, modeRecord)

        # Computes the adiabatic quantities `rhoR` and `rhoI` for the mode `k` at time `x`.
        rhoR = EnergyDensityR(x, s, fieldRecord, modeRecord, metricRecord)
        rhoI = EnergyDensityI(x, s, fieldRecord, modeRecord, metricRecord)

        # Incorporates all contributions from the mode of wavenumber `k` to the the field energy density.
        rhoRI += (q**2 * ( rhoR*R + rhoI*I ))/(4. * s**3 * np.pi**2)

    # Returns the total field energy density due to the interference terms `R` and `I` at `x`.
    return rhoRI


def EnergyDensityVacuum(x, s, fieldRecord, modeRecord, metricRecord):
    """ Computes the contribution to the the energy density due to the adiabatic vacuum-like term of a given field mode.

    Args:
        x: a float which characterizes the `x` time coordinate on the Real axis.
        s: a float specifying the scale factor `s(x)` at the time `x`.
        fieldRecord: a dictionary containing a dictionary entry for each field mode.
        modeRecord: the dictionary entry for a particular field mode in `fieldRecord`.
        metricRecord: a metric dictionary containing NumPy arrays which detail the evolution of the spacetime metric.

    Returns:
        The portion of the field energy density due to the adiabatic vacuum-like term of a given field mode.

    WARNING:
        This function is valid only for massive conformally coupled fields. A more general version still needs to be implemented.
    """

    # Extracts `Xi` and `K` from `fieldRecord` and `metricRecord`.
    Xi = fieldRecord['Xi']
    K = metricRecord['K']

    # Extracts `WSample`, `VSample`, and `mfSample` from `modeRecord`.
    WSample = modeRecord['phaseIntegralW']
    VSample = modeRecord['phaseIntegralV']
    mfSample = modeRecord['minkowskiFrequency']

    # Extracts `hSample` from `metricRecord`.
    hSample = metricRecord['s1']/metricRecord['s']

    # Creates a linear interpolation of the functions `W(x)`, `V(x)`, `mf(x)`, and `h(x)` over the Real axis.
    W = interp1d(fieldRecord['x'], WSample, fill_value="extrapolate")
    V = interp1d(fieldRecord['x'], VSample, fill_value="extrapolate")
    mf = interp1d(fieldRecord['x'], mfSample, fill_value="extrapolate")
    h = interp1d(metricRecord['x'], hSample, fill_value="extrapolate")

    # Constructs the function `rhoVacuum(x)`.
    rhoVacuum = ( (W(x) - mf(x))**2 + 0.25*(V(x) - h(x))**2 )/W(x) - 0.25*h(x)**2/mf(x)**5

    # Returns `rhoVacuum(x)` computed at the time `x`.
    return rhoVacuum


def EnergyDensityVac(x, s, fieldRecord, metricRecord):
    """ Computes the contribution to the field energy density due the adiabatic vacuum-like term.

    Args:
        x: a float which characterizes the `x` time coordinate on the Real axis.
        s: a float specifying the scale factor `s(x)` at the time `x`.
        fieldRecord: a dictionary containing a dictionary entry for each field mode.

    Returns:
        The total contribution to field energy density at the time `x` due to the adiabatic vacuum-like term.

    WARNING:
        This function is valid only for massive conformally coupled fields. A more general version still needs to be implemented.
    """

    # Initializes the field energy density due to adiabatic particles.
    rhoVac = 0.

    # Loops over all field modes being tracked by `fieldRecord`.
    for k in fieldRecord.keys()[:-5]:

        # Creates a copy of the `modeRecord` associated with the mode of wavenumber `k`.
        modeRecord = copy.deepcopy(fieldRecord[k])

        # Extracts the mode wavenumber from `modeRecord`.
        q = modeRecord['k']

        # Computes the adiabatic quantities `R`, and `I` for the mode `k` at time `x`.
        R = AdiabaticR(x, fieldRecord, modeRecord)
        I = AdiabaticI(x, fieldRecord, modeRecord)

        # Computes the adiabatic quantities `rhoN`, `rhoR`, and `rhoI` for the mode `k` at time `x`.
        rhoVacuum = EnergyDensityVacuum(x, s, fieldRecord, modeRecord, metricRecord)

        # Incorporates all contributions from the mode of wavenumber `k` to the the field energy density.
        rhoVac += (q**2 * ( rhoVacuum ))/(8. * s**3 * np.pi**2)

    # Returns the total field energy density at `x`.
    return rhoVac


def Friedmann(s, x, rho, Lambda, K, planckRatio):
    """ Defines the Friedmann equation.

    Args:
        s: a float specifying the scale factor `s(x)` at the time `x`.
        x: a float specifying the time coordinate of a point `x` along the Real axis.
        rho: a callable function `rho(x)` which outputs a the field energy density as a function of the time variable `x`.
        Lambda: a float characterizing the cosmological constant.
        K: an integer equal to `-1`, `0`, or `1` indicating the background curvature.
        planckRatio: a float specifying the ratio between the field mass and the Planck mass.

    Returns:
        A float which characterizes `s1(x)`, i.e., the first derivative of the scale factor `s(x)` with respect to the time coordinate computed at the point `x`.
    """

    # Defines the Planck mass scale in the Friedmann equation.
    massScale = (8. * np.pi * planckRatio**2)/3.

    # Defines the Friedmann equation, which is the differential equation satisfied by the scale factor `s(x)` describing the Universe's expansion.
    s1 = -np.sqrt((massScale * rho(x) * s**2) + (Lambda * s**2)/3. - K)

    # Returns the left-hand side of the differential equations.
    return s1
