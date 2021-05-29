# NumPy
import numpy as np

# Control
from control import *


def CreateBox(x, xLength, yLength, xPoints, yPoints):
    """ Creates a box bisected by the Real axis which delimits a rectangular region of the complex plane over which other functions will operate.
        
    Args:
        x: a float specifying the reference Real axis coordinate.
        xLength: a float specifying the box's length along the `x` direction.
        yLength: a float specifying the box's length along the `y` direction.
        xPoints: the integer number of points along the `x` direction.
        yPoints: the integer number of points along the `y` direction.
        
    Returns:
        A tuple of the form `(xRange, yRange)`, where the entries is are 1-dimensional NumPy arrays containing `xPoints` and `yPoints` elements, respectively. The array `xRange` spams the interval `[x, x + xLength]` and `yRange` spams the interval `[-yLength/2., yLength/2.]`.
        
    Examples:
        >>> InitBox(x=1., xLength=1., yLength=4., xPoints=2, yPoints=5)
        >>> (np.array([1., 2.]), np.array([-2., -1, 0., 1., 2.]))
        
    """
    
    # The number of points `yPoints` along the `y` direction is incremented by `1` if its input value is an even number. This is necessary in order to make sure the Real axis is included as a subset of the output `(xRange, yRange)`.
    yPoints = yPoints - yPoints%2 + 1
    
    # Generates the arrays spamming the `x` and `y` directions.
    xRange = np.linspace(x, x + xLength, xPoints)
    yRange = np.linspace(-yLength/2., yLength/2., yPoints)
    
    # Returns the box.
    return (xRange, yRange)


def ScaleFactor(xRange, H, K):
    """ Sets the initial conditions for scale factor evolution within the initial range `xRange`.
        
    Args:
        xRange: an (m,)-dimensional NumPy array covering the initial patch over the Real axis.
        H: a float specifying the scale factor's initial expansion rate.
        K: an integer equal to `-1`, `0`, or `1` indicating the background curvature.
        
    Returns:
        An (m,)-dimensional NumPy array storing the initial scale factor evolution within `Box`.
        
    """

    # Closed Universe.
    if K == 1:
        return np.cosh(H*xRange)/H

    # Flat Universe.
    if K == 0:
        return np.exp(H*xRange)

    # Open Universe.
    if K == -1:
        return np.sinh(H*xRange)/H


def HubbleParameter(xRange, H, K):
    """ Sets the initial conditions for Hubble parameter evolution within the initial range `xRange`.
        
    Args:
        xRange: an (m,)-dimensional NumPy array covering the initial patch over the Real axis.
        H: a float specifying the scale factor's initial expansion rate.
        K: an integer equal to `-1`, `0`, or `1` indicating the background curvature.
        
    Returns:
        An (m,)-dimensional NumPy array storing the initial Hubble parameter evolution within `Box`.
        
    """
    
    # Closed Universe.
    if K == 1:
        return H*np.tanh(H*xRange)
    
    # Flat Universe.
    if K == 0:
        return H
    
    # Open Universe.
    if K == -1:
        return H/np.tanh(H*xRange)


def Hubble1Parameter(xRange, H, K):
    """ Sets the initial conditions for the evolution of the first derivative of the Hubble parameter within the initial range `xRange`.
        
    Args:
        xRange: an (m,)-dimensional NumPy array covering the initial patch over the Real axis.
        H: a float specifying the scale factor's initial expansion rate.
        K: an integer equal to `-1`, `0`, or `1` indicating the background curvature.
        
    Returns:
        An (m,)-dimensional NumPy array storing the first derivative of the initial Hubble parameter evolution within `Box`.
        
    """
    
    # Closed Universe.
    if K == 1:
        return H**2/np.cosh(H*xRange)**2
    
    # Flat Universe.
    if K == 0:
        return 0.0
    
    # Open Universe.
    if K == -1:
        return -H**2/(np.cosh(H*xRange)*np.tanh(H*xRange))**2


def RicciScalar(xRange, H, K):
    """ Sets the initial conditions for Ricci scalar evolution within the initial range `xRange`.
        
    Args:
        xRange: an (m,)-dimensional NumPy array covering the initial patch over the Real axis.
        H: a float specifying the scale factor's initial expansion rate.
        K: an integer equal to `-1`, `0`, or `1` indicating the background curvature.
        
    Returns:
        An (m,)-dimensional NumPy array storing the initial Ricci scalar evolution within `Box`.
        
    """
    
    # Computes the initial scale factor evolution within the initial range `xRange`.
    s = ScaleFactor(xRange, H, K)
    
    # Computes the initial Hubble parameter evolution within the initial range `xRange`.
    h = HubbleParameter(xRange, H, K)
    
    # Computes the derivative of the initial Hubble parameter evolution within the initial range `xRange`.
    h1 = Hubble1Parameter(xRange, H, K)

    # Calculates the Ricci scalar within the initial range `xRange`.
    Ricci = 6.*(h1 + 2.*h**2 + K/s**2)

    # Returns the Ricci scalar within `xRange`.
    return Ricci


def FrequencySquared(x, y, k, m, Xi, H, K):
    """ Computes the initial frequency function `R(x)` at the point `z = x + 1j*y` of the complex plane for a given field mode of wavenumber `k`.
        
    Args:
        x: a float specifying the Real coordinate of the point where the mode frequency function is to be evaluated.
        y: a float specifying the Imaginary coordinate of the point where the mode frequency function is to be evaluated.
        k: a float indicating the wavenumber of the field mode under consideration.
        m: a float indicating the field mass.
        Xi: a float indicating the field coupling.
        H: a float specifying the scale factor's initial expansion rate.
        K: an integer equal to `-1`, `0`, or `1` indicating the background curvature.
        
    Returns:
        An (m,n)-dimensional NumPy array storing the frequency function `R(x)` of the field mode of wavenumber `k` at the point `z = x + 1j*y` on the complex plane.
        
    """

    # Constructs the point in the complex plane from its input Real `x` and Imaginary `y` coordinates.
    z = x + 1j*y

    # Evaluates the initial scale factor at the point `z = x + 1j*y` of the complex plane.
    s = ScaleFactor(z, H, K)
    
    # Computes the initial Hubble parameter at the point `z = x + 1j*y` of the complex plane.
    h = HubbleParameter(z, H, K)
    
    # Calculates the derivative of the initial Hubble parameter at the point `z = x + 1j*y` of the complex plane.
    h1 = Hubble1Parameter(z, H, K)
    
    # Evaluates the initial Ricci scalar at the point `z = x + 1j*y` of the complex plane.
    ricci = RicciScalar(z, H, K)
    
    # Defines the conformal factor.
    cFactor = Xi - 1.0/6.0
    
    # Computes the frequency function sample at `z = x + 1j*y`.
    frequencySquared = (k/s)**2 + m**2 + cFactor*ricci - 0.5*(h1 + h**2/2.0)
    
    # Reuturns a sample of the frequency function squared for the field mode of wavenumber `k`.
    return frequencySquared


def MinkowskiFrequency(x, k, m, H, K):
    """ Computes the initial Minkowski frequency function `w(x)` at the point `x` for a given field mode of wavenumber `k`.
        
    Args:
        x: a float specifying the Real coordinate of the point where the mode frequency function is to be evaluated.
        k: a float indicating the wavenumber of the field mode under consideration.
        m: a float indicating the field mass.
        H: a float specifying the scale factor's initial expansion rate.
        K: an integer equal to `-1`, `0`, or `1` indicating the background curvature.
        
    Returns:
        An (n,)-dimensional NumPy array storing the Minkowski frequency function `w(x)` of the field mode of wavenumber `k` at Real axis point `x`.
        
    """
    
    # Evaluates the initial scale factor at the Real axis point `x`.
    s = ScaleFactor(x, H, K)
    
    # Computes the Minkowski frequency function sample at `x`.
    minkowskiFrequency = np.sqrt((k/s)**2 + m**2)
    
    # Reuturns a sample of the Minkowski frequency function for the field mode of wavenumber `k`.
    return minkowskiFrequency


def InitField(m, planckRatio, Xi, waveNumbers, N, f, Box, H, K):
    """ Creates a dictionary containing the initial conditions for each field mode. Each mode has a corresponding entry which is also dictionary. The initial conditions are set by `f`.
        
    Args:
        m: a float indicating the field mass.
        Xi: a float indicating the field coupling.
        waveNumbers: a (n,)-dimensional NumPy array containing the field wave numbers `k` to be considered.
        N: a callable function `N(k)` which outputs the initial adiabatic particle distribution associated with the field.
        f: a callable function `f(x, y, k, m, Xi, H, K)` which outputs a mode's frequency squared defined over the complex plane.
        Box: a tuple of the form `(xRange, yRange)`, where each entry is a (m,)-dimensional NumPy array spamming one of the box dimensions.
        H: a float specifying the scale factor's initial expansion rate.
        K: an integer equal to `-1`, `0`, or `1` indicating the background curvature.
        
    Returns:
        A dictionary containing a dictionary entry for each field mode, as well as entries for the complex plane variables `x` and `y`, the field mass `m`, and coupling `Xi`.
        
    """

    # Creates a box with the appropriate dimensions for setting the field initial conditions. The box spams the interval `[Box[0].min(), CropBox(Box)[0].max()]` in the `x` direction, and spams the same interval as `Box` along the `y` direction.
    xRange = Box[0][Box[0] <= CropBox(Box)[0].max()]
    yRange = Box[1]
    
    # Creates a mesh grid over the complex plane delimited by `xRange` and `yRange`.
    y, x = np.meshgrid(yRange, xRange, indexing='xy')

    # Constructs a dictionary template for an individual field mode.
    modeTemplate = {'k':0, 'frequency':np.zeros(y.shape, dtype='complex'), 'root':[], 'eventTime':[], 'eventAmplitude':[], 'eventPhase':[], 'eventWidth':[], 'initialN':0, 'minkowskiFrequency':np.zeros(xRange.shape), 'phaseIntegralW':np.zeros(xRange.shape), 'phaseIntegralV':np.zeros(xRange.shape)}

    # Creates an empty dictionary which will carry a `modeTemplate` for each field mode.
    fieldRecord = {}

    # Initializes each mode's frequency record.
    for i, k in enumerate(waveNumbers):
        fieldRecord[i] = modeTemplate.copy()
        fieldRecord[i]['k'] = k
        fieldRecord[i]['initialN'] = N(k)
        fieldRecord[i]['frequency'] = f(x, y, k, m, Xi, H, K)
        fieldRecord[i]['minkowskiFrequency'] = MinkowskiFrequency(xRange, k, m, H, K)
        fieldRecord[i]['phaseIntegralW'] = PhaseIntegralW(xRange, np.sqrt(np.real(f(xRange, 0., k, m, Xi, H, K))), xRange, Order=2, N=-1)
        fieldRecord[i]['phaseIntegralV'] = PhaseIntegralV(xRange, np.sqrt(np.real(f(xRange, 0., k, m, Xi, H, K))), xRange, Order=2, N=-1)
    
    # Stores the `x` and `y` coordinates over which the modes are initialized.
    fieldRecord['x'] = xRange
    fieldRecord['y'] = yRange

    # Stores the field mass `m`, coupling `Xi`, and `planckRatio` as dictionary entries.
    fieldRecord['m'] = m
    fieldRecord['Xi'] = Xi
    fieldRecord['planckRatio'] = planckRatio

    # Returns an initialized dictionary containing a dictionary record for each field mode.
    return fieldRecord


def InitMetric(s, h, H, K, Lambda, Box):
    """ Creates a dictionary containing the initial conditions for the spacetime metric within `Box`.
        
    Args:
        s: a callable function `s(x)` which outputs a the evolution of the scale factor within `Box`.
        s1: a callable function `s1(x)` which outputs a the evolution of the first derivative of the scale factor within `Box`.
        K: an integer equal to `-1`, `0`, or `1` indicating the background curvature.
        Lambda: a float characterizing the cosmological constant.
        Box: a tuple of the form `(xRange, yRange)`, where each entry is a (m,)-dimensional NumPy array spamming one of the box dimensions.
        
    Returns:
        A dictionary with entries 'x', 's', and 's1' containing the initial conditions for the `x` coordinate, the scale factor, and the first derivative of the scale factor, respectively.
        
    """

    # Extracts the `xRange` array from `Box`.
    xRange = Box[0]

    # Constructs the dictionary template.
    metricRecord = {'x':0, 's':0, 's1':0, 'K':0, 'Lambda':0}

    # Initializes each dictionary entry.
    metricRecord['x'] = xRange
    metricRecord['s'] = s(xRange, H, K)
    metricRecord['s1'] = s(xRange, H, K)*h(xRange, H, K)
    metricRecord['K'] = K
    metricRecord['Lambda'] = Lambda

    # Returns the initialized dictionary containing
    return metricRecord
