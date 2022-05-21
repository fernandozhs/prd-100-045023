# SciPy B-spline Interpolation:
from scipy.interpolate import interp1d
from scipy.interpolate import splev
from scipy.interpolate import splrep
from scipy.interpolate import RectBivariateSpline

# SciPy Integrator:
from scipy.integrate import quad

# SciPy Pade Approximants:
from scipy.misc import pade

# SciPy Root Finder:
from scipy.optimize import fsolve

# SciPy Special Functions:
from scipy.special import factorial

# PyGSL
from pygsl import bspline, multifit

# NumPy
import numpy as np


def Downsample(Array, N):
    """ Uniformly downsamples a (n,)-dimensional NumPy array.

    Args:
        Array: the input (n,)-dimensional NumPy array to be downsampled.
        N: the integer number of elements in the resulting downsampled array.

    Returns:
        If `N < len(Array)`, returns a (N,)-dimensional NumPy array.
        The entries correspond to `N` equispaced elements of `Array`.

        If `N >= len(Array)`, no downsampling is performed and `Array` is returned.

    Examples:
        >>> A = np.array([0., 1., 2., 3., 4., 5., 6.])
        >>> Downsample(A, 3)
        >>> array([0., 3., 6.])

        >>> Downsample(A, 10)
        >>> array([0., 1., 2., 3., 4., 5., 6.])
    """

    # Determines the necessary index spacing for the uniform sampling of `N` elements of `Array` which includes both its first and last elements.
    indexSampling = np.linspace(0, len(Array) - 1, N, dtype=int)

    # Notice that if `N >= len(Array)` the instruction above will return indices of `Array` with repetitions. In that case the repeated indices are excluded.
    indexSampling = np.unique(indexSampling)

    # Returns the a uniformly downsampled version of `Array`.
    return Array[indexSampling]


def Crop(Array, Mask):
    """ Crops an n-dimensional NumPy array along its 0-th axis.

    Args:
        Array: the input (n, m)-dimensional NumPy array to be cropped.
        Mask: a (n,)-dimensional NumPy Boolean array containing `False` or `True` as elements.

    Returns:
        A copy of `Array` with all elements with 0-th (`n`) axis entry labeled as `False` in `Mask` removed.

    Examples:
        >>> A = np.array([0., 1., 2., 3.])
        >>> M = np.array([False, False, False, True])
        >>> Crop(A, M)
        >>> array([3.])

        >>> B = np.array([[0., 1.], [2., 3.], [4., 5.], [6., 7.]])
        >>> Crop(B, A > 1.)
        >>> array([[4., 5.], [6., 7.]])
    """

    # Returns a copy of `Array` containing only those elements with 0-th (`n`) axis entries labeled as `True` in `Mask`.
    return Array[Mask]


def Distance(Array, Point):
    """ Computes the distance between `Point` and each entry in `Array`.

    Args:
        Array: the input (n, 2)-dimensional NumPy array storing a collection of `n` points `(x, y)`.
        Point: a (2,)-dimensional NumPy array of the form `np.array([x0, y0])` labeling a reference point.

    Returns:
        An (n,)-dimensional NumPy array with entries corresponding to the distances from each element of `Array` to `Point`.

    Example:
        >>> A = np.array([[0., 0.], [1., 2.]])
        >>> P = np.array([5., 5.])
        >>> Distance(A, P)
        >>> array([7.07106781, 5.])
    """

    # Computes the component-wise difference squared between `Array` and `Point`, and performs their sum to obtain the distances squared.
    distanceSquared = np.sum((Array - Point)**2, axis=1)

    # Returns the (n,)-dimensional Array with entries corresponding to the distances from each element of `Array` to `Point`.
    return np.sqrt(distanceSquared)


def BSplineParameters(Order):
    """ Generates the adequate parameters for the B-spline differentiation of a function sample.

    Args:
        Order: an integer which determines the differentiation order.

    Returns:
        A tuple containing the integer B-spline order most appropriate for the intended order of differentiation and the integer number of points which should be used in the B-spline interpolation.

    Example:
        >>> BSplineParameters(1.)
        >>> (3, 11)

    WARNING:
        The parameters returned by this function have been tunned for the problem at hand and might require revision in order to be used in a different context.
    """

    # Identifies the B-spline order which suits the intended order of differentiation.
    if (Order < 2):
        BSplineOrder = 3

    elif (Order >= 2 and Order < 6):
        BSplineOrder = 2*Order + 1

    elif (Order >= 6 and Order < 8):
        BSplineOrder = 2*Order - 1

    elif (Order >= 8):
        BSplineOrder = 14

    # Identifies the number of points `N` which should be used in the B-spline interpolation.
    if BSplineOrder < 10:
        N = 11
    else:
        N = BSplineOrder + 1

    # Returns a tuple containing the B-spline order and the number of points.
    return (BSplineOrder, N)


def BSpline(fSample, xSample, Order, N):
    """ Generates a B-spline interpolation from a function `f = f(x)` sample.

    Args:
        fSample: a (n,)-dimensional NumPy array sampling the function `f` to be interpolated.
        xSample: a (n,)-dimensional NumPy array sampling the independent variable `x`.
        Order: an integer which determines the order of the B-spline interpolation.
        N: the integer number of points used in the interpolation.

    Clearly the first two arguments are related through `f` by `fSample = f(xSample)`.

    Returns:
        A tuple containing the B-spline internal knots, coefficients, and order.
        This tuple can be used as the input of SciPy's `splev` and `splrep` routines.
    """

    # The B-spline is constructed using PyGSL routines. The spline object is initialized with uniformly distributed knots.
    BSpline = bspline.bspline(Order + 1, N - Order + 1)
    BSpline.knots_uniform(xSample[0], xSample[-1])

    # The optimal spline coefficients are obtained through a PyGSL least-squares routine.
    Y = np.zeros((N, N))
    Y = BSpline.eval_vector(xSample)
    Coefficients, Covariances, Chi_Squared = multifit.linear(Y, fSample)

    # Returns a tuple containing the B-spline internal knots, coefficients, and order.
    return (BSpline.get_internal_knots(), Coefficients, Order)


def Derivative(x, fSample, xSample, Order, N=-1):
    """ Computes the derivative of a function `f = f(x)` from its numerical sample.

    Args:
        x: a float or (m,)-dimensional NumPy array specifying the point(s) at which the derivative should be computed.
        fSample: a (n,)-dimensional NumPy array sampling the function `f` to be differentiated.
        xSample: a (n,)-dimensional NumPy array sampling the independent variable `x`.
        Order: an integer which determines the differentiation order.
        N: the integer number of points used in the computation.

    Clearly the second and third arguments are related through `f` by `fSample = f(xSample)`.

    Returns:
        A (m,)-dimensional NumPy array sampling the derivative of order `Order` of `f` at the point(s) specified by `x`.

    Example:
        >>> xSample = np.linspace(-1., 1., 300)
        >>> fSample = np.exp(xSample)
        >>> Derivative(0., fSample, xSample, 1, N=len(xSample))
        >>> array(1.000000000002948)
    """

    # If the number of points to be used in the computation of the derivative has not been specified (i.e., `N == -1`), that information is obtained through the function `BSplineParameters`, which also returns the adequate B-spline order for interpolation. Otherwise, only the B-spline order is obtained through `BSplineParameters`.
    if N == -1:
        BSplineOrder, N = BSplineParameters(Order)
    else:
        BSplineOrder, __ = BSplineParameters(Order)

    # `fSample` and `xSample` are downsampled to `N` points each prior to the B-spline interpolation.
    fSample = Downsample(fSample, N)
    xSample = Downsample(xSample, N)

    # The object containing the B-Spline construction is built using PyGSL routines.
    BSplineObj = BSpline(fSample, xSample, BSplineOrder, N)

    # Returns the value of the B-spline derivative of order `Order` at the point `x`:
    return splev(x, BSplineObj, der=Order)


def TaylorCoefficient(x, fSample, xSample, Order, N=-1):
    """ Computes the coefficients of a function's `f = f(x)` Taylor series representation from the function's numerical sample.

    Args:
        x: a float specifying the point around which the Taylor series is being computed.
        fSample: a (n,)-dimensional NumPy array sampling the function `f` being expanded.
        xSample: a (n,)-dimensional NumPy array sampling the independent variable `x`.
        Order: an integer which determines the order of the Taylor series coefficient.
        N: the integer number of points used in the computation.

    Clearly the second and third arguments are related through `f` by `fSample = f(xSample)`.

    Returns:
        A float which corresponds to the coefficient of order `Order` of `f`'s Taylor expasion around `x`.

    Example:
        >>> xSample = np.linspace(-1., 1., 300)
        >>> fSample = np.exp(xSample)
        >>> TaylorCoefficient(0., fSample, xSample, N=len(xSample))
        >>> 1.0
    """

    # Computes `f`'s derivative of order `Order` at the point `x` from `fSample` and `xSample`.
    fDerivative = Derivative(x, fSample, xSample, Order, N)

    # Returns the Taylor coefficient associated with term of order 'Order' in the series.
    return fDerivative/factorial(Order)


def PadeApproximant(x, y, fSample, xSample, Order, N=-1):
    """ Computes the Pade approximant of a function's `f = f(x)` truncated Taylor series at the point `z = x + 1j*y` of the complex plane.

    Args:
        x: a float specifying the Real coordinate of the point where the Pade approximant is to be evaluated.
        y: a float of (m,)-dimensional NumPy array specifying the Imaginary coordinate of the point(s) where the Pade approximant is to be evaluated.
        fSample: a (n,)-dimensional NumPy array sampling the function `f` along the Real axis.
        xSample: a (n,)-dimensional NumPy array sampling the independent variable `x` over the Real axis.
        Order: an integer which determines the order of the Pade approximant.
        N: the integer number of points used in the computation.

    Clearly the third and fourth arguments are related through `f` by `fSample = f(xSample)`.

    Returns:
        An (m,)-dimensional NumPy array containing the Pade approximant of `f` at the point(s) specified by `x + 1j*y`.

    Example:
        >>> xSample = np.linspace(-1., 1., 300)
        >>> fSample = np.exp(xSample)
        >>> PadeApproximant(0., np.array([0., 1.]), fSample, xSample, Order=6, N=len(xSample))
        >>> array([ 1.00000000+0.j, 0.54030231+0.84147097j])
    """

    # Computes `f`'s Taylor coefficients at `z = x + 1j*0` up to order `10`:
    fTaylor = np.array([TaylorCoefficient(x, fSample, xSample, 0, N),
                        TaylorCoefficient(x, fSample, xSample, 1, N),
                        TaylorCoefficient(x, fSample, xSample, 2, N),
                        TaylorCoefficient(x, fSample, xSample, 3, N),
                        TaylorCoefficient(x, fSample, xSample, 4, N),
                        TaylorCoefficient(x, fSample, xSample, 5, N),
                        TaylorCoefficient(x, fSample, xSample, 6, N),
                        TaylorCoefficient(x, fSample, xSample, 7, N),
                        TaylorCoefficient(x, fSample, xSample, 8, N),
                        TaylorCoefficient(x, fSample, xSample, 9, N),
                        TaylorCoefficient(x, fSample, xSample, 10, N)])

    # Computes the Pade approximant to `f`'s truncated Taylor series:
    P, Q = pade(fTaylor, Order)

    # Returns the Pade approximant evaluated at the point `z = x + 1j*y` of the complex plane.
    return P(0. + 1j*y)/Q(0. + 1j*y)


def AnalyticContinuation(fSample, Box, Order, N=-1):
    """ Analitically continues a function `f = f(x)` over the complex plane from its numerical sample along an interval of the Real axis.

    Args:
        fSample: a (n,)-dimensional NumPy array sampling the function `f` along the Real axis.
        Box: a tuple of the form `(xRange, yRange)`, where the entries are an (n,)-dimensional and (m,)-dimensional NumPy arrays spamming the box dimensions.
        Order: an integer which determines the order of the Pade approximant used in the analytical continuation.
        N: the integer number of points over the Real axis used in the computation.

    Clearly the first, second and third arguments are related through `f` by `fSample = f(xSample + 1j*ySample)`.

    Returns:
        An (n, m)-dimensional array containing the analytical continuation of `f` over the mesh grid defined by `Box`.

    Example:
        >>> xSample = np.linspace(-1., 1., 300)
        >>> fSample = np.exp(xSample)
        >>> PadeApproximant(0., np.array([0., 1.]), fSample, xSample, Order=6, N=len(xSample))
        >>> array([ 1.00000000+0.j, 0.54030231+0.84147097j])
    """

    # Extracts the `xSample` and `ySample` arrays from `Box`.
    xSample = Box[0]
    ySample = Box[1]

    # Creates a mesh grid over the complex plane delimited by `xSample` and `ySample`.
    y, x = np.meshgrid(ySample, xSample, indexing='xy')

    # Initializes the NumPy array which will store the analytic continuation of `f`.
    fContinuation = np.zeros(y.shape, dtype='complex')

    # The function `f` is analitically extended to `z = x + 1j*y` through its truncated Taylor representation around `z = x + 1j*0` and the use of Pade approximants, which accelerate the convergence of the truncated series.
    for i in range(0, len(x)):
        fContinuation[i] = PadeApproximant(x[i][0], y[i], fSample, xSample, Order, N)

    # Returns the function's analytic continuation `fContinuation` over the region defined by `Box` on the complex plane.
    return fContinuation


def FunctionInterpolation(fSample, Box):
    """ Use bivariate B-splines to interpolate a function `f` from its sample values over a region of the complex plane.

    Args:
        fSample: a (n, m)-dimensional NumPy array sampling the function `f` over the region defined by `Box` in the complex plane.
        Box: a tuple of the form `(xRange, yRange)`, where the entries are an (n,)-dimensional and (m,)-dimensional NumPy arrays spamming the box dimensions.

    Clearly the first and second arguments are related through `fSample = f(Box)`.

    Returns:
        Two callable functions. The first `f(x, y)` interpolates the input sample of `f`, while the second `fSqrt(x, y)` interpolates the square root of the sampled `f`.
    """

    # Extracts the `xRange` and `yRange` arrays from `Box`.
    xRange = Box[0]
    yRange = Box[1]

    # Separates the Real and Imaginary parts stored in `fSample`.
    fReSample = np.real(fSample)
    fImSample = np.imag(fSample)

    # Interpolates each function sample using a bivariate B-spline.
    fRe = RectBivariateSpline(xRange, yRange, fReSample, s=0)
    fIm = RectBivariateSpline(xRange, yRange, fImSample, s=0)

    # Puts Real and Imaginary parts together to return `f(x,y) = fRe(x,y) + 1j*fIm(x,y)` and its square root.
    f = lambda x, y: fRe(x, y)[0][0] + 1j*fIm(x, y)[0][0]
    fSqrt = lambda x, y: np.sqrt(fRe(x, y)[0][0] + 1j*fIm(x, y)[0][0])

    # Returns the interpolated functions.
    return f, fSqrt


def StokesLine(x, y, f):
    """ Defines the differential equation satisfied by the Stokes lines of a complex function `f`.

    Args:
        x: a float specifying the Real axis coordinate of a point `z = x + 1j*y` on the complex plane.
        y: a float specifying the Imaginary axis coordinate of a point `z = x + 1j*y` on the complex plane. This argument also parametrizes the Stokes line through `x = x(y)`.
        f: a callable function with two arguments `f(x, y)` defined over the complex plane.

    Returns:
        A float which characterizes `x1(z)`, i.e., the first derivative of the Stokes line coordinates `x(y)` with respect to the coordinate `y` computed at the point `z`.
    """

    # Defines the differential equations satisfied by the coordinates `x(y)` describing a Stokes line of `f`:
    x1 = np.imag(f(x, y))/np.real(f(x, y))

    # Returns the left-hand side of the differential equation:
    return x1


def RootFinder(f, Box):
    """ Finds the root of the complex function `f` which is closest to the Real axis and within a specified box.

    Args:
        f: a callable function with two arguments `f(x, y)` defined over the complex plane.
        Box: a tuple of the form `(xRange, yRange)`, where the entries are an (n,)-dimensional and (m,)-dimensional NumPy arrays spamming the box dimensions.

    Returns:
        If a root has been found, the function returns `Root`, a (2,)-dimensional NumPy array of the form `np.array([x, y])`.
        If the search failed, the function returns `False`.
    """

    # Recasts `f(x, y)` to a form which is a better suited input to the root-finding algorithm.
    fSquared = lambda z: np.array([np.real(f(z[0], z[1])), np.imag(f(z[0], z[1]))])

    # The starting points for root-searching are distributed along a line that vertically bisects the input `Box`:
    x = np.average(Box[0])
    yRange = np.linspace(0., Box[1].max(), 10)

    # All roots are temporarily stored in `tempRoots`.
    tempRoots = np.zeros((10, 2))

    # Search for roots using as starting guesses all the points with coordinates `(x, yRange)`.
    for i, y in enumerate(yRange):
        candidateRoot = fsolve(fSquared, x0=np.array([x, y]), full_output=True)

        # If a given search converges, the root is stored in `tempRoots`.
        if (candidateRoot[3] == 'The solution converged.'):
            tempRoots[i,:] = candidateRoot[0]

        # If the search fails, the loop proceeds to the next starting guess.
        else:
            continue

    # Pre-selects roots found within `Box` with non-zero Imaginary part:
    withinBox = ( (tempRoots[:,1] <> 0.)
                 *(tempRoots[:,0] > Box[0].min())
                 *(tempRoots[:,0] < Box[0].max())
                 *(tempRoots[:,1] > Box[1].min())
                 *(tempRoots[:,1] < Box[1].max()) )

    # If there are no roots remaining once the above filter is applied, the function ends returning `False`.
    if np.all(withinBox <> True):
        return False

    # Otherwise the root which is closest to the Real axis is selected, and `Root` is returned.
    else:
        closestRe = np.abs(tempRoots[withinBox, 1] - 0.).argmin()
        Root = tempRoots[withinBox,:][closestRe]

        # Here we ensure the Imaginary part of `Root` is positive. (This is a matter of convention, and is adopted here since for particle production all roots in the complex plane come in symmetric pairs with respect to the Real axis.)
        Root = np.array([Root[0], np.abs(Root[1])])

        # Returns remaining root.
        return Root


def Epsilon0(x, QSample, xSample, N=-1):
    """ Computes the phase-integral function `e0(x)` from a sample of the phase-integral base function `Q(x)`. The base function is chosen to be the 2nd order O.D.E. frequency function `f(x)`.

    Args:
        x: a float or (n,)-dimensional NumPy array specifying the point(s) at which the function should be computed.
        QSample: a (n,)-dimensional NumPy array sampling the function `Q` over the Real axis.
        xSample: a (n,)-dimensional NumPy array sampling the independent variable `x` over the Real axis.
        N: the integer number of points over the Real axis used in the computation.

    Clearly the second and third arguments are related through `QSample = Q(xSample)`.

    Returns:
        The phase-integral function `e0(x)` computed at `x`.
    """

    # Generates a sample of the function `g = 1/Q` from `QSample`.
    gSample = 1./QSample

    # Computes the second derivative of `g(x)` over the range defined by `xSample`.
    gDerivative = Derivative(x, gSample**0.5, xSample, 2, N)

    # Calculates `e0(x)`.
    e0Sample = gDerivative * gSample**1.5

    # Returns the function `e0(x)` sampled at `x`.
    return e0Sample


def Epsilon1(x, QSample, xSample, N=-1):
    """ Computes the phase-integral function `e1(x)` from a sample of the phase-integral base function `Q(x)`. The base function is chosen to be the 2nd order O.D.E. frequency function `f(x)`.

    Args:
        x: a float or (n,)-dimensional NumPy array specifying the point(s) at which the function should be computed.
        QSample: a (n,)-dimensional NumPy array sampling the function `Q` over the Real axis.
        xSample: a (n,)-dimensional NumPy array sampling the independent variable `x` over the Real axis.
        N: the integer number of points over the Real axis used in the computation.

    Clearly the second and third arguments are related through `QSample = Q(xSample)`.

    Returns:
        The phase-integral function `e1(x)` computed at `x`.
    """

    # Generates a sample of the function `g = 1/Q` from `QSample`.
    gSample = 1./QSample

    # Computes the many different factors involving derivatives of `g(x)`.
    factor1 = Derivative(x, gSample**1.5, xSample, 1, N)*Derivative(x, gSample**0.5, xSample, 2, N)*gSample
    factor2 = Derivative(x, gSample**0.5, xSample, 3, N)*gSample**2.5

    # Calculates `e1(x)`.
    e1Sample = factor1 + factor2

    # Returns the function `e0(x)` sampled at `x`.
    return e1Sample


def Epsilon2(x, QSample, xSample, N=-1):
    """ Computes the phase-integral function `e2(x)` from a sample of the phase-integral base function `Q(x)`. The base function is chosen to be the 2nd order O.D.E. frequency function `f(x)`.

    Args:
        x: a float or (n,)-dimensional NumPy array specifying the point(s) at which the function should be computed.
        QSample: a (n,)-dimensional NumPy array sampling the function `Q` over the Real axis.
        xSample: a (n,)-dimensional NumPy array sampling the independent variable `x` over the Real axis.
        N: the integer number of points over the Real axis used in the computation.

    Clearly the second and third arguments are related through `QSample = Q(xSample)`.

    Returns:
        The phase-integral function `e2(x)` computed at `x`.
    """

    # Generates a sample of the function `g = 1/Q` from `QSample`.
    gSample = 1./QSample

    # Computes the many different factors involving derivatives of `g(x)`.
    factor1 = Derivative(x, gSample, xSample, 1, N)*Derivative(x, gSample**1.5, xSample, 1, N)*Derivative(x, gSample**0.5, xSample, 2, N)*gSample
    factor2 = Derivative(x, gSample, xSample, 1, N)*Derivative(x, gSample**0.5, xSample, 3, N)*gSample**2.5
    factor3 = Derivative(x, gSample**1.5, xSample, 2, N)*Derivative(x, gSample**0.5, xSample, 2, N)*gSample**2
    factor4 = 2.*Derivative(x, gSample**1.5, xSample, 1, N)*Derivative(x, gSample**0.5, xSample, 3, N)*gSample**2
    factor5 = Derivative(x, gSample**0.5, xSample, 4, N)*gSample**3.5

    # Calculates `e2(x)`.
    e2Sample = factor1 + factor2 + factor3 + factor4 + factor5

    # Returns the function `e0(x)` sampled at `x`.
    return e2Sample


def PhaseIntegralW(x, QSample, xSample, Order=2, N=-1):
    """ Computes the phase-integral corrections `W` up to the input order `Order`. The base phase-integral base function `Q(x)` is chosen to be the 2nd order O.D.E. frequency function `f(x)`.

    Args:
        x: a (n,)-dimensional NumPy array specifying the points at which the function should be computed.
        QSample: a (n,)-dimensional NumPy array sampling the function `Q` over the Real axis.
        xSample: a (n,)-dimensional NumPy array sampling the independent variable `x` over the Real axis.
        Order: an integer equal to `0`, `2`, or `4` which determines up to which order the phase-integral correction will be computed.
        N: the integer number of points over the Real axis used in the computation.

    Clearly the second and third arguments are related through `QSample = Q(xSample)`.

    Returns:
        The phase-integral correction `W(x)` computed at all points of `x` up to order `Order`.
    """

    # Selects the appropriate expression according to the input top `Order`.
    if Order == 0:
        # Returns a sample of 0th order function `W0(x)`.
        return QSample

    elif Order == 2:
        # Computes a sample of `e0(x)`.
        e0Sample = Epsilon0(x, QSample, xSample, N)

        # Returns a sample of 2nd order function `W2(x)`.
        return QSample*(1. + 0.5*e0Sample)

    elif Order == 4:
        # Computes a sample of `e0(x)`.
        e0Sample = Epsilon0(x, QSample, xSample, N)

        # Computes a sample of `e2(x)`.
        e2Sample = Epsilon2(x, QSample, xSample, N)

        # Returns a sample of 2nd order function `W4(x)`.
        return QSample*(1. + 0.5*e0Sample - 0.125*(e0Sample**2 + e2Sample))

    # If `Order` does not match any of the above possibilities, an error message is raised.
    else:
        raise ValueError("Order must be equal to `1`, `2`, or `4`. Higher orders have not been implemented.")


def PhaseIntegralV(x, QSample, xSample, Order=2, N=-1):
    """ Computes the phase-integral corrections `V` up to the input order `Order`. The base phase-integral base function `Q(x)` is chosen to be the 2nd order O.D.E. frequency function `f(x)`.

    Args:
        x: a (n,)-dimensional NumPy array specifying the points at which the function should be computed.
        QSample: a (n,)-dimensional NumPy array sampling the function `Q` over the Real axis.
        xSample: a (n,)-dimensional NumPy array sampling the independent variable `x` over the Real axis.
        Order: an integer equal to `0` or `2` which determines up to which order the phase-integral correction will be computed.
        N: the integer number of points over the Real axis used in the computation.

    Clearly the second and third arguments are related through `QSample = Q(xSample)`.

    Returns:
        The phase-integral correction `V(x)` computed at all points of `x` up to order `Order`.
    """

    # Selects the appropriate expression according to the input top `Order`.
    if Order == 0:
        # Computes the first order derivative of `W0(x)`.
        W01Sample = Derivative(x, QSample, xSample, 1, N)

        # Returns a sample of 0th order function `V0(x)`.
        return -1.*(W01Sample/QSample)

    elif Order == 2:
        # Computes the first order derivative of `W0(x)`.
        W01Sample = Derivative(x, QSample, xSample, 1, N)

        # Computes a sample of `e0(x)`.
        e0Sample = Epsilon0(x, QSample, xSample, N)

        # Computes a sample of `e1(x)`.
        e1Sample = Epsilon1(x, QSample, xSample, N)

        # Returns a sample of 2nd order function `V2(x)`.
        return -1.*(W01Sample/QSample + 0.5*e1Sample/(1. + 0.5*e0Sample))

    # If `Order` does not match any of the above possibilities, an error message is raised.
    else:
        raise ValueError("Order must be equal to `1` or `2`. Higher orders have not been implemented.")
