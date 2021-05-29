# Copy
import copy

# Pickle
import cPickle

# SciPy B-spline Interpolation:
from scipy.interpolate import splrep

# NumPy
import numpy as np

# Tools, Phys, and Init
from tools import *
from phys import *
from init import *


def CropBox(Box, xRatio=20.):
    """ Equally crops the left and right borders of a rectanguar box leaving a central band bisected by the Real axis. The width of the cropped box is smaller by a factor of `xRatio`.
        
    Args:
        Box: a tuple of the form `(xRange, yRange)`, where the entries are an (n,)-dimensional and (m,)-dimensional NumPy arrays spamming the box dimensions.
        xRatio: an integer specifying the ratio between the width of the input box and the resulting cropped box.
        
    Returns:
        A tuple of the form `(croppedRange, yRange)` representing a box with left and right borders cropped by equal amounts.
        
    Examples:
        >>> Box = (np.array([1., 2., 3., 4., 5., 6., 7.]), np.array([-2., -1, 0., 1., 2.]))
        >>> SqueezeBox(Box, xFactor=3)
        >>> (np.array([3., 4., 5.]), np.array([-2., -1, 0., 1., 2.]))
        
    """
    
    # Retrieves the arrays spamming the `x` and `y` directions from the input `Box`.
    xRange = Box[0]
    yRange = Box[1]
    
    # Finds the central index in the input box's `x` direction.
    indexCenter = int(len(xRange)/2.)

    # Computes the integer number of indices that needed to be added or subtracted from `indexCenter` to create a central region corresponding to a box which is smaller by a factor of `xRatio` in the `x` direction.
    indexLength = int(len(xRange)/(2.*xRatio))
    
    # Creates the indices bounding the squeezed box.
    indexLeft = indexCenter - indexLength
    indexRight = indexCenter + indexLength + 1
    
    # Returns a box with left and right borders cropped by equal amounts.
    return (xRange[indexLeft:indexRight], yRange)


def ShiftBox(Box, xShift):
    """ Shifts a rectanguar box along the Real axis by `xShift`.
        
    Args:
        Box: a tuple of the form `(xRange, yRange)`, where the entries are an (n,)-dimensional and (m,)-dimensional NumPy arrays spamming the box dimensions.
        xShift: a float specifying the amount by which the Box is to be shifted along the `x` direction.
        
    Returns:
        A tuple of the form `(xRange + xShift, yRange)`.
        
    Examples:
        >>> Box = (np.array([1., 2.]), np.array([-2., -1, 0., 1., 2.]))
        >>> ShiftBox(Box, xShift=0.1)
        >>> (np.array([1.1, 2.1]), np.array([-2., -1, 0., 1., 2.]))
        
    """
    
    # Retrieves the arrays spamming the `x` and `y` directions from the input `Box`.
    xRange = Box[0]
    yRange = Box[1]
    
    # Returns the updated box shifted by `xShift` along the `x` direction.
    return (xRange + xShift, yRange)


def UpdateMetric(metricRecord, Box, s, h, H, K):
    """ Updates the `metricRecord` to incorporate the most recent values of its entries 'x', 's', and 's1'.
        
    Args:
        metricRecord: a metric dictionary containing NumPy arrays which detail the evolution of the spacetime metric.
        Box: a tuple of the form `(xRange, yRange)`, where the entries are an (n,)-dimensional and (m,)-dimensional NumPy arrays spamming the box dimensions.
        s: a callable function `s(x)` which outputs a the value of the scale factor at `x`.
        h: a callable function `h(x)` which outputs a the value of the Hubble parameter at `x`.
        
    Returns:
        An updated version of `metricRecord` which includes the most recent values of its entries 'x', 's', and 's1'.
    
    WARNING:
        This function needs to reformulated in order to update `metricRecord` according to the Friedman equations. At the moment the update merely imposes a pre-determined evolution which does not include backreaction effects.
        
    """

    # Extracts the latest `x` value.
    x = Box[0].max()
    
    # Updates the entries 'x', 's', and 's1' of `metricRecord`.
    metricRecord['x'] = np.append(metricRecord['x'], x)
    metricRecord['s'] = np.append(metricRecord['s'], s(x, H, K))
    metricRecord['s1'] = np.append(metricRecord['s1'], s(x, H, K)*h(x, H, K))

    # Returns the updated `metricRecord`
    return metricRecord


def UpdateMetricNew(metricRecord, fieldRecord, Box, rho):
    """ Updates the `metricRecord` to incorporate the most recent values of its entries 'x', 's', and 's1'.
        
    Args:
        metricRecord: a metric dictionary containing NumPy arrays which detail the evolution of the spacetime metric.
        fieldRecord: a dictionary containing a dictionary entry for each field mode.
        Box: a tuple of the form `(xRange, yRange)`, where the entries are an (n,)-dimensional and (m,)-dimensional NumPy arrays spamming the box dimensions.
        rho: a callable function `rho(x, s)` which outputs a the field energy density as a function of the time variable `x` and the scale factor `s`.
        
    Returns:
        An updated version of `metricRecord` which includes the most recent values of its entries 'x', 's', and 's1'.
        
    """
    
    # Extracts the value of the scale factor `s(x)` added in the previous update of `metricRecord`.
    s = metricRecord['s'][-1]
    
    # Extracts the `x` value to be added.
    xRange = Box[0][-2:]
    
    # Extracts the curvature `K` and the cosmological constant `Lambda` from `metricRecord`.
    K = metricRecord['K']
    Lambda = metricRecord['Lambda']
    
    # Extracts the ratio between the field mass and the Planck mass from `fieldRecord`.
    planckRatio = fieldRecord['planckRatio']
    
    # Encapsulates the field enegy density into the callable function `rho(x, s)`.
    #rho = lambda x, s: 0 #To turn on the back-reaction, set the function to return: EnergyDensity(x, s, fieldRecord)
    
    # Solves for the scale factor `s(x)`.
    sNew = odeint(Friedmann, s, xRange, args=(rho, Lambda, K, planckRatio))
    
    # Determine the derivative of the scale factor `s1(x)`.
    s1New = Friedmann(sNew[-1], xRange[-1], rho, Lambda, K, planckRatio)
    
    # Updates the entries 'x', 's', and 's1' of `metricRecord`.
    metricRecord['x'] = np.append(metricRecord['x'], xRange[-1])
    metricRecord['s'] = np.append(metricRecord['s'], sNew[-1])
    metricRecord['s1'] = np.append(metricRecord['s1'], s1New)
    
    # Returns the updated `metricRecord`
    return metricRecord


def RootCheck(modeRecord, Root, absoluteError=0.1):
    """ Checks whether `Root` has already been recorded to `modeRecord['root']` within a given absolute error.
        
    Args:
        Root: a (2,)-dimensional NumPy array of the form `np.array([x, y])` specifying the root coordinates on the complex plane.
        modeRecord: a field mode dictionary containing a Python list of roots under the 'root' entry.
        absoluteError: a float specifying the amount by which `Root` must deviate from the roots stored in `modeRecord` in order to be considered a new root.
        
    Returns:
        If `Root` has already been recorded to `modeRecord` up to variations of the order of `absoluteError`, `False` is returned.
        
        If `Root` has not been previously recorded to `modeRecord` up to variations of the order of `absoluteError`, `True` is returned.
        
    """
    
    # Converts the list of roots in `modeRecord` to an array.
    rootRecord = np.array(modeRecord['root'])
    
    # If the `rootRecord` array is empty, `Root` is the first root to be found and therefore a new entry. In this case `True` is returned.
    if len(rootRecord) == 0:
        return True
    
    # When the above condition is not met, the distance between `Root` and all roots stored in `modeRecord` is computed.
    distances = Distance(rootRecord, Root)
    
    # If any of the distances computed above is smaller than `absoluteError`, `Root` is considered a repeated entry and `False` is returned.
    if np.any(distances <= absoluteError):
        return False
    
    # If none of the distances computed above is smaller than `absoluteError`, `Root` is considered a new entry and `True` is returned.
    else:
        return True


def AmplitudeCheck(modeRecord, Root, f, thresholdAmplitude=1e-5):
    """ Checks whether the amplitude of the particle production event sourced by `Root` is negligible compared to previously recorded events stored in `modeRecord`.
        
    Args:
        Root: a (2,)-dimensional NumPy array of the form `np.array([x, y])` specifying the root coordinates on the complex plane.
        modeRecord: a field mode dictionary containing a Python list of roots under the 'root' entry.
        f: a callable function with two arguments `f(x, y)` represeting the frequency function over the complex plane for the field mode of interest.
        
    Returns:
        If the ratio between the amplitude of the event sourced by `Root` and the amplitude of all previous events recorded in `modeRecord` is smaller or equal to `thresholdAmplitude`, `False` is returned.
        
        If the ratio between the amplitude of the event sourced by `Root` and the amplitude of any previous events recorded in `modeRecord` is larger than `thresholdAmplitude`, `True` is returned.
        
    """
    
    # If the `modeRecord['eventAmplitude']` list is empty, the event is the first to be found and therefore no comparison is needed. In this case `True` is returned.
    if len(modeRecord['eventAmplitude']) == 0:
        return True
    
    # Calculates the amplitude of the particle production event sourced by `Root`.
    eventAmplitude = EventAmplitude(Root, f)
    
    # Computes the ratio between `eventAmplitude` and all entries of `modeRecord['eventAmplitude']`.
    amplitudeRatio = eventAmplitude/np.array([modeRecord['eventAmplitude']])
    
    # If all ratios stored in `amplitudeRatio` are smaller or equal to `thresholdAmplitude`, `False` is returned.
    if np.all(amplitudeRatio <= thresholdAmplitude):
        return False
    
    # Otherwise, if any of the ratios stored in `amplitudeRatio` is larger than `thresholdAmplitude`, `True` is returned.
    else:
        return True


def UpdateFrequencies(fieldRecord, modeRecord, metricRecord, Box):
    """ Updates `modeRecord['frequency']`, `modeRecord['minkowskiFrequency']`, `modeRecord['phaseIntegralW']`, and `modeRecord['phaseIntegralV']` using the input information obtained from `metricRecord` within the region delimited by `Box`.
        
    Args:
        fieldRecord: a dictionary containing a dictionary entry for each field mode.
        modeRecord: the dictionary entry for a particular field mode in `fieldRecord`.
        metricRecord: a metric dictionary containing NumPy arrays which detail the evolution of the spacetime metric.
        Box: a tuple of the form `(xRange, yRange)`, where the entries are an (n,)-dimensional and (m,)-dimensional NumPy arrays spamming the box dimensions.
        
    Returns:
        An updated version of `modeRecord` which includes the phase-integral functions `mf`, `W` and `V`, and the analytic continuation of the mode frequency function squared sampled over an area of the complex plane which had not been previously covered.

    """

    # Extracts `xRange` from `Box`.
    xRange = Box[0]
    
    # Creates a mask to select the metric records which lie within `Box`.
    withinBox = (metricRecord['x'] >= xRange.min())*(metricRecord['x'] <= xRange.max())

    # Restricts the `metricRecord` arrays which track the evolution of the metric to those values which lie within `Box`.
    xSample = metricRecord['x'][withinBox]
    sSample = metricRecord['s'][withinBox]
    s1Sample = metricRecord['s1'][withinBox]
    K = metricRecord['K']

    # Retrieves the field and mode properties from `modeRecord` and `fieldRecord`.
    k = modeRecord['k']
    m = fieldRecord['m']
    Xi = fieldRecord['Xi']
    
    # Generates a sample of the frequency function squared for the field mode of wavenumber `k`.
    fSample = FrequencySample(k, m, Xi, xSample, sSample, s1Sample, K)**2
    
    # Generates a sample of the Minkowski frequency function for the field mode of wavenumber `k`.
    mfSample = MinkowskiFrequencySample(k, m, sSample)

    # Performs the analytic continuation of the mode frequency function squared over `Box`.
    fContinuation = AnalyticContinuation(fSample, Box, Order=6, N=-1)
    
    # Computes samples of the phase-integral functions `W` and `V`.
    WSample = PhaseIntegralW(xSample, np.sqrt(fSample), xSample, Order=2, N=-1) #Before: N=len(xSample)
    VSample = PhaseIntegralV(xSample, np.sqrt(fSample), xSample, Order=2, N=-1) #Before: N=len(xSample)

    # Creates a mask to select the entries of `fContinuation` to be recorded.
    newEntries = (xRange > fieldRecord['x'].max())*(xRange <= CropBox(Box)[0].max())

    # Extracts those values of `fContinuation`, `mfSample`, `WSample`, and `VSample` which are new and must be included in `modeRecord`.
    fNew = fContinuation[newEntries]
    mfNew = mfSample[newEntries]
    WNew = WSample[newEntries]
    VNew = VSample[newEntries]

    # Updates `modeRecord` with the new entries obtained above.
    modeRecord['frequency'] = np.concatenate((modeRecord['frequency'], fNew))
    modeRecord['minkowskiFrequency'] = np.concatenate((modeRecord['minkowskiFrequency'], mfNew))
    modeRecord['phaseIntegralW'] = np.concatenate((modeRecord['phaseIntegralW'], WNew))
    modeRecord['phaseIntegralV'] = np.concatenate((modeRecord['phaseIntegralV'], VNew))

    # Returns the updated `modeRecord`.
    return modeRecord


def JumpForward(fieldRecord, modeRecord, metricRecord, Box, Root, f, rho):#, s, h, H, K):
    """ Temporarily steps the mode frequency function `f` forward in time.
        
    Args:
        fieldRecord: a dictionary containing a dictionary entry for each field mode.
        modeRecord: the dictionary entry for a particular field mode in `fieldRecord`.
        metricRecord: a metric dictionary containing NumPy arrays which detail the evolution of the spacetime metric.
        Box: a tuple of the form `(xRange, yRange)`, where the entries are an (n,)-dimensional and (m,)-dimensional NumPy arrays spamming the box dimensions.
        Root: a (2,)-dimensional NumPy array of the form `np.array([x, y])` specifying the root coordinates on the complex plane.
        f: a callable function with two arguments `f(x, y)` represeting the frequency function over the complex plane for the field mode of interest.
        rho: a callable function `rho(x)` which outputs a the field energy density as a function of the time variable `x`.
        
    Returns:
        A callable function with two arguments `f(x, y)` represeting the interpolated mode frequency function over the complex plane stepped forward in time. This function can be used to compute the parameters characterizing the particle production event sourced by `Root`.
        
    """

    # Creates temporary copies of `fieldRecord`, `modeRecord, and `metricRecord`.
    fieldTemp = fieldRecord.copy()
    modeTemp = modeRecord.copy()
    metricTemp = metricRecord.copy()

    # Initializes the variable `eventTime` which stores the most recent estimate of the time of particle production sourced by `Root`.
    eventTime = EventTime(Root, f)
        
    # Establishes the required time increment for stepping `fieldTemp`, `modeTemp`, and `metricTemp` forward in time.
    Increment = np.linspace(Box[0].min(), Box[0].max(), len(Box[0]), retstep=True)[1]

    # Evolves `fieldTemp`, `modeTemp`, and `metricTemp` forward in time until the time of particle production is smaller than `Box[0].mean()`.
    while (eventTime > Box[0].mean()):
    
        # Creates a mask to select the entries of the cropped `Box` which will be recorded in `fieldTemp['x']`.
        newEntries = (CropBox(Box)[0] > fieldTemp['x'].max())
        
        # Interpolates the mode frequency function over the entire region of the complex plane which has been covered.
        _, f = FunctionInterpolation(modeTemp['frequency'], (np.append(fieldTemp['x'], CropBox(Box)[0][newEntries]), fieldTemp['y']))
            
        # Recomputes the time of particle production.
        eventTime = EventTime(Root, f)
        
        # Updates `fieldTemp['x']` with the new `x` coordinate entries.
        fieldTemp['x'] = np.append(fieldTemp['x'], CropBox(Box)[0][newEntries])
            
        # Shifts `tempBox` by one time step.
        Box = ShiftBox(Box, Increment)
          
        # Updates `metricTemp`.
        metricTemp = UpdateMetricNew(metricTemp, fieldRecord, Box, rho)
           
        # Updates the `modeTemp['frequency'].
        modeTemp = UpdateFrequencies(fieldTemp, modeTemp, metricTemp, Box)

    # Once `eventTime <= Box[0].mean()`, the function `f` which stores the interpolated mode frequency stepped forward in time is returned. This function can be used to compute the parameters characterizing the particle production event sourced by `Root`.
    return f


def UpdateEvent(fieldRecord, modeRecord, metricRecord, Box, Root, f, rho):#, s, h, H, K):
    """ Updates `modeRecord['root']`, `modeRecord['eventTime']`, `modeRecord['eventAmplitude']`, `modeRecord['eventWith']`, and `modeRecord['eventPhase']` which describe a particle production event sourced by a `Root` of `f`.
        
    Args:
        fieldRecord: a dictionary containing a dictionary entry for each field mode.
        modeRecord: the dictionary entry for a particular field mode in `fieldRecord`.
        metricRecord: a metric dictionary containing NumPy arrays which detail the evolution of the spacetime metric.
        Box: a tuple of the form `(xRange, yRange)`, where the entries are an (n,)-dimensional and (m,)-dimensional NumPy arrays spamming the box dimensions.
        Root: a (2,)-dimensional NumPy array of the form `np.array([x, y])` specifying the root coordinates on the complex plane.
        f: a callable function with two arguments `f(x, y)` represeting the frequency function over the complex plane for the field mode of interest.
        rho: a callable function `rho(x)` which outputs a the field energy density as a function of the time variable `x`.
        
    Returns:
        An updated version of `modeRecord` which includes all quantities characterizing a new particle creation event.
        
    """

    # Adds `Root` to `modeRecord`.
    modeRecord['root'].append(Root)

    # If the time at which particles are produced happened in the past, all quantities associated with the particle production event are computed and `modeRecord` is updated to incorporate this information.
    if (EventTime(Root, f) <= Box[0].mean()):
    
        # Updates each event entry.
        modeRecord['eventTime'].append(EventTime(Root, f))
        modeRecord['eventAmplitude'].append(EventAmplitude(Root, f))
        modeRecord['eventWidth'].append(EventWidth(Root, modeRecord['eventTime'][-1], f))
        # Old reference point not appropriate for backreaction: modeRecord['eventPhase'].append(EventPhase(modeRecord['eventTime'][-1], modeRecord['eventTime'][0], f))
        modeRecord['eventPhase'].append(EventPhase(modeRecord['eventTime'][-1], fieldRecord['x'][0], f))
    
    # If the time at which particles are produced happens in the future, the metric and the mode under consideration are temporarily stepped forward in time so that all quantities associated with the particle production event can be computed.
    else:

        # Temporarily steps the mode frequency function `f` forward in time in order to compute all quantities associated with the particle production event.
        f = JumpForward(fieldRecord, modeRecord, metricRecord, Box, Root, f, rho)#, s, h, H, K)

        # Once `f` has been appropriately stepped forward in time, each entry of `modeRecord` associated with the particle production event is updated.
        modeRecord['eventTime'].append(EventTime(Root, f))
        modeRecord['eventAmplitude'].append(EventAmplitude(Root, f))
        modeRecord['eventWidth'].append(EventWidth(Root, modeRecord['eventTime'][-1], f))
        # Old reference point not appropriate for backreaction: modeRecord['eventPhase'].append(EventPhase(modeRecord['eventTime'][-1], modeRecord['eventTime'][0], f))
        modeRecord['eventPhase'].append(EventPhase(modeRecord['eventTime'][-1], fieldRecord['x'][0], f))

    # Returns the updated `modeRecord`.
    return modeRecord


def UpdateField(fieldRecord, metricRecord, Box, rho):#, s, h, H, K):
    """ Updates each `modeRecord` in `fieldRecord` to include the analytical continuation of mode's frequency function over `Box`, and information regarding possible new particle production events.
        
    Args:
        fieldRecord: a dictionary containing a dictionary entry for each field mode.
        metricRecord: a metric dictionary containing NumPy arrays which detail the evolution of the spacetime metric.
        Box: a tuple of the form `(xRange, yRange)`, where the entries are an (n,)-dimensional and (m,)-dimensional NumPy arrays spamming the box dimensions.
        rho: a callable function `rho(x)` which outputs a the field energy density as a function of the time variable `x`.
        
    Returns:
        An updated version of `fieldRecord` which includes the analytical continuation of mode frequency functions over `Box` and quantities characterizing new particle production events.
        
    """
    
    # Creates a mask to select the entries of the cropped `Box` which will be recorded in `fieldRecord['x']`. This update is delayed until all modes have been examined.
    newEntries = (CropBox(Box)[0] > fieldRecord['x'].max())
    
    # Loops over all field modes being tracked by `fieldRecord`.
    for k in fieldRecord.keys()[:-5]:
        
        # Creates a copy of the `modeRecord` associated with the mode of wavenumber `k`.
        modeRecord = copy.deepcopy(fieldRecord[k])
    
        # Updates the `modeRecord['frequency'].
        modeRecord = UpdateFrequencies(fieldRecord, modeRecord, metricRecord, Box)
        
        # Interpolates (both squared and non-squared versions of) the mode frequency function over the entire region of the complex plane which has been covered.
        fSquared, f = FunctionInterpolation(modeRecord['frequency'], (np.append(fieldRecord['x'], CropBox(Box)[0][newEntries]), fieldRecord['y']))
        
        # Searches for new roots.
        Root = RootFinder(fSquared, CropBox(Box))
    
        # If a new root has been found which sources a non-negligible particle production amplitude, `modeRecord['root']` is updated to include `Root` and all quantities characterizing the new particle creation event sourced by `Root`.
        if np.any(Root) and RootCheck(modeRecord, Root) and AmplitudeCheck(modeRecord, Root, f):
            modeRecord = UpdateEvent(fieldRecord, modeRecord, metricRecord, Box, Root, f, rho)#, s, h, H, K)

        # Consolidates all `modeRecord` updates into `fieldRecord`.
        fieldRecord[k] = copy.deepcopy(modeRecord)

    # Updates `fieldRecord['x']` with the new `x` coordinate entries.
    fieldRecord['x'] = np.append(fieldRecord['x'], CropBox(Box)[0][newEntries])

    # Returns the updated `fieldRecord`.
    return fieldRecord


def EvolveSystem(metricRecord, fieldRecord, Box, xEnd, Increment, Save=True, Verbose=True, Backreaction=True):
    """ Iteratively updates `metricRecord` and `fieldRecord` according to the semi-classical Friedmann equations.
        
    Args:
        metricRecord: a metric dictionary containing NumPy arrays which detail the evolution of the spacetime metric.
        fieldRecord: a dictionary containing a dictionary entry for each field mode.
        Box: a tuple of the form `(xRange, yRange)`, where the entries are an (n,)-dimensional and (m,)-dimensional NumPy arrays spamming the box dimensions.
        xEnd: a float specifying at which value of the time variable `x` the iteration is to be interrupted.
        Increment: a float specifying the spacing between the interation time steps.
        Save: a boolean parameters which, if set to `True`, continuously saves the intermediate evolution steps in the Pickle files `fieldRecordTemp.p` and `metricRecordTemp.p`.
        Backreaction: a boolean parameter which, if set to `True`, takes into account the backreaction of the produced particles. If set to `False`, the backreaction is turned off.
        
    Returns:
        The evolved `fieldRecord` and `metricRecord` obtained through the backreacting solution of the semi-classical Friednmann equations.
        
    """

    # If `Backreaction == True`, `fieldRecord` and `metricRecord` are updated iteratively:
    if Backreaction == True:
        
        # Initializes the `xIteration` variable which contains the value of the time variable `x` at which a given iteration must stop. Each iteration adds `Increment` to the `xIteration` of the previous iteration.
        xIteration = Box[0].max() + Increment
        
        # Initializes the variables `xSample` and `rhoSample` which encapsulate the enegy density `rho(x)`. These variables are related through `rhoSample = rho(xSample)`. Otherwise, the energy density is set to zero, i.e., `rho(x) = 0`.
        xSample = metricRecord['x']
        rhoSample = EnergyDensity(metricRecord['x'], metricRecord['s'], fieldRecord, metricRecord)
    
        # Updates `modeRecord` and `fieldRecord` through multiple iterations starting from their initial input states. Each new iteration takes into consideration additional particle production events found in the preceeding iteration.
        while xIteration < xEnd:
        
            # Outputs a summary of the last iteration if `Verbose == True`.
            if Verbose == True:
                print "Current iteration time: ", xIteration
        
            # Initializes the iterated dictionaries `metricRecordIteration` and `fieldRecordIteration`.
            metricRecordIteration = copy.deepcopy(metricRecord)
            fieldRecordIteration = copy.deepcopy(fieldRecord)
        
            # Initializes the iterated box `BoxIteration`.
            BoxIteration = Box
        
            # Updates `metricRecordIteration` and `fieldRecordIteration` until the time variable reaches `xIteration`.
            while CropBox(BoxIteration)[0].max() < xIteration:
    
                # Shifts `BoxIteration` by an `Increment`.
                BoxIteration = ShiftBox(BoxIteration, Increment)
            
                # Creates a filter which picks the entries of `xSample` which are within `BoxIteration`.
                withinBox = (xSample >= BoxIteration[0].min())*(xSample <= BoxIteration[0].max())

                # Defines the field energy density function from the last iterations of `xSample` and `rhoSample` which are within `BoxIteration`. Currently `s=0` so that no smoothing is applied to the interpolation.
                rho = lambda x: splev(x, splrep(xSample[withinBox], rhoSample[withinBox], s=0)) # Old weight: w = np.exp(-xSample[withinBox] + xSample[withinBox][0])

                # Updates `metricRecordIteration` given the field changes found in the previous iteration.
                metricRecordIteration = UpdateMetricNew(metricRecordIteration, fieldRecordIteration, BoxIteration, rho)
    
                # Updates `fieldRecordIteration` given the metric changes found in the previous iteration.
                fieldRecordIteration = UpdateField(fieldRecordIteration, metricRecordIteration, BoxIteration, rho)

            # Updates `xSample` and `rhoSample` to reflect possible particle production events found in the previous iteration.
            xSample = metricRecordIteration['x']
            rhoSample = EnergyDensity(metricRecordIteration['x'], metricRecordIteration['s'], fieldRecordIteration, metricRecordIteration)

            # Updates `xIteration` by adding `Increment` to its value.
            xIteration += Increment

            # Saves the intermediary iteration steps if `Save == True`.
            if Save == True:
                cPickle.dump(fieldRecordIteration, open("fieldRecordTemp.p", "wb"))
                cPickle.dump(metricRecordIteration, open("metricRecordTemp.p", "wb"))


    # If `Backreaction == False`, `fieldRecord` and `metricRecord` are updated independently:
    else:
            
        # Initializes the iterated dictionaries `metricRecordIteration` and `fieldRecordIteration`.
        metricRecordIteration = copy.deepcopy(metricRecord)
        fieldRecordIteration = copy.deepcopy(fieldRecord)
        
        # Initializes the iterated box `BoxIteration`.
        BoxIteration = Box
        
        # Sets the field energy density function to zero, i.e., `rho(x) = 0`.
        rho = lambda x: 0
    
        # Updates `metricRecord` and `fieldRecord` until the time variable reaches `xEnd`.
        while CropBox(BoxIteration)[0].max() < xEnd:
            
            # Outputs a summary of the last iteration if `Verbose == True`.
            if Verbose == True:
                print "Current time: ", CropBox(BoxIteration)[0].max()

            # Shifts `Box` by an `Increment`.
            BoxIteration = ShiftBox(BoxIteration, Increment)
    
            # Updates `metricRecord`.
            metricRecordIteration = UpdateMetricNew(metricRecordIteration, fieldRecordIteration, BoxIteration, rho)
    
            # Updates `fieldRecord`.
            fieldRecordIteration = UpdateField(fieldRecordIteration, metricRecordIteration, BoxIteration, rho)

            # Saves the intermediary iteration steps if `Save == True`.
            if Save == True:
                cPickle.dump(fieldRecordIteration, open("fieldRecordTemp.p", "wb"))
                cPickle.dump(metricRecordIteration, open("metricRecordTemp.p", "wb"))


    # Once the iteration process is completed, a copy of the final iterations of `metricRecordIteration` and `fieldRecordIteration` are stored in `metricRecord` and `fieldRecord`.
    metricRecord = copy.deepcopy(metricRecordIteration)
    fieldRecord = copy.deepcopy(fieldRecordIteration)
    
    # Returns `metricRecord` and `fieldRecord`.
    return metricRecord, fieldRecord
