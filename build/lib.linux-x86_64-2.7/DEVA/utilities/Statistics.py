#!/usr/bin/env python
# -*- coding: utf-8 -*-
######################################################################
## Filename:      Statistics.py
## Version:       $Revision: 1.41 $
## Description:   Statistic utilities.
## Author:        Yannick Copin <yannick@ipnl.in2p3.fr>
## Created at:    2008-12-20 21:07:39
## Author:        $Author: ycopin $
## $Id: Statistics.py,v 1.41 2013/01/22 17:47:46 ycopin Exp $
######################################################################

r"""
.. _Statistics:

ToolBox.Statistics - Statistical utilities
==========================================

Statistical tests
-----------------

* `R one- and two-sample tests
  <http://cran.r-project.org/doc/manuals/R-intro.html#One_002d-and-two_002dsample-tests>`_
* `Quantitative Techniques
  <http://www.itl.nist.gov/div898/handbook/eda/section3/eda35.htm>`_
  in `NIST/SEMATECH e-Handbook of Statistical Methods
  <http://www.itl.nist.gov/div898/handbook/index.htm>`_
* `Statistical tests
  <http://www.fao.org/docrep/W7295E/w7295e08.htm#6.4%20statistical%20tests>`_
  in `Guidelines for Quality Management in Soil and Plant Laboratories
  <http://www.fao.org/docrep/W7295E/w7295e00.htm>`_ (FAO Soils
  Bulletin - 74)

Robust statistics
-----------------

* `IDL-astro Robust Statistics Procedures
  <http://idlastro.gsfc.nasa.gov/ftp/pro/robust/aaareadme.txt>`_

Weighted mean and variance
--------------------------

For the inverse-variance weighting :math:`w_i = 1/\sigma_i^2`, the
weighted mean :math:`\mu_w = \sum_i w_i x_i / \sum_i w_i` is the
minimal-:math:`\chi^2` estimator, and its associated error is:

.. math:: \sigma_{\mu_w}^2 = 1/\sum_i w_i.

.. Attention:: :math:`\sigma_i^2` is there supposed to be the *total*
   variance of normally-distributed :math:`x_i`, including both
   measurement error :math:`\epsilon_i` *and* intrinsic distribution
   dispersion :math:`\sigma` if any: :math:`\sigma_i^2 = \epsilon_i^2
   + \sigma^2`.

The reduced :math:`\chi^2`, defined as:

.. math:: \chi_\nu^2 = 1/(N-1) × \sum_i (x_i - \mu_w)^2 / \sigma_i^2,

should equal 1 for :math:`\sigma_i` compatible to the observed
dispersion. If this is not the case, :math:`\sigma_i` has probably
been misestimated, either because it ignored an intrinsic dispersion
:math:`\sigma` or because the measurement errors are misestimated.

- **If the intrinsic dispersion is null** (:math:`\sigma = 0` and
  :math:`\sigma_i^2 = \epsilon_i^2`), the error on the weighted mean
  can be corrected for under- or over-dispersion with respect to the
  estimated measurement errors:

  .. math:: \sigma_{\mu_w}^{*2} = \chi^2_\nu × \sigma_{\mu_w}^2.

  This is formally equivalent to a prior rescaling of the estimated
  errors :math:`\epsilon_i` by a factor :math:`\eta =
  \sqrt{\chi^2_\nu}` to reach a final :math:`\chi^2_\nu` of one.

- When the intrinsic dispersion :math:`\sigma` is **not** null, one
  should not just rescale the total error :math:`\sigma_i` (this would
  over-weight the measurements with least errors), but add the
  intrinsic scatter needed to retrieve a :math:`\chi^2_\nu` of one.

  In practice, one should:

  * *simultaneously* solve :math:`\chi^2_{\nu}(\eta,\sigma) = 1` *and*
    :math:`\sigma_w^2(\eta,\sigma) = \sigma^2` (see weighted variance
    below) for measurement error rescaling :math:`\eta` and intrinsic
    distribution scatter :math:`\sigma`,
  * use properly estimated variances :math:`\sigma_{i}^{2} =
    (\eta\epsilon_{i})^{2} + \sigma^{2}` for weighted mean and its
    associated error.

The *unbiased* weighted variance (wrt. weighted :math:`\mu_w`) is
given by (no matter the chosen weighting):

.. math::

   \sigma_w^2 = \frac{\sum w_i}{\left(\sum w_i\right)^2 - \sum w_i^2} ×
   \sum_i w_i (x_i - \mu_w)^2.

.. Note:: the inverse-variance weighting :math:`w_i = 1/\sigma_i^2`
   probably still corresponds to the minimal-:math:`\chi^2` estimator
   of :math:`\sigma^2` (not demonstrated though).

**Reference:** `Weighted mean
<http://en.wikipedia.org/wiki/Weighted_mean>`_

Weighted RMS
............

The weighted RMS around :math:`\mu_w` is usually defined as:

.. math:: \textrm{WRMS}^2 = \sum_i w_i (x_i - \mu_w)^2 / \sum_i w_i.

.. Attention:: this corresponds to the *biased* weighted variance.

When using inverse-variance weighting :math:`w_i = 1/\sigma_i^2`, one
has:

.. math:: \textrm{WRMS}^2 = \chi^2 × \sigma_{\mu_w}^2.

**Under the assumption that the intrinsic scatter is null**, the
weighted RMS can be used as a proxy for the *rescaled* error on the
weighted mean:

.. math::

   \sigma_{\mu_w}^{*2} &= \chi^2_\nu × \sigma_{\mu_w}^2 \\
                       &= \textrm{WRMS}^2/(N-1).

.. Todo::

   - uniformization of weighted sample variance normalization
   - convergence of :func:`basic_stats`, :func:`weighted_mean` and
     :func:`sample_mean`.
"""

__author__ = "Yannick Copin <y.copin@ipnl.in2p3.fr>"
__version__ = "$Id: Statistics.py,v 1.41 2013/01/22 17:47:46 ycopin Exp $"

import numpy as N
import scipy.stats as S

if not hasattr(N, "percentile"):        # From Numpy 1.4 only
    # Low-key percentile (no support for axis keyword)
    N.percentile = lambda arr,percents: \
                   N.array([ S.scoreatpercentile(arr, p)
                             for p in N.atleast_1d(percents) ]).squeeze()

# R interaction ==============================

def r2py(robject):
    """
    Transform an R object (rpy2 output) into a python dictionary
    """

    import rpy2.robjects as RO

    out = {}
    for n, v in zip(robject.names, robject):
        if isinstance(v, RO.Vector) and v.names.r_repr() != 'NULL':
            #going recursively
            out[n] = r2py(v)
        elif isinstance(v, (int, float)):
            out[n] = v
        else:
            #out[n] = list(v) if len(list(v)) > 1 else list(v)[0] # >= 2.5
            out[n] = len(list(v))>1 and list(v) or list(v)[0]

    return out


# Statistic functions ==============================

def proportion_CI(success, sample, cl=0.6827, relative=True):
    """Compute binomial population proportion -- observed *success*
    count/*sample* size -- and its confidence interval *cl* from a
    Bayesian approach.

    Reference: Cameron, 2011PASA...28..128C

    >>> [ round(t,5) for t in proportion_CI(6, 27, cl=0.6827, relative=True) ]
    [0.22222, -0.05948, 0.09806]
    >>> [ round(t,5) for t in proportion_CI(60, 270, cl=0.6827, relative=False) ]
    [0.22222, 0.199, 0.24953]
    """

    assert 0 < cl < 1, "Confidence level should be in ]0,1["

    p = float(success)/sample

    bppf = S.distributions.beta.ppf
    p_lower = bppf(  (1-cl)/2.,success+1,sample-success+1)
    p_upper = bppf(1-(1-cl)/2.,success+1,sample-success+1)

    if relative:
        # Confidence interval relative to estimated proportion,
        # e.g. (0.22,-0.06,+0.10) for 22_{-6}^{+10}%
        return (p,p_lower-p,p_upper-p)
    else:
        # Absolute confidence interval
        # e.g. (0.22,0.16,0.32) for 22_{-6}^{+10}%
        return (p,p_lower,p_upper)


def pull(x, dx):
    r"""Compute the pull from *x*, *dx*.

    * Input data: *x* = :math:`x_i`, associated error *dx* =
      :math:`s_i`
    * Optimal (variance-weighted) mean:
      :math:`E(x) = (\sum_i x_i/s_i^2) / (\sum_i 1/s_i^2)`
    * Variance on weighted mean: :math:`\sigma_E^2 = 1/(\sum_i 1/s_i^2)`
    * Pull: :math:`p_i = (x_i - E_i(x))/\sqrt{\sigma_{E_i}^2 + s_i^2}`
      where :math:`E_i` and :math:`\sigma_{E_i}^2` are computed
      *without* point i.

    If errors *dx* are correct, the pull distribution is centered on 0
    with standard deviation of 1.

    >>> import numpy as N; N.random.seed(0); 
    >>> pulls = pull(N.random.normal(size=100),N.ones(100));
    >>> print pulls.mean(),pulls.std(ddof=1)
    -1.23685783837e-17 1.01806287345
    """

    assert x.ndim==dx.ndim==1, "pull implemented on 1D-arrays only"
    assert len(x)==len(dx), "x and dx should be of the same length"

    n = len(x)

    i = N.resize(N.arange(n), n*n) # 0,...,n-1,0,...n-1,..., n times (n²,)
    i[::n+1] = -1                  # Mark successively 0,1,2,...,n-1
    j = i[i>=0].reshape((n,n-1))   # Remove marked indices & reshape (n,n-1)

    v = dx**2                      # Variance
    w = 1/v                        # Variance (optimal) weighting

    # Weighted mean, and variance on weighted mean (n,)
    Ei,wsum = N.average(x[j], weights=w[j], axis=1, returned=True)
    vEi = 1./wsum

    p = (x - Ei)/N.sqrt(vEi + v)   # Pull (n,)

    return p


# Sample statistics ==============================

def weighted_mean(y, dy, weighted=True, addIntrinsic=False):
    """Compute [variance-weighted] mean and standard deviation of mean
    for input array *y*.

    Regarding weighted (intrinsic) standard deviation, GSL advocates a
    normalization by `(w.sum()**2 - (w**2).sum())/w.sum()`, which
    reverts to (N-1)/N for unity-weighting. However, in agreement with
    `N.std` and `N.var` using N-normalization, we uses plain
    `w.sum()`.

    .. Note:: When adding intrinsic standard deviation, aren't we
              surestimating dy by ~sqrt(2)? It's normal to measure an
              intrinsic variance when measures have errors dy, which
              is accounted for in variance on mean.
    """

    if addIntrinsic:
        import warnings
        warnings.warn("weighted_mean: addIntrinsic keyword is deprecated",
                      DeprecationWarning)

    bad = N.isnan(y) | N.isnan(dy)      # Discard NaNs if any
    l  = N.asarray(y)[~bad]
    dl = N.asarray(dy)[~bad]

    v = dl**2                           # Variance

    if weighted:                        # 1/Var optimal weighting
        w = 1/v
        my,wsum = N.average(l, weights=w, returned=True) # Weighted mean
        # Variance on mean: ((w*dl)**2).sum() / w.sum()**2 = 1/w.sum()
        var = 1/wsum
        if addIntrinsic:
            wdelt = w*(l-my)
            # var += (wdelt**2).sum() / wsum # Intrinsic variance
            var += N.dot(wdelt,wdelt) / wsum # Faster (~x5)
    else:                               # No weighting
        my = l.mean()                   # Plain mean
        var = v.sum() / len(l)**2       # Variance on mean
        if addIntrinsic:
            var += l.var()              # Intrinsic variance

    dmy = N.sqrt(var)                   # Standard error on [weighted] mean

    return my,dmy


def basic_stats(a, weights=None, axis=None):
    """Compute [weighted] mean and standard deviation of array *a*
    along *axis*.
    
    Reference: `GSL Weighted samples
    <http://www.gnu.org/software/gsl/manual/html_node/Weighted-Samples.html>`_
    and `Weighted sample variance
    <http://en.wikipedia.org/wiki/Weighted_mean#Weighted_sample_variance>`_.
    """

    if weights is not None:
        assert (weights > 0).all(), "weights are not all strictly positive"

    # [Weighted] mean
    mean = N.average(a, axis=axis, weights=weights)
    if axis is None:
        umean = mean                    # Scalar
        n = N.size(a)                   # Nb of elements
    else:
        umean = N.expand_dims(mean, axis) # Same ndim as a
        n = N.shape(a)[axis]
    var = N.average((a - umean)**2, axis=axis, weights=weights)
    if weights is None:                 # No weighting, use (N-1) normalization
        if n>1:
            var *= n/(n-1.)
        else:
            var = N.nan
    else:                               # Use unbiased normalisation
        sumw  = N.sum(weights,    axis=axis)
        sumw2 = N.sum(weights**2, axis=axis)
        denom = sumw**2 - sumw2
        if denom > 0:
            var *= sumw**2 / denom
        else:
            var = N.nan

    return mean,N.sqrt(var)             # Mean and std deviation


def sample_mean(x, dx=None, axis=None, corrected=False):
    """Compute arithmetic mean and standard error on mean of *x* along
    given *axis*. If *dx*, the computations are inverse
    variance-weighted (`w=1/dx**2`), and potentially *corrected* for
    intrinsic dispersion ('intrinsic') and error scaling ('scaling').

    Reference: `Weighted mean
    <http://en.wikipedia.org/wiki/Weighted_mean>`_
    """

    if dx is None:                        # Plain average
        mean = N.average(x, axis=axis)    # Mean
        var = N.var(x, axis=axis, ddof=1) # Variance on mean
        if axis is None:
            var /= N.size(x)
        else:
            var /= N.shape(x)[axis]
    else:                                 # Variance-weighted average
        assert (dx > 0).all(), "errors are not all strictly positive"
        if not corrected:
            w = 1./dx**2                  # Inverse-variance weighting
        else:
            raise NotImplementedError 
        mean,wsum = N.average(x, axis=axis, weights=w, returned=True)
        var = 1./wsum                   # Variance on mean

    return mean,N.sqrt(var)             # Mean and associated error


def sample_median(x, dx=None, axis=None, scale=1.2533):
    """Compute median and (normal) error on median of *x* along given
    *axis*.

    The standard error of the median for **large samples and normal
    distributions** is :math:`\sqrt{\pi/2}×s/\sqrt{n}`, where
    dispersion `s` will be estimated from :func:`nMAD`.

    Source: `Sampling Distribution of Median
    <http://davidmlane.com/hyperstat/A106993.html>`_
    """

    if dx is not None:
        if axis is not None:
            raise NotImplementedError("implemented on 1D-arrays only")
        weights = N.where(dx > 0, 1/dx**2, 0)
    else:
        weights = None

    # Median and nMAD (default scale = 1.2533 = sqrt(pi/2))
    med,err = median_stats(x, weights=weights, axis=axis)
    if axis is None:
        err *= scale/N.sqrt(N.size(x))
    else:
        err *= scale/N.sqrt(N.shape(x)[axis])

    return med,err


# Robust statistics ==============================

def wpercentile(a, q, weights=None):
    """Compute weighted percentiles *q* [%] of input 1D-array *a*."""

    a = N.asarray(a)
    if a.ndim > 1:
        raise NotImplementedError("implemented on 1D-arrays only")

    if weights is None:
        weights = N.ones_like(a)
    else:
        assert len(weights)==len(a), "incompatible weight and input arrays"
        assert (weights>0).all(), "weights are not always strictly positive"

    isorted = N.argsort(a)
    sa = a[isorted]
    sw = weights[isorted]
    sumw = N.cumsum(sw)                        # Strictly increasing
    scores = 1e2*(sumw - 0.5*sw)/sumw[-1]      # 0-100 score at center of bins

    def interpolate(q):
        i = scores.searchsorted(q)
        if i==0:                        # Below 1st score
            val = sa[0]
        elif i==len(a):                 # Above last score
            val = sa[-1]
        else:                           # Linear score interpolation
            val = (sa[i-1]*(scores[i] - q) + sa[i]*(q - scores[i-1])) / \
                  (scores[i] - scores[i-1])
        return val

    out = N.array([ interpolate(qq) for qq in N.atleast_1d(q) ])

    return out.reshape(N.shape(q))      # Same shape as input q


def median_stats(a, weights=None, axis=None, scale=1.4826):
    """Compute [weighted] median and :func:`nMAD` of array *a* along
    *axis*. Weighted computation is implemented for *axis* = None
    only."""

    if weights is not None:
        if axis is not None:
            raise NotImplementedError("implemented on 1D-arrays only")
        else:
            med  = wpercentile(a, 50., weights=weights)
            nmad = wpercentile(N.absolute(a - med), 50., weights=weights)*scale
    else:
        med = N.median(a, axis=axis)
        if axis is None:
            umed = med                      # Scalar
        else:
            umed = N.expand_dims(med, axis) # Same ndim as a
        nmad = N.median(N.absolute(a - umed), axis=axis) * scale

    return med,nmad


def nMAD(a, weights=None, axis=None, scale=1.4826):
    """Normalized Median Absolute Deviation (nMAD) along given *axis*
    of array *a*::

      median(abs(a - median(a))) * scale

    For normally distributed data, *scale* should be set to::

      1/scipy.stats.norm.ppf(0.75) = 1.4826022185056018...

    Source: `Median Absolute Deviation
    <http://en.wikipedia.org/wiki/Median_absolute_deviation>`_
    """

    med,nmad = median_stats(a, axis=axis, weights=weights, scale=scale)
    return nmad


def medcov(a, scale=1.4826):
    """Median covariance of *a*, representing *ncols* realisations of
    *nrows* variables."""

    # Remove median: (nr,nc) - (nr,1) → (nr,nc)
    arr = N.array(a, ndmin=2)
    arr -= N.median(arr, axis=1).reshape(-1,1)
    # Compute cross-product: (nr,1,nc) * (nr,nc) → (nr,nr,nc)
    arr = N.expand_dims(arr, 1) * arr

    # Compute nMAD along last axis → (nr,nr)
    return N.median(arr, axis=-1) * scale**2


def iqr(a, axis=None, percentile=0.25):
    """Return inter-quartile range [25%,75%] of array *a* along
    *axis*.

    For normally distributed data, *iqr* is related to standard
    deviation through::

      std = iqr / ( 2*scipy.stats.norm.isf(percentile) )
    """

    l,h = N.percentile(a, [percentile*100,(1-percentile)*100], axis=axis)

    return h-l


def trim(arr, axis=None, trimming=(0.25,0.75), winsorized=False, whiskers=0):
    """Return a trimmed array (as a masked array) along
    *axis*. *Trimming* is expressed as fraction (in [0,1]). If
    *winsorized*, replace trimmed value by largest/smallest remaining
    value. If *whiskers* > 0, extend acceptable range by
    `whiskers*IQR`."""

    assert 0 <= trimming[0] <= trimming[1] <= 1, \
        "Trimming values should be expressed as sorted [0,1] fractional values"

    # Low and high values
    l,h = N.percentile(arr, [trimming[0]*100,trimming[1]*100], axis=axis)
    if axis is not None:
        l = N.expand_dims(l, axis)
        h = N.expand_dims(h, axis)

    # Whiskers
    if whiskers:
        l -= whiskers*(h-l)
        h += whiskers*(h-l)

    # Trimmed array (trimmed values are masked)
    tarr = N.ma.masked_where((arr<l)|(arr>h),arr)

    if winsorized:              # Replace trimmed values by remeaning min/max
        # Mins and maxs of trimmed array
        if axis is None:        # mins and maxs are scalars
            mins = tarr.min()
            maxs = tarr.max()
        else:                   # mins and maxs are masked arrays
            mins = N.expand_dims(tarr.min(axis=axis).filled(), axis)
            maxs = N.expand_dims(tarr.max(axis=axis).filled(), axis)
        # Returned array is a masked array for consistency, even
        # though there's no masked values
        tarr = N.ma.array(arr)
        tarr = N.where(tarr<l, mins, tarr) # Replace trimmed values
        tarr = N.where(tarr>h, maxs, tarr)

    if axis is None:
        tarr = tarr.ravel()

    return tarr


def robust_mean(arr, axis=None, trimming=(0.25,0.75), whiskers=1.5):
    """Mean computed over the [1st quartile - whiskers×IQR, 3rd
    quartile + whiskers×IQR] range."""

    tarr = trim(arr, axis=axis, trimming=trimming, whiskers=whiskers)
    if axis is None:
        return tarr.mean()      # Scalar
    else:
        return tarr.mean(axis=axis).filled() # Array


def robust_std(arr, axis=None, trimming=(0.25,0.75), whiskers=1.5):
    """Standard deviation computed over the [1st quartile -
    whiskers×IQR, 3rd quartile + whiskers×IQR] range."""

    tarr = trim(arr, axis=axis, trimming=trimming, whiskers=whiskers)
    if axis is None:
        return tarr.std(ddof=1)      # Scalar
    else:
        return tarr.std(axis=axis, ddof=1).filled() # Array


# Distributions ==============================

def logNormal(mu, sig):
    """Return log-normal distribution with given *mu* and *sigma*.

    Source: `Log-normal <http://en.wikipedia.org/wiki/Log-normal>`_
    """

    # sqrt(2*pi) = 2.5066282746310002
    return lambda x: \
           N.exp(-0.5*((N.log(x)-mu)/sig)**2) / (x*sig*2.5066282746310002)


def logNormalMeanStd(mean, std):
    """Return log-normal distribution with given *mean* and *std*."""

    sig2 = N.log(1 + (std/mean)**2)
    mu = N.log(mean) - 0.5*sig2
    sig = N.sqrt(sig2)

    print "Log-normal distribution with mean=%.2f, std=%.2f: " \
          "mu=%.2f, sigma=%.2f" % (mean, std, mu, sig)

    return logNormal(mu, sig)


# Statistical tests ==============================

def t_test(A, B=None, alt='<>', cl=0.95, equalVar=False,
           mu=0, paired=False, verbose=False):
    """
    Test null-hypothesis mean(*A*) = mean(*B*)+*mu* against
    alternative hypothesis mean(*A*) *alt* mean(*B*) + *mu* where
    *alt* = '<', '>' or '<>' (default), at (fractional) confidence
    level *cl*, with or without assuming the equality of variances.

    Returns *t*, *df*, *pt*, *tint* (statistic, DoF, p-value and
    confidence interval). The alternative hypothesis is rejected at
    the requested confidence level if tint[0] < t < tint[1].

    See equivalent `t.test
    <https://svn.r-project.org/R/trunk/src/library/stats/R/t.test.R>`_
    in R with: t.pdf=dt, t.cdf=pt, t.sf=pt(lower.tail=FALSE),
    t.ppf=qt, t.isf(x)=qt(1-x)
    """

    if B is None:
        raise NotImplementedError
    if paired:
        raise NotImplementedError

    n1 = len(A)
    m1 = N.mean(A)
    v1 = N.var(A, ddof=1)       # uses N-1
    se12 = v1/n1

    if B is not None:
        n2 = len(B)
        m2 = N.mean(B)
        v2 = N.var(B, ddof=1)   # uses N-1
        se22 = v2/n2

    if equalVar:
        title = " t-test for equal means (equal variance) "
        df = n1 + n2 - 2
        stderr = N.sqrt( ((n1-1)*v1 + (n2-1)*v2)/df * (1./n1 + 1./n2) )
    else:
        title = " Welch t-test for equal means (unequal variance) "
        # Welch-Satterthwaite equation (non-equal variance)
        stderr = N.sqrt( se12 + se22 )
        df = (se12 + se22)**2 / ( se12**2/(n1-1) + se22**2/(n2-1) )
    t = (m1 - m2 - mu) / stderr

    if alt == '<>':
        # 2-sided p-value: >0 or <0.
        # S.t.cdf(-t,df) + S.t.sf(t,df) = 2*S.t.sf(t,df)
        pt = 2*S.t.cdf(-abs(t),df)
        dt = S.t.isf((1-cl)/2,df)
        tint = N.array([t-dt, t+dt])*stderr + mu
    elif alt == '>':
        # greater p-value: >0
        pt = S.t.sf(t,df)
        tint = N.array([t - S.t.ppf(alpha,df), N.inf])*stderr + mu
    elif alt == '<':
        # less p-value: <0, also = S.t.sf(-t,df)
        pt = S.t.cdf(t,df) # = 1 - sf
        tint = N.array([-N.inf, t + S.t.ppf(alpha,df)])*stderr + mu
    else:
        raise ValueError("Unknown alternative '%s'" % alt)

    if verbose:
        print title.center(70,'=')
        print "A: %d elements, mean=%f, stdErr(mean)=%f" % \
            (n1,m1,N.sqrt(se12))
        if B is not None:
            print "B: %d elements, mean=%f, stdErr(mean)=%f" % \
                (n2,m2,N.sqrt(se22))
        print "Observed difference (B-A)=%f, stdErr(diff)=%f" % \
            (m2-m1,stderr)
        print "Null hypothesis: true difference of means = %s" % mu
        print "Alternative hypothesis: true difference of means %s %s" % \
            (alt, mu)
        print "Statistic: t=%f, df=%f, p-value=%f" % (t,df,pt)
        print "Interval (%d%%): %f,%f" % (cl*100,tint[0],tint[1])
        print "Alternative hypothesis [%s]: %s at the %d%% confidence level" % \
            (alt, (tint[0]<t<tint[1]) and "REJECTED" or "ACCEPTED", cl*100)

    return t,df,pt,tint


def F_test(A, B, alt='<>', cl=0.95, ratio=1, verbose=False, R=False):
    """
    Test null-hypothesis var(*A*) = var(*B*) × *ratio* against
    alternative hypothesis var(*A*) *alt* var(*B*) × *ratio* where
    *alt* = '<', '>' or '<>' (default), at (fractional) confidence
    level *cl*.

    Returns *F*, (*numdf*, *denumdf*), *pf*, *fint* (statistic, DoFs
    of numerator and denominator, p-value and confidence
    interval). The alternative hypothesis is rejected at the requested
    confidence level if fint[0] < F < Fint[1].

    See equivalent `var.test
    <https://svn.r-project.org/R/trunk/src/library/stats/R/var.test.R>`_
    in R with: f.pdf=df, f.cdf=pf, f.sf=pf(lower.tail=FALSE),
    f.ppf=qf, f.isf(x)=qf(1-x)
    """

    n1 = len(A)
    m1 = N.mean(A)
    v1 = N.var(A, ddof=1)       # uses N-1

    n2 = len(B)
    m2 = N.mean(B)
    v2 = N.var(B, ddof=1)       # uses N-1

    if R:                       # Use R
        import rpy2.robjects as RO
        #necessary to make rpy2 accept numpy arrays
        import rpy2.robjects.numpy2ri
        altdict = {'>': 'g', '<': 'l', '<>': 't'}
        # R output is a python dictionary
        Rout = r2py(RO.r('var.test')(A,B, ratio=ratio,
                                     alternative=altdict[alt], conf_level=cl))

        F = Rout['statistic']['F']
        ndf, ddf =  Rout['parameter']['num df'], Rout['parameter']['denom df']
        pf = Rout['p.value']
        fint = Rout['conf.int']

    else:
        F = v1/v2 / ratio
        ndf = n1 - 1
        ddf = n2 - 1

        if alt == '<>':
            # two-sided p-value: != ratio
            # S.f.cdf(1/F,n2-1,n1-1) + S.f.sf(F,n1-1,n2-1) =
            # 2 * S.f.sf(F,n1-1,n2-1) = 2 * min( pf(greater),pf(less) )
            pf = S.f.sf(F,ndf,ddf)
            pf = 2*min(pf,1-pf)
            fint = [F/S.f.isf((1-cl)/2,ndf,ddf), F/S.f.ppf((1-cl)/2,ndf,ddf)]
        elif alt == '>':
            # greater p-value: > ratio
            pf = S.f.sf(F,ndf,ddf) # = S.fprob(n1-1,n2-1,F)
            fint = [F/S.f.ppf(cl,ndf,ddf), N.inf]
        elif alt == '<':
            # less p-value: < ratio, also = S.f.sf(1/F,n2-1,n1-1)
            pf = S.f.cdf(F,ndf,ddf) # = 1 - sf
            fint = [0, F/S.f.isf(1-cl,ndf,ddf)]
        else:
            raise ValueError("Unknown alternative '%s'" % alt)

    if verbose:
        print " F-test for equal variances ".center(70,'=')
        print "A: %d elements, mean=%f, stdErr=%f" % (n1,m1,N.sqrt(v1))
        print "B: %d elements, mean=%f, stdErr=%f" % (n2,m2,N.sqrt(v2))
        print "Null hypothesis: true ratio of variances = %s" % ratio
        print "Alternative hypothesis: true ratio of variances %s %s" % \
            (alt, ratio)
        print "Statistic: F=%f, df=%d/%d, p-value=%f" % (F,ndf,ddf,pf)
        print "Interval (%d%%): %f,%f" % (cl*100,fint[0],fint[1])
        print "Alternative hypothesis [%s]: %s at the %d%% confidence level" % \
            (alt, (fint[0]<F<fint[1]) and "REJECTED" or "ACCEPTED", cl*100)

    return F,(ndf,ddf),pf,fint


def KruskalWallis_test(samples, sl=0.05, verbose=False):
    """
    Perform a Kruskal-Wallis test for comparing *samples* with unknown
    distributions at significance level *sl*. Null hypothesis is: 'all
    samples come from the same population'.

    Return (*htest*, *cval*). Null hypothesis is rejected at the
    requested significance level if htest>cval.

    Example:

    >>> KruskalWallis_test([[4.2,4.6,3.9,4.0,4.5],
    ...                     [3.3,2.4,2.6,3.8,2.8],
    ...                     [1.9,2.4,2.1,2.7,1.8],
    ...                     [3.5,3.1,3.7,4.1,4.4]], sl=0.05)
    (14.768571428571434, 7.8147279032511792)

    Since 14.8>7.8, null hypothesis -- all samples come from the same
    population -- is *rejected* at the 5%-significance level
    (i.e. null hypothesis is rejected with a risk of 5% that it was
    actually true). Conversely, alternative hypothesis -- samples come
    from different populations -- is *accepted* at the 95%-confidence
    level (i.e. it is accepted with 95% of chances that it is really
    true).

    Source: `Kruskal-Wallis test
    <http://www.itl.nist.gov/div898/handbook/prc/section4/prc41.htm>`_
    """

    # Test individual sample size
    for i,sample in enumerate(samples):
        if len(sample)<=4:
            print "WARNING: Sample #%d is too small (%d observations)" % \
                  (i+1,len(sample))

    nsamp = len(samples)                # Nb of samples
    all = N.concatenate(samples)
    nobs = len(all)                     # Total nb of observations

    ranks = S.rankdata(all)             # Ranks (w/ handling of ties)
    # Boundaries in 'all' such that samples = split(all,bounds)
    bounds = N.cumsum([ len(s) for s in samples[:-1] ])

    # Test statistic
    htest = 12./(nobs*(nobs+1)) * \
            sum([ r.sum()**2/len(r) for r in N.split(ranks,bounds) ]) \
            - 3*(nobs+1)
    # Critical value: chi2(1-sl,dof=nsamp-1)
    cval = S.chi2.ppf(1-sl,nsamp-1)

    if verbose:                         # Sample summary
        print " Kruskal-Wallis Test ".center(70,'=')
        print "%d observations in %d samples" % (nobs,nsamp)
        for i,(sample,r) in enumerate(zip(samples, N.split(ranks,bounds))):
            print "Sample #%d: %d observations, min=%.2f, max=%.2f" % \
                  (i+1, len(sample),min(sample),max(sample))
            print "   Mean=%.2f, stddev=%.2f" % \
                  (N.mean(sample), N.std(sample, ddof=1))
            print "   Median=%.2f, nMAD=%.2f" % \
                  (median_stats(sample))
            print "   Sum of ranks=%.1f" % (r.sum())
            #print "Sample:",sample
            #print "Ranks:",r
        print "Null hypothesis: all samples come from the same population"
        print "Alt. hypothesis: samples come from different populations"
        print "Test statistics: H=%f" % htest
        print "Critical value: chi2(SL=%.2f, dof=%d)=%f" % \
              (sl,nsamp-1,cval)
        print "Null hypothesis %s at the %.1f%% significance level" % \
              (htest>cval and "REJECTED" or "ACCEPTED", sl*100)
        print "Alt. hypothesis %s at the %.1f%% confidence level" % \
              (htest<cval and "REJECTED" or "ACCEPTED", (1-sl)*100)

    return htest,cval


def pvalue2sigma(p):
    """Express the input one-sided *p*-value as a sigma equivalent
    significance from a normal distribution (the so-called z-value).

    =====  =======  =================
    sigma  p-value  terminology
    =====  =======  =================
    1      0.1587
    1.64   0.05     significant
    2      0.0228
    2.33   0.01     higly significant
    3      0.0013   evidence
    3.09   0.001
    5      2.9e-7   discovery
    =====  =======  =================

    >>> pvalue2sigma(1e-3) # p=0.1% corresponds to a ~3-sigma significance
    3.0902323061678132
    """

    return S.distributions.norm.isf(p)  # isf = ppf(1 - p)


def sigma2pvalue(s):
    """Express the input *s*-sigma significance from a normal
    distribution as a one-sided p-value (two-sided p-value is twice
    the one-sided value).

    >>> sigma2pvalue(3)       # one-sided p-value at 3-sigma
    0.0013498980316301035
    >>> 1 - 2*sigma2pvalue(3) # two-sided confidence level at 3-sigma
    0.99730020393673979
    """

    return S.distributions.norm.sf(s)   # sf = 1 - cdf

# Chi2 ==============================

def comp_chi2(y, dy, m=None, axis=None):
    r"""Compute chi2 of observations (*y*, *dy*) with respect to to model *m*:

    .. math:: \chi^2 = \sum ((y-m)/dy)^2

    where the summation is done over *axis*. The model *m* defaults to
    inverse-variance weighted mean of *y* if *axis* is `None` (no
    default otherwise).
    """

    if m is None:                       # Default to weighted average of y
        if axis is None:
            m = N.average(y, weights=1/dy**2)
        else:
            # It is not possible to know which axis to average over
            raise NotImplementedError

    return (((y - m)/dy)**2).sum(axis=axis)

def prob_chi2(chi2, dof, sigma=False):
    """Compute one-tailed p-value of *chi2*, i.e. the probability that
    x > *chi2* when x follows a chi2-distribution with *dof* degrees
    of freedom. If *sigma*, express the result as a sigma equivalent
    significance from a normal distribution (:func:`pvalue2sigma`).
    """

    p = S.distributions.chi2.sf(chi2, dof) # Survival function
    if sigma:
        p = pvalue2sigma(p)

    return p

# Correlation coefficients ==============================

def correlation_CI(rho, n, cl=0.95):
    """Compute Pearson's correlation coefficient confidence interval,
    for (fractional) confidence level *cl*.

    cl=0.6827 corresponds to a 1-sigma error, 0.9973 for a 3-sigma
    error (`2*scipy.stats.norm.cdf(n)-1` for a n-sigma error).

    Sources: `Confidence Interval of rho
    <http://vassarstats.net/rho.html>`_, `Correlation CI
    <http://onlinestatbook.com/chapter8/correlation_ci.html>`_
    """

    assert -1 < rho < 1, "Correlation coefficient should be in ]-1,1["
    assert n >= 6, "Insufficient sample size"
    assert 0 < cl < 1, "Confidence level should be in ]0,1["

    z = N.arctanh(rho)                  # Fisher's transformation
    # z is normally distributed with std error = 1/sqrt(N-3)
    zsig = S.distributions.norm.ppf(0.5*(cl+1)) / N.sqrt(n-3)
    # Confidence interval on z is [z-zsig,z+zsig]

    return N.tanh([z-zsig,z+zsig])      # Confidence interval on rho


def correlation_significance(rho, n, directional=True, sigma=False):
    """Significance of (Pearson's or Spearman's) correlation
    coefficient *rho*, given the size *n* of the sample.

    If non-*directional*, this is the (two-sided) probability `p =
    Prob_N(|r| >= |rho|)` to find such an extreme correlation
    coefficient (no matter the sign of the correlation) from a purely
    uncorrelated population. If *directional*, this is the (one-sided)
    probability `p = Prob_N(r > rho)` for rho>0 (resp. `Prob_N(r <
    rho)` for rho<0). The directional (one-sided) probability is just
    half the non-directional (two-sided) one.

    If *sigma*, express the result as a sigma equivalent significance
    from a normal distribution (:func:`pvalue2sigma`).

    Sources: *Introduction to Error Analysis* (Taylor, 1997),
    `Significance of a Correlation Coefficient
    <http://vassarstats.net/rsig.html>`_
    """

    assert -1 < rho < 1, "Correlation coefficient should be in ]-1,1["
    assert n >= 6, "Insufficient sample size"

    # t is distributed as Student's T distribution with DoF=n-2
    t = rho * N.sqrt((n - 2.)/(1 - rho**2))
    p = S.distributions.t.sf(t, n-2)    # directional (one-sided) p-value
    if sigma:
        p = pvalue2sigma(p)
    elif not directional:               # non-directional (two-sided) p-value
        p *= 2

    return p


def correlation(x, y, method='pearson', error=False, confidence=0.6827):
    """Compute Pearson/Spearman (unweighted) coefficient correlation
    rho between *x* and *y*.

    If `error=True`, returns non-symmetric errors on rho, for a given
    confidence, see :func:`correlation_CI` (only implemented for
    Pearson).
    """

    assert len(x)==len(y), "Incompatible input arrays x and y"

    if method.lower()=='pearson':
        rho,p = S.pearsonr(x,y)
    elif method.lower()=='spearman':
        rho,p = S.spearmanr(x,y)
    else:
        raise ValueError("Unknown correlation method '%s'" % method)

    if not error:
        return rho

    # Compute error on correlation coefficient
    if method.lower() != 'pearson':
        raise NotImplementedError("Error on correlation coefficient is " \
                                  "implemented for Pearson's correlation only.")

    rho_dn,rho_up = correlation_CI(rho, n=len(x), cl=confidence)

    return rho,rho-rho_dn,rho_up-rho    # Assymetric errors


def correlation_weighted(x, y, wx=None, wy=None, axis=None):
    """Compute (weighted) Pearson correlation coefficient between *x*
    and *y*.

    .. Note:: Checked against `scipy.stats.pearsonr`.
    """

    x = N.asarray(x)
    y = N.asarray(y)
    assert x.shape==y.shape, "Incompatible input arrays x and y"

    if x.ndim != 1 or axis is not None:
        raise NotImplementedError("implemented for 1D-arrays only")

    # Weights
    if wx is None:
        wx = N.ones_like(x)
    else:
        wx = N.where(N.isfinite(wx),wx,0) # Discard NaN's and Inf's
        assert wx.shape==x.shape, "Incompatible arrays x and wx"
    if wy is None:
        wy = N.ones_like(y)
    else:
        wy = N.where(N.isfinite(wy),wy,0)
        assert wy.shape==y.shape, "Incompatible arrays y and wy"

    n = len(x)
    # Weighted means
    mx = N.average(x, weights=wx, axis=axis)
    my = N.average(y, weights=wy, axis=axis)
    xm = x - mx
    ym = y - my
    # Weighted covariance
    cov,sumwxy = N.average(xm*ym, weights=wx*wy, axis=axis, returned=True)
    cov *= (1 - 1/sumwxy)               # = cov(ddof=0) * (1 - 1/sumwxy)**2
    # Weighted variance
    vx,sumwx = N.average(xm**2, weights=wx, axis=axis, returned=True)
    vx *= (1 - 1/sumwx)                 # = var(x) * (1-1/sumwx)**2
    vy,sumwy = N.average(ym**2, weights=wy, axis=axis, returned=True)
    vy *= (1 - 1/sumwy)

    # Weighted correlation
    r = cov/N.sqrt(vx*vy)
    assert -1 <= r <= +1, "Invalid correlation coefficient: r = %f" % r

    return r


def weighted_correlation(*args, **kwargs):

    import warnings
    warnings.warn("Deprecated: use `correlation weighted`", DeprecationWarning)

    return correlation_weighted(*args, **kwargs)


def correlation_MRRtest(rx1y, rx2y, rx1x2, n):
    """Meng-Rosenthal-Rubin method comparing two correlated
    correlation coefficients.

    :param float rx1y: the correlation coefficient between variables x1 and y
    :param float rx2y: the correlation coefficient between variables x2 and y
    :param float rx1x2: the correlation coefficient between variables x1 and x2
    :param int n: sample size
    :return: two-tailed p-value

    Reference: Meng et al. (1992), Psychol. Bull. 111, 172,
    """

    assert (-1<=rx1y<=1) and (-1<=rx2y<=1) and (-1<=rx1x2<=1), \
           "Correlation coefficients should be in ]-1,1["

    from math import sqrt, atanh

    z1 = atanh(rx1y)                    # Fisher's Z-transformation
    z2 = atanh(rx2y)
    rs = (rx1y**2 + rx2y**2)/2.
    f = min((1 - rx1x2)/(2*(1 - rs)), 1)
    h = (1 - f*rs) / (1 - rs)
    dz = (z1 - z2) * sqrt((n - 3)/(2 * (1 - rx1x2) * h))
    p = 2 * S.distributions.norm.sf(abs(dz)) # Two-tailed probability

    return p


# Histogram utilities ==============================

def get_range(x, range=None, log=False, percentiles=False):
    """Get range from *x* and *range*=(*min*,*max*) or `None`. If
    *min* (resp. *max*) is `None`, *xmin* is set to min of *x*, or of
    strictly positive *x* if *log* (resp. *xmax* is set to the max of
    *x*). If *percentiles*, *range* is actually expressed in
    percentiles (in percents).

    >>> import numpy as N
    >>> get_range(N.linspace(0,10,101), range=(5,95), percentiles=True)
    (0.5, 9.5)
    """

    if range is not None:               # Specified range
        if percentiles:                 # Range in percentiles
            vmin,vmax = range
            if vmin is None:            # Take the min
                vmin = 0
            if vmax is None:            # Take the max
                vmax = 100
            xmin,xmax = N.percentile(x, (vmin,vmax)) # Range in values
        else:
            xmin,xmax = range
    else:                               # Full range
        xmin,xmax = None,None

    xx = N.ravel(x)
    if xmin is None:                    # Automatic xmin = min(x)
        if log:                         # xmin = min(x>0)
            xmin = xx[xx>0].min()       # Might raise ValueError
        else:
            xmin = xx.min()
    if xmax is None:                    # Automatic xmax = max(x)
        xmax = xx.max()

    return xmin,xmax


def hist_binwidth(x, choice='FD', range=None):
    """Optimal histogram binwidth. Choices are:

    - 'S': Scott's choice
    - 'FD': Freedman & Diaconis (1981), fast, fair if single-peaked [default]
    - 'SS': Shimazaki and Shinomoto (2007), slow, best choice if double-peaked
    - 'BR': Birgé and Rozenholc (2006), slow

    Analysis is restricted to *range*=(*min*,*max*) if not `None`
    (full range by default).

    References:

    - `Histogram <http://en.wikipedia.org/wiki/Histogram>`_
    - `Histogram Bin-width Optimization
      <http://176.32.89.45/~hideaki/res/histogram.html>`_
    """

    xx = N.ravel(x)
    xmin,xmax = get_range(xx, range)
    xx = xx[(xx>=xmin) & (xx<=xmax)]

    if choice=='FD':                     # Freedman and Diaconis (1981)
        l,h = N.percentile(xx, [25.,75.])
        h = 2*(h-l)/len(xx)**(1./3.)
    elif choice=='S':                    # Scott's choice
        h = 3.49*N.std(xx, ddof=1)/len(xx)**(1./3.)
    elif choice=='BR':                   # Birgé and Rozenholc (2006)
        def penalty(nbin):
            return nbin - 1 + N.log(nbin)**2.5
        def likelihood(nbin):
            hist,bins = N.histogram(xx,bins=nbin)
            return (hist*N.log(nbin*N.maximum(hist,1)/float(len(xx)))).sum()
        nbins = N.arange(2,round(len(xx)/N.log(len(xx)))+1, dtype='i')
        nbin = nbins[N.argmax([ likelihood(n)-penalty(n) for n in nbins ])]
        h = (xmax-xmin)/nbin
    elif choice=='SS':                   # Shimazaki and Shinomoto (2007)
        # http://web.mit.edu/hshimaza/www//res/histogram.html
        def objf(nbin):
            hist,bins = N.histogram(xx,bins=nbin)
            delta = bins[1]-bins[0]
            k = hist.mean()
            v = hist.var(ddof=0)
            #print "nbin",nbin,delta,k,v,(2*k - v)/delta**2
            return (2*k - v)/delta**2
        nbins = N.arange(2,round(len(xx)/N.log(len(xx)))+1, dtype='i')
        nbin = nbins[N.argmin([ objf(n) for n in nbins ])]
        h = (xmax-xmin)/nbin
    else:
        raise ValueError("Unknow histogram binwidth's choice '%s'" % choice)

    return h


def hist_nbin(x, choice='FD', range=None):
    """Optimal number of bins. See :func:`hist_binwidth` for details."""

    xmin,xmax = get_range(x, range=range)

    return int( N.ceil((xmax-xmin) /
                       hist_binwidth(x, choice=choice, range=range)) )


def hist_bins(x, choice='FD', range=None, log=False):
    """Optimal binning. See :func:`hist_binwidth` for details."""

    xmin,xmax = get_range(x, range=range, log=log)
    if log:
        from math import log10
        lxmin,lxmax = log10(xmin),log10(xmax)
        xx = N.ravel(x)
        xx = xx[xx>=xmin]
        return N.logspace(lxmin,lxmax,
                          hist_nbin(N.log10(xx), choice=choice,
                                    range=(lxmin,lxmax)))
    else:
        return N.linspace(xmin, xmax, hist_nbin(x, choice=choice, range=range))


# ############################################################

if __name__ == '__main__':

    import matplotlib.pyplot as P

    if False:
        # R test-case
        A= [79.98, 80.04, 80.02, 80.04, 80.03, 80.03, 80.04,
            79.97, 80.05, 80.03, 80.02, 80.00, 80.02]
        B = [80.02, 79.94, 79.98, 79.97, 79.97, 80.03, 79.95, 79.97]
    else:
        # AUTO83B.DAT test-case from
        # http://www.itl.nist.gov/div898/handbook/eda/section3/eda353.htm
        A = [18,15,18,16,17,15,14,14,14,15,15,14,15,14,22,18,21,21,10,10,11,
             9,28,25,19,16,17,19,18,14,14,14,14,12,13,13,18,22,19,18,23,26,25,
             20,21,13,14,15,14,17,11,13,12,13,15,13,13,14,22,28,13,14,13,14,15,
             12,13,13,14,13,12,13,18,16,18,18,23,11,12,13,12,18,21,19,21,15,16,
             15,11,20,21,19,15,26,25,16,16,18,16,13,14,14,14,28,19,18,15,15,16,
             15,16,14,17,16,15,18,21,20,13,23,20,23,18,19,25,26,18,16,16,15,22,
             22,24,23,29,25,20,18,19,18,27,13,17,13,13,13,30,26,18,17,16,15,18,
             21,19,19,16,16,16,16,25,26,31,34,36,20,19,20,19,21,20,25,21,19,21,
             21,19,18,19,18,18,18,30,31,23,24,22,20,22,20,21,17,18,17,18,17,16,
             19,19,36,27,23,24,34,35,28,29,27,34,32,28,26,24,19,28,24,27,27,26,
             24,30,39,35,34,30,22,27,20,18,28,27,34,31,29,27,24,23,38,36,25,38,
             26,22,36,27,27,32,28,31]
        B = [24,27,27,25,31,35,24,19,28,23,27,20,22,18,20,31,32,31,32,24,26,29,
             24,24,33,33,32,28,19,32,34,26,30,22,22,33,39,36,28,27,21,24,30,34,
             32,38,37,30,31,37,32,47,41,45,34,33,24,32,39,35,32,37,38,34,34,32,
             33,32,25,24,37,31,36,36,34,38,32,38,32]

    alpha=0.95                              # Confidence level

    # Welsh t-test, unequal variance
    t,df,pt,tint = t_test(A,B, cl=alpha, verbose=True)
    # F-test
    F,(dfnum,dfden),pf,fint = F_test(A,B, cl=alpha, verbose=True)
    # t-test, equal variance
    t2,df2,pt2,tint2 = t_test(A,B, cl=alpha, equalVar=True, verbose=True)

    fig = P.figure()
    fig.subplots_adjust(hspace=0.3)

    ax1 = fig.add_subplot(2,1,1,
                          xlabel='t',
                          title="t-test (df(unequal)=%.2f, df(equal)=%.2f)" % \
                              (df,df2))
    xt = N.linspace(-4,4,101)
    ax1.plot(xt, S.t.pdf(xt,df), 'b-', label='pdf')
    ax1.plot(xt, S.t.cdf(xt,df), 'g-', label='cdf')
    ax1.plot(xt, S.t.sf(xt,df), 'r-', label='sf=1-cdf')
    ax1.axvline(t, color='k', label='t(unequal)')
    ax1.axvspan(tint[0],tint[1],fc='0.8', label='%.0f%%-CI' % (alpha*100),
                alpha=0.5)
    ax1.plot(xt, S.t.pdf(xt,df2), 'b--', label='_')
    ax1.plot(xt, S.t.cdf(xt,df2), 'g--', label='_')
    ax1.plot(xt, S.t.sf(xt,df2), 'r--', label='_')
    ax1.axvline(t2, color='k', ls='--', label='t(equal)')
    ax1.axvspan(tint2[0],tint2[1], label='_',
                linestyle='dashed', fc='0.8', alpha=0.5)
    ax1.legend(loc='best')
    ax1.grid(True)

    ax2 = fig.add_subplot(2,1,2,
                          xlabel='F',
                          title="F-test (df=%d/%d)" % (dfnum,dfden))
    xf = N.linspace(0,4,100)
    ax2.plot(xf, S.f.pdf(xf,dfnum,dfden), label='pdf')
    ax2.plot(xf, S.f.cdf(xf,dfnum,dfden), label='cdf')
    ax2.plot(xf, S.f.sf(xf,dfnum,dfden), label='sf=1-cdf')
    ax2.axvline(F, color='k', label='_')
    ax2.axvspan(fint[0],fint[1], label='%.0f%%-CI' % (alpha*100),
                fc='0.8', alpha=0.5)
    ax2.legend(loc='best')
    ax2.grid(True)

    P.show()
