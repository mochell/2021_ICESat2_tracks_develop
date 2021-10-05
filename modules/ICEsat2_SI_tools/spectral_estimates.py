import numpy as np


# basic functions
def create_chunk_boundaries(L, dsize, ov= None,  iter_flag=True):
    """
    returns all need chunk boudaries and center position given L, and ov
    inputs:
    L desired length of window,
    dsize size of the data

    if ov is None, = L/2
    if iter_flag True returns iter else it returns an ndarray

    """
    ov=int(np.round(L/2)) if ov is None else ov


    xleft = np.arange(0,dsize-int(L-ov),int(L-ov))
    xright = np.arange(int(L-ov)*2,dsize+1,int(L-ov))
    xcenter_pos = np.arange(int(L-ov),dsize-int(L-ov)+1,int(L-ov))
    max_size = min([xleft.size , xcenter_pos.size, xright.size])
    # if xright[max_size-1] < dsize:
    #     print('left out last ' + str(dsize- xright[max_size-1]) + ' data points' )
    #print([xleft[0:max_size], xcenter_pos[0:max_size], xright[0:max_size]])
    position_stancil = np.vstack([xleft[0:max_size], xcenter_pos[0:max_size], xright[0:max_size]])

    if iter_flag is True:
        return iter(position_stancil.T.tolist())
    else:
        return position_stancil


def create_chunk_boundaries_unit_lengths(L_unit, data_limits, ov= None,  iter_flag=True):
    """
    returns all need chunk boudaries and center position given L, and ov
    inputs:
    L desired length of window in units of the x axis of the data,
    data_limits (x_min, x_max) tuple with the beginning and end the the derived window stancils

    if ov is None, = L/2
    if iter_flag True returns iter else it returns an ndarray

    """
    L= L_unit
    ov=np.round(L/2) if ov is None else ov
    #print(ov)
    dl = (L-ov)
    xleft   = np.arange(data_limits[0]           ,   data_limits[1]-dl,   dl )
    xcenter_pos = np.arange(data_limits[0]+ L/2   ,   data_limits[1]-dl+1, dl )
    xright  = np.arange(data_limits[0] + L    ,   data_limits[1]+1,           dl )


    max_size = min([xleft.size , xcenter_pos.size, xright.size])
    # if xright[max_size-1] < data_limits[1]:
    #     print('left out last ' + str(data_limits[1]- xright[max_size-1]) + ' data points' )
    #print([xleft[0:max_size], xcenter_pos[0:max_size], xright[0:max_size]])
    position_stancil = np.vstack([xleft[0:max_size], xcenter_pos[0:max_size], xright[0:max_size]])

    if iter_flag is True:
        return iter(position_stancil.T.tolist())
    else:
        return position_stancil

# 2nd cal spectra
def calc_spectrum_fft(phi, df, N):
    """ compute the 1d spectrum of a field phi
    inputs:

    df      frequency / or wavenumber step
    N       length of data (= L)
    neven   bool, if True
    """
    neven = True if (N%2) else False

    phih = np.fft.rfft(phi)

    # the factor of 2 comes from the symmetry of the Fourier coeffs
    spec = 2.*(phih*phih.conj()).real / df /N**2

    # the zeroth frequency should be counted only once
    spec[0] = spec[0]/2.
    if neven:
        spec[-1] = spec[-1]/2.

    return spec


def LS_power_to_PSD( ls_power, L , dff):
    """
    returns Power spectral density (unit^2/dfreq)
    ls_power    output of astropy.timeseries.LombScargle.power with normalization='psd'
    """
    return 2 * ls_power / L /dff

def calc_spectrum_LS( x, y, k,  LS= None, dk =None):
    """
    returns Power spectral density of y given postitions x, for wanumbers k
    LS is a Lomb-scargel object. if None its initlized again

    """
    if LS is None:
        from astropy.timeseries import LombScargle
        LS = LombScargle(x , y, fit_mean=True)
    else:
        LS.t = x
        LS.y = y
    ls_power    = LS.power(k, normalization='psd', assume_regular_frequency='True')

    dk          = np.diff(k).mean() if dk is None else dk

    return  2 * ls_power / y.size / dk

def calc_freq_fft(x_grid, N):
    """ calculate array of spectral variable (frequency or
            wavenumber) in cycles per unit of L """

    neven = True if (N%2) else False
    dx=np.diff(x_grid).mean()
    df = 1./((N-1)*dx)

    if neven:
        f = df*np.arange(N/2+1)
    else:
        f = df*np.arange( (N-1)/2.  + 1 )
    return f,df

def calc_freq_fft_given_dx(dx, N):
    """
    calculate array of spectral variable (frequency or
            wavenumber) in cycles per unit of L
    dx      given resolution of the data
    N       number of datapoints used in window
    """

    neven = True if (N%2) else False
    df = 1./((N-1)*dx)

    if neven:
        f = df*np.arange(N/2+1)
    else:
        f = df*np.arange( (N-1)/2.  + 1 )
    return f,df


def calc_freq_LS(x, N, method='fftX2', dx=None, minimum_frequency=None, maximum_frequency=None, samples_per_peak=0.01):
    """
    calculate array of spectral variable (frequency or
    wavenumber) in cycles per unit of N (window length in number of data points)
    x can be unevenly spaced
    method:
        "fftX2"     defined the frequencyt grid as for FFT, but with double its resolution
        "LS_auto"   using LS algorithm with samples_per_peak=0.1

        minimum_frequency, maximum_frequency only used for LS_auto
    """

    if method is 'fftX2':
        neven = True if (N%2) else False
        dx = np.diff(x).mean() if dx is None else dx
        df = 1./((N-1)*dx) /2
        if neven:
            f = df*np.arange(df, N+1)
        else:
            f = df* np.arange(df, (N-1)  + 1 )

    elif method is 'fft':
        neven = True if (N%2) else False
        dx = np.diff(x).mean() if dx is None else dx
        df = 1./((N-1)*dx)
        if neven:
            f = df*np.arange(N/2+1)
        else:
            f = df*np.arange( (N-1)/2.  + 1 )

    elif method is 'LS_auto':
        from astropy.timeseries import LombScargle
        f = LombScargle(x , np.random.randn(len(x)), fit_mean=True).autofrequency(minimum_frequency=minimum_frequency, maximum_frequency=maximum_frequency, samples_per_peak=samples_per_peak)##0.1)

        df = np.diff(f).mean()
        df = np.round(df, 5)

    elif method is 'fixed_ratio':

        neven = True if (N%2) else False
        dx = np.diff(x).mean() if dx is None else dx
        df = dx / 50
        if neven:
            f = df*np.arange(df, N +1)
        else:
            f = df* np.arange(df, N )

    return f ,df

def create_window(L, window=None):
    """
    define window function and weight it to conserve variance
    if window is not None it show have a length of N
    """
    if window is None:
        win=np.hanning(L)
    else:
        win=window

    factor=np.sqrt(L/(win**2).sum())
    win*=factor
    return win

def spec_error(E,sn,ci=.95):

    """ Computes confidence interval for one-dimensional spectral
        estimate E (the power spectra).

        Parameters
        ===========
        - sn is the number of spectral realizations;
                it can be either an scalar or an array of size(E)
        - ci = .95 for 95 % confidence interval

        Output
        ==========
        lower (El) and upper (Eu) bounds on E """



    def yNlu(sn,yN,ci):
        """ compute yN[l] yN[u], that is, the lower and
                    upper limit of yN """

        from scipy.special import gammainc
        # cdf of chi^2 dist. with 2*sn DOF
        cdf = gammainc(sn,sn*yN)

        # indices that delimit the wedge of the conf. interval
        fl = np.abs(cdf - ci).argmin()
        fu = np.abs(cdf - 1. + ci).argmin()

        return yN[fl],yN[fu]

    dbin = .005
    yN = np.arange(0,2.+dbin,dbin)

    El, Eu = np.empty_like(E), np.empty_like(E)

    try:
        n = sn.size
    except AttributeError:
        n = 0

    if n:

        assert n == E.size, " *** sn has different size than E "

        for i in range(n):
            yNl,yNu = yNlu(sn[i],yN=yN,ci=ci)
            El[i] = E[i]/yNl
            Eu[i] = E[i]/yNu

    else:
        yNl,yNu = yNlu(sn,yN=yN,ci=ci)
        El = E/yNl
        Eu = E/yNu

    return El, Eu




class wavenumber_spectrogram(object):
    def __init__(self, x_grid, data, L, ov=None, window=None):
        """
        returns a wavenumber spectrogram with the resolution L-ov
        this uses standard fft and assumes equally gridded data

        inputs:
        data    data
        x       grid the fft is taken on
        L       window length in number of grid points
        ov      (default=None) number of grid points the windows should overlab
        window  (default=np.hanning) numpy window function

        returns:
        xr.Dataset with x, k as cooridates of the spectrogram and the mean error
            other arributes are in the .attr dict.
        """

        self.L      = L
        self.ov = int(L/2) if ov is None else ov #when not defined in create_chunk_boundaries then L/2

        self.data   = data

        # create subsample k
        self.k, self.dk  = calc_freq_fft(x_grid, L)
        # create window
        self.win    = create_window(L)

    def cal_spectrogram(self, data=None, name=None):

        """
        defines apply function and calculated all sub-sample sprectra using map
        """
        import xarray as xr

        DATA = self.data if data is None else data
        L, dk = self.L, self.dk
        win =self.win

        def calc_spectrum_apply(stancil):

            "returns spectrum per stencil, detrends and windows the data"
            from scipy.signal import detrend

            idata = DATA[stancil[0]:stancil[-1]]
            idata = detrend(idata) * win

            return stancil[1], calc_spectrum_fft(idata , dk, L)

        # def test_func(i_stancil):
        #     return i_stancil[1], yy[i_stancil[0]:i_stancil[-1]].shape

        # %% derive L2 stancil
        stancil_iter = create_chunk_boundaries(L, DATA.size, ov= self.ov)
        # apply func to all stancils
        D_specs = dict(map(calc_spectrum_apply,stancil_iter))

        chunk_positions = np.array(list(D_specs.keys()))
        self.N_stancils = len(chunk_positions) # number of spectal relazations

        # repack data, create xarray
        self.spec_name = 'power_spec' if name is None else name
        G =dict()
        for xi,I in D_specs.items():
            G[xi] = xr.DataArray(I,  dims=['k'], coords={'k': self.k, 'x': xi } , name=self.spec_name)

        self.G = xr.concat(G.values(), dim='x').T#.to_dataset()
        if self.G.k[0] == 0:
            self.G = self.G[1:, :]

        self.G.attrs['ov'] = self.ov
        self.G.attrs['L'] = self.L

        return self.G

    # cal variance
    def calc_var(self):
        """ Compute total variance from spectragram """
        return self.dk*self.G.mean('x').sum().data  # do not consider zeroth frequency

    def mean_spectral_error(self, confidence = 0.95):
        "retrurns spetral error for the x-mean spectral estimate and stores it as coordindate in the dataarray"
        #  make error estimate
        El_of_mean, Eu_of_mean = spec_error(self.G.mean('x'), self.N_stancils, confidence )
        El_of_mean.name = 'El_mean'
        Eu_of_mean.name = 'Eu_mean'

        self.G.coords['mean_El'] = (('k'), El_of_mean)
        self.G.coords['mean_Eu'] = (('k'), Eu_of_mean)

    def parceval(self, add_attrs=True ):
        "test Parceval theorem"
        DATA = self.data
        L = self.L


        # derive mean variances of stancils
        stancil_iter = create_chunk_boundaries(L, DATA.size)

        def get_stancil_var_apply(stancil):
            from scipy.signal import detrend
            "returns the variance of yy for stancil"
            idata = DATA[stancil[0]:stancil[-1]]
            idata = detrend(idata)# * win
            return stancil[1], idata.var()

        D_vars = dict(map(get_stancil_var_apply,stancil_iter))

        stancil_vars =list()
        for I in D_vars.values():
            stancil_vars.append(I)

        print('Parcevals Theorem:')
        print('variance of unweighted timeseries: ',DATA.var())
        print('mean variance of detrended chunks: ', np.array(stancil_vars).mean())
        #print('variance of weighted timeseries: ',self.phi.var() )
        #self.calc_var(self)
        print('variance of the pwelch Spectrum: ', self.calc_var())

        if add_attrs:
            self.G.attrs['variance_unweighted_data']        = DATA.var()
            self.G.attrs['mean_variance_detrended_chunks']  = np.array(stancil_vars).mean()
            self.G.attrs['mean_variance_pwelch_spectrum']   = self.calc_var()



class wavenumber_spectrogram_LS_even(object):
    def __init__(self, x, data, L, waven_method = 'fftX2' , dy=None ,  ov=None, window=None, kjumps=1):
        """
        returns a wavenumber spectrogram with the resolution L-ov
        this uses Lombscargle

        inputs:
        data    data
        x       grid the fft is taken on
        dy      passed to LS object "error or sequence of observational errors associated with times t"
        waven_method     ('auto' (default), or )
        L       window length in number of grid points
        ov      (default=None) number of grid points the windows should overlab
        window  (default=np.hanning) numpy window function

        returns:
        xr.Dataset with x, k as cooridates of the spectrogram and the mean error
            other arributes are in the .attr dict.
        """
        from astropy.timeseries import LombScargle
        self.L      = L
        self.ov = int(L/2) if ov is None else ov #when not defined in create_chunk_boundaries then L/2

        self.x      = x
        self.data   = data
        self.dy     = dy


        # create subsample k
        #print(waven_method)
        if type(waven_method) is str:
            self.k, self.dk  = calc_freq_LS(x, L, method = waven_method )
        elif type(waven_method) is np.ndarray:
            self.k, self.dk  = waven_method, np.diff(waven_method).mean()
        else:
            raise ValueError('waven_method is neither string nor an array')

        self.k, self.dk  = self.k[::kjumps], self.dk*kjumps
        # create window
        self.win    = None #create_window(L)

    def cal_spectrogram(self, x = None, data=None, name=None, dx=1):

        """
        defines apply function and calculated all sub-sample sprectra using map
        dx      nominal resolution of the data resolutionif not set, dx= 1
        """
        from astropy.timeseries import LombScargle
        import xarray as xr

        X       = self.x if x is None else x # all x positions
        DATA    = self.data if data is None else data # all data points
        L, dk   = self.L, self.dk
        win     = self.win
        self.dx = dx
        # init Lomb scargle object with noise as nummy data ()
        #dy_fake= np.random.randn(len(dy))*0.001 if self.dy is not None else None
        self.LS = LombScargle(X[0:L] , np.random.randn(L)*0.001, fit_mean=True)


        def calc_spectrum_apply(stancil):

            "returns spectrum per stencil, detrends and windows the data"
            from scipy.signal import detrend

            x = X[stancil[0]:stancil[-1]]
            #x_mask= (stancil[0] < X) & (X <= stancil[-1])
            #x = X[x_mask]
            idata = DATA[stancil[0]:stancil[-1]]
            y = detrend(idata)# * win

            return stancil[1], calc_spectrum_LS( x, y, self.k,  LS= self.LS, dk =self.dk)

        # % derive L2 stancil
        stancil_iter = create_chunk_boundaries(L, DATA.size, ov= self.ov)
        # apply func to all stancils
        D_specs = dict(map(calc_spectrum_apply,stancil_iter))

        chunk_positions = np.array(list(D_specs.keys()))
        self.N_stancils = len(chunk_positions) # number of spectal relazations

        # repack data, create xarray
        self.spec_name = 'power_spec' if name is None else name
        G =dict()
        for xi,I in D_specs.items():
            G[xi] = xr.DataArray(I,  dims=['k'], coords={'k': self.k, 'x': xi * self.dx } , name=self.spec_name)

        self.G = xr.concat(G.values(), dim='x').T#.to_dataset()
        if self.G.k[0] == 0:
            self.G = self.G[1:, :]

        self.G.attrs['ov'] = self.ov
        self.G.attrs['L'] = self.L

        return self.G

    def calc_var(self):
        return wavenumber_spectrogram.calc_var(self)

    def parceval(self, add_attrs=True ):
        return wavenumber_spectrogram.parceval(self, add_attrs= add_attrs )

    def mean_spectral_error(self, confidence = 0.95 ):
        return wavenumber_spectrogram.mean_spectral_error(self, confidence= confidence )


class wavenumber_spectrogram_LS(object):
    def __init__(self, x, data, L, waven_method = 'fftX2' ,  ov=None, window=None, kjumps=1):
        """
        returns a wavenumber spectrogram with the resolution L-ov
        this uses Lombscargle

        inputs:
        data    data
        x       x-positions of where the data is taken
        dy      passed to LS object "error or sequence of observational errors associated with times t"
        waven_method     ('auto' (default), or something else if "waven_method" is an array these wavenumbers are used )

        L       window length in units of x
        ov      (default=None) number of grid points the windows should overlab
        window  (default=np.hanning) numpy window function

        returns:
        xr.Dataset with x, k as cooridates of the spectrogram and the mean error
            other arributes are in the .attr dict.
        """
        from astropy.timeseries import LombScargle
        self.L      = L
        self.ov = int(L/2) if ov is None else ov #when not defined in create_chunk_boundaries then L/2

        self.x      = x
        self.dx     = np.diff(x).mean()
        self.data   = data
        self.Lpoints= int(self.L/self.dx)


        # create subsample k
        #print(waven_method)
        if type(waven_method) is str:
            self.k, self.dk  = calc_freq_LS(x, self.Lpoints, method = waven_method )
        elif type(waven_method) is np.ndarray:
            self.k, self.dk  = waven_method, np.diff(waven_method).mean()
        else:
            raise ValueError('waven_method is neither string nor an array')

        self.k, self.dk  = self.k[::kjumps], self.dk*kjumps
        # create window
        self.win    = None #create_window(L)

    def cal_spectrogram(self, x = None, data=None, name=None, dx=1, xlims =None):

        """
        defines apply function and calculated all sub-sample sprectra using map
        """
        from astropy.timeseries import LombScargle
        import xarray as xr
        import copy

        X       = self.x if x is None else x # all x positions
        DATA    = self.data if data is None else data # all data points
        L, dk   = self.L, self.dk
        win     = self.win
        self.dx = dx
        self.xlims = ( np.round(X.min()), X.max() ) if xlims is None else xlims

        # init Lomb scargle object with noise as nummy data ()
        #dy_fake= np.random.randn(len(dy))*0.001 if self.dy is not None else None
        #self.LS = LombScargle(X[0:L] , np.random.randn(L)*0.001, fit_mean=True)


        def calc_spectrum_apply(stancil):

            """
            windows the data accoding to stencil and applies LS spectrogram
            returns: stancil center, spectrum for this stencil, number of datapoints in stancil
            """
            from scipy.signal import detrend

            #x = X[stancil[0]:stancil[-1]]
            x_mask= (stancil[0] < X) & (X <= stancil[-1])
            x = X[x_mask]
            if x.size < 100: # if there are not enough photos set results to nan
                return stancil[1], self.k*np.nan, x.size

            y = DATA[x_mask]

            #print(x.shape, y.shape, self.k,  self.LS)
            return stancil[1], calc_spectrum_LS( x, y, self.k,  LS= None, dk =self.dk), x.size

        # % derive L2 stancil
        self.stancil_iter = create_chunk_boundaries_unit_lengths(L, self.xlims, ov= self.ov, iter_flag=True)
        #stancil_iter = create_chunk_boundaries_unit_lengths(L, ( np.round(X.min()), X.max() ), ov= self.ov, iter_flag=True)

        # apply func to all stancils
        # Spec_returns=list()
        # for ss in stancil_iter:
        #     print(ss)
        #     Spec_returns.append( calc_spectrum_apply(ss) )

        Spec_returns = list(map( calc_spectrum_apply, copy.copy(self.stancil_iter)   ))

        # unpack resutls of the mapping:
        D_specs = dict()
        N_per_stancil = list()
        for I in Spec_returns:
            D_specs[I[0]] = I[1]
            N_per_stancil.append(I[2])

        self.N_per_stancil = N_per_stancil
        chunk_positions = np.array(list(D_specs.keys()))
        self.N_stancils    = len(chunk_positions) # number of spectral realizatiobs

        # repack data, create xarray
        self.spec_name = 'power_spec' if name is None else name
        G =dict()
        for xi,I in D_specs.items():
            G[xi] = xr.DataArray(I,  dims=['k'], coords={'k': self.k, 'x': xi * self.dx } , name=self.spec_name)

        self.G = xr.concat(G.values(), dim='x').T#.to_dataset()
        if self.G.k[0] == 0:
            self.G = self.G[1:, :]

        self.G.attrs['ov'] = self.ov
        self.G.attrs['L'] = self.L
        self.G.coords['N_per_stancil'] = ( ('x'), N_per_stancil)

        return self.G

    def calc_var(self):

        Gmean = np.nanmean(self.G, 1)
        infmask = np.isinf(Gmean)

        return self.dk * Gmean[~infmask].sum().data

    # def parceval(self, add_attrs=True ):
    #     return wavenumber_spectrogram.parceval(self, add_attrs= add_attrs )

    def parceval(self, add_attrs=True ):
        "test Parceval theorem"
        import copy
        DATA = self.data
        L = self.L
        X = self.x

        # derive mean variances of stancils
        #stancil_iter = create_chunk_boundaries_unit_lengths(L, self.xlims, ov= self.ov )

        def get_stancil_var_apply(stancil):
            from scipy.signal import detrend
            "returns the variance of yy for stancil"
            x_mask= (stancil[0] < X) & (X <= stancil[-1])
            idata = DATA[x_mask]
            if len(idata) < 1:
                return stancil[1], np.nan, len(idata)
            idata = detrend(idata)# * win
            return stancil[1], idata.var(), len(idata)

        D_vars = list(map(get_stancil_var_apply, copy.copy(self.stancil_iter)  ))

        stancil_vars, Nphotons =list(), 0
        for I in D_vars:
            stancil_vars.append(I[1]  * I[2])
            Nphotons  += I[2]

        stancil_weighted_variance = np.nansum(np.array(stancil_vars))/Nphotons

        print('Parcevals Theorem:')
        print('variance of timeseries: ', DATA.var())
        print('mean variance of stancils: ', stancil_weighted_variance )
        #print('variance of weighted timeseries: ',self.phi.var() )
        #self.calc_var(self)
        print('variance of the pwelch LS Spectrum: ', self.calc_var())

        if add_attrs:
            self.G.attrs['variance_unweighted_data']        = DATA.var()
            self.G.attrs['mean_variance_stancils']  = np.nanmean(np.array(stancil_vars) )
            self.G.attrs['mean_variance_LS_pwelch_spectrum']   = self.calc_var()


    def mean_spectral_error(self, confidence = 0.95 ):
        return wavenumber_spectrogram.mean_spectral_error(self, confidence= confidence )



# class for getting standard Pwelch spectrum. old version, deprechiate
class wavenumber_pwelch(object):
    def __init__(self,data, x, L, ov=None, window=None, save_chunks=False, plot_chunks=False):
        """
        returns a wavenumber spectrum using the pwelch method

        inputs:
        data    data
        x       grid the fft is taken on
        L       window length in number of grid points
        ov      (default=None) number of grid points the windows should overlab
        window  (default=np.hanning) numpy window function
        save_chunks     if True, self.chunks contains all compute chunks
        plot_chunks     if True, it plots all chunks

        returns:
        self.spec_est   mean power spectrum
        self.n_spec
        self.n
        self.dx
        self.n_spec
        """
        from scipy import signal


        self.data       =   data      # field to be analyzed
        self.dx         =   np.diff(x)[0] # sampling interval
        self.save_chunks=   save_chunks
        dsize           =   data.size

        ov=int(np.round(L/2)) if ov is None else ov

        self.n = L
        if window is None:
            win=np.hanning(self.n)
        else:
            win=window

        factor=np.sqrt(self.n/(win**2).sum())
        win*=factor

        # test if n is even
        if (self.n%2):
            self.neven = False
        else:
            self.neven = True
        #calculate freq
        self.k  = self.calc_freq()
        #del(self.df)

        #print(data.size, L, ov, int(L-ov) )
        nbin=int(np.floor(dsize/(L-ov)))
        #print(nbin)

        if save_chunks:
            chunks=np.empty([int(nbin),int(L)])

        self.specs=np.empty([int(nbin),self.k.size])
        #print(chunks.shape)
        #result_array = np.empty((0, 100))
        #if plot_chunks:
            #M.figure_axis_xy()
        last_k=0
        k=0
        #print('iter range', np.arange(0,data.size,int(L-ov)))
        for i in np.arange(0,dsize-int(L-ov)+1,int(L-ov)):

            if (plot_chunks) and (i >= dsize-6*int(L-ov)):
                M.figure_axis_xy()

            self.phi=data[int(i):int(i+L)]

            #self.ii=np.append(self.ii,[i,i+L])
            #print(self.phi.max())

            #print(self.phi.mean())
            #print(self.phi.shape)
            #print('i',int(i), int(i+L))
            #print(chunk.size, l)
            if int(i+L) <= data.size-1:
                if save_chunks:
                    chunks[k,:]=self.phi


                self.phi=signal.detrend(self.phi)*win
                if plot_chunks:
                    #MT.stats_format(self.phi, 'chunk '+str(i))
                    plt.plot(self.phi)

                self.specs[k,:]= self.calc_spectrum()
                last_k=k
                last_used_TS=int(i+L)
                #if plot_chunks:
                #    MT.stats_format(self.spec, 'spec '+str(i))

            else:
                if plot_chunks:
                    print('end of TS is reached')
                    print('last spec No: '+str(last_k))
                    print('spec container: '+str(specs.shape))
                    print('last used Timestep: '+str(last_used_TS))
                    print('length of TS '+ str(dsize) +'ms')

            k+=1


        if save_chunks:
            self.chunks=chunks
            #del(chunks)

        self.spec_est=self.specs.mean(axis=0)
        # if prewhite is None:
        #     self.specs=specs[:last_k,:]
        #     self.spec_est=self.specs.mean(axis=0)
        # elif prewhite ==1:
        #     self.specs=specs[:last_k,:]*(2*np.pi*self.f)
        #     self.spec_est=self.specs.mean(axis=0)
        # elif prewhite ==2:
        #     self.specs=specs[:last_k,:]*(2*np.pi*self.f)**2
        #     self.spec_est=self.specs.mean(axis=0)


        self.n_spec,_=self.specs.shape
        self.calc_var()
        #self.phi=self.data
        #self.phi*=win*np.sqrt(factor)

    def calc_freq(self):
        """ calculate array of spectral variable (frequency or
                wavenumber) in cycles per unit of L """

        self.df = 1./((self.n-1)*self.dx)

        if self.neven:
            f = self.df*np.arange(self.n/2+1)
        else:
            f = self.df*np.arange( (self.n-1)/2.  + 1 )
        return f

    def calc_spectrum(self):
        """ compute the 1d spectrum of a field phi """

        self.phih = np.fft.rfft(self.phi)

        # the factor of 2 comes from the symmetry of the Fourier coeffs
        spec = 2.*(self.phih*self.phih.conj()).real / self.df /self.n**2

        # the zeroth frequency should be counted only once
        spec[0] = spec[0]/2.
        if self.neven:
            spec[-1] = spec[-1]/2.

        return spec


    def error(self, ci=0.95):
        self.El, self.Eu =spec_error(self.spec_est,self.n_spec,ci=ci)

    def parceval(self):
        print('Parcevals Theorem:')
        print('variance of unweighted timeseries: ',self.data.var())
        print('mean variance of timeseries chunks: ',self.chunks.var(axis=1).mean() if self.save_chunks is True else 'data not saved')
        #print('variance of weighted timeseries: ',self.phi.var() )
        #self.calc_var(self)
        print('variance of the pwelch Spectrum: ',self.var)

    def calc_var(self):
        """ Compute total variance from spectrum """
        self.var = self.df*self.specs[1:].mean(axis=0).sum()  # do not consider zeroth frequency
