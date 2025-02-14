
import numpy as np
from numpy import pi
from scipy.special import gammainc
from scipy import signal
import matplotlib.pyplot as plt
try:
    import mkl
    np.use_fastnumpy = True
except ImportError:
    pass

def calc_freq(self):
    """ calculate array of spectral variable (frequency or
            wavenumber) in cycles per unit of L """

    self.df = 1./((self.n-1)*self.dt)

    if self.neven:
        self.f = self.df*np.arange(self.n/2+1)
    else:
        self.f = self.df*np.arange( (self.n-1)/2.  + 1 )

def calc_spectrum(self):
    """ compute the 1d spectrum of a field phi """

    self.phih = np.fft.rfft(self.phi)

    # the factor of 2 comes from the symmetry of the Fourier coeffs
    self.spec = 2.*(self.phih*self.phih.conj()).real / self.df /self.n**2

    # the zeroth frequency should be counted only once
    self.spec[0] = self.spec[0]/2.
    if self.neven:
        self.spec[-1] = self.spec[-1]/2.

def calc_var(self):
        """ Compute total variance from spectrum """
        self.var = self.df*self.spec[1:].sum()  # do not consider zeroth frequency

def create_timeaxis_collection(time_as_datetime64):
    import datetime as DT
    from matplotlib import dates
    T=dict()
    T['sec']=      time_as_datetime64.astype('M8[s]').astype('float')
    T['day']=      time_as_datetime64.astype('M8[D]').astype('float')
    T['datetime']= time_as_datetime64.astype(DT.datetime)

    T['float_plot']=dates.date2num(T['datetime'])
    T['dt64']=      time_as_datetime64
    return T

def spicke_remover(data, nstd=20.0, spreed=500.0, max_loops=10.0 , verbose=False):
    """
    This function removes spickes from timeseries based on its std and maximum values
    inputs:
    data     timeseries as1d arrays
    nstd     number of standard deviations that have to be exceeded by the maximum data value
    spreed   width of data between wich will be linear interpolated around the spike. width is in units of data points (dt)
    max_loops number of maximum possible loobs until the spike removal is stopped, even it the first creteria is not true
    """
    datastd=np.nanstd(data)
    data2=np.copy(data)

    peak_remove=True
    looper_count=0
    act_flag=False
    while peak_remove is True:
        if  nstd* data.std() < np.max(np.abs(data2)):
            act_flag=True
            if verbose:
                print('true: '+ str(nstd* datastd) +' < '+str( np.max(np.abs(data)) ) )
            data2=M.spickes_to_mean(data2, nloop=0, spreed=spreed, gaussian=False)
            looper_count+=1
        else:
            if verbose:
                print('False: '+ str(nstd* datastd) +' > '+str( np.max(np.abs(data)) ) )
            peak_remove=False

        if looper_count > max_loops:
            peak_remove=False
            if verbose:
                print('stoped by max#')


    if verbose:
        plt.plot(data, 'r')
        plt.plot(data2, 'b')

    return data2 , act_flag

class Spectrum(object):
    """ A class that represents a single realization of
            the one-dimensional spectrum  of a given field phi """

    def __init__(self,data,dt, win_flag=1,pwelch=False, verbose=False):

        import m_tools_ph3 as MT
        self.data = data     # field to be analyzed
        self.dt = dt        # sampling interval

        self.n = data.size

        self.hist='Spectorgram'
        win =  np.hanning(self.n)
        factor=np.sqrt(self.n/(win**2).sum())

        if win_flag:
            if verbose:
                print('window')
            MT.write_log(self.hist, 'window')
            self.win_flag=win_flag
            self.phi=np.copy(self.data[:])
            self.phi*=win*np.sqrt(factor)
        else:
            if verbose:
                print('no window')
            MT.write_log(self.hist, 'no window')
            self.win_flag=0
            self.phi =np.copy(self.data[:])
        # test if n is even
        if (self.n%2):
            self.neven = False
        else:
            self.neven = True

        # calculate frequencies
        calc_freq(self)

        # calculate spectrum
        calc_spectrum(self)

        # calculate total var
        calc_var(self)

#    def calc_var(self):
#        """ Compute total variance from spectrum """
#        self.var = self.df*self.spec[1:].sum()  # do not consider zeroth frequency

    def parceval(self):

        print('Parcevals Theorem:')
        print('variance of unweighted timeseries: ',self.data.var())
        print('variance of weighted timeseries: ',self.phi.var() if self.win_flag is 1 else 'data not windowed')
        print('variance of weighted timeseries: ',self.phi.var() )
        print('variance of the Spectrum: ',self.var)



class moments(object):
    def __init__(self,data_org,dt, L=None, ov=None,window=None, save_chunks=False, plot_chunks=False, prewhite=None):
        """
        This function calculates the spectral moments from a station (buoy, GPS, or seismic station) that measures
        displacement in 3 direction in a right hand coordinate system. The directions should have the same untis.

        data        list of data in the Directions [Z, E, N]
        prewhite    None(default)
                    1 = 1 timederivative
                    2=  2 timederivatives
        commment: the 0 freq is left out.
        """
        self.data = np.array(data_org)
        if prewhite is None:
            data = np.array(data_org)    # field to be analyzed
        elif prewhite ==1:
            print('prewhite =1')
            data=np.gradient(np.array(data_org)  , axis=1)
        elif prewhite ==2:
            data=np.gradient(np.gradient(np.array(data_org), axis=1), axis=1)

        #self.data = data      # field to be analyzed
        self.dt = dt        # sampling interval
        self.save_chunks=save_chunks

        data_size=np.shape(data)[1]
        data_dim=np.shape(data)[0]

        L=int(np.round(data_size/10)) if L is None else L
        if type(L) != int:
            M.echo_dt(L)
            L=L.item().total_seconds()
        #else:
            #print('unknown L type')

        ov=int(np.round(L/2)) if ov is None else ov

        self.n = L
        if window is None:
            win=np.hanning(self.n)
        else:
            win=window

        #weigthing factor for preserved spectra
        factor=np.sqrt(self.n/(win**2).sum())
        win*=factor

        # test if n is even
        if (self.n%2):
            self.neven = False
        else:
            self.neven = True
        #calculate freq
        calc_freq(self)

        # exclude 0 freq.
        self.f=self.f[1:]
        #print('fsize', self.f.size)
        #print(data.size, L, ov, int(L-ov) )
        nbin=int(np.floor(data_size/(L-ov)))-1
        self.nbin=nbin
        if save_chunks:
            chunks=np.empty([data_dim, int(nbin),int(L)])

        #container for spectra
        specs=np.empty([data_dim, int(nbin),self.f.size])
        self.mom_list=list()
        k=0
        #print('iter range', np.arange(0,data.size,int(L-ov)))
        for i in np.arange(0,data_size-int(L-ov)+1,int(L-ov)):

            if (plot_chunks) and (i >= data_size-6*int(L-ov)):
                M.figure_axis_xy()

            self.phi=data[:,int(i):int(i+L)]

            #self.ii=np.append(self.ii,[i,i+L])
            #print(self.phi.max())

            #print(self.phi.mean())
            #print(self.phi.shape)
            #print('i',int(i), int(i+L))
            #print(chunk.size, l)
            if int(i+L) <= data_size-1:
                if save_chunks:
                    chunks[:, k,:]=self.phi

                #normalizing
                #print('normalizing')
                #self.phi=signal.detrend(self.phi,axis=1)*win
                #self.phi=(self.phi.T/np.std(self.phi, axis=1)).T

                #non normalizing
                #de trending and prewhitening
                self.phi=signal.detrend(self.phi,axis=1)*win

                if plot_chunks:
                    #MT.stats_format(self.phi, 'chunk '+str(i))
                    plt.plot(self.phi)

                X=self.calc_spectral_density()
                #print(np.shape(X))
                self.mom=self.calc_moments(X)

                self.mom_list.append(self.mom)

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
                    print('length of TS '+ str(data.size))

            k+=1

        #print('last k', last_k)
        #print('len mom_list', len(self.mom_list)  )
        #print('nbin', nbin)

        if save_chunks:
            self.chunks=chunks
            #del(chunks)

        self.moments_stack=dict()
        self.moments_est=dict()
        for k in self.mom.keys():
            #print(k)
            stack=np.empty([nbin,self.f.size])
            for i in range(len(self.mom_list)):
                #stack=np.vstack((stack,I[k])) #stack them
                #print(i, k)
                stack[i,:]=self.mom_list[i][k] #stack them
                #print(stack.shape)
            #mom_est[k]=stack.mean(axis=0)
            #mean them and decide prewhitheing or nots
            if prewhite is None:
                self.moments_stack[k]=stack
                self.moments_est[k]=np.nanmean(stack,axis=0)
            elif prewhite ==1:
                self.moments_stack[k]=stack*(2*np.pi*self.f)
                self.moments_est[k]=np.nanmean(stack, axis=0)*(2*np.pi*self.f)
            elif prewhite ==2:
                self.moments_stack[k]=stack*(2*np.pi*self.f)**2
                self.moments_est[k]=np.nanmean(stack, axis=0)*(2*np.pi*self.f)**2


        self.moments_unit='[data]^2'
        self.n_spec=len(self.mom_list)
        #self.calc_var()
        #self.phi=self.data
        #self.phi*=win*np.sqrt(factor)

    def calc_spectral_density(self):
        """ compute the 1d spectraln density of a field phi """

        # the factor of 2 comes from the symmetry of the Fourier coeffs (rfft does only calculate the f>=0 part.
        # to preserve variance, we multiply by 2)
        phih = 2 * np.fft.rfft(self.phi)/self.n**2
        # the 1/n**2 is debatalbe. Sarahs notes say, 1/n should be enough.

        # the zeroth frequency should be counted only once
        #phih[0] = phih[0]/2.

        if self.neven:
            phih[-1] = phih[-1]/2.

        return phih

    def calc_moments(self, X):
        mom=dict()
        # 1: for killing 0 freq.
        mom['P11']=(X[0,1:]*X[0,1:].conj()).real
        mom['P22']=(X[1,1:]*X[1,1:].conj()).real
        mom['P33']=(X[2,1:]*X[2,1:].conj()).real

        mom['Q12']=(X[0,1:]*X[1,1:].conj()).imag
        mom['Q13']=(X[0,1:]*X[2,1:].conj()).imag
        mom['Q23']=(X[1,1:]*X[2,1:].conj()).imag

        mom['P13']=(X[0,1:]*X[2,1:].conj()).real
        mom['P12']=(X[0,1:]*X[1,1:].conj()).real
        mom['P23']=(X[1,1:]*X[2,1:].conj()).real

        return mom

    def plot_moments(self):
        num=1
        F=plt.figure(figsize=(9,10) )

        for k,I in self.moments_est.iteritems():
            plt.subplot(len(self.moments_est),1,num )
            plt.plot(self.f, I/self.f, c='k')
            num+=1
            plt.xlim(0.001,.1 )
            plt.ylabel(self.moments_unit)
            plt.title(k)


    def cal_MEM(self, theta=None, flim=(0.01, .5)):
        """ wrapping method that calculated Maximum entrophy method for self.moments_est

        RETURNS
        self.MEM (dict)     D(\theta, f) normalized directional specturm
                            S(f)         Power spectrum ([data]^2/Hz)
                            E(\theta, f) Directional Power Spectrum ([data]^2/Hz)
                            freq         frequency axis
                            theta        theta exis
                            unit_power       unit of S
                            unit_dir         unit of E

        """
        #theta=np.arange(1,361) if theta is None else theta

        MEM=dict()
        MEM['D'], MEM['S'], MEM['E'], MEM['freq'], MEM['theta'] = MEM_cal(self.moments_est,self.f,theta=theta ,flim=flim)
        MEM['fmin']=flim[0]
        MEM['fmax']=flim[1]
        MEM['unit_power']=self.moments_unit + '/Hz'
        MEM['unit_dir']=self.moments_unit + '/deg/Hz'

        self.MEM =MEM

    #def error(self, ci=0.95):
    #    self.El, self.Eu =spec_error(self.spec_est,self.n_spec,ci=ci)

    #def parceval(self):
        ## method not checked jet
    #    print('Parcevals Theorem:')
    #    print('variance of unweighted timeseries: ',self.data.var())
    #    print('mean variance of timeseries chunks: ',self.chunks.var(axis=1).mean() if self.save_chunks is True else 'data not saved')
        #print('variance of weighted timeseries: ',self.phi.var() )
        #self.calc_var(self)
    #    print('variance of the pwelch Spectrum: ',self.var)

    #def calc_var(self):

    #    """ Compute total variance from spectral density """
    #    self.var = self.specs[1:].mean(axis=0).sum()  # do not consider zeroth frequency



class pwelch(object):
    def __init__(self,data,dt,L=None, ov=None,window=None, save_chunks=False, plot_chunks=False, periodogram=False, prewhite=None):
        """
        prewhite    None(default)
                    1 = 1 timederivative
                    2=  2 timederivatives

                    L is the number of timesteps to overlab. must be an int.
        """

        from m_general_ph3 import echo_dt
        if prewhite is None:
            self.data = data      # field to be analyzed
        elif prewhite ==1:
            print('prewhite =1')
            self.data=np.gradient(data)
        elif prewhite ==2:
            self.data=np.gradient(np.gradient(data))

        self.data       =   data      # field to be analyzed
        self.dt         =   dt        # sampling interval
        self.save_chunks=   save_chunks
        dsize           =   dt * data.size

        L=int(np.round(data.size/10)) if L is None else L
        if type(L) != int:
            echo_dt(L)
            L=L.item().total_seconds()
        #else:
            #print('unknown L type')

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
        calc_freq(self)
        #print(data.size, L, ov, int(L-ov) )
        nbin=int(np.floor(dsize/(L-ov)))
        #print(nbin)
        if periodogram:
            self.nbin=nbin
            self.dt_periodogram=L-ov

        if save_chunks:
            chunks=np.empty([int(nbin),int(L)])

        specs=np.empty([int(nbin),self.f.size])
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

                calc_spectrum(self)
                specs[k,:]=self.spec
                last_k=k
                last_used_TS=int(i+L)
                #if plot_chunks:
                #    MT.stats_format(self.spec, 'spec '+str(i))
                del(self.spec)
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

        if prewhite is None:
            self.specs=specs[:last_k,:]
            self.spec_est=self.specs.mean(axis=0)
        elif prewhite ==1:
            self.specs=specs[:last_k,:]*(2*np.pi*self.f)
            self.spec_est=self.specs.mean(axis=0)
        elif prewhite ==2:
            self.specs=specs[:last_k,:]*(2*np.pi*self.f)**2
            self.spec_est=self.specs.mean(axis=0)


        self.n_spec,_=self.specs.shape
        self.calc_var()
        #self.phi=self.data
        #self.phi*=win*np.sqrt(factor)

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

class Spectogram_subsample(pwelch):
    def __init__(self, data, dt, dt_unit=None, timestamp=None, L=None, ov=None, subL=None, window=None, save_chunks=False, plot_chunks=False, ci=0.95 ):

        import m_general_ph3 as M
        self.hist='Subsampled Spectogram'
        L=int(np.round(data.size/10)) if L is None else L
        if type(L) != int:
            #print('Length = ')
            M.echo_dt(L)
            self.write_log('Length = '+ M.echo_dt(L, as_string=True))
            L=int(L.item().total_seconds())
        else:
            print('unknown L type')
            self.write_log('Length = '+ 'unknown L type')

        subL=int(np.round(L/10)) if subL is None else subL
        if type(subL) != int:
            #print('Subsample Length= ')
            #M.echo_dt(subL)
            self.write_log('Subsample Length= '+ M.echo_dt(subL, as_string=True))
            subL=int(subL.item().total_seconds())
        else:
            #print('unknown subL type')
            self.write_log('Length = '+ 'unknown subL type')

        ov=int(np.round(L/2.0)) if ov is None else ov
        if type(ov) != int:
            #print('overlab Length= ')
            #M.echo_dt(ov)
            ov=int(ov.item().total_seconds())
        else:
            pass
            #print('ov type:', type(ov))

        self.n = subL
        if window is None:
            win=np.hanning(self.n)
        else:
            win=window


        self.dt = float(dt)
        # test if n is even
        if (self.n%2):
            self.neven = False
        else:
            self.neven = True

        calc_freq(self)

        #dt_inms             = np.timedelta64(dt, dt_unit).astype('m8[ms]').astype(float)
        data_size_adjust    = dt * data.size # datasize in seconds

        nbin=np.floor(data_size_adjust  /(L-ov))
        self.nbin=nbin
        self.dt_periodogram=L-ov
        specs=np.empty([int(nbin),self.f.size])
        error_El=[]#np.empty([int(nbin),self.f.size])
        error_Ey=[]#np.empty([int(nbin),self.f.size])
        dbin = .005
        yN = np.arange(0,2.+dbin,dbin)

        n_specs=[]
        k=0
        print('subL', subL)
        print('L', L)
        print(data_size_adjust)
        #print(specs.shape, np.arange(0,data.size,int(L-ov)), np.arange(0,data.size-L,int(L-ov)).shape, )
        #print(specs.shape,np.arange(0,data_size_adjust-int(L-ov)+1,int(L-ov)) )
        for i in np.arange(0,data_size_adjust-int(L-ov)+1,int(L-ov)):

            phi=data[int(i):int(i+L)]
            #self.ii=np.append(self.ii,[i,i+L])
            #print(self.phi.max())

            #print(self.phi.mean())
            #print(phi.shape)
            #print('i',i, i+L, data.size, i+L <= data.size)
            #print(chunk.size, l)
            if i+L <= data.size:

                if plot_chunks:
                    plt.plot(phi)
                if save_chunks:
                    chunks[k,:]=phi

                phi=signal.detrend(phi)#*win

                pwelch.__init__(self,phi,dt,L=subL,window=window, save_chunks=False, plot_chunks=False, periodogram=False)

                specs[k,:]=self.spec_est
                sn,_=self.specs.shape
                n_specs=np.append(n_specs, sn)
                El,Ey=yNlu(sn,yN,ci=.95)
                error_El=np.append(error_El, El)
                error_Ey=np.append(error_Ey, Ey)
                #error_El[k,:]=El
                #error_Ey[k,:]=Ey

                del(self.spec_est)
            elif i+L > data.size:
                phi=data[-L:]
                #print('phi shape',phi.shape )
                if plot_chunks:
                    plt.plot(phi)
                if save_chunks:
                    chunks[k,:]=phi

                phi=signal.detrend(phi)#*win

                pwelch.__init__(self,phi,dt,L=subL,window=window, save_chunks=False, plot_chunks=False, periodogram=False)

                specs[k,:]=self.spec_est
                sn,_=self.specs.shape
                n_specs=np.append(n_specs, sn)
                El,Ey=yNlu(sn,yN,ci=.95)
                error_El=np.append(error_El, El)
                error_Ey=np.append(error_Ey, Ey)
                #error_El[k,:]=El
                #error_Ey[k,:]=Ey

                del(self.spec_est)

            #print(k)
            #print('------')
            k+=1


        self.n_specs=n_specs
        # assign output values
        self.data=specs # main data matrix
        self.specs =specs
        self.error_El=error_El
        self.error_Ey=error_Ey
        self.n_spec=n_specs
        self.calc_var()



        if dt_unit is not None and type(dt_unit) == str:
            dt_timedelta=np.timedelta64(int(dt),dt_unit)

        elif dt_unit is None:
            dt_unit='s'
            dt_timedelta=np.timedelta64(dt,dt_unit)

        #print('sample resolution:')
        #M.echo_dt(dt_timedelta)
        self.write_log('basic sample resolution:'+ M.echo_dt(dt_timedelta, as_string=True))

        timeres=np.timedelta64(int(self.dt_periodogram), 's')
        #print('time resolution:')
        #M.echo_dt(timeres)
        self.write_log('Spectogram time res:'+ M.echo_dt(timeres, as_string=True))

        if timestamp is None:
            start_time=np.datetime64(0,dt_unit)
            end_time=start_time+timeres*self.n_spec
        else:
            start_time=timestamp[0]
            end_time=timestamp[-1]

        #print('Periodogram starttime and endtime:')
        #print(start_time)
        #print(end_time)
        self.write_log('Spectogram starttme:'+ str(start_time))
        self.write_log('            endtime:'+ str(end_time))

        #print(type(start_time), type(L), type(end_time), type(dt_timedelta), print(timeres) )
        time = np.arange(start_time+int(L/2), end_time+dt_timedelta, timeres)
        #print(time)
        if time.shape[0] != self.data.shape[0]:
            time=time[0:self.data.shape[0]]
        self.time=time
        self.timeres=timeres
        self.dt_unit=dt_unit
        #print(self.time)
        self.time_dict=create_timeaxis_collection(self.time)
        self.write_log('created timestamp collection')

        self.log()

    def write_log(self, s, verbose=False):
        from m_tools_ph3 import write_log
        self.hist=write_log(self.hist, s, verbose=verbose)
    def log(self):
        print('.hist variable')
        print(self.hist)
    def power_anomalie(self, clim=None):
        dd=10*np.log10(self.data[:,:])
        #print(dd)
        #if anomalie is True:
        #    dd_tmp=dd.mean(axis=0).repeat(self.time.size-1)
        #        dd=dd- dd_tmp.reshape(self.fs.size,self.time.size-1).T

        self.data_power_mean=np.nanmedian(dd, axis=0) if clim is None else 10*np.log10(clim)
        dd_tmp=self.data_power_mean.repeat(self.time.size)
        print(self.data_power_mean.shape)
        print(self.f.size,self.time.size)
        self.data_power_ano=dd- dd_tmp.reshape(self.f.size,self.time.size).T

    def anomalie(self, clim=None):
        #dd=
        #print(dd)
        #if anomalie is True:
        #    dd_tmp=dd.mean(axis=0).repeat(self.time.size-1)
        #        dd=dd- dd_tmp.reshape(self.fs.size,self.time.size-1).T

        self.data_mean=np.nanmedian(self.data, axis=0) if clim is None else clim
        dd_tmp=self.data_mean.repeat(self.time.size)
        self.data_ano=self.data- dd_tmp.reshape(self.f.size,self.time.size).T


class Periodogram(pwelch):
    def __init__(self,data,dt,dt_unit=None, timestamp=None,L=None, ov=None,window=None, save_chunks=False, plot_chunks=False):

        pwelch.__init__(self,data,dt,L, ov,window, save_chunks, plot_chunks, periodogram=True)
        self.data=self.specs
        if dt_unit is not None and type(dt_unit) == str:
            dt_timedelta=np.timedelta64(int(dt),dt_unit)

        elif dt_unit is None:
            dt_unit='s'
            dt_timedelta=np.timedelta64(dt,dt_unit)

        print('sample resolution:')
        M.echo_dt(dt_timedelta)
        timeres=np.timedelta64(int(self.dt_periodogram), dt_unit)
        print('time resolution:')
        M.echo_dt(timeres)

        if timestamp is None:
            start_time=np.datetime64(0,dt_unit)
            end_time=start_time+timeres*self.n_spec
        else:
            start_time=timestamp[0]
            end_time=timestamp[-1]

        print('Periodogram starttime and endtime:')
        print(start_time)
        print(end_time)


        time = np.arange(start_time+L/2, end_time+dt_timedelta, timeres)
        #print(time)
        if time.shape[0] != self.data.shape[0]:
            time=time[0:self.data.shape[0]]
        self.time=time
        self.timeres=timeres
        self.dt_unit=dt_unit
        #G.dt_data=dt_data

        #self.time=np.arange(,nin*)
    def save_data(self, path=None, S=None):
        P=save_data_periodogram(self , S=S)
        print(P.meta)
        print('constructed class for saving')
        save_file(P, path)


class save_data_periodogram(object):
    def __init__(self,P, S=None):
        self.meta=S.meta if S is not None else ''
        self.data_unit=S.unit if S is not None else ''
        self.time=P.time
        self.f=P.f
        self.data=P.data
        self.dt_unit=P.dt_unit
        self.time=P.time
        #self.y_lim=y_lim
        self.timeres=P.timeres


def MEM_cal(moments_est,freq,theta=None ,flim=None):
    import getmem as getmem

    """
    computes Maximum Entrophy deinsity form calculated Moments

    Parameters:
    moments_est         dictionary with spectral moment estimates (units in  [data]^2)
                        should have the co-spectra P_nm and quadrature spectra Q_nm.

    freq (float)        frequency axis that corresponds to momment_est

    theta (np.array)    vector in degree (0, 360) that determines resolution of the directional spectrum
                        if not set, its 1 degree on 0 to 361.

    flim (float)        boolen of upper and lower limit for MEM calculation. If None flim is set to (freq[0], freq[-1])

    Returns:
    D,S, E, freq_sel,theta

    D(theta, freq_sel)      normalized spectrum with dimensions [360, N],
                        where N is the number of frequency bands
    S                   Power spectra ([data]^2/Hz)
    E(theta, freq_sel)  Directional Spectral Energy (([data]^2/deg/Hz))

    freq_sel            cutted freqency vector corresponding to flim
    theta               vector in degree that determines resolution of the directional spectrum


    """

    theta=np.arange(1,361) if theta is None else theta
    flim=(freq[0], freq[-1]) if flim is None else flim

    freq_sel_bool=M.cut_nparray(freq,flim[0], flim[1])
    freq_sel=freq[freq_sel_bool]


    N_sel=dict()
    for k,I in moments_est.iteritems():
        N_sel[k]=moments_est[k][freq_sel_bool]

    d1=N_sel['Q12'] / np.sqrt(N_sel['P11'] *(N_sel['P22'] + N_sel['P33'] ))
    # Lygre and Krongstad 1986 have here sqrt(N_sel['P11'] *(N_sel['P22'] + N_sel['P33']). I guess its a typo.
    d2=N_sel['Q13'] / np.sqrt(N_sel['P11'] *(N_sel['P22'] + N_sel['P33']) )

    d3=(N_sel['P22'] - N_sel['P33']) / (N_sel['P22'] + N_sel['P33'])
    d4=2*N_sel['P23'] / (N_sel['P22'] + N_sel['P33'])

    S=N_sel['P11']/freq_sel # define power spectrum

    D= GetMem(d1, d2, d3, d3, theta )
    E=np.tile(S.T, (theta.size,1)).T*D

    return D, S, E, freq_sel, theta


def GetMem(a1,b1,a2,b2, theta):
    import numpy as np
    from numpy import pi, sin, cos

    """
    function from Bia.
    10/25/17: Slidely alterted, theta IO

    Compute the directional spectrum using the Maximum Entropy
    Method (MEM) from the four angular moments (a1, b1, a2, b2)
    based on Lygre and Krogstad (1986).

    To generate  a directional spectrum with energy
    density of m^2/deg-Hz, multiply the N normalized distributions
    by the wave energy at each frequency.

    This script is based on O'Reilly's Matlab code
    Parameters
    --------------------------------------
    (a1, b1, a2, b2): float
        1st and 2nd moment normalized directional fourier coefficients.

    theta (np.array)    vector in degree (0, 360) that determines resolution of the directional spectrum
    Returns
    ----------------------------------------
    norm_mem: array
        normalized spectrum with dimensions [360, N],
        where N is the number of frequency bands
    """

    d1 = np.reshape(a1,[len(a1),1])
    d2 = np.reshape(b1,[len(b1),1])
    d3 = np.reshape(a2,[len(a2),1])
    d4 = np.reshape(b2,[len(b2),1])

    c1 = (1.+ 0j)*d1 + (0 + 1.j)*d2
    c2 = (1.+ 0j)*d3 + (0 + 1.j)*d4

    p1 = (c1 - c2*c1.conjugate())/(1-abs(c1)**2)
    p2 = c2 - c1*p1
    x = 1.-p1*c1.conjugate() - p2*c2.conjugate()

    # calculate MEM using 1 deg resolution

    a = theta.T
    #a = np.arange(1,361).T
    a = a*pi/180.
    a = np.reshape(a,[1,len(a)])

    e1 = (1. + 0j)*cos(a) - (0 + 1.j)*sin(a)
    e2 = (1. + 0j)*cos(2*a) - (0 + 1.j)*sin(2*a)
    y = abs((1.+ 0j) - np.dot(p1,e1) - np.dot(p2,e2))**2

    mem = abs((x*np.ones([1,a.size]))/y)

    summing = np.sum(mem, axis=1)
    sum_mem = 1./np.reshape(summing, [len(summing),1])
    norm_mem = np.dot(sum_mem, np.ones([1,a.size]))*mem

    return norm_mem


def peak_angle(data, freq,fpos,  max_jump=5, smooth_l=5):
    """
    This method computes the frequency of D mean between the given frequency limits  fpos (fpos[0], fpos[-1])
    it identifies the maxima and save them as degree from true north

    the maxima are detected using a smoothing and minimum distance between them

    Return:
    peak angles index orded buy the peak height
    """
    #max_jump=5
    #smooth_l=5

    if data.shape[0] != freq.size:
        raise Warning('shape of data in dim=0 must be equal to freq')

    f_bool=M.cut_nparray(freq, fpos[0], fpos[-1])
    dregree_series=data[f_bool,:].mean(0)
    #plt.plot(dregree_series)
    #print(dregree_series)

    a=M.find_max_ts(dregree_series,smooth=True, spreed=smooth_l, plot=False, jump=max_jump,verbose=False)
    if len(a[0]) ==0:
        a=M.find_max_ts(dregree_series,smooth=False, spreed=1, plot=False, jump=None,verbose=False)

    angles=np.flipud(np.sort(a[0]))

    return angles

def save_file( data, path):
    import pickle
    outfile=path
    f = open(outfile, 'wb')
    pickle.dump(data,f, pickle.HIGHEST_PROTOCOL)
    print('saved to:',outfile)
    f.close()

## ceasars funcitons
def yNlu(sn,yN,ci):
    """ compute yN[l] yN[u], that is, the lower and
                upper limit of yN """

    # cdf of chi^2 dist. with 2*sn DOF
    cdf = gammainc(sn,sn*yN)

    # indices that delimit the wedge of the conf. interval
    fl = np.abs(cdf - ci).argmin()
    fu = np.abs(cdf - 1. + ci).argmin()

    return yN[fl],yN[fu]


def spec_error(E,sn,ci=.95):

    """ Computes confidence interval for one-dimensional spectral
        estimate E.

        Parameters
        ===========
        - sn is the number of spectral realizations;
                it can be either an scalar or an array of size(E)
        - ci = .95 for 95 % confidence interval

        Output
        ==========
        lower (El) and upper (Eu) bounds on E """

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

# for llc output only; this is temporary
def spec_est(A,d,axis=2,window=True,detrend=True, prewhiten=False):

    l1,l2,l3,l4 = A.shape

    if axis==2:
        l = l3
        if prewhiten:
            if l3>1:
                _,A,_ = np.gradient(A,d,d,1.)
            else:
                _,A = np.gradient(A.squeeze(),d,d)
                A = A[...,np.newaxis]
        if detrend:
            A = signal.detrend(A,axis=axis,type='linear')
        if window:
            win = np.hanning(l)
            win = (l/(win**2).sum())*win
            win = win[np.newaxis,np.newaxis,:,np.newaxis]
        else:
            win = np.ones(l)[np.newaxis,np.newaxis,:,np.newaxis]
    elif axis==1:
        l = l2
        if prewhiten:
            if l3 >1:
                A,_,_ = np.gradient(A,d,d,1.)
            else:
                A,_ = np.gradient(A.squeeze(),d,d)
                A = A[...,np.newaxis]
        if detrend:
            A = signal.detrend(A,axis=1,type='linear')
        if window:
            win = np.hanning(l)
            win = (l/(win**2).sum())*win
            win = win[np.newaxis,...,np.newaxis,np.newaxis]
        else:
            win = np.ones(l)[np.newaxis,...,np.newaxis,np.newaxis]


    df = 1./(d*l)
    f = np.arange(0,l/2+1)*df

    Ahat = np.fft.rfft(win*A,axis=axis)
    Aabs = 2 * (Ahat*Ahat.conjugate()) / l

    if prewhiten:

        if axis==1:
            fd = 2*np.pi*f[np.newaxis,:, np.newaxis]
        else:
            fd = 2*np.pi*f[...,np.newaxis,np.newaxis]

        Aabs = Aabs/(fd**2)
        Aabs[0,0] = 0.

    return Aabs,f
