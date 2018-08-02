# Read the header of a Scan file

from array import array
import os
import math
import tkFileDialog
#from PlotHistogram import PlotHistogram
#from ExpPdf import ExpPdf
#from Simplex import Simplex
#from SimplexDC import SimplexDC


def readHeader (fname,verb=1):
    'Reads header of a scan file'

    header = {}                            #          a blank list to put each header item into
    floats = array ('f')                #        make dummy arrays to read floats, integers (LONG in C) and short integers from the header
    ints = array('i')
    shorts = array ('h')
    doubles = array('d')
    longs = array('l')
    
    try:
        f=open(fname, 'r')                    #        open the .ssd file as read only
    
        if verb: 
            print "Full header from SCAN file %s is as follows:" %fname
        else:
            print "Highlights of SCAN header from %s" %fname
        
        ints.fromfile(f,1)                #        read a short into the dummy array
        header['iscanver'] = ints.pop()                #        pop the short from the array into a variable
        # new scan files- version 104, 103 (simulated) and -103
        #if verbose: print 'version', iscanver
        
        ints.fromfile(f,1)
        header['ioffset'] = ints.pop()
        #if verbose: print 'ioffset', ioffset
        
        ints.fromfile(f,1)
        header['nint'] = ints.pop()
        #if verbose: print 'nint', nint
        
        header['title'] = f.read(70)
        
        header['date'] = f.read(11)
        
        header['defname'] = f.read(6)
        
        header['tapeID'] = f.read(24)
        
        ints.fromfile(f,1)
        header['ipatch'] = ints.pop()
        
        ints.fromfile(f,1)
        header['npatch'] = ints.pop()
        
        floats.fromfile(f,1)
        header['Emem'] = floats.pop()
        
        floats.fromfile(f,1)
        header['temper'] = floats.pop()
        
        header['adcfil'] = f.read(30)
        
        header['qfile1'] = f.read(35)
        
        # logical; true if data from CJUMP file
        ints.fromfile(f,1)
        header['cjump'] = ints.pop()
        
        ints.fromfile(f,1)
        header['nfits'] = ints.pop()
        
        ints.fromfile(f,1)
        header['ntmax'] = ints.pop()
        
        ints.fromfile(f,1)
        header['nfmax'] = ints.pop()

        # Number of data points read into memory at each disk read -bigger
        # the better (max depends on how much RAM you have).
        # nbuf=131072		!=1024*128
        ints.fromfile(f,1)
        header['nbuf'] = ints.pop()
        
        # Number of extra points read in, at each end of data section to
        # allow display of transitions on section boundaries; 2048 is OK usually.
        ints.fromfile(f,1)
        header['novlap'] = ints.pop()
        
        # Sample rate (Hz)
        floats.fromfile(f,1)
        header['srate'] = floats.pop()
        
        # finter = microsec between data points; finter=1.e6/srate
        floats.fromfile(f,1)
        header['finter'] = floats.pop()
        
        # TSECT=time (microsec) from first point of one section to first point of next
        # tsect=float(nbuf)*finter
        floats.fromfile(f,1)
        header['tsect'] = floats.pop()
        
        # The first data point in data file, idata(1) starts at byte (record #) ioff+1
        ints.fromfile(f,1)
        header['ioff'] = ints.pop()
        
        ints.fromfile(f,1)
        header['ndat'] = ints.pop()
        
        # calc nsec etc here, in case default nbuf altered 
        # if(ndat.lt.nbuf) nbuf=ndat    !allocate smaller array
        # nsec= 1 + (ndat-1)/nbuf  !number of sections
        ints.fromfile(f,1)
        header['nsec'] = ints.pop()
        
        # nrlast=ndat - (nsec-1)*nbuf  !number of idata in last section
        ints.fromfile(f,1)
        header['nrlast'] = ints.pop()
        
        floats.fromfile(f,1)
        header['avtot'] = floats.pop()
        
        ints.fromfile(f,1)
        header['navamp'] = ints.pop()
        
        floats.fromfile(f,1)
        header['avamp'] = floats.pop()
        
        floats.fromfile(f,1)
        header['rms'] = floats.pop()
        
        # Data will be written to disk at (approx) every nth transition, so
        # analysis can be restarted by using the ''restart'' option when SCAN reentered.
        ints.fromfile(f,1)
        header['nwrit'] = ints.pop()

        # nwsav=0		!used for auto disc write
        ints.fromfile(f,1)
        header['nwsav'] = ints.pop()
        
        # logical
        ints.fromfile(f,1)
        header['newpar'] = ints.pop()

        # logical
        ints.fromfile(f,1)
        header['opendown'] = ints.pop()

        # logical; Invert trace (openings must be downwards)
        ints.fromfile(f,1)
        header['invert'] = ints.pop()
        
        # logical; usepots=.false.
        ints.fromfile(f,1)
        header['usepots'] = ints.pop()
        
        # in SCAN: Display only (no step-response function)
        ints.fromfile(f,1)
        header['disp'] = ints.pop()
        
        # if(iscrit.eq.1): Percentage of full amplitude for critical level
        # (Scrit) beyond which transition is deemed to occur.
        # if(iscrit.eq.2): Multiple of RMS noise to define critical level
        # (Scrit) beyond which transition is deemed to occur.
        # if(iscrit.eq.1): smult=0.14 !scrit=0.14*avamp
        # if(iscrit.eq.2): smult=5. scrit=5.0*rms
        floats.fromfile(f,1)
        header['smult'] = floats.pop()
        
        ints.fromfile(f,1)
        header['scrit'] = ints.pop()
        
        ints.fromfile(f,1)
        header['vary'] = ints.pop()

        # Number of consecutive points beyond Scrit for a transition to be
        # deemed to have occurred. Default: ntrig=2
        ints.fromfile(f,1)
        header['ntrig'] = ints.pop()
        
        # navtest=ntrig-1 ; if(navtest.le.0) navtest=1
        # navtest=number averaged before average curlev is used, rather than
        # input curlev in FINDTRANS (NB must be less than ntrig, or, for
        # example, if input baseline is not close to current baseline
        # (ie baseline has drifted since last time) then will get a 'trigger'
        # straight away!
        ints.fromfile(f,1)
        header['navtest'] = ints.pop()

        # Trace will be amplified by this factor before display (but better to
        # amplify correctly BEFORE sampling). DGAIN=1.0
        floats.fromfile(f,1)
        header['dgain'] = floats.pop()

        # IBOFF=0		!BASELINE OFFSET FOR DISPLAY (ADC)
        ints.fromfile(f,1)
        header['iboff'] = ints.pop()
        
        # Factor by which trace is expanded when ''expand'' is first hit.
        # expfac=2.
        floats.fromfile(f,1)
        header['expfac'] = floats.pop()

        # Position of baseline on screen is offset to this level after initial
        # ''get piece of baseline on screen'' is completed.
        # bdisp=0.75 if openings downwards; bdisp=0.25 if openings upwards
        floats.fromfile(f,1)
        header['bdisp'] = floats.pop()
        
        ints.fromfile(f,1)
        header['ibflag'] = ints.pop()
        
        # Auto-fit to avoid sublevels if possible. In case of doubt fit brief
        # open-shut-open rather than fitting a sublevel.
        ints.fromfile(f,1)
        header['iautosub'] = ints.pop()

        # When opening crosses the red trigger line display stops with the
        # opening transition at this point on the x-axis of display.
        # xtrig=0.2: trigger at 20% of X axis on screen
        floats.fromfile(f,1)
        header['xtrig'] = floats.pop()
        
        # ndev='C:'; disk partition for Windows
        header['ndev'] = f.read(2)
        
        header['cdate'] = f.read(11)
        print header['cdate']
        
        header['adctime'] = f.read(8)

        ints.fromfile(f, 1)
        header['nsetup'] = ints.pop()

        header['filtfile'] = f.read(20)

        # Low pass filter (Hz, -3dB)
        # later needs to be converted to kHz
        floats.fromfile(f, 1)
        header['ffilt'] = floats.pop()

        # npfilt=number of points to jump forward after a transition, to start
        # search for next transition
        # npfilt1= number of data points for filter to go from 1% to 99%
        # npfilt1=ifixr((tf99-tf1)/finter)
        # npfilt=ifixr(float(npfilt1)*facjump)
        ints.fromfile(f, 1)
        header['npfilt'] = ints.pop()

        # sfac1=(yd2-yd1)/65536.
        # sfac1=sfac1*dgain			!true scal fac for ADC to pixel units
        floats.fromfile(f, 1)
        header['sfac1'] = floats.pop()

        # nscale=1 + ifix(alog(4096./(yd2-yd1))/alog(2.))
        # sfac2=sfac1*float(2**nscale)	!Converts ADC units to intermed units
        floats.fromfile(f, 1)
        header['sfac2'] = floats.pop()

        # sfac3=1.0/float(2**nscale) 	!converts intermed to pixel units
        floats.fromfile(f, 1)
        header['sfac3'] = floats.pop()

        ints.fromfile(f, 1)
        header['nscale'] = ints.pop()

        # Calibration factor (pA per ADC unit)
        floats.fromfile(f, 1)
        header['calfac'] = floats.pop()

        # calfac1=calfac/sfac1		!converts pixel display units to pA
        floats.fromfile(f, 1)
        header['calfac1'] = floats.pop()

        # calfac2=calfac/sfac2		!converts intermed units to pA
        floats.fromfile(f, 1)
        header['calfac2'] = floats.pop()

        # iyoff=ifixr(yd1 + 0.5*(yd2-yd1))	!zero in centre until baseline done
        # (NB iyoff is in pixel units)
        ints.fromfile(f, 1)
        header['iyoff'] = ints.pop()

        ints.fromfile(f, 1)
        header['ioff1'] = ints.pop()

        # Show position of guessed transition points on screen as purple line
        # + blue line to mark end of transition.
        ints.fromfile(f, 1)
        header['disptran'] = ints.pop()

        # When first derivative used to identify two closely-spaced
        # transitions, display it below the trace.
        ints.fromfile(f, 1)
        header['dispderiv'] = ints.pop()

        # dispguess=.true.
        ints.fromfile(f, 1)
        header['dispguess'] = ints.pop()

        # Amplitude difference (as fraction of full amp) below which openings
        # are deemed to have ''same'' amplitude: for (a) elim of short gaps
        # (b) setting guesses. ampfac=0.05
        floats.fromfile(f, 1)
        header['ampfac'] = floats.pop()

        # Length of fitted event (microsec) below which refit, omitting short
        # events, is offered.  Events guessed to be shorter than this are
        # rejected before fitting. tmin=15
        floats.fromfile(f, 1)
        header['tmin'] = floats.pop()

        # Length (multiple of risetime) of event below which its amplitude is
        # fixed (also length guessed from peak amplitude).
        # tsfac=2.0		!tsfac*trise=tshort
        floats.fromfile(f, 1)
        header['tsfac'] = floats.pop()

        # Length (multiple of risetime) of event above which amplitude is
        # ''well-defined'' so usable to fix length of an adjacent brief opening.
        # tlfac=3.0		!tlfac*trise=tlong
        floats.fromfile(f, 1)
        header['tlfac'] = floats.pop()

        # sdone=.false.	!no baseline SD yet -NO -set BEFORE INSCAN
        ints.fromfile(f, 1)
        header['sdone'] = ints.pop()

        # real*8 in fortran;double precizion value of finter
        # Need double prec versions of finter if time of transition from 1st point
        # in CONSAM to be recorded accurately (see FITSUB).  At present, cannot
        # have smaller finter than 0.25 microsec (with 4 MHz clock) so want to get
        # rid of non sig figs when DBLE(finter) is calc
        doubles.fromfile(f, 1)
        header['dfinter'] = doubles.pop()

        doubles.fromfile(f, 1)
        header['tlast'] = doubles.pop()

        ints.fromfile(f, 1)
        header['shut'] = ints.pop()

        ints.fromfile(f, 1)
        header['shutprev'] = ints.pop()

        ints.fromfile(f, 1)
        header['backward'] = ints.pop()

        ints.fromfile(f, 1)
        header['prevlevel'] = ints.pop()

        floats.fromfile(f, 1)
        header['t0sav'] = floats.pop()

        floats.fromfile(f, 1)
        header['y0sav'] = floats.pop()

        floats.fromfile(f, 1)
        header['vard'] = floats.pop()

        # Number of points before first, and after last, transition to be
        # fitted in auto mode. Number of shut points to be fitted at ends.
        # Default: nshutfit=50;
        ints.fromfile(f, 1)
        header['nshutfit'] = ints.pop()

        ints.fromfile(f, 1)
        header['infit'] = ints.pop()

        ints.fromfile(f, 1)
        header['infirst'] = ints.pop()

        ints.fromfile(f, 1)
        header['ixfprev'] = ints.pop()

        ints.fromfile(f, 1)
        header['idiskq'] = ints.pop()

        ints.fromfile(f, 1)
        header['ifirst'] = ints.pop()

        ints.fromfile(f, 1)
        header['base'] = ints.pop()

        ints.fromfile(f, 1)
        header['basevga'] = ints.pop()

        ints.fromfile(f, 1)
        header['ibasevga'] = ints.pop()

        ints.fromfile(f, 1)
        header['itrig'] = ints.pop()

        ints.fromfile(f, 1)
        header['itrigvga'] = ints.pop()

        ints.fromfile(f, 1)
        header['itriglev'] = ints.pop()

        ints.fromfile(f, 1)
        header['inc'] = ints.pop()

        ints.fromfile(f, 1)
        header['incabs'] = ints.pop()

        ints.fromfile(f, 1)
        header['indfast'] = ints.pop()

        ints.fromfile(f, 1)
        header['isdfst'] = ints.pop()

        ints.fromfile(f, 1)
        header['isec'] = ints.pop()

        ints.fromfile(f, 1)
        header['ndisp'] = ints.pop()

        ints.fromfile(f, 1)
        header['ndisp1'] = ints.pop()

        #? idatyp=0 for usual CONSAM, idatyp=1 for pdpdata and 
        #? idatyp=2 for Axon data
        ints.fromfile(f, 1)
        header['idatyp1'] = ints.pop()

        header['cdate1'] = f.read(11)

        # number of channels in patch
        ints.fromfile(f, 1)
        header['nchan'] = ints.pop()

        # Length (multiple of risetime) of interval between two transitions
        # (in same direction) below which an attempt is made to fit brief
        # events rather than sublevel.
        # tcfac=4.		!tcfac*trise=tclose
        # ''Close'' transitions (multiple of risetime)'	!tclose
        floats.fromfile(f, 1)
        header['tcfac'] = floats.pop()

        # Fraction of step-response (filter) length (1-99%) allowed after a
        # transition before search for next transition is started. facjump=0.6
        floats.fromfile(f, 1)
        header['facjump'] = floats.pop()

        ints.fromfile(f, 1)
        header['shutsav'] = ints.pop()

        ints.fromfile(f, 1)
        header['goback'] = ints.pop()

        ints.fromfile(f, 1)
        header['imin'] = ints.pop()

        ints.fromfile(f, 1)
        header['imax'] = ints.pop()

        # Factor by which initial guess must be reduced before Simplex
        # converges; (e.g. 0.01=low precision; 0.0001=high precision)
        # errfac=0.005
        floats.fromfile(f, 1)
        header['errfac'] = floats.pop()

        # Multiple of SD of 1st deriv used to find inflections; small value
        # e.g. 2.0 makes it more likely that multiple transitions fitted,
        # rather than sublevel. Sensitivity for multiple trans (vs sublevel).
        # derivfac=3.
        floats.fromfile(f, 1)
        header['derivfac'] = floats.pop()

        # Controls how fast the simplex contracts around a putative minimum.
        # Usually 0.5; smaller value (down to 0.2) gives faster convergence
        # but fit may be worse.
        floats.fromfile(f, 1)
        header['confac'] = floats.pop()

        ints.fromfile(f, 1)
        header['nsweep'] = ints.pop()

        # njdim=nsamp/njump
        ints.fromfile(f, 1)
        header['njdim'] = ints.pop()

        doubles.fromfile(f, 1)
        header['tzerod'] = doubles.pop()

        floats.fromfile(f, 1)
        header['intzero'] = floats.pop()

        # real*8 in fortran
        doubles.fromfile(f, 1)
        header['tsample'] = doubles.pop()

        floats.fromfile(f, 1)
        header['ktjump'] = floats.pop()

        ints.fromfile(f, 1)
        header['njfit'] = ints.pop()

        ints.fromfile(f, 1)
        header['njump'] = ints.pop()

        ints.fromfile(f, 1)
        header['nnull'] = ints.pop()

        floats.fromfile(f, 1)
        header['ktlast'] = floats.pop()

        #Default zoom factor (must be an integer power of 2 (i.e. 1,2,4,8...etc)
        ints.fromfile(f, 1)
        header['izoom'] = ints.pop()

        # Filter cut-off (Hz, -3 dB) when zoomed (default = same as input data,
        # i.e. no extra filtering when zoomed)
        floats.fromfile(f, 1)
        header['fcz'] = floats.pop()

        # fczoom=actual fc applied to achieve fcz (in kHz)
        # fczoom=1./sqrt((1.0/(fcz*fcz)) - (1.0/(ffilt*ffilt)))
        floats.fromfile(f, 1)
        header['fczoom'] = floats.pop()

        # ''Full amplitude'' to be used to identify transitions when zoomed
        # (approx size of small channels).
        # avampz=alternative 'full amplitude' to use for event detection
        # when zoomed/filtered (NB use zoomfac=1 or fczoom=filtf to get 
        # filtering only or zoom only). avampz in intermed units 
        # (same as avamp); ampz=same thing pA
        floats.fromfile(f, 1)
        header['ampz'] = floats.pop()

        floats.fromfile(f, 1)
        header['avampsav'] = floats.pop()

        # No. of iterations done by SIMPLEX before swapping to DFPMIN
        # (Davidon-Fletcher-Powell minimisation) when both used.
        # itsimp=1200
        ints.fromfile(f, 1)
        header['itsimp'] = ints.pop()

        # 1=SIMPLEX only (as before). 2=SIMPLEX for fixed number of iterations,
        # then swap to DFPMIN. 3=Davidon-Fletcher-Powell (DFPMIN) only.
        ints.fromfile(f, 1)
        header['minmeth'] = ints.pop()

        # Minimum # of points at the shut level within a fit required for this
        # shut level to be used to estimate the amplitude of nearby openings.
        # nbasemin=10
        ints.fromfile(f, 1)
        header['nbasemin'] = ints.pop()

        # Definition of threshold for detecting transitions: 1=fraction of
        # mean full amp; 2=multiple of rms noise. (Level is set as
        # parameter 3, page 1).
        # iscrit=1 to use scrit=smult*avamp (recent versions of scan)
        # iscrit=2 to use scrit=smult*rms (as in original version of scan)
        # latter could be better for records with wide range of amps
        ints.fromfile(f, 1)
        header['iscrit'] = ints.pop()

        # Use reduced critical amplitude for detecting transitions when all
        # amplitudes in the fit are smaller than the full amplitude.
        # Use variable critical level (see HELP)
        # default scritvar=.false.
        # if scritvar=true, smaller scrit is used for fittings that have only
        # small amplitudes. In this case value of scrit changed after region
        # to be fitted has been defined, with minimum smult=smultmin
        ints.fromfile(f, 1)
        header['scritvar'] = ints.pop()

        # Minimum crit level when using lower threshold for detecting
        # transitions when all amplitudes in the fit are smaller than the
        # full amplitude (as mult of RMS).
        # smultmin=2.5		!2.5*rms1 in intermed units
        # Minimum critical level (see HELP): multiple of RMS
        floats.fromfile(f, 1)
        header['smultmin'] = floats.pop()

        floats.fromfile(f, 1)
        header['stpfac'] = floats.pop()

        ints.fromfile(f, 1)
        header['nlig'] = ints.pop()

        floats.fromfile(f, 1)
        header['conc1'] = floats.pop()



        f.close() #    close the file
    
    except:
        print "Scan file header not read correctly. If you are lucky you can input what is needed later"
        
    return header

def erf(z):
    'from: http://www.cs.princeton.edu/introcs/21function/ErrorFunction.java.html \
    Implements the Gauss error function. \
    erf(z) = 2 / sqrt(pi) * integral(exp(-t*t), t = 0..z) \
    fractional error in math formula less than 1.2 * 10 ^ -7. \
    although subject to catastrophic cancellation when z in very close to 0 \
    from Chebyshev fitting formula for erf(z) from Numerical Recipes, 6.2.'

    t = 1.0 / (1.0 + 0.5 * abs(z))
    # use Horner's method
    ans = 1 - t * math.exp( -z*z -  1.26551223 +
                                            t * ( 1.00002368 +
                                            t * ( 0.37409196 +
                                            t * ( 0.09678418 +
                                            t * (-0.18628806 +
                                            t * ( 0.27886807 +
                                            t * (-1.13520398 +
                                            t * ( 1.48851587 +
                                            t * (-0.82215223 +
                                            t * ( 0.17087277))))))))))
    if z >= 0.0:
            return ans
    else:
            return -ans

def readData(fname, ioffset, nint, calfac2):
    'Reads data- intervals, amplitudes, options- from SCAN file.'

    # make dummy arrays to read floats, integers (LONG in C) and short integers from the header
    floats = array ('f')    # 4 byte float
    ints = array('i')     # 4 byte integer
    shorts = array ('h')    # 2 byte integer
    ops = array('b')    # 1 byte integer
    # Data=
    # real*4 tint(1...nint) 	   4nint bytes
    # integer*2 iampl(1..nint)   2nint bytes
    # integer*1 iprops(1..nint)  nint  bytes
    # Total storage needed=7*nint bytes

    tint = []
    iampl = []
    iprops = []
    'integer*1 iprops(i) holds properties of ith duration and amp \
    (integer*1 has range -128 to +127 (bit 7 set gives -128; can use bits 0-6) \
    0 = all OK; \
    1 = amplitude dubious = bit 0; \
    2 = amplitude fixed = bit 1; \
    4 = amplitude of opening constrained (see fixamp) = bit 2; \
    8 = duration unusable = bit 3; etc \
    and keep sum of values of more than one property is true.'

    f=open(fname, 'rb')
    f.seek(ioffset-1)
    for i in range(0, nint):
        floats.fromfile(f, 1)
        tint.append(floats.pop())
    for i in range(0, nint):
        shorts.fromfile(f, 1)
        iampl.append(shorts.pop()*calfac2)
    for i in range(0, nint):
        ops.fromfile(f, 1)
        iprops.append(ops.pop())

    if tint[-1] == 0:
        if iprops[-1] != 8:
            iprops[-1] = 8
            print 'Last interval in file was shut. It was set as unusable.'
    else:
        tint.append(-1.0)
        iampl.append(0.0)
        iprops.append(8)
        print 'Last interval in file was open. An unusable shut time was \
        \n inserted at the end. Total number of intervals increased by one.'

    f.close()
    return tint, iampl, iprops

def falseEvents(tres,fc,rms,avamp):
    'Version for EKDIST/new SCAN (avamp, rms already in pA). \
    To calc false event rate (per sec) in EKDIST (from RESINT.) \
    First calc threshold as amp attained by pulse of length=tres (in ms)'

    w = tres * 1.e-3    # sec
    u = erf(2668. * fc * w)
    amp = avamp    # full amp (pA)
    phi = u * amp    #'threshold' (pA)
    var = (rms) ** 2    # noise variance (pA)**2

    # Calc rate from C & Sigworth eq. 9, with k=1
    frate = 1000. * fc * math.exp(-(phi * phi) / (2. * var))

    return frate    # false event rate

def imposeResolution(tres, itint, iampl, iprops):
    ''

    # if last interval is shut, then set it as unusable
    if iampl[-1] == 0: iprops[-1] = 8

    # check for negative intervals and set them unusable. Skip those already set unusable
    for i in range(0, len(itint)):
        if (itint[i] < 0) and (iprops[i] != 8):
            iprops[i] = 8
            print '\n interval %d set unusable.'

    # IF THE FIRST INTERVAL IS BOTH USABLE AND RESOLVABLE
    # THIS IS STARTING POINT. IF NOT LOOK FOR FIRST INTERVAL THAT
    # IS BOTH, AND IS PRECEDED BY AN RESOLVABLE INTERVAL TOO (OTHERWISE
    # ITS START WILL BE DEFINED ONLY BY THE POSITION OF THE PRECEDING
    # UNRESOLVABLE INTERVAL AND SO WILL BE UNRELIABLE)

    n = 0
    firstResolved = 0
    if (itint[n] > tres) and (iprops[n] != 8):
        firstResolved = 1    # first interval is usable and resolvable
        print 'first usable interval is ', n

    while not firstResolved:
        n = n + 1
        if (itint[n] > tres) and (iprops[n] != 8) and (itint[n-1] > tres)and (iprops[n-1] != 8):
            firstResolved = 1    # first interval is usable and resolvable
            print 'first usable interval is ', n
        else: n = n + 1

    # NOW START TO LOOK FOR UNRESOLVABLE INTERVALS

    # (1) A concantenated shut period starts with a good, resolvable shutting
    #     and ends when first good resolvable opening found.
    #     Length of concat shut period=sum of all durations before the resol opening
    #     Amplitude of concat shut period=0
    # (2) A concantenated open period starts with a good, resolvable opening
    #     and ends when first good resolvable interval is found that
    #     has a different amplitude (either shut, or open but diff amplitude).
    #     Length of concat open period=sum of all concatenated durations
    #     Amplitude of concat open period weighted mean amp of all concat intervals
    # Simpler to have separate code for shut groups and for open groups. If first
    # interval of group is shut then set shutint=true.

    # First interval of any concat group must be good and resolvable so
    # insert warning to check this

    # First interval in each concatenated group must be resolvable, but may
    # be bad (in which case next group will be bad). !!! This in DC's prog does not work!!!

    setdub = 0    # false.     !initialise setdub

    otint = []
    oampl = []
    oprops = []

    nc = 1    # number of intervals concat in each group
    # if nc=1 transfer properties directly from input to output
    ni = 0    # counts intervals in output lists

    otint.append(itint[n])    # tint[ni)=tint0(n)    !start new concat group
    oampl.append(iampl[n])
    oprops.append(iprops[n])
    ttemp = otint[-1]
    aavtemp = oampl[-1] * otint[-1]
    atemp = iampl[n]
    n = n + 1

    n_sa = 0
    n_fx = 0

    while n < len(itint):

        if itint[n] < tres:
            otint[-1] = otint[-1] + itint[n]
            if (oampl[-1] != 0) and (iampl[n] != 0):
                aavtemp = aavtemp + iampl[n]*itint[n]
                ttemp = ttemp + itint[n]
                oampl[-1] = aavtemp / ttemp
            if iprops[n] == 8:
                oprops[-1] = 8
            n = n + 1
        else:
            if iprops[n] == 4 and oampl[-1] != 0 and oampl[-1] == iampl[n]:
                otint[-1] = otint[-1] + itint[n]
                #aavtemp = aavtemp + iampl[n]*itint[n]
                #ttemp = ttemp + itint[n]
                n = n + 1
                n_fx = n_fx + 1
            elif iprops[n] == 6 and oampl[-1] != 0 and oampl[-1] == iampl[n]:
                otint[-1] = otint[-1] + itint[n]
                #aavtemp = aavtemp + iampl[n]*itint[n]
                #ttemp = ttemp + itint[n]
                n = n + 1
                n_fx = n_fx + 1

            elif iampl[n-1] == iampl[n]:    # elif
                otint[-1] = otint[-1] + itint[n]
                aavtemp = aavtemp + iampl[n]*itint[n]
                ttemp = ttemp + itint[n]
                oampl[-1] = aavtemp / ttemp
                n = n + 1
                n_sa = n_sa + 1
            elif iampl[n] == 0 and oampl[-1] == 0:
                otint[-1] = otint[-1] + itint[n]
                n = n + 1
            else:
                otint.append(itint[n])
                oampl.append(iampl[n])
                oprops.append(iprops[n])
                ttemp = otint[-1]
                aavtemp = oampl[-1] * otint[-1]
                atemp = iampl[n]
                n = n + 1

    print n_sa, 'pairs of same amplitude openings.'
    print n_fx, 'removed openings of fixed amplitude '

    return otint, oampl, oprops

def getOpenShutPeriods(otint, oampl, oprops):
    '\n There may be many small amplitude transitions during one \
    \n ''opening'', each of which will count as an individual \
    \n opening, so generally better to look at ''open periods'''

    ' Options for distribution type:'
    ' (1) Duration of individual apparent openings'
    ' (2) Duration of contiguous open times (open PERIODS)'
    '   and mean length of open periods that are adjacent to shut'
    '   times or open periods with duration in specified ranges'
    ' (3) Distribution of lengths of individual app. openings that'
    '    are adjacent to a gap with duration in a specified range'
    ' (4) Distribution of lengths of apparent open PERIODS that'
    '    are adjacent to a gap with duration in a specified range'
    ' (5) Distribution of lengths of apparent open PERIODS that'
    '    bounded on BOTH sides by shut times in a specified range'

    'Specify amplitude range for openings'
    amplo = -10000.0    # default is to include openings of any amplitude
    amphi = 10000.0
    '  Low amp, high amp (pA with sign) = '

    ' Exclude openings with ''dubious'' amps [Y] ? '
    ' Exclude openings with ''assumed'' amps [Y] ? '
    exass = 0    # logical exclude assumed or dubious = 1
    ' (i.e. ANY dubious amp excludes whole open period)'

    #  Look for start of a group of openings i.e. any opening that has
    #  defined duration (i.e. usable).  A single unusable opening in a group
    #  makes its length undefined so it is excluded.
    # NEW VERSION -ENSURES EACH OPEN PERIOD STARTS WITH SHUT-OPEN TRANSITION
    # Find start of a group (open period) -valid start must have a good shut
    # time followed by a good opening -if a bad opening is found as first (or
    # any later) opening then the open period is abandoned altogether, and the
    # next good shut time sought as start for next open period, but for the
    # purposes of identifying the nth open period, rejected ones must be counted
    # as an open period even though their length is undefined.

    opint = []
    oppro = []
    opamp = []
    shint = []
    shpro = []
    n = 0
    tint = 0
    prop = 0
    ampl = 0
    first = 0
    while n < len(otint):
        if oampl[n] != 0:
            tint = tint + otint[n]
            ampl = ampl + oampl[n]*otint[n]
            if oprops[n] == 8: prop = 8
            avamp = ampl / tint
            n = n + 1
            first = 1
        else:
            shint.append(otint[n])
            shpro.append(oprops[n])
            n = n + 1
            if first:
                opamp.append(avamp)
                avamp = 0
                opint.append(tint)
                tint = 0
                oppro.append(prop)
                prop = 0

    return opint, opamp, oppro, shint, shpro

def fitHistogram(tres):
    'Fits dwell time histogram.'

    xLow = tres
    xHigh = 10.0**xmax    # here have to decide if max observed or max displayed

	

if __name__ == '__main__':
    
    ftypes=(('SCN file', '*.scn'), ('scn file', '*.SCN'))
    scnFile = tkFileDialog.askopenfilename(filetypes=ftypes)

    ioffset = 0
    nint = 0

    if os.path.exists(scnFile):

        header = readHeader(scnFile)
        for eachKey in sorted(header):
            print eachKey, '=', header[eachKey]

        ioffset = header['ioffset']
        nint = header['nint']
        calfac2 = header['calfac2']
        #print ioffset
        #print nint
        itint, iampl, iprops = readData(scnFile, ioffset, nint, calfac2)
        for i in range(0, 10):
            print i+1, 'int=', itint[i], '; ampl=', iampl[i], '; op= ', iprops[i]

        print 'last: int=', itint[-1], '; ampl=', iampl[-1], 'op=', iprops[-1]

        ffilt = header['ffilt']    # kHz
        rms = header['rms']    # pA
        avamp = header['avamp'] * calfac2    # pA

        trise = 1000 * 0.3321 / (ffilt)    # in mikrosec
        
        #tres = float(raw_input('\n Resolution for open and shut times (microsec)= '))
        tres = 20    # get resolution in microsec
        tres1 = 0.001 * tres    # convert resolution to ms
        
        flsEvRate = falseEvents(tres1, ffilt, rms, avamp)
        zo = 1000. * tres1 / trise    # tres in ms, trise in mus
        aamaxo = erf(0.88604 * zo)
        print '\n At resolution %f microsec false event rate (per sec) for openings and shuttings ' %tres
        print ' = %f ( = %f risetimes, A/Amax = %f' %(flsEvRate, zo, aamaxo)
        
        otint, oampl, oprops = imposeResolution(tres1, itint, iampl, iprops)
        print '\n After imposing the resolution tres=', tres, 'of original', len(itint)
        print 'intervals were left', len(oampl)

        for i in range(0, 10):
            print i+1, 'int=', otint[i], '; ampl=', oampl[i], '; op= ', oprops[i]

        opint, opamp, oppro, shint, shpro = getOpenShutPeriods(otint, oampl, oprops)

        shintav = 0
        opintav = 0
        for i in range(0, len(opint)):
            shintav = shintav + shint[i]
            opintav = opintav + opint[i]
        shintav = shintav / len(opint)
        opintav = opintav / len(opint)

        print '\n number of open periods =', len(opint)
        print 'Average= ', opintav
        print 'Range:', min(opint), 'to', max(opint)
        print '\n number of shut times =', len(shint)
        print 'Average= ', shintav
        print 'Range:', min(shint), 'to', max(shint)

        #a = PlotHistogram(opint, tres1, 0)
        #b = PlotHistogram(shint, tres1, 1)

        nset = 1
        ncomp = [3]    # number of components in a set
        tau1 = [0.01, 0.1, 1]
        tau = []
        tau.append(tau1)    # mean of a component of a set
        tfix1 = [0, 0, 0]
        tfix = []
        tfix.append(tfix1)
        area1 = [0.8, 0.1, 0.1]
        area = []
        area.append(area1)
        afix1 = [0, 0, 0]
        afix = []
        afix.append(afix1)
        nfit = [len(shint)]     # number of points to be fitted
        ylow = [min(shint)]     # fitted range; lowest value
        yhigh = [max(shint)]    # fitted range; highest value
        yval = shint

        stpfac = 0.1
        errfac = 0.0001
        nevmax = 100

        #expf = ExpPdf(nset, ncomp, tau, tfix, area, afix, nfit, ylow, yhigh, yval, tres1)
        #simplex = Simplex(expf.theta, stpfac, errfac, nevmax, expf)
        #simplex = SimplexDC(expf.theta, expf)
        #simplex = SimplexNM(expf.theta, expf)

        #x, y = expf.calculatePlot()
        #print 'y', y
        #print 'x', x
        #b = PlotHistogram(shint, tres1, 1, x, y)

    else:
        print 'file not found'

    if os.path.exists(scnFile):
        pass
    else:
        print 'file not found'
