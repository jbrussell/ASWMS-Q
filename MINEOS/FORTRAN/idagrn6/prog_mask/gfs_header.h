      integer*4 num_pts,data_time(5),event_time(5),num_poles,num_zeros
      integer*4 djulian, ejulian, flag(10)
      real*4 dig, elat, elon, slat, slon, selev, dsec, esec
      real*4 poles(60), zeros(40), scal(20)
      real*4 moment(6), cal_const, magnit, delta, depth
      character statn*4, stype*4, comp*4, catch*8, model*12, source*24, 
     & synid*12
      character comment*80
      real*4 hdr(200)
      common /head$/ num_pts, statn, stype, comp, catch, data_time,
     &               dsec, djulian, dig, elat, elon, depth, magnit,
     &               delta, slat, slon, selev, event_time, esec,
     &               ejulian, cal_const, num_poles, num_zeros, poles,
     &               moment, flag, model, source, synid, zeros, scal, 
     &               comment
      equivalence (hdr(1), num_pts)
c
c num_pts     1,  1   integer*4        number of points in timeseries
c statn       2,  2   char*4           station name, such as KONO
c stype       3,  3   char*4           station type, such as ASRO or DWWS
c comp        4,  4   char*4           component, such as LPZ or SPE
c catch       5,  6   char*8           'catch phrase' descriptor, such as 'raw data' 
c data_time   7, 11   integer*4        integer array with data start time
c                                      data_time(1) - year
c                                      data_time(2) - month
c                                      data_time(3) - day 
c                                      data_time(4) - hour
c                                      data_time(5) - minute
c dsec        12,12   real*4           secs of data start time
c djulian     13,13   integer*4        day of year for data start time
c dig         14,14   real*4           !!! sampling rate in sec !!!!
c elat        15,15   real*4           event latitude (degrees)
c elon        16,16   real*4           event longitude (degrees)
c depth       17,17   real*4           event depth (km)
c magnit      18,18   real*4           event magnitude
c delta       19,19   real*4           event-station distance in degrees
c slat        20,20   real*4           station latitude (degrees) 
c slon        21,21   real*4           station longitude (degrees)
c selev       22,22   real*4           station elevation (meters)
c event_time  23,27   integer*4        integer array with event origin time
c                                      event_time(1) - year
c                                      event_time(2) - month
c                                      event_time(3) - day 
c                                      event_time(4) - hour
c                                      event_time(5) - minute
c esec        28,28   real*4           secs of origin time
c ejulian     29,29   integer*4        day of year for origin time
c cal_const   30,30   real*4           a0 in the gdsn terminology (a0*scal(1) = scale for complex instrument response)
c num_poles   31,31   integer*4        number of poles in the complex instrument response
c num_zeros   32,32   integer*4        number of zeros in the complex instrument response
c poles       33,92   real*4           real array storing complex poles (real poles(60) = complex cpoles(30))
c moment      93,98   real*4           six components of the event moment tensor (Lind assumes 27 dyne-cm scaling)
c                                      moment(1) = Mrr  
c                                      moment(2) = Mtt
c                                      moment(3) = Mpp
c                                      moment(4) = Mrt
c                                      moment(5) = Mrp
c                                      moment(6) = Mtp
c flag        99,108  integer*4        array for keeping track of data flags (0 = ok)
c                                      flag(1) - units of instrument response (as stored) (0 = du/m; 1 = du/m/s; 2 = du/m/s2)
c                                      flag(2) - units of "raw" response (same as above)
c                                      flag(3) - units of current data (0 - 
c                                      flag(4) - units of displacement (0=m, 3=mm, 6=microns) 
c                                                can be used as a conversion factor -- 
c						 10**flag(4) units/m  -- JBG 2/93
c                                      flag(5-10) - unused at this time
c model       109,111 char*12          character string describing model used for synthetics
c source      112,117 char*24          character string describing sourrce region (from CMT catalog)
c synid       118,120 char*12          character string describing synthetic type, such as 'mode'
c zeros       121,160 real*4           real array storing complex zeros (real zeros(40) = complex czeros(20))
c scal        161,180 real*4           real array storing useful stuff
c                                      scal(1) - avg digital sensitivity (for instrument response)
c                                      scal(2) - instrument dip - lsg
c                                      scal(3) - instrument azimuth - lsg
c                                      scal(13) - (formerly scal(2) -- component's digital sensitivity
c                                      scal(14) - frequency of measurement for scal(13)
c                                      Lind Gee specific values:
c                                      scal(4-7) - filter parameters 
c                                      scal(8-10) - ray synthetic values (ray parameter, q, travel time of particular ray synthetic)
c                                      scal(11) - azimuth from source to receiver - lsg addition 7/14/92
c                                      scal(12) - source half duration (in sec) - lsg addition 10/20/92
c                                      scal(15-20) - unused at this time
c comment     181,200 char*80          character string used to describe the path of data/synthetic

