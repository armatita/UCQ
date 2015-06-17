# UCQ
UCQ (upscale by convolution quality) is a simple alternative algorithm (very experimental) to upscale well data considering how well it will relate with the seismic data behind it. User needs to input well data (original and upscaled), seismic data (numpy binary), and wavelet (2 columns text file) and the algorithm will return the upscale.

## How does it work
It works by building bins for every depth step (let's say 4 ms) and iteration by iteration switch a value from the original well (in the range of 4 ms) to the upscaled. If the wavelet convolution gives a better result the switch is accepted otherwise moves on. In the end you should have a considerable better upscaling using the seismic convergence criteria.

## What is loaded
You have a few functions inside it: initial_stats() and load_stuff() that are loading the necessary data to run the algorithm. The code is a mess at this point but should be simple enough for python programmer (among others) to understand.

## What are LAS files loaded
LAS is a file format that generally looks like this:

```
~VERSION INFORMATION
VERS.                  2.0 : CWLS LOG ASCII STANDARD VERSION 2.0
WRAP.                  NO  : DATA IS NOT WRAPPED
~Well Information Section
#MNEM.UNIT          VALUE/NAME                 DESCRIPTION
#----------------------------------------------------------------------
STRT.m               1475.48681 : START DEPTH
STOP.m               5259.37288 : STOP DEPTH
STEP.m              0.152400002 : STEP
NULL.                -999.25 : NULL VALUE
COMP.                                             : COMPANY
WELL.    ABCDEFGH                                 : WELL
FLD .                                             : FIELD
LOC .         xxxxxxx,      xxxxxxx               : LOCATION
CNTY.                                             : COUNTY
STAT.                                             : STATE
CTRY.                                             : COUNTRY
SRVC.    Fugro-Jason                              : SERVICE
DATE.                                             : DATE
API .                                             : API
~Curve Information Section
#MNEM.UNIT          API CODE                   DESCRIPTION
#----------------------------------------------------------------------
DEPTH.m                          : * MD (TRACK/INDEX)
TVD  .m                          : * TVD (TRACK/INDEX)
TVDSS.m                          : * TVDSS (TRACK/INDEX)
DEVX .m                          : TRACK X-COORD
DEVY .m                          : TRACK Y-COORD
TIME .ms                         : TIME 
PHIE .none                       : POROSITY EFFECTIVE
DTS  .us/ft                      : S-SONIC
IS   .g/cm^3*m/s                 : S-IMPEDANCE
DT   .us/ft                      : P-SONIC
AI   .g/cm^3*m/s                 : P-IMPEDANCE
GR   .gAPI                       : GAMMA RAY
RHOB .g/cm^3                     : DENSITY
~Parameter Information Section
#MNEM.UNIT          API CODE                   DESCRIPTION
#----------------------------------------------------------------------
XCRD .m                 1234567              : WELL POSITION X-COORD
YCRD .m                 1234567              : WELL POSITION Y-COORD
DATUM.m                      75.000          : WELL DATUM KB
#----------------------------------------------------------------------
~A DEPTH        TVD            TVDSS          DEVX           DEVY           TIME           PHIE           DTS            IS             DT             AI             GR             RHOB           
1075.4868      1075.4868      1050.4868      6122148.2429   7133475.3362   2815.4884      -999.25        -999.25        -999.25        -999.25        -999.25        -999.25        -999.25        
1075.6392      1075.6392      1050.6392      6122148.2429   7133475.3362   2815.5563      -999.25        -999.25        -999.25        -999.25        -999.25        -999.25        -999.25    
```

The algorithm loads this as original files (I would upload the originals to github for you to give it a try but they're confidential information).

## What's a wavelet file
This is an example of a wavelet file:
```
-4	-150410
-3	-188835
-2	-109970
-1	 121088
 0	 401148
 1	 545489
 2	 451922
 3	 191310
 4	-54597.1
```
 First column is just the step (ignore it but put it there), the second the actual wavelet signal. The wavelet size is actually the number of points until you get to 0 in the wavelet step. In our example is 4.
 
 ## What is a well log file
 A well log file is a file with simple X,Y,Z,AI columned information such as this:
 
      326.66667       144.66667         0.25000      9278.43000
      326.66667       144.66667         1.25000      9278.24000
      326.66667       144.66667         2.25000      9278.05000
      326.66667       144.66667         3.25000      9277.86000
      326.66667       144.66667         4.25000      9277.66000
      326.33333       144.66667         5.25000      9277.47000
      326.33333       144.66667         6.25000      9277.28000
      326.33333       144.66667         7.25000      9277.09000
      326.33333       144.66667         8.25000      9276.89000
      326.33333       144.66667         9.25000     13694.50000
      326.33333       144.66667        10.25000     14945.32000
      
Please notice that the well should come in regular GRID coordinates so x=326 is actually cell row=326 in the seismic. Not cell size or first coordinate is being considered for this version of the algorithm. All of this is important because you'll probably have to transform your data to use this algorithm.

## What is a seismic file
The seismic file is binary and I could build an analogue in Python using the following code:
```Python
import numpy as np
g = np.zeros((100,100,30),dtype='float32') # So building a 3D matrix of type float32 with number of nodes "x,y,x" of (100,100,30).
np.save('zero_matrix.npy,g)                # I'm just saving a 3D cube with zeros inside it to a binary file. You should actually populate it with your seismic values.
```

## Case study
I've run this upscale algorithm in a real well log and compared it with the typical upscale made by the client (and used as base information by this algorithm. The quasi-correlation went from 64 % to 97 % which is an amazing result. The comparison can be seen in the following figure.

![alt tag](https://lh4.googleusercontent.com/PhFDkdPGJoqsIG2DFcQ4BahkkfhOC56pmMkosC8ZrgOBSTD0jNVsHs08aSLeKxhA3V-qXVEUWapT9l4=w1117-h645)
