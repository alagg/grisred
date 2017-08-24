V6 changes
==========

gris_v6 mainly handles some errors that have appeared with the derotator positioning angle, required for the polarimetric calibration and telescope correction.

The parameters of the mirrors have been updated and a keyword has been added:

1.- The parameters of the mirrors at the beginning of the 2017 campaign were pretty similar to the old values. The difference in the reduced data is hardly noticeable using 2017 or the previous set of parameters. Nevertheless, the main gris routine has been modified to use automatically the old or the new parameter set, depending on the date.

We have three new routines: parameter_mirrors_2016.pro, parameter_mirrors_2017.pro and gris_v6

2. The acquisition program (polar) stores in the header two values of the derotator angle: at the beginning and at the end of each map. This way it can be checked whether the derotator was tracking or not, depending on whether the values are different or not.

There were some time series for which the last value was missing and the reduction routine did not know whether tracking was on or off. It happened that the routine "decided" that in this case the default was "derotator stopped", and this was not true.  To solve this, a new keyword has been introduced in gris_v6: /rotator_track, to force tracking of the derotator during the reduction.

Routine(s) affected: gris_v6

3. Routine gregormod (which calculates the Mueller matrix of the trasnfer optics) can now handle any angle of the rotator (except for -999). Previosly it only handled angles below 360 degrees 

Routine(s) affected: gregormod_v6

With all this, there are four new routines:

- gris_v6
- gregormod_v6.pro
- parameter_mirrors_2016.pro
- parameter_mirrors_2017.pro

The main routine (gris_v6) can be used for the old data, just by changing gris_v5 by gris_v6 in the cal***.pro routine. There should be no effects.
