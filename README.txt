This is version 4 of all routines required to reduce GRIS data.

Changes with respect to version 3 are:

1.- Added keywords to main routine "gris": /dust and /maxshift

These keywords have been included to handle the image drift along the slit that happened during the first half of 2015. 

/dust takes into account that dust on the slit moves with the images, while pixel gain factor stays with the pixel.

/maxshift is the maximum shift along the slit allowed between the most-separated images. Default is 5 pixels

2.- Added keywords to main routine "gris": /checksync

This keyword has been introduced to check the synchronisation between the modulators and image acquisition. This syncronisation was lost during one campaign in 2015 because of a bad electrical conection in the camera electronics. 

3.- Positive-Q in reduced data points in the solar N-S direction. Previosly it pointed parallel to the elevation axis, which varies during the day in a fixed solar reference frame. Routines gris and gregormod modified.

4.- Routine get_flat has been modified for better fringe and bad pixel correction

5.- Teide observatory coordinates corrected and unified in the routines that use them (gregormod and r_frame_asp)

6.- Added some new bad pixels in routine badpixels2006

Routines modified:

- badpixels2006.pro 
- get_flat.pro
- gregormod.pro
- gris.pro
- r_frame_asp.pro

