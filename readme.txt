Mail of Manolo, May 23 2014:

Dear all,

during the last days/weeks I have been working on the improvement of the reduction of the GRIS data.

I attach the latest version of the routines in the tar file. For those who are not aware of it, the main routine now is called "gris", not "acum2..." anymore. Its use follows the same philosophy as in the past.

There is now a routine "get_teles", which derives the polarization parameters of the mirrors M4, M5, M6 and M7 using several telescope calibration files.

I attach the reduction script for May 3rd, so that you can see how to use it. The output parameter xtau includes these polarization parameters of mirrors M4, M5, M6 and M7. I have checked some combinations of calibration files, using files of the same day or from different days, and I get very similar parameters (see the script). We will need to check how stable these numbers are with time and see whether can skip this calculation from one day to the next. If we are confident on the xtau values, we can just give the numbers (i.e., uncommment the line "xtau=....." , and comment the "get_teles" command) before doing the reduction itself.

I have reduced the data from May 3rd and I am satisfied with the results. We can apply it to the rest of the days.

I have not checked this yet with the 10830 data. I will do it in the coming days. It might be that some "fine tuning" is required for this wavelength.

Cheers,
Manolo 
