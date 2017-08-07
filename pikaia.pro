;
;  IDL program pikaia:
;
;  Designed to follow as closely as possible to Charbonneau and
;  Knapp.   For that reason, most of the comment blocks are 
;  intact, and Fortran definitions of variables remain in the form of
;  additional comments.  
;
;
;  The biggest change made is that the function being
;  evaluated is no longer external with its name ("ff" in pikaia.f)
;  passed throughout, but is an IDL standard function called "func"
;  in this IDL program.  The function "two_d" from  pikaia.f thus appears 
;  here as function "func", and the user can simply remove this 
;  example and substitute any other function desired, naming it "func".
;
;  Also, because IDL indexes from  0 -> n-1 rather than 1 -> n,
;  arrays have been changed accordingly.
;
;  The order of functions has been changed somewhat so that 
;  they are always above the function/procedure that calls them.
;
;  July, 1996
; 


;*********************************************************************
      FUNCTION ran2
;
;     Common block to make iseed visible to urand_init (and to save
;     it between calls)

      common share1,iseed
      common share2,iv,iy,idum2
;
;      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
;      REAL ran2,AM,EPS,RNMX
;      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
;   *IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,
;     *NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
;      INTEGER idum2,j,k,iv(NTAB),iy
;      SAVE iv,iy,idum2
;      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
;
     IM1=2147483563 
     IM2=2147483399 
     AM=1./IM1 
     IMM1=IM1-1
     IA1=40014 
     IA2=40692 
     IQ1=53668 
     IQ2=52774 
     IR1=12211 
     IR2=3791
     NTAB=32 
     NDIV=1+IMM1/NTAB 
     EPS=1.2e-7 
     RNMX=1.-EPS
;
     idum = iseed
;
      if (idum le 0) then begin
        idum= -idum > 1
        idum2=idum
        for j=NTAB+8,1,-1 do begin
          k=idum/IQ1
          idum=(IA1*(idum-k*IQ1)-k*IR1)
          if (idum lt 0) then idum=(idum+IM1)
          if (j le NTAB) then iv(j-1)=idum
	endfor
        iy=iv(0)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum lt 0) then idum=idum+IM1 
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2 lt 0) then idum2=idum2+IM2 
      j=1+iy/NDIV
      iy=iv(j-1)-idum2
      iv(j-1)=idum
      if(iy lt 1) then iy=iy+IMM1 
      ran2= AM*iy < RNMX
;
      iseed = (idum)
;
      return,ran2
      END
;  (C) Copr. 1986-92 Numerical Recipes Software .@)1.<BcL&1i4-130Rk#.
;   (adapted to IDL)
;*********************************************************************
function func_pikaia,n,x        ;a. lagg, renamed function from
                                ;'func'  to 'func_pikaia'
;
;   
;*************
;     previously called two_d
;
;***************
;     Compute sample fitness function (2-d landscape)
;
;     Input:
;     integer n
;     real x(n)
;
;     Output:
;     real two_d
;
;

      if  x(0) gt 1  or x(1) gt 1 then stop
      rr=sqrt( (0.5-x(0))^2+ (0.5-x(1))^2)
      ff=cos(rr*9.*!pi)
      ff=16.*x(0)*x(1)*(1.-x(0))*(1.-x(1))*ff*ff

      func = ff

      return,func
      end
;********************************************************************
      function urand

;     Common block to make iseed visible to urand_init (and to save
;     it between calls)

      common share1,iseed

;
;     Return the next pseudo-random deviate from a sequence which is
;     uniformly distributed in the interval [0,1]
;
;     Uses the function ran2, from Press et al, Numerical Recipes,
;     2nd ed., Cambridge Univ. Press, 1992.
;
;     Input - none
;
;     Output
;     real    urand
;
;     Local
;      integer  iseed
;      real     ran2
;      external ran2
;
 
;
;      can either use numerical recipes random number generator ran2, 
;      or idl random number generator urand
;
;      urand = randomu( iseed )
       urand = ran2()
      return,urand
      end

;******************************************************************
      pro indexx,n,arr,indx
      
      M=11
      NSTACK=50
;
      istack = lonarr(NSTACK)
;
      arr1 = fltarr(n+1)
      for i = 0,n-1 do arr1(i+1) = arr(i)
      indx1 = indgen(n+1)
      indx = lonarr(n)
;
      jstack=0
      l=1
      ir=n
;
      label:
;
      if ir-l lt M then begin
	for j = l+1,ir do begin
          indxt=indx1(j)
          a=arr1(indxt)
          for i=j-1,l,-1 do begin
            if arr1(indx1(i)) le a then goto,label2
            indx1(i+1)=indx1(i)
	  endfor
          i=l-1
	  label2:
          indx1(i+1)=indxt
	endfor
        if(jstack eq 0) then begin
	   for ll = 0,n-1 do indx(ll) = indx1(ll+1)
           return
	endif
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      endif else begin
        k=(l+ir)/2
        itemp=indx1(k)
        indx1(k)=indx1(l+1)
        indx1(l+1)=itemp
        if(arr1(indx1(l+1)) gt arr1(indx1(ir))) then begin
          itemp=indx1(l+1)
          indx1(l+1)=indx1(ir)
          indx1(ir)=itemp
        endif
        if(arr1(indx1(l)) gt arr1(indx1(ir))) then begin
          itemp=indx1(l)
          indx1(l)=indx1(ir)
          indx1(ir)=itemp
        endif
        if(arr1(indx1(l+1)) gt arr1(indx1(l))) then begin
          itemp=indx1(l+1)
          indx1(l+1)=indx1(l)
          indx1(l)=itemp
        endif
        i=l+1
        j=ir
        indxt=indx1(l)
        a=arr1(indxt)
        label3:
          i=i+1
        if(arr1(indx1(i)) lt a) then goto,label3
	label4:
          j=j-1
        if(arr1(indx1(j)) gt a) then goto,label4
        if(j lt i) then goto,label5
        itemp=indx1(i)
        indx1(i)=indx1(j)
        indx1(j)=itemp
        goto,label3
	label5:       
	indx1(l)=indx1(j)
        indx1(j)=indxt
        jstack=jstack+2
        if(jstack gt NSTACK) then print,'NSTACK too small in indexx'
        if(ir-i+1 ge j-1) then begin
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
	endif else begin
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endelse
      endelse
      goto,label
      END
;  (C) Copr. 1986-92 Numerical Recipes Software .@)1.<BcL&1i4-130Rk#.
;********************************************************************
       pro xpikaia

      common share2,iv,iy,idum2
;
;     Sample driver program for pikaia.f
;
     NTAB=32 
      idum2=123456789
      iv=lonarr(NTAB) 
      iv(*)= 0
      iy=0
      n=2
;
;     (two_d ("func" here) is an example fitness function, a 
;     smooth 2-d landscape)
;
;     First, initialize the random-number generator
;
      label:
      read, "Random number seed (I*4)? ",seed
      seed = -abs(long(seed))
      urand_init,seed
;
;     Set control variables (use defaults)
;
      ctrl = fltarr(12)
      ctrl(*) = -1
;
;     Now call pikaia
;
;**********************************************************************
;     note, I've had to change the call of all the procedures to not
;     include the function name, and the function throughout
;     the code is referred to as func. 
;**********************************************************************
;
        pikaia,n,ctrl,x,f,status
;
;     Print the results
      print," status: ",status
      print,"      x: ",x
      print,"      f: ",f
      print,"      ctrl: ",ctrl
;
      goto,label
      end
;*******************************************************************

      pro urand_init,urand_seed 

;     Common block to communicate with urand

      common share1,iseed
;
;     Initialize random number generator urand with given seed
;
;
;     Input
;     integer urand_seed
;
;     Output - none
;
;     Local
;     integer iseed
;
 
;     Set the seed value
      iseed = urand_seed

      return

      end
;
;********************************************************************
      pro setctl,ctrl,n,np,ngen,nd,pcross,pmutmn,pmutmx,pmut,imut,fdif,irep,ielite,ivrb,status


;     Set control variables and flags from input and defaults
;
;
;     Input
;     integer  n
;
;     Input/Output
;     real     ctrl(12)
;
;     Output
;     integer  np, ngen, nd, imut, irep, ielite, ivrb, status
;     real     pcross, pmutmn, pmutmx, pmut, fdif
;
;     Local
;     integer  i
;     real     DFAULT(12)
;     save     DFAULT
;     data     DFAULT /100,500,5,.85,2,.005,.0005,.25,1,1,1,0/

      DFAULT = [100,20,4,.85,2,.005,.0005,.25,1,1,1,1]
;
      test = where(ctrl lt 0)
      if min(test) gt -1 then ctrl(test) = DFAULT(test)
; 
      np = long(ctrl(0))
      ngen = long(ctrl(1))
      nd = long(ctrl(2))
      pcross = ctrl(3)
      imut = long(ctrl(4))
      pmut = ctrl(5)
      pmutmn = ctrl(6)
      pmutmx = ctrl(7)
      fdif = ctrl(8)
      irep = long(ctrl(9))
      ielite = long(ctrl(10))
      ivrb = long(ctrl(11))
      status = 0
;
;     Print a header
;
      if (ivrb gt 0) then begin
 
         f2 = '(/1x,60("*"),/'
         f2 = f2+'" *",13x,"PIKAIA Genetic Algorithm Report ",13x,"*",/'
	 f2 = f2 + '1x,60("*"),//,"   Number of Generations evolving: ",i4,/'
         f2 = f2 + '"       Individuals per generation: ",i4,/'
         f2 = f2 +'"    Number of Chromosome segments: ",i4,/'
	 f2 = f2 + '"    Length of Chromosome segments: ",i4,/'
	 f2 = f2 + '"            Crossover probability: ",f9.4,/'
	 f2 = f2 + '"            Initial mutation rate: ",f9.4,/'
	 f2 = f2 + '"            Minimum mutation rate: ",f9.4,/'
	 f2 = f2 + '"            Maximum mutation rate: ",f9.4,/'
	 f2 = f2 + '"    Relative fitness differential: ",f9.4)'
	print,format=f2,ngen,np,n,nd,pcross,pmut,pmutmn,pmutmx,fdif
;
        f3='("                    Mutation Mode: ",A)'
         if (imut eq 1) then print,format=f3, 'Constant'
         if (imut eq 2) then print,format=f3, 'Variable'
;
        f4 = '("                Reproduction Plan: ",A)'
         if (irep eq 1) then print,format=f4, 'Full generational replacement'
         if (irep eq 2) then print,format=f4, 'Steady-state-replace-random'
         if (irep eq 3) then print,format=f4, 'Steady-state-replace-worst'
;
      endif
 
;     Check some control values
      f10  = '" ERROR: illegal value for imut (ctrl(5))"'
      if (imut ne 1 and imut ne 2) then begin
         print,f10 
         status = 5
      endif
 
      f11 = '" ERROR: illegal value for fdif (ctrl(9))"'
      if (fdif gt 1) then begin
         print,f11
         status = 9
      endif
 
      f12 = '" ERROR: illegal value for irep (ctrl(10))"'
      if (irep ne 1  and  irep ne 2  and  irep ne 3) then begin
	 print,f12
         status = 10
      endif
 
      f13 = ' ERROR: illegal value for pcross (ctrl(4))'
      if (pcross gt 1.0  or  pcross lt 0 ) then begin
	 print,f13
         status = 4
      endif
 
      f14 = ' ERROR: illegal value for pcross (ctrl(4))'
      if (ielite ne 0  and  ielite ne 1) then begin
 	 print,f14
         status = 11
      endif
 
      f15 = '" WARNING: dangerously high value for pmut (ctrl(6));",/" (Should enforce elitism with ctrl(11)=1.)"'
      if (irep eq 1  and  imut eq 1  and  pmut gt 0.5  and ielite eq 0) then begin
	 print,f15
      endif
 
      f16 = '" WARNING: dangerously high value for pmutmx (ctrl(8));",/" (Should enforce elitism with ctrl(11)=1.)"'
      if (irep eq 1 and imut eq 2 and pmutmx gt 0.5 and ielite eq 0) then begin
	 print,f16
      endif
 
      f17 = ' WARNING: dangerously low value of fdif (ctrl(9))'
      if (fdif lt 0.33) then begin
	 print,f17
      endif
 
      return
      end
;********************************************************************
      pro report,ivrb,ndim,n,np,nd,oldph,fitns,ifit,pmut,ig,nnew

	common datstore,bestft,pmutpv
 
;     Write generation report to standard output
;
;     Input:
;     integer ifit(np),ivrb,ndim,n,np,nd,ig,nnew
;     real oldph(ndim,np),fitns(np),pmut
;
;     Output: none
;
;     Local
;     real bestft,pmutpv
;     save bestft,pmutpv
;     integer ndpwr,k
;     logical rpt
;     data bestft,pmutpv /0,0/
;
      rpt=0.
 
      if (pmut ne pmutpv) then begin
         pmutpv=pmut
         rpt=1.
      endif
 
      if (fitns(ifit(np-1)) ne bestft) then begin
         bestft=fitns(ifit(np-1))
         rpt=1.
      endif
 
      if (rpt ne 0.  or  ivrb ge 1) then begin
 
;        Power of 10 to make integer genotypes for display

         ndpwr = round(10.^nd)
         
                                ;changed verbosity (A. Lagg)
         if ivrb ge 2 then $
           print,format='(i6,i6,f12.6,4f12.6)', ig+1,nnew,pmut,fitns(ifit(np-1)), fitns(ifit(np-2)), fitns(ifit((np+1)/2 -1))
;          if ivrb eq 2 then begin
;             for k=0,n-1 do begin
;               print,format='(22x,3i10)',round(ndpwr*oldph(k,ifit(np-1))),round(ndpwr*oldph(k,ifit(np-2))),round(ndpwr*oldph(k,ifit((np+1)/2 -1)))
;             endfor
;          endif
 
      endif
      end

;**********************************************************************
;                         GENETICS MODULE
;**********************************************************************
;
;     ENCODE:    encodes phenotype into genotype
;                called by: PIKAIA
;
;     DECODE:    decodes genotype into phenotype
;                called by: PIKAIA
;
;     CROSS:     Breeds two offspring from two parents
;                called by: PIKAIA
;
;     MUTATE:    Introduces random mutation in a genotype
;                called by: PIKAIA
;
;     ADJMUT:    Implements variable mutation rate
;                called by: PIKAIA
;
;**********************************************************************
      pro encode,n,nd,ph,gn
;======================================================================
;     encode phenotype parameters into integer genotype
;     ph(k) are x,y coordinates [ 0 < x,y < 1 ]
;======================================================================
;
;
;     Inputs:
;     integer   n, nd
;     real      ph(n)
;
;     Output:
;     integer   gn(n*nd)
;
;     Local:
;     integer   ip, i, j, ii
;     real      z
;
      z=10.^nd
      ii=0
      for i=0,n-1 do begin
         ip = long(ph(i)*z)
         for  j=nd-1,0,-1 do begin
            gn(ii+j)=ip mod 10
            ip=ip/10
	endfor
        ii=ii+nd
       endfor	
 
      return
      end
 
;**********************************************************************
      pro decode,n,nd,gn,ph
;======================================================================
;     decode genotype into phenotype parameters
;     ph(k) are x,y coordinates [ 0 < x,y < 1 ]
;======================================================================
;
;
;     Inputs:
;     integer   n, nd, gn(n*nd)
;
;     Output:
;     real      ph(n)
;
;     Local:
;     integer   ip, i, j, ii
;     real      z
;
      z=10.^(-nd)
      ii=0
      for i=0,n-1 do begin
         ip=0
         for j=0,nd-1 do begin	
            ip=10*ip+gn(ii+j)
	 endfor
         ph(i)=ip*z
         ii=ii+nd
      endfor
 
      return
      end
 
;**********************************************************************
      pro cross,n,nd,pcross,gn1,gn2
;======================================================================
;     breeds two parent chromosomes into two offspring chromosomes
;     breeding occurs through crossover starting at position ispl
;======================================================================
;
;     Inputs:
;     integer        n, nd
;     real           pcross
;
;     Input/Output:
;     integer        gn1(n*nd), gn2(n*nd)
;
;     Local:
;     integer        i, ispl, t
;
;     Function
;     real           urand
;     external       urand
;
 
;     Use crossover probability to decide whether a crossover occurs

      if (urand() lt pcross) then begin
 
;        Compute crossover point

         ispl=long(urand()*n*nd)+1
 
;        Swap genes at ispl and above

         for i=ispl-1,n*nd-1 do begin
            t=gn2(i)
            gn2(i)=gn1(i)
            gn1(i)=t
	 endfor
      endif
 
      return
      end
 
;**********************************************************************
      pro mutate,n,nd,pmut,gn
;======================================================================
;     Mutations occur at rate pmut at all gene loci
;======================================================================
;
;
;     Input:
;     integer        n, nd
;     real           pmut
;
;     Input/Output:
;     integer        gn(n*nd)
;
;     Local:
;     integer        i
;
;     Function:
;     real           urand
;     external       urand
;
;     Subject each locus to mutation at the rate pmut

      for i=0,n*nd-1 do begin
         if urand() lt pmut then begin
            gn(i)=long(urand()*10.)
         endif
      endfor
 
      return
      end
 
;**********************************************************************
      pro adjmut,np,fitns,ifit,pmutmn,pmutmx,pmut
;======================================================================
;     dynamical adjustment of mutation rate; criterion is relative
;     difference in absolute fitnesses of best and median individuals
;======================================================================
;
;     implicit none
;
;     Input:
;     integer        np, ifit(np)
;     real           fitns(np), pmutmn, pmutmx
;
;     Input/Output:
;     real           pmut
;
;     Local:
;     real           rdif, rdiflo, rdifhi, delta
;     parameter      (rdiflo=0.05, rdifhi=0.25, delta=1.5)
           
      rdiflo=0.05
      rdifhi=0.25  
      delta=1.5
 
      rdif=abs(fitns(ifit(np-1))-fitns(ifit(np/2-1)))/(fitns(ifit(np-1))+fitns(ifit(np/2-1)))
      if rdif le rdiflo then pmut=pmutmx<pmut*delta
      if rdif ge rdifhi then pmut=pmutmn>pmut/delta
 
      return
      end

;**********************************************************************
;                       REPRODUCTION MODULE
;**********************************************************************
;
;     SELECT:   Parent selection by roulette wheel algorithm
;               called by: PIKAIA
;
;     RNKPOP:   Ranks initial population
;               called by: PIKAIA, NEWPOP
;
;     GENREP:   Inserts offspring into population, for full
;               generational replacement
;               called by: PIKAIA
;
;     STDREP:   Inserts offspring into population, for steady-state
;               reproduction
;               called by: PIKAIA
;               calls:     FF
;
;     NEWPOP:   Replaces old generation with new generation
;               called by: PIKAIA
;               calls:     FF, RNKPOP
;
;**********************************************************************
      pro select,np,jfit,fdif,idad
;======================================================================
;     Selects a parent from the population, using roulette wheel
;     algorithm with the relative fitnesses of the phenotypes as
;     the "hit" probabilities [see Davis 1991, chap. 1].
;======================================================================
;
;     implicit none
;
;     Input:
;     integer        np, jfit(np)
;     real           fdif
;
;     Output:
;     integer        idad
;
;     Local:
;     integer        np1, i
;     real           dice, rtfit
;
;     Function:
;     real           urand
;     external       urand
;
;
      np1 = np+1
      dice = urand()*np*np1
      rtfit = 0.
      for i=0,np-1 do begin
         rtfit = rtfit+np1+fdif*(np1-2*(jfit(i)+1))
         if (rtfit ge dice) then begin
            idad=i
            goto,label22
         endif
      endfor
;     Assert: loop will never exit by falling through
 
    label22: 
    return
      end
 
;**********************************************************************
      pro rnkpop,n,arrin,indx,rank
;======================================================================
;     Calls external sort routine to produce key index and rank order
;     of input array arrin (which is not altered).
;======================================================================
;
;     implicit none
;
;     Input
;     integer    n
;     real       arrin(n)
;
;     Output
;     integer    indx(n),rank(n)
      indx = lonarr(n)
      rank = lonarr(n)
;
;     Local
;     integer    i
;
;     Numerical Recipes subroutine
;     external indexx
;
;
;     Compute the key index
      indexx,n,arrin,indx
;
      indx = indx-1
;     ...and the rank order
      for  i=0,n-1 do begin
         rank(indx(i)) = n-i-1
      endfor
      return
      end
 
;***********************************************************************
      pro genrep,ndim,n,np,ip,ph,newph
;=======================================================================
;     full generational replacement: accumulate offspring into new
;     population array
;=======================================================================
;
;     implicit none
;
;     Input:
;     integer        ndim, n, np, ip
;     real           ph(ndim,2)
;
;     Output:
;     real           newph(ndim,np)
;
;     Local:
;     integer        i1, i2, k
;
;
;     Insert one offspring pair into new population
      i1=2*ip
      i2=i1+1
      for  k=0,n-1 do begin
         newph(k,i1)=ph(k,0)
         newph(k,i2)=ph(k,1)
      endfor
      return
      end
 
;**********************************************************************
      pro stdrep,ndim,n,np,irep,ielite,ph,oldph,fitns,ifit,jfit,nnew
;======================================================================
;     steady-state reproduction: insert offspring pair into population
;     only if they are fit enough (replace-random if irep=1 or
;     replace-worst if irep=2).
;======================================================================
;
;     implicit none
;
;     Input:
;     integer        ndim, n, np, irep, ielite
;     real           ff, ph(ndim,2)
;     external       ff
;
;     Input/Output:
;     real           oldph(ndim,np), fitns(np)
;     integer        ifit(np), jfit(np)
;
;     Output:
;     integer        nnew
;
;     Local:
;     integer        i, j, k, i1, if1
;     real           fit
;
;     External function
;     real           urand
;     external       urand
;
;
      nnew = 0
      for j=0,1 do begin
 
;        1. compute offspring fitness (with caller's fitness function)
         fit=func(n,ph(0:n-1,j))
 
;        2. if fit enough, insert in population
         for i=np-1,0,-1 do begin
            if (fit gt fitns(ifit(i))) then begin
 
;              make sure the phenotype is not already in the population
               if (i lt np-1) then begin
                  for  k=0,n-1 do begin
                     if (oldph(k,ifit(i+1)) ne ph(k,j)) then goto,label6
		  endfor
                  goto,label11
		label6:
               endif
 
;              offspring is fit enough for insertion, and is unique
 
;              (i) insert phenotype at appropriate place in population
               if (irep eq 3) then i1=0 else begin
		 if (ielite eq 0 or i eq np-1) then begin
			 i1=long(urand()*np)
		 endif else begin
			 i1=long(urand()*(np-1))
		 endelse
	       endelse
               if1 = ifit(i1)
               fitns(if1)=fit
	       oldph(*,if1) = ph(*,j)
;
;              (ii) shift and update ranking arrays
               if (i lt i1) then begin
 
;                 shift up
                  jfit(if1)=np-i-2
                  for k=i1-1,i+1,-1 do begin
                     jfit(ifit(k))=jfit(ifit(k))-1
                     ifit(k+1)=ifit(k)
		  endfor
                  ifit(i+1)=if1
	       endif else begin
 
;                 shift down
                  jfit(if1)=np-i-1
                  for k=i1+1,i do begin
                     jfit(ifit(k))=jfit(ifit(k))+1
                     ifit(k-1)=ifit(k)
		  endfor
                  ifit(i)=if1
               endelse
               nnew = nnew+1
               goto,label11
            endif
	endfor
 
    label11: 
    endfor
 
    return
    end
 
;**********************************************************************
      pro newpop,ielite,ndim,n,np,oldph,newph,ifit,jfit,fitns,nnew
;======================================================================
;     replaces old population by new; recomputes fitnesses & ranks
;======================================================================
;
;     implicit none
;
;     Input:
;     integer        ndim, np, n, ielite
;     real           ff
;     external       ff
;
;     Input/Output:
;     real           oldph(ndim,np), newph(ndim,np)
;
;     Output:
;     integer        ifit(np), jfit(np), nnew
;      real           fitns(np)
;
;     Local:
;     integer        i, k
;
;
      nnew = np
 
;     if using elitism, introduce in new population fittest of old
;     population (if greater than fitness of the individual it is
;     to replace)
      if (ielite eq 1 and func(n,newph(0:n-1,0)) lt fitns(ifit(np-1))) then begin
         for k=0,n-1 do begin
            newph(k,0)=oldph(k,ifit(np-1))
	 endfor
         nnew = nnew-1
      endif
 
;     replace population
;     and get fitness using caller's fitness function
;     but be careful of zeros
      for i = 0,np-1 do begin
	 if newph(0,i) ne 0. then begin
		oldph(*,i) = newph(*,i)
		fitns(i)=func(n,oldph(0:n-1,i))
	 endif
      endfor
 
;     compute new population fitness rank order
      rnkpop,np,fitns,ifit,jfit
 
      return
      end
      
;A. Lagg: put pikaia routine to the end of the file (Compilation!))
;*************************************************************************
       pro pikaia,n,ctrl,x,xall,f,status,init=init
;
      common share1,iseed
      common datstore,bestft,pmutpv
;
;     Optimization (maximization) of user-supplied "fitness" function
;     ff  over n-dimensional parameter space  x  using a basic genetic
;     algorithm method.
;
;********************
;     note, ff is func at the top
;
;********************
;
;     Paul Charbonneau & Barry Knapp
;     High Altitude Observatory
;     National Center for Atmospheric Research
;     Boulder CO 80307-3000
;     <paulchar@hao.ucar.edu>
;     <knapp@hao.ucar.edu>
;
;     Version of 1995 April 13
;
;	(IDL version 1996 July - comments kept intact)
;
;     Genetic algorithms are heuristic search techniques that
;     incorporate in a computational setting, the biological notion
;     of evolution by means of natural selection.  This subroutine
;     implements the three basic operations of selection, crossover,
;     and mutation, operating on "genotypes" encoded as strings.
;
;     References:
;
;        Charbonneau, Paul.  "Genetic Algorithms in Astronomy and
;           Astrophysics."  Astrophysical J. (Supplement), vol 101,
;           in press (December 1995).
;
;        Goldberg, David E.  Genetic Algorithms in Search, Optimization,
;           & Machine Learning.  Addison-Wesley, 1989.
;
;        Davis, Lawrence, ed.  Handbook of Genetic Algorithms.
;           Van Nostrand Reinhold, 1991.
;
;
;     implicit none
;
;     Input:
;     integer   n
;     real      ff
;     external  ff
;
;      o Integer  n  is the parameter space dimension, i.e., the number
;        of adjustable parameters.  (Also the number of chromosomes.)
;
;      o Function  ff  is a user-supplied scalar function of n vari-
;        ables, which must have the calling sequence f = ff(n,x), where
;        x is a real parameter array of length n.  This function must
;        be written so as to bound all parameters to the interval [0,1];
;        that is, the user must determine a priori bounds for the para-
;        meter space, and ff must use these bounds to perform the appro-
;        priate scalings to recover true parameter values in the
;        a priori ranges.
;
;        By convention, ff should return higher values for more optimal
;        parameter values (i.e., individuals which are more "fit").
;        For example, in fitting a function through data points, ff
;        could return the inverse of chi**2.
;
;        In most cases initialization code will have to be written
;        (either in a driver or in a separate subroutine) which loads
;        in data values and communicates with ff via one or more labeled
;        common blocks.  An example exercise driver and fitness function
;        are provided in the accompanying file, xpikaia.f.
;
;
;      Input/Output:
;      real ctrl(12)
;
;      o Array  ctrl  is an array of control flags and parameters, to
;        control the genetic behavior of the algorithm, and also printed
;        output.  A default value will be used for any control variable
;        which is supplied with a value less than zero.  On exit, ctrl
;        contains the actual values used as control variables.  The
;        elements of ctrl and their defaults are:
;
;           ctrl( 1) - number of individuals in a population (default
;                      is 100)
;           ctrl( 2) - number of generations over which solution is
;                      to evolve (default is 500)
;           ctrl( 3) - number of significant digits (i.e., number of
;                      genes) retained in chromosomal encoding (default
;                      is 6)  (Note: This number is limited by the
;                      machine floating point precision.  Most 32-bit
;                      floating point representations have only 6 full
;                      digits of precision.  To achieve greater preci-
;                      sion this routine could be converted to double
;                      precision, but note that this would also require
;                      a double precision random number generator, which
;                      likely would not have more than 9 digits of
;                      precision if it used 4-byte integers internally.)
;           ctrl( 4) - crossover probability; must be  <= 1.0 (default
;                      is 0.85)
;           ctrl( 5) - mutation mode; 1/2=steady/variable (default is 2)
;           ctrl( 6) - initial mutation rate; should be small (default
;                      is 0.005) (Note: the mutation rate is the proba-
;                      bility that any one gene locus will mutate in
;                      any one generation.)
;           ctrl( 7) - minimum mutation rate; must be >= 0.0 (default
;                      is 0.0005)
;           ctrl( 8) - maximum mutation rate; must be <= 1.0 (default
;                      is 0.25)
;           ctrl( 9) - relative fitness differential; range from 0
;                      (none) to 1 (maximum).  (default is 1.)
;           ctrl(10) - reproduction plan; 1/2/3=Full generational
;                      replacement/Steady-state-replace-random/Steady-
;                      state-replace-worst (default is 3)
;           ctrl(11) - elitism flag; 0/1=off/on (default is 0)
;                      (Applies only to reproduction plans 1 and 2)
;           ctrl(12) - printed output 0/1/2=None/Minimal/Verbose
;                      (default is 0)
;
;
;     Output:
;     real      x(n), f
      
                                ;just do the compilation of the routines
      if keyword_set(init) then return
      
	x = fltarr(n)

;     integer   status
;
;      o Array  x(1:n)  is the "fittest" (optimal) solution found,
;         i.e., the solution which maximizes fitness function ff
;
;      o Scalar  f  is the value of the fitness function at x
;
;      o Integer  status  is an indicator of the success or failure
;         of the call to pikaia (0=success; non-zero=failure)
;
;
;     Constants
;     integer   NMAX, PMAX, DMAX
;

      ;; NMAX = 96
      ;; PMAX = 512
      ;; DMAX = 8
      NMAX = 96*2
      PMAX = 512*2
      DMAX = 8
;
;      o NMAX is the maximum number of adjustable parameters
;        (n <= NMAX)
;
;      o PMAX is the maximum population (ctrl(1) <= PMAX)
;
;      o DMAX is the maximum number of Genes (digits) per Chromosome
;        segement (parameter) (ctrl(3) <= DMAX)
;
;
;     Local variables
;     integer        np, nd, ngen, imut, irep, ielite, ivrb, k, ip, ig,
;    +               ip1, ip2, new, newtot
;     real           pcross, pmut, pmutmn, pmutmx, fdif
;
;     real           ph(NMAX,2), oldph(NMAX,PMAX), newph(NMAX,PMAX)

	ph = fltarr(nmax,2)
	oldph = fltarr(nmax,pmax)
	newph = fltarr(nmax,pmax)
	xall = fltarr(nmax,pmax)

;
;     integer        gn1(NMAX*DMAX), gn2(NMAX*DMAX)
;     integer        ifit(PMAX), jfit(PMAX)

	gn1 = lonarr(nmax*dmax)
	gn2 = lonarr(nmax*dmax)
	ifit = lonarr(pmax)
	jfit = lonarr(pmax)

;     real           fitns(PMAX)

	fitns = fltarr(pmax)

;
;
;     User-supplied uniform random number generator
;     real           urand
;     external       urand
;
;     Function urand should not take any arguments.  If the user wishes
;     to be able to initialize urand, so that the same sequence of
;     random numbers can be repeated, this capability could be imple-
;     mented with a separate subroutine, and called from the user's
;     driver program.  An example urand function (and initialization
;     subroutine) which uses the function ran2 from Press, et al,
;     Numerical Recipes, 2nd ed., Cambridge Univ Press, 1992, is
;     provided in the accompanying file, xpikaia.f.
;
;
;  rather than do a data block, I will zero these out now,
;  and keep them in a common block shared between pro pikaia
;  and pro report
;
      pmutpv = 0.
      bestft = 0.
;
;     Set control variables from input and defaults


      setctl,ctrl,n,np,ngen,nd,pcross,pmutmn,pmutmx,pmut,imut, fdif,irep,ielite,ivrb,status

      if (status  ne  0) then begin
          print," Control vector (ctrl) argument(s) invalid"
         return
      endif
 
;     Make sure locally-dimensioned arrays are big enough

      if (n gt NMAX  or  np gt PMAX  or  nd gt DMAX) then begin
         print," Number of parameters, population, or genes too large"
         status = -1
         return
      endif
 
;     Compute initial (random but bounded) phenotypes

      for ip=0,np-1 do begin
         for k=0,n-1 do begin
            oldph(k,ip)=urand()
	 endfor
         fitns(ip) = func(n,oldph(0:n-1,ip))
      endfor	
 
;     Rank initial population by fitness order

      rnkpop,np,fitns,ifit,jfit
 
 
;     Main Generation Loop

      for ig=0l,ngen-1 do begin
 
;        Main Population Loop

         newtot=0
         for ip=0,np/2-1 do begin
 
;           1. pick two parents

            select,np,jfit,fdif,ip1
	    label21:
            select,np,jfit,fdif,ip2
            if (ip1 eq ip2) then goto,label21
 
;           2. encode parent phenotypes

            encode,n,nd,oldph(0:n-1,ip1),gn1
            encode,n,nd,oldph(0:n-1,ip2),gn2
 
;           3. breed

            cross,n,nd,pcross,gn1,gn2
            mutate,n,nd,pmut,gn1
            mutate,n,nd,pmut,gn2
 
;           4. decode offspring genotypes

	    phsub = fltarr(n)
	    phsub = ph(0:n-1,0)
            decode,n,nd,gn1,phsub
	    ph(0:n-1,0) = phsub(*)
	    phsub = ph(0:n-1,1)
            decode,n,nd,gn2,phsub
	    ph(0:n-1,1) = phsub(*)

;           5. insert into population

            if (irep eq 1) then begin
               genrep,NMAX,n,np,ip,ph,newph
	    endif else begin
               stdrep,NMAX,n,np,irep,ielite,ph,oldph,fitns,ifit,jfit,new
               newtot = newtot+new
            endelse
 
;        End of Main Population Loop
 	 endfor 
 
;        if running full generational replacement: swap populations

         if (irep eq 1) then newpop,ielite,NMAX,n,np,oldph,newph,ifit,jfit,fitns,newtot
 
;        adjust mutation rate?

         if (imut eq 2) then adjmut,np,fitns,ifit,pmutmn,pmutmx,pmut
;
;        print generation report to standard output?

         if (ivrb gt 0) then report,ivrb,NMAX,n,np,nd,oldph,fitns,ifit,pmut,ig,newtot
 
;     End of Main Generation Loop
   endfor
;
;     Return best phenotype and its fitness

      x(*) = oldph(0:n-1,ifit(np-1))
      xall=oldph
      f = fitns(ifit(np-1))
;
    end
