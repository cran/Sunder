***********************************************************
*     Main Fortran subroutine

      subroutine tangle (nsite,nloc,nal,nalmax,g,DG,DE,theta,thetacur,
     &     thetaprop,nit,thinning,fcur,fprop,xcur,xprop,ycur,yprop,zcur,
     &     zprop,Sigmacur,Sigmaprop,Uchol,am,bdm,bem,gm,dm,zeta,fsav,ud,
     &     zsav,swi)
      implicit none
*     computing options
      integer nit,thinning, swi
*     data 
      integer nsite,nloc,nal,nalmax,g
      double precision s,DG,DE
      double precision am, bdm, bem,gm,dm
      dimension  s(2,nsite),g(nsite,nloc,nalmax),
     &     nal(nloc),DG(nsite,nsite),DE(nsite,nsite)
*     parameters/variables to be simulated
      double precision theta,thetacur,thetaprop,zeta,fsav,zsav
      dimension theta(5,nit/thinning),thetacur(5),thetaprop(5),
     & zeta(5,nit/thinning), fsav(nit/thinning,nsite,nloc,nalmax),
     & zsav(nit/thinning,nsite,nloc,nalmax)
*     hidden variables
      double precision fcur,fprop,xcur,xprop,ycur,yprop,zcur,zprop,
     &     Sigmacur,Sigmaprop,Uchol
      dimension fcur(nsite,nloc,nalmax),fprop(nsite,nloc,nalmax),
     &     xcur(nsite,nloc,nalmax),xprop(nsite,nloc,nalmax),
     &     ycur(nsite,nloc,nalmax),yprop(nsite,nloc,nalmax),
     &     zcur(nsite,nloc,nalmax),zprop(nsite,nloc,nalmax),
     &     Sigmacur(nsite,nsite), Sigmaprop(nsite,nsite),
     &     Uchol(nsite,nsite)
*     local variables
      integer iit,fracnit,i,j,wit, isite, iloc, ial,ud
      double precision ka, to
      dimension ud(5),ka(5)
!     call intpr('************************************',-1,0,0)
!     call intpr('***  Entering Fortran program    ***',-1,0,0)
!     call intpr('************************************',-1,0,0)
      call rndstart()


      to=0
         do i = 1,5
         if ( ud(i)==1 ) then
            to=to+1
            ka(i)=to
         else
            ka(i)=0
         endif
      enddo

      fracnit = nit/thinning
      wit = 1

            

*     initialize
      call init(nsite,nloc,nal,nalmax,s,
     &     DG,DE,thetacur,nit,thinning, 
     &     fcur,xcur,ycur,zcur,Sigmacur,Uchol,swi)





      do iit = 1,nit

*     storing outputof current iteration into MCMC storage array
         if(mod(iit,thinning) .eq. 0) then 
            do i=1,5
               theta (i,wit) = thetacur(i)
               zeta(i,wit) = zcur(1,i,1)
            enddo
            do  isite = 1, nsite
               do iloc = 1, nloc
                  do ial = 1, nalmax
                     fsav(wit,isite,iloc,ial)=fcur(isite,iloc,ial)
                     zsav(wit,isite,iloc,ial)=zcur(isite,iloc,ial)
                  enddo
               enddo
            enddo
            wit = wit + 1
         endif

         call mhcov(nsite,nloc,nal,nalmax,g,s,DG,DE,thetacur,thetaprop,
     &        fcur,fprop,xcur,xprop,ycur,yprop,zcur,iit,fracnit,
     &        Sigmacur,Sigmaprop,Uchol,am,bdm,bem,gm,dm,ud,ka,to,swi)

         call mhzed(nsite,nloc,nal,nalmax,g,s,DG,DE,thetacur,thetaprop,
     &        fcur,fprop,xcur,xprop,ycur,yprop,zcur,zprop,iit,fracnit,
     &        Sigmacur,Uchol,am, bdm,bem,gm,dm)

      enddo







      call rndend()  
!     call intpr('************************************',-1,0,0)
!     call intpr('***  Leaving Fortran program     ***',-1,0,0)
!     call intpr('************************************',-1,0,0)
      end subroutine tangle




************************************************************************
*     propose change in parameters of covariance function
*     updates Choleski decomposition 
*     upates variables f
      subroutine mhcov(nsite,nloc,nal,nalmax,g,s,DG,DE,thetacur,
     &     thetaprop,fcur,fprop,xcur,xprop,ycur,yprop,zcur,iit,fracnit,
     &     Sigmacur,Sigmaprop,Uchol,am,bdm,bem,gm,dm,ud,ka,to,swi
     &     )
      implicit none
*     data
      integer nsite,nloc,nal,nalmax,g,iit,fracnit
      double precision s,DG,DE
      double precision am, bdm,bem,gm,dm
      dimension  s(2,nsite),g(nsite,nloc,nalmax),
     &     nal(nloc),DG(nsite,nsite),DE(nsite,nsite)
*     parameters/variables to be simulated
      double precision thetacur,thetaprop
      dimension thetacur(5),thetaprop(5)
*     computing options
      integer nit,thinning,swi
*     hidden variables
      double precision fcur,fprop,xcur,xprop,ycur,yprop,zcur,
     &     Sigmacur,Sigmaprop,Uchol,dice
      dimension fcur(nsite,nloc,nalmax),fprop(nsite,nloc,nalmax),
     &     xcur(nsite,nloc,nalmax),xprop(nsite,nloc,nalmax),
     &     ycur(nsite,nloc,nalmax),yprop(nsite,nloc,nalmax),
     &     zcur(nsite,nloc,nalmax),Sigmacur(nsite,nsite),
     &     Sigmaprop(nsite,nsite),Uchol(nsite,nsite)
*     local variables
      integer i,j,info,iloc,ial,isite,jsite,ud
      double precision sum,ggrunif,alpha,accept,r,lr,p,anam,yy,ka,to
      dimension ud(5), ka(5)


      


      do i=1,5
         thetaprop(i)=thetacur(i)
      enddo




*     propose new values
      dice = ggrunif(dble(0),dble(1))


      if (ud(1)==1 .and. (ka(1)-1)/to < dice .and. 
     &     dice < ka(1)/to) then
         thetaprop(1) = thetacur(1) +  ggrunif(dble(-.01),dble(.01))
      elseif (ud(2)==1 .and. (ka(2)-1)/to < dice .and. 
     &        dice < ka(2)/to) then
         thetaprop(2) = thetacur(2) +  ggrunif(dble(-.01),dble(.01))
      elseif (ud(3)==1 .and. (ka(3)-1)/to < dice .and. 
     &        dice < ka(3)/to) then
         thetaprop(3) = thetacur(3) +  ggrunif(dble(-.01),dble(.01))
      elseif (ud(4)==1 .and. (ka(4)-1)/to < dice .and.
     &        dice < ka(4)/to) then
         thetaprop(4) = thetacur(4) +  ggrunif(dble(-.01),dble(.01))
      elseif (ud(5)==1 .and. (ka(5)-1)/to < dice .and.
     &        dice < ka(5)/to)  then
         thetaprop(5) = thetacur(5) +  ggrunif(dble(-.01),dble(.01))
      endif
      
      if ( (((thetaprop(1) .gt. 0  .and.  thetaprop(1) .le. am) .and.
     &     (thetaprop(2) .gt. 0  .and.  thetaprop(2) .le. bdm )).and.
     &     ((thetaprop(3) .gt. 0  .and.  thetaprop(3) .le. bem) .and.
     &     (thetaprop(4) .gt. 0  .and.  thetaprop(4) .le. gm ))) .and.
     &     (thetaprop(5) .ge. 0 .and.  thetaprop(5) .le. dm ) ) then

*     computes covariance matrix
         call buildcov(Sigmaprop,DG,DE,nsite,thetaprop,swi)

*     computes proposed Choleski decomposition
         call dpofa(Sigmaprop,nsite,nsite,info)
*     now Sigma contains U
         if(info .ne. 0) then
            call intpr('non-0 exit of dpofa in mhcov',-1,0,0)
         endif 
         
*     computes proposed y (y^t = z^t * U)
         do iloc=1,nloc
            do ial=1,nal(iloc)
               yprop(1,iloc,ial) = Sigmaprop(1,1)*zcur(1,iloc,ial)
               if(nsite .gt. 1) then
                  do j = 2,nsite
                     yprop(j,iloc,ial) = 0
                     do i = 1,j
                        yprop(j,iloc,ial) = yprop(j,iloc,ial) + 
     &                       Sigmaprop(i,j) * zcur(i,iloc,ial)
                     enddo
                  enddo
               endif
            enddo
         enddo
         
*     transform into Gamma distributed values
         do iloc=1,nloc
            do ial=1,nal(iloc)
               do isite=1,nsite
                  xprop(isite,iloc,ial) = 
     &                 anam(yprop(isite,iloc,ial),thetaprop(1))
               enddo
            enddo
         enddo
         
*     divide f(isite,iloc,.) by sum of f(isite,iloc,.)  
*     so as to have sum = 1
         do isite=1,nsite
            do iloc=1,nloc
               sum = 0
               do ial=1,nal(iloc)
                  sum = sum + xprop(isite,iloc,ial)
               enddo
               if(sum .gt. 0) then 
                  do ial=1,nal(iloc)
                     fprop(isite,iloc,ial) = 
     &                    xprop(isite,iloc,ial)/sum
                  enddo
               endif
               if(sum .eq. 0) then
               endif
            enddo
         enddo
         
*     MH ratio
         lr = 0
         do isite=1,nsite
            do iloc=1,nloc
               do ial=1,nal(iloc)
                  lr = lr + g(isite,iloc,ial)*
     &                 (dlog(fprop(isite,iloc,ial)/
     &                 fcur(isite,iloc,ial)))
               enddo
            enddo
         enddo
         p = dmin1(dexp(lr),1.d0)
         accept = ggrunif(dble(0),dble(1))
         if(accept .le. p) then
            do i=1,5
               thetacur(i) = thetaprop(i) 

            enddo
            do isite=1,nsite
               do iloc=1,nloc
                  do ial=1,nal(iloc)
                     fcur(isite,iloc,ial) = fprop(isite,iloc,ial)
                  enddo
               enddo
            enddo

            do isite=1,nsite
               do jsite=1,nsite
                  Sigmacur(isite,jsite) = Sigmaprop(isite,jsite)
               enddo
            enddo
            

         endif

      endif

      end subroutine mhcov
      

      subroutine mhzed(nsite,nloc,nal,nalmax,g,s,DG,DE,
     &     thetacur,thetaprop,fcur,fprop,xcur,xprop,ycur,yprop,
     &     zcur,zprop,iit,fracnit,
     &     Sigmacur,Uchol,
     &     am, bdm,bem,gm,dm
     &     )
      implicit none
*     data
      integer nsite,nloc,nal,nalmax,g,iit,fracnit
      double precision s,DG,DE
      double precision am, bdm,bem,gm,dm
      dimension  s(2,nsite),g(nsite,nloc,nalmax),
     &     nal(nloc),DG(nsite,nsite),DE(nsite,nsite)
*     parameters/variables to be simulated
      double precision thetacur,thetaprop
      dimension thetacur(5),thetaprop(5)
*     computing options
      integer nit,thinning
*     hidden variables
      double precision fcur,fprop,xcur,xprop,ycur,yprop,zcur,zprop,
     &     Sigmacur,Uchol,
     &     dice
      dimension fcur(nsite,nloc,nalmax),fprop(nsite,nloc,nalmax),
     &     xcur(nsite,nloc,nalmax),xprop(nsite,nloc,nalmax),
     &     ycur(nsite,nloc,nalmax),yprop(nsite,nloc,nalmax),
     &     zcur(nsite,nloc,nalmax),zprop(nsite,nloc,nalmax),
     &     Sigmacur(nsite,nsite),Uchol(nsite,nsite)
*     local variables
      integer i,j,info,iloc,ial,isite
      double precision sum,ggrunif,ggrnorm,alpha,accept,r,lr,p,anam,yy

      do iloc = 1,nloc

*     Proposes a Z different only for the "slice" iloc
         do isite=1,nsite
            do ial=1,nal(iloc)
               zprop(isite,iloc,ial) = zcur(isite,iloc,ial) +              
     &              ggrnorm(dble(.0),dble(.01))
            enddo
         enddo
         

*     computes "slice" iloc of y (y^t = z^t * U)
         do ial=1,nal(iloc)


            yprop(1,iloc,ial) = Sigmacur(1,1)*zprop(1,iloc,ial)
            if(nsite .gt. 1) then

               do j = 2,nsite
                  yprop(j,iloc,ial) = 0
                  do i = 1,j
                     yprop(j,iloc,ial) = yprop(j,iloc,ial) + 
     &                    Sigmacur(i,j) * zprop(i,iloc,ial)
                  enddo
               enddo

            endif



         enddo
         
         
*     transforms into Gamma distributed values
         do ial=1,nal(iloc)
            do isite=1,nsite
               xprop(isite,iloc,ial) = 
     &              anam(yprop(isite,iloc,ial),thetacur(1))
            enddo
         enddo
         
*     divides f(isite,iloc,.) by sum of f(isite,iloc,.)  
*     so as to have sum = 1
         do isite=1,nsite
            sum = 0
            do ial=1,nal(iloc)
               sum = sum + xprop(isite,iloc,ial)
            enddo
            if(sum .gt. 0) then 
               do ial=1,nal(iloc)
                  fprop(isite,iloc,ial) = 
     &                 xprop(isite,iloc,ial)/sum
                 
               enddo
            endif
            if(sum .eq. 0) then
               call intpr('AAA all freq = 0',-1,0,0) 
            endif
         enddo
         

*     computes MH ratio
         lr = 0
         do isite=1,nsite
            do ial=1,nal(iloc)
                
               lr = lr + g(isite,iloc,ial)*
     &              (dlog(fprop(isite,iloc,ial)/
     &              fcur(isite,iloc,ial))) -
     &              (zprop(isite,iloc,ial)**2 -
     &              zcur(isite,iloc,ial)**2)/2

            enddo
         enddo

         p = dmin1(dexp(lr),1.d0)
         accept = ggrunif(dble(0),dble(1))

         if(accept .le. p) then

            do isite=1,nsite
               do ial=1,nal(iloc)
                  zcur (isite,iloc,ial) = zprop (isite,iloc,ial)
               enddo
            enddo
            
            do isite=1,nsite
               do ial=1,nal(iloc)
                  fcur(isite,iloc,ial) = fprop(isite,iloc,ial)
               enddo
            enddo

         endif
      enddo
      
      end subroutine mhzed






*************************************************
*     build covariance matrix for MVN variables y
*     upper part only
      subroutine buildcov(Sigma,DG,DE,nsite,thetacur,swi)
      implicit none
      integer nsite
      double precision Sigma,DG,DE,thetacur
      dimension Sigma(nsite,nsite),DG(nsite,nsite),DE(nsite,nsite),
     &     thetacur(5)
      integer isite1,isite2,swi
*     thetacur = (alpha,beta_G,beta_E,gamma,delta)
      do isite1=1,nsite
         do isite2=isite1,nsite

            select case(swi)
            case(1)
            Sigma(isite1,isite2) = (1-thetacur(5))*
     &           dexp(-(DG(isite1,isite2)/thetacur(2)+
     &       DE(isite1,isite2)/thetacur(3))**thetacur(4))

         case(2)
            Sigma(isite1,isite2) = (1-thetacur(5))*
     &           dexp(-(DG(isite1,isite2)/thetacur(2))**thetacur(4))

         case(3)
            Sigma(isite1,isite2) = (1-thetacur(5))*
     &           dexp(-(DE(isite1,isite2)/thetacur(3))**thetacur(4))


            end select


            if(isite2 .eq. isite1) then
               Sigma(isite1,isite2) = Sigma(isite1,isite2) + 
     &              thetacur(5)
            endif
         enddo
      enddo
      end subroutine buildcov


*************************************************
*     transform  Gaussian N(0,1) 
*     into Gamma  G(shape=alpha,rate=1) 
      double precision function anam(y,alpha)
      implicit none 
      double precision y,p,alpha,ggpnorm,ggqgam
      p = ggpnorm(y,0.d0,1.d0,1,0)
      anam = ggqgam(p,alpha,1.d0,1,0)
      end function 



* Give initial values for the M-H algorithm
***********************************
      subroutine init(nsite,nloc,nal,nalmax,s,
     &     DG,DE,thetacur,nit,thinning,
     &     fcur,xcur,ycur,zcur,Sigmacur,Uchol,swi)
      implicit none
*     data 
      integer nsite,nloc,nal,nalmax,g
      double precision s,DG,DE
      dimension  s(2,nsite),
     &     nal(nloc),DG(nsite,nsite),DE(nsite,nsite)

*     parameters/variables to be simulated
      double precision thetacur
      dimension thetacur(5)

*     computing options
      integer nit,thinning,swi

*     hidden variables
      double precision fcur,fprop,xcur,xprop,ycur,yprop,zcur,
     &     Sigmacur, Uchol
      dimension fcur(nsite,nloc,nalmax),
     &     xcur(nsite,nloc,nalmax),
     &     ycur(nsite,nloc,nalmax),
     &     zcur(nsite,nloc,nalmax),
     &     Sigmacur(nsite,nsite),Uchol(nsite,nsite)
      integer iloc,isite,ial,i,j,info
      double precision sum,anam

*     computes covariance matrix
      call buildcov(Sigmacur,DG,DE,nsite,thetacur,swi)
         
*     computes proposed Choleski decomposition
         call dpofa(Sigmacur,nsite,nsite,info)
*     now Sigma contains U
         if(info .ne. 0) then
            call intpr('non-0 exit of dpofa in mhcov',-1,0,0)
         endif 
         


         
*     computes proposed y (y^t = z^t * U)
         do iloc=1,nloc
            do ial=1,nal(iloc)
               ycur(1,iloc,ial) = Sigmacur(1,1)*zcur(1,iloc,ial)
               if(nsite .gt. 1) then
                  do j = 2,nsite
                     ycur(j,iloc,ial) = 0
                     do i = 1,j
                        ycur(j,iloc,ial) = ycur(j,iloc,ial) + 
     &                       Sigmacur(i,j) * zcur(i,iloc,ial)
                     enddo
                  enddo
               endif
            enddo
         enddo

               


*     transform into Gamma distributed values
         do iloc=1,nloc
            do ial=1,nal(iloc)
               do isite=1,nsite
                  xcur(isite,iloc,ial) = 
     &                 anam(ycur(isite,iloc,ial),thetacur(1))
               enddo
            enddo
         enddo
         



*     divide f(isite,iloc,.) by sum of f(isite,iloc,.)  
*     so as to have sum = 1
         do isite=1,nsite
            do iloc=1,nloc
               sum = 0
               do ial=1,nal(iloc)
                  sum = sum + xcur(isite,iloc,ial)
               enddo
               if(sum .gt. 0) then 
                  do ial=1,nal(iloc)
                     fcur(isite,iloc,ial) = 
     &                    xcur(isite,iloc,ial)/sum
                  enddo
               endif
               if(sum .eq. 0) then
                  call intpr('BBB all freq = 0',-1,0,0) 
c init
               endif
            enddo
         enddo


 


      end subroutine init
