program cytokinesis

use declare
use cdk
USE IEEE_ARITHMETIC
implicit none
integer omp_get_thread_num,omp_get_num_threads,tid


   CALL system_clock(count_rate=cr)
   rate = REAL(cr)
   call system_clock(timestart)



!-----------------------------------------
!   Select job:
   call setjob(nthreads,jobid,mbrad,mwid,rrad,rwid,rthick,falen,fanum,myonum,crlknum1,crlknum2, &
        ntether,kmemb,tension1,tension2,crlkorient,j_dep,cofilin,jmyoturn,jxldel,jxlturnover,ltether,halftime)


!$omp parallel &
!$omp default(none) &
!$omp private(tid) &
!$omp shared(nthreads)

tid=omp_get_thread_num()
if(tid==0)then
nthreads=omp_get_num_threads()
end if

!$omp end parallel


   write(*,*)'================================================'
   write(*,*)'Cytokinesis modeling version 1'

!-----------------------------------------------------

!  common parameters:
   pi=3.141592653589793239D0
   delta=0.000001d0
   invdelta=1.0d0/delta
   beta=0.00000000000001D0

!  Spring constant and relaxed length for actin:

   k_a=100.0d0 ! unit = 1e-20 J/nm**2
   l_a=5.4d0    ! unit length of 2 G-actin segment
   invl_a=1.0d0/l_a

!  F-actin persistence length L_p ~ 10 micron
!  Bending rigidity of F-actin calculated as k_B*T*L_p/l_a at room temp = 295 K

   ka_thet=240.d0 ! unit = 1e-20 J
   thet_fa=pi


!  Myosin head-body angle in two conformations:


   kmh_thet=50.0d0
   kmb_thet=50.0d0

   thet_mh1=pi/3
   thet_mh2=2*pi/3
   thet_mb=pi


!  Myosin head length:

   l_mh=10.0d0

!  Distance between myosin beads on body:

   l_mb=10.0d0

!  Assuming a four-head molecule is 150 nm long
!  the number of myosin beads/molecule not including the four heads:

   mybodylen=16

!  Modeling crosslinks after alpha-actinin

   lklen=3

   l_lk1=11.0d0

   l_lk2=5.0d0

!  assume linear rigidity of crosslinkers are the same with those of actin

   k_lk=k_a*l_a/l_lk1

!  assume bending rigidity of crosslinkers 10 times smaller than that of actin

   klk_thet=ka_thet*l_lk1/l_a/10

   thet_lk=pi


   ktether=k_a*l_a/ltether


!  not all beads should be tethered, but just a part of them, e.g. 10%

   p_tether=1.0d0*ntether/fanum/falen!0.05d0

!  bead size for membrane beads:

   l_mem=10.0d0

   k_mem=10.0d0


!  cell wall stiffness:

   wthick=20.0d0

   lsqueez=5.0d0


   kwall=0.005d0 !fturgor/wallgap



!-----------------------------------------------------
!  Generate random seed:
   CALL INIT_RANDOM_SEED()

!======================================================
!  Direction of the job:

   if(jobid==1)then
      print*,'to generate the fellowship of the ring ...'
   elseif(jobid==2)then
      goto 29
   else
      print*,'suitable jobid was not found'
      stop
   end if

!--------------------------------------------------

!  This section is to create the ring.


!  First create a cell wall

!  define cell wall radius = membrane radius + wall thickness

   walrad=mbrad+wthick

!  bead size:

   wallbead=10.0d0

!  number of cell wall bead seprated at 25 nm on one hoop of radius =1.0*walrad, also number of solid angle:

   nphisol=2*pi*walrad/25

!  the angle between these beads is also the small solid angle component:
   dphisol=2*pi/nphisol


!  number of hoops separated xdel= 20 nm within membrane width, also component of solid angle:

   dxsol=20.0d0

!  number of cell wall hoops:
   nxsol=(mwid+1.0d0)/dxsol

   if(mod(nxsol,2)==0)then
      nxsol=nxsol+1
   end if

   allocate(xwall(nphisol,nxsol),ywall(nphisol,nxsol),zwall(nphisol,nxsol))

   allocate(xnorwall(nphisol,nxsol),ynorwall(nphisol,nxsol),znorwall(nphisol,nxsol))

!   rwall=walrad


!  the total number of wall beads:
   nwall=nxsol*nphisol


   nw=0

   xmin=-(nxsol-1)/2*dxsol

   xmax=-xmin


   xboundmin=xmin
   xboundmax=xmax


   do jx=1,nxsol

      x0=xmin+(jx-1)*dxsol

      xwall(1:nphisol,jx)=x0

      xnorwall(1:nphisol,jx)=0.0d0


      do jp=1,nphisol

         nw=nw+1

         ywall(jp,jx)=walrad*cos(jp*dphisol)

         zwall(jp,jx)=walrad*sin(jp*dphisol)


         ynorwall(jp,jx)=-ywall(jp,jx)/walrad
         znorwall(jp,jx)=-zwall(jp,jx)/walrad


      end do

!     adding bonds between hoops:


   end do

!  to visualize cell wall growth, start with the surface:

   nwalsurf_p=2*pi*walrad/wallbead

   nwalsurf_x=(mwid+1.0d0)/wallbead

   nwalsurf=nwalsurf_p*nwalsurf_x

   allocate(xwalsurf(nwalsurf),ywalsurf(nwalsurf),zwalsurf(nwalsurf))

   dphi=2*pi/nwalsurf_p

!  REMEMBER, for convenience, surface coordinates are scaled down by 10:

   xwalsurf(1:nwalsurf_p)=xmin/10

   do n=1,nwalsurf_p

      phi=(n-1)*dphi

      ywalsurf(n)=0.1*walrad*cos(phi)

      zwalsurf(n)=0.1*walrad*sin(phi)

   end do

   nwalsurf=nwalsurf_p

   do n=2,nwalsurf_x

      xwalsurf(nwalsurf+1:nwalsurf+nwalsurf_p)=(xmin+(n-1)*wallbead)/10

      ywalsurf(nwalsurf+1:nwalsurf+nwalsurf_p)=ywalsurf(1:nwalsurf_p)

      zwalsurf(nwalsurf+1:nwalsurf+nwalsurf_p)=zwalsurf(1:nwalsurf_p)

      nwalsurf=nwalsurf+nwalsurf_p

   end do

!  to visualize radial beads:


   nlayer=walrad/wallbead

   nwalrad=2*nlayer*nwalsurf_p

   allocate(xwalrep(nwalrad),ywalrep(nwalrad),zwalrep(nwalrad),wposid(2,nwalrad))

   nwalrad=0

   do n=1,nlayer

      rad=walrad-n*wallbead

      if(rad<50.0)then
         exit
      end if

      nphi=2*pi*rad/wallbead

      dphi=2*pi/nphi

      do j=1,nphi

         phi=j*dphi

         nwalrad=nwalrad+1

         ywalrep(nwalrad)=rad*cos(phi)

         zwalrep(nwalrad)=rad*sin(phi)

         wposid(1,nwalrad)=(phi+dphisol/2)/dphisol

         if(wposid(1,nwalrad)==0)then
            wposid(1,nwalrad)=nphisol
         end if

         if(wposid(1,nwalrad)>nphisol)then
            wposid(1,nwalrad)=nphisol
         end if

      end do

   end do

   xwalrep(1:nwalrad)=xmin

   wposid(2,1:nwalrad)=1

   xwalrep(nwalrad+1:2*nwalrad)=-xmin

   ywalrep(nwalrad+1:2*nwalrad)=ywalrep(1:nwalrad)

   zwalrep(nwalrad+1:2*nwalrad)=zwalrep(1:nwalrad)

   wposid(1,nwalrad+1:2*nwalrad)=wposid(1,1:nwalrad)

   wposid(2,nwalrad+1:2*nwalrad)=nxsol

   nwalrad=nwalrad*2

!  at beginning, store radial beads on the surface:
!  REMEMBER, scale up coordinates by 10:

   allocate(xwalrad(nwalrad),ywalrad(nwalrad),zwalrad(nwalrad))

   n=0

   do j=1,nwalrad

      n=n+1

      xwalrad(j)=10*xwalsurf(n)
      ywalrad(j)=10*ywalsurf(n)
      zwalrad(j)=10*zwalsurf(n)

      if(n==nwalsurf)then
         n=0
      end if

   end do

!------------------------------------------------------

!  Create a membrane composed of a single-layer

!  if the bead size is l_mem then number of bead per hoop:

   nmeho=2*pi*mbrad/l_mem

   dphi=2*pi/nmeho

!  number of hoops:


   nhoop=(xmax-xmin+1.0d0)/l_mem

   if(mod(nhoop,2)==0)then
      nhoop=nhoop+1
   end if

!  total number of beads:

   nmemb=nhoop*nmeho

   allocate(xmemb(nmemb),ymemb(nmemb),zmemb(nmemb))


   nm=0

   do jh=1,nhoop

      x0=xmin+(jh-1)*l_mem

      xmemb(nm+1:nm+nmeho)=x0


      do jm=1,nmeho

         nm=nm+1

         ymemb(nm)=mbrad*cos(dphi*jm)

         zmemb(nm)=mbrad*sin(dphi*jm)



      end do

   end do

!  coarse representation of membrane surface:

   allocate(xsurf(nphisol,nxsol),ysurf(nphisol,nxsol),zsurf(nphisol,nxsol))

   xsurf=xwall

!  but:

   ysurf=ywall*mbrad/walrad

   zsurf=zwall*mbrad/walrad



!------------------------------------------------------

!  Now generate F-actin

   nfa=2*falen*fanum

!  coordinates:
   allocate(xfa(nfa),yfa(nfa),zfa(nfa))

!  filament identification:
   allocate(filid(nfa),apos(nfa))


!  membrane tethering:
   allocate(a2mem(nfa))
   a2mem=0
   allocate(mark(nmemb))
   mark=0


!  filament configuration:
   allocate(astart(fanum),alen(fanum))

   astart0=1
   length=0

   alen=0

   nfa=0
   nfa1=0

   ntethertem=0

   fanum1=fanum/2

   fanum2=fanum-fanum1

   do n=1,fanum

!     start index on filament:

      astart(n)=astart0+length


      astart0=astart(n)

!     pick the F-actin length randomly from falen/2 to 3falen/2

      call random_number(r)

      length=falen/2+falen*r

!     filament id:
      filid(nfa+1:nfa+length)=n

      alen(n)=length

!     pick the x position of F-actin:

      call random_number(r)

      x0=-rwid/2+rwid*r

!     pick the angle of first bead:

      call random_number(r)

      phi0=2*pi*r



!     pick the "radius" of the F-actin hoop

      call random_number(r)

      rad=rrad-rthick*r

!     increment in angle between beads:

      dphi=l_a/rad


      if(n<=fanum1)then


         nfa1=nfa1+length

         jdir=1

      else


         jdir=-1

      end if

!     now assign coordinates:

      xfa(nfa+1:nfa+length)=x0

      do jf=1,length

         phi=phi0+(jf-1)*dphi*jdir

         nfa=nfa+1

         apos(nfa)=jf

         yfa(nfa)=rad*cos(phi)

         zfa(nfa)=rad*sin(phi)


      end do

!     assign tethering on the beads:

      if(mbrad-rad<ltether)then

         do nf=nfa-length+1,nfa

            call random_number(r)

            if(p_tether>r.and.ntethertem<ntether)then

               jh=(x0-xmin)/l_mem+1

               jm0=(jh-1)*nmeho

               dist2=(ltether+l_mem)**2

               npick=0

               do jm=jm0+1,jm0+2*nmeho

                  if(mark(jm)==1)then
                     cycle
                  end if

                  dx=xfa(nf)-xmemb(jm)
                  dy=yfa(nf)-ymemb(jm)
                  dz=zfa(nf)-zmemb(jm)

                  d2=dx*dx+dy*dy+dz*dz

                  if(d2<dist2)then

                     dist2=d2

                     npick=jm

                  end if

               end do

               if(npick==0)then
                  print*,'could not find npick'

               else

                  mark(npick)=1

               end if

               a2mem(nf)=npick

               ntethertem=ntethertem+1

            end if

         end do

      end if


   end do

   nfa2=nfa-nfa1


!---------------------------------------------------

!  generating a set of myosin molecules ...

   nmyo=myonum*(4+mybodylen)

!  coordinates:
   allocate(xmyo(nmyo),ymyo(nmyo),zmyo(nmyo))


!  head and body addresses:
   allocate(myhead(4,myonum),mybody(mybodylen,myonum))

!  bonds:
   allocate(bondmyo(2,nmyo))

   nmyo=0
   nb=0

   do nm=1,myonum

!     pick the x position of myosin:

      call random_number(r)

      x0=-rwid/2+rwid*r

!     pick the angle of first body bead:

      call random_number(r)

      phi0=2*pi*r


!     pick the "radius" of the myosin hoop

      call random_number(r)

      rad=rrad-rthick*r

!     increment in angle between body beads:

      dphi=l_mb/rad

!     pick the elongating direction

      call random_number(r)

      if(r>0.5d0)then
         jdir=1
      else
         jdir=-1
      end if


!     now assigning coordinates for body beads ...

      xmyo(nmyo+1:nmyo+mybodylen)=x0

      do jm=1,mybodylen

         phi=phi0+(jm-1)*dphi*jdir

         nmyo=nmyo+1

         mybody(jm,nm)=nmyo

         ymyo(nmyo)=rad*cos(phi)

         zmyo(nmyo)=rad*sin(phi)

         if(jm<mybodylen)then

            nb=nb+1

            bondmyo(1,nb)=nmyo

            bondmyo(2,nb)=nmyo+1

         end if

      end do

!     assigning coordinates for head beads ...
!     the first two heads are attached to the first body bead

      n1=mybody(1,nm)

21    call random_number(r)

      dx=r-0.5d0

      call random_number(r)

      dy=r-0.5d0

      call random_number(r)

      dz=r-0.5d0

      dist=sqrt(dx*dx+dy*dy+dz*dz)

      dx=dx*l_mh/dist
      dy=dy*l_mh/dist
      dz=dz*l_mh/dist

      xh=xmyo(n1)+dx
      yh=ymyo(n1)+dy
      zh=zmyo(n1)+dz

      if(yh*yh+zh*zh>mbrad*mbrad)then
         goto 21
      end if


      nmyo=nmyo+1

      xmyo(nmyo)=xh
      ymyo(nmyo)=yh
      zmyo(nmyo)=zh

      myhead(1,nm)=nmyo

      nb=nb+1

      bondmyo(1,nb)=n1
      bondmyo(2,nb)=nmyo

!     do similar steps for head bead 2:

22    call random_number(r)

      dx=r-0.5d0

      call random_number(r)

      dy=r-0.5d0

      call random_number(r)

      dz=r-0.5d0

      dist=sqrt(dx*dx+dy*dy+dz*dz)

      dx=dx*l_mh/dist
      dy=dy*l_mh/dist
      dz=dz*l_mh/dist

      xh=xmyo(n1)+dx
      yh=ymyo(n1)+dy
      zh=zmyo(n1)+dz

      if(yh*yh+zh*zh>mbrad*mbrad)then
         goto 22
      end if


      nmyo=nmyo+1

      xmyo(nmyo)=xh
      ymyo(nmyo)=yh
      zmyo(nmyo)=zh


      myhead(2,nm)=nmyo

      nb=nb+1

      bondmyo(1,nb)=n1
      bondmyo(2,nb)=nmyo

!     for head beads 3 and 4, repeat similar steps:

      n1=mybody(mybodylen,nm)

23    call random_number(r)

      dx=r-0.5d0

      call random_number(r)

      dy=r-0.5d0

      call random_number(r)

      dz=r-0.5d0

      dist=sqrt(dx*dx+dy*dy+dz*dz)

      dx=dx*l_mh/dist
      dy=dy*l_mh/dist
      dz=dz*l_mh/dist

      xh=xmyo(n1)+dx
      yh=ymyo(n1)+dy
      zh=zmyo(n1)+dz

      if(yh*yh+zh*zh>mbrad*mbrad)then
         goto 23
      end if


      nmyo=nmyo+1

      xmyo(nmyo)=xh
      ymyo(nmyo)=yh
      zmyo(nmyo)=zh


      myhead(3,nm)=nmyo

      nb=nb+1

      bondmyo(1,nb)=n1
      bondmyo(2,nb)=nmyo

!     head 4:

24    call random_number(r)

      dx=r-0.5d0

      call random_number(r)

      dy=r-0.5d0

      call random_number(r)

      dz=r-0.5d0

      dist=sqrt(dx*dx+dy*dy+dz*dz)

      dx=dx*l_mh/dist
      dy=dy*l_mh/dist
      dz=dz*l_mh/dist

      xh=xmyo(n1)+dx
      yh=ymyo(n1)+dy
      zh=zmyo(n1)+dz

      if(yh*yh+zh*zh>mbrad*mbrad)then
         goto 24
      end if


      nmyo=nmyo+1

      xmyo(nmyo)=xh
      ymyo(nmyo)=yh
      zmyo(nmyo)=zh


      myhead(4,nm)=nmyo

      nb=nb+1

      bondmyo(1,nb)=n1
      bondmyo(2,nb)=nmyo

   end do

   nbondmyo=nb

!------------------------------

!  Generating crosslinkers in the ring:

   crlknum=crlknum1+crlknum2

   nlk=crlknum*lklen

   allocate(xlk(nlk),ylk(nlk),zlk(nlk))

   allocate(bondlk(2,nlk))

   allocate(lkstart(crlknum))


   nlk=0

   nbondlk=0

   do n=1,crlknum

      if(n<=crlknum1)then
         l_lk=l_lk1
      else
         l_lk=l_lk2
      end if

!     pick the x position of crosslinker center:

      call random_number(r)

      x0=-rwid/2+rwid*r

!     pick the angular position of the center:

      call random_number(r)

      phi0=2*pi*r



!     pick the "radius" of the center

      call random_number(r)

      rad=rrad-rthick*r

      y0=rad*cos(phi0)

      z0=rad*sin(phi0)

!     pick the direction of the crosslinker

      call random_number(r)

      dx=r-0.5d0

      call random_number(r)

      dy=r-0.5d0

      call random_number(r)

      dz=r-0.5d0

      dist=sqrt(dx*dx+dy*dy+dz*dz)

      dx=dx/dist
      dy=dy/dist
      dz=dz/dist

25    if((y0-(lklen-1)*l_lk/2*dy)**2+(z0-(lklen-1)*l_lk/2*dz)**2>mbrad**2)then
         y0=y0*(mbrad-5.0)/mbrad
         z0=z0*(mbrad-5.0)/mbrad
         goto 25
      end if

26    if((y0+(lklen-1)*l_lk/2*dy)**2+(z0+(lklen-1)*l_lk/2*dz)**2>mbrad**2)then
         y0=y0*(mbrad-5.0)/mbrad
         z0=z0*(mbrad-5.0)/mbrad
         goto 26
      end if

      lkstart(n)=nlk+1

      do j=1,lklen

         nlk=nlk+1

         xlk(nlk)=x0+(j-1)*l_lk*dx-(lklen-1)*l_lk/2*dx
         ylk(nlk)=y0+(j-1)*l_lk*dy-(lklen-1)*l_lk/2*dy
         zlk(nlk)=z0+(j-1)*l_lk*dz-(lklen-1)*l_lk/2*dz


         if(j>1)then
            nbondlk=nbondlk+1

            bondlk(1,nbondlk)=nlk-1
            bondlk(2,nbondlk)=nlk
         end if

      end do

   end do

!------------------------------
!  write out the ring:

   call makering(nwall,nxsol,nphisol,nwalsurf,nwalrad,nmemb,fanum1,fanum2,nfa1,nfa2, &
               ntether,myonum,mybodylen,nmyo,nbondmyo,crlknum1,crlknum2,lklen,nlk,nbondlk, &
                astart,alen,filid,a2mem,apos,lkstart,mybody,myhead,bondmyo,bondlk,wposid, &
                 xboundmin,xboundmax,dxsol,dphisol,xwall,ywall,zwall,xnorwall,ynorwall,znorwall, &
                  xwalsurf,ywalsurf,zwalsurf,xwalrad,ywalrad,zwalrad,xwalrep,ywalrep,zwalrep, &
                   xmemb,ymemb,zmemb,xsurf,ysurf,zsurf,xfa,yfa,zfa,xmyo,ymyo,zmyo,xlk,ylk,zlk)



   stop

!=====================================================
29 print*,'Constriction of the ring ...'

!  getting the initial ring:

   call getinfo(nstart,jfile,time0,runtime,natom,natom_a1,natom_a2,natom_a2m, &
              nwall,nxsol,nphisol,nwalsurf,nwalrad,nmemb,nact1,nact2,nfamax1,nfamax2, &
              nfa1,nfa2,fanum1,fanum2,ntether,nmyo,myonum,nlk,crlknum1,crlknum2,crlknummax,crlknumactive,nmyoturn)


!  maximum possible number of F-actin:

   fanummax1=nact1/20+1

   fanummax1=max(fanummax1,fanum1)

   fanummax2=nact2/20+1

   fanummax2=max(fanummax2,fanum2)

   fanummax=fanummax1+fanummax2

   nact=nact1+nact2

   nfa=nfa1+nfa2

   fanum=fanum1+fanum2

!  total crosslinkers:

   crlknum=crlknum1+crlknum2


!  setting parameters:

!  when binding, the relaxed distance between actin and myosin is:
   lbind=5.0d0  ! unit = nm

!  binding strength represented with a spring constant:
   kbind=10.0d0


!  probability for myosin to hydrolyze ATP, to change mytyp from 1 --> 2
   p1_hydr0=0.00005d0/2 ! unit = inserve micro second
   p1_hydr=0.0d0

!  this is used to calculate activation probability:
   invmylen2=(mybodylen-1)*l_mb
   invmylen2=invmylen2*invmylen2
   invmylen2=1.0d0/invmylen2

!  probability for ADP + Pi myosin to bind actin, to change mytyp from 2 --> 3
   p2_bind=0.0001d0/2

!  probability to release Pi, to change mytyp from 3 --> 4
   p3_pi_rele=0.00005d0/2

!  probability to release ADP, to change mytyp from 4 --> 5
   p4_adp_rele=0.00005d0/2

!  probability to bind ATP and unbind actin, to change mytyp from 5 --> 1
   p5_ubind=0.0003d0/2

!  for myosin head status, use mytyp
!  mytyp = 1 for ATP-bound, non-actin-binding, and angle = thet_mh1
!  mytyp = 2 for ADP+Pi, non-actin-binding, and angle = thet_mh2
!  mytyp = 3 for ADP+Pi, actin-binding, and angle = thet_mh2
!  mytyp = 4 for ADP, actin-binding, and angle = thet_mh1
!  mytyp = 5 for nucleotide-free, actin-binding, and angle = thet_mh1

!  myosin turnover rate (14 seconds):

!-------------------------------------------------------------

!  skiping a number of steps when calculating time step:
   jskip=10

!  skiping a big number of step when checking turnover:

   nturnover=100000

!  myosin turnover rate (Mithilesh insisted ~ 14 seconds):

   pmyturn=1.0d0/halftime/1000000*nturnover/jskip ! unit = inverse micro second

!  depolymerization of F-actin once a second:
   p_dep=j_dep*1.0d0/1000000*nturnover/jskip ! unit = inverse micro second



!  To model binding and unbinding of crosslinkers to actin:

   plk_bind=1.0d0/10000 ! unit = inverse micro second
   plk_ubind1=3.3d0/1000000
   plk_ubind2=0.05d0/1000000

!  removing crosslinkers once every 1000 seconds:

   plk_remove=1.0d0/1000.0/1000000*nturnover/jskip ! unit = inverse micro second

!  crosslinker turnover:
   plk_turn=1.0d0/20.0d0/1000000*nturnover/jskip


!--------------------------------------------

   allocate(apar(2,nact),mytyp(4,myonum))

!  to tell status of a f-actin bead, use apar (see 2015jan08A)

   allocate(fa1stbound(fanummax))
!  fa1stbound tell the lowest index on a filament that is bound to a crosslinker.


!  configuration of the system:

   allocate(xwall(nphisol,nxsol),ywall(nphisol,nxsol),zwall(nphisol,nxsol))
   allocate(xwalrad(nwalrad),ywalrad(nwalrad),zwalrad(nwalrad))
   allocate(xwalrep(nwalrad),ywalrep(nwalrad),zwalrep(nwalrad))
   allocate(xwalsurf(nwalsurf),ywalsurf(nwalsurf),zwalsurf(nwalsurf))
   allocate(waltyp(nwalrad),wposid(2,nwalrad))
   allocate(xnorwall(nphisol,nxsol),ynorwall(nphisol,nxsol),znorwall(nphisol,nxsol))

   allocate(xmemb(nmemb),ymemb(nmemb),zmemb(nmemb),jmbsol(2,nmemb))
   allocate(xsurf(nphisol,nxsol),ysurf(nphisol,nxsol),zsurf(nphisol,nxsol),nsurf(nphisol,nxsol))
   allocate(jsursol(2,nmemb))
   allocate(xnorsurf(nphisol,nxsol),ynorsurf(nphisol,nxsol),znorsurf(nphisol,nxsol))

   allocate(xfa(nact),yfa(nact),zfa(nact))
   allocate(astart(fanummax),alen(fanummax),jfasol(2,nact),fadist(nact))
   allocate(xmyo(nmyo),ymyo(nmyo),zmyo(nmyo),myhead(4,myonum),mybody(mybodylen,myonum),jmysol(2,nmyo),mydist(nmyo))
   allocate(xlk(nlk),ylk(nlk),zlk(nlk),lkstart(crlknum),jlksol(2,nlk),lkdist(nlk))

!  to tell if a bead is close to the membrane for constraint to apply

   fadist=1
   mydist=1
   lkdist=1

!  filament identification:
   allocate(filid(nact),apos(nact))


!  modeling binding of myosin to actin with fa2myo(jh,nb) which points to an actin bead
   allocate(fa2myo(4,myonum))

!  binding of crosslinkers to actin with fa2lk(jh,nb) which points to an actin bead
   allocate(fa2lk(2,crlknum))

!  tethering of actin to membrane:
   allocate(a2mem(nact))
   a2mem=0

!  crosslinker type to tell if it is deleted:

   allocate(lktyp(crlknum))

   call ringin(nstart,nxsol,nphisol,nmemb,fanum,nfa,myonum,nmyo,nwalrad,nwalsurf,nwall, &
               mybodylen,crlknum,nlk, &
                astart,alen,filid,a2mem,apos,fa1stbound,waltyp,lkstart,lktyp, &
                 mybody,myhead,mytyp,apar,wposid,fa2myo,fa2lk,xboundmin,xboundmax,dxsol,dphisol, &
                  xwall,ywall,zwall,xnorwall,ynorwall,znorwall,xwalsurf,ywalsurf,zwalsurf, &
                   xwalrad,ywalrad,zwalrad,xwalrep,ywalrep,zwalrep,xfa,yfa,zfa,xmyo,ymyo,zmyo, &
                    xlk,ylk,zlk,xmemb,ymemb,zmemb,xsurf,ysurf,zsurf,xnorsurf,ynorsurf,znorsurf)


!------------------------------------------------------------

!  forces:

   allocate(fxmemb(nmemb),fymemb(nmemb),fzmemb(nmemb))

   allocate(fxfa(nact),fyfa(nact),fzfa(nact))
   allocate(fxfarep(nact),fyfarep(nact),fzfarep(nact))

   allocate(fxmyo(nmyo),fymyo(nmyo),fzmyo(nmyo))
   allocate(fxmyorep(nmyo),fymyorep(nmyo),fzmyorep(nmyo))
   allocate(fxanglmyo(nmyo),fyanglmyo(nmyo),fzanglmyo(nmyo))

   allocate(fxlk(nlk),fylk(nlk),fzlk(nlk))
   allocate(fxlkrep(nlk),fylkrep(nlk),fzlkrep(nlk))
   allocate(fxangllk(nlk),fyangllk(nlk),fzangllk(nlk))




!  to add random forces:

   nrforce=1028

   allocate(rxmemb(nmemb,nrforce),rymemb(nmemb,nrforce),rzmemb(nmemb,nrforce))
   allocate(rxfa(nact,nrforce),ryfa(nact,nrforce),rzfa(nact,nrforce))
   allocate(rxmyo(nmyo,nrforce),rymyo(nmyo,nrforce),rzmyo(nmyo,nrforce))
   allocate(rxlk(nlk,nrforce),rylk(nlk,nrforce),rzlk(nlk,nrforce))


   jrforce=nrforce

   k_scale=0.5d0
   k_scale_lk=4*k_scale

!  to constrain filaments and myosin within the membrane:


!  set of random numbers for convenience:

   nrand=10000000!*jreset

   allocate(rands(nrand))
   jrand=nrand

!  actin neighbors of myo-heads and crosslinkers:
   neinum=100
   allocate(mynei(neinum,4,myonum),lknei(neinum,2,crlknum))
   neinum5=500
   allocate(mynei5(neinum5,4,myonum),lknei5(neinum5,2,crlknum))

!  to add volume exclusion effect:

   nnei=500

!  between membrane beads:
   allocate(pair5_mb2(2,nmemb*nnei))

!  between different components of the ring:
   allocate(pair5_myac(2,max(nact,nmyo)*nnei),pair5_lkac(2,max(nact,nlk)*nnei),pair5_mylk(2,max(nmyo,nlk)*nnei))

!  between same components:
   allocate(pair5_ac2(2,nact*nnei),pair5_my2(2,nmyo*nnei),pair5_lk2(2,nlk*nnei))



!  between membrane beads:
   allocate(pair_mb2(2,nmemb*100))


!  between different components of the ring:
   allocate(pair_myac(2,max(nact,nmyo)*100),pair_lkac(2,max(nact,nlk)*100),pair_mylk(2,max(nmyo,nlk)*100))

!  between same components:
   allocate(pair_ac2(2,nact*100),pair_my2(2,nmyo*100),pair_lk2(2,nlk*100))


!  list of tetrahedron used to preserve membrane layer:

   allocate(pairpart(2,nmemb*100))

!  boundary problem:

   allocate(boundtyp(nmemb))


!   allocate(boundpairpart(2,nboundpair))

!  exclusion parameters for the ring components:

   r_off=lbind

   rvdw=r_off-1.0d0

   r_on=rvdw+(r_off-rvdw)/10

   r_on2=r_on**2

   r_off2=r_off**2

   kvdw=kbind

   fvdwmax=kvdw*(r_off-r_on)**2/(r_on-rvdw)**2


!  parameters for attraction between membrane beads:

   kpair=2.0d0

   l_pair=2.0d0*l_mem


   thetpair=2.5d0*l_mem/mbrad

   cos_t_2=(1.0d0-(2*thetpair)**2/2)**2

   thet2by2=thetpair*thetpair/2


!  radial constraint:
   krad=0.001d0



   gap=0.1d0



!------------------------------------------------------------

!  To write out trajectory




   call random_number(r)
   junit=80*r+11
   junit1=junit+1
   junit2=junit+2

   jfile=jfile+1
   call dcdheader(junit,jfile,natom)

   nframe=0


         call writedcd(junit,nframe,natom,nxsol,nphisol,nwalsurf,nwalrad,nmemb,fanum1,fanum2, &
              natom_a1,natom_a2,nfa,natom_a2m,nmyo,nlk,alen,a2mem,astart,waltyp,wposid, &
              xwalsurf,ywalsurf,zwalsurf,xwalrep,ywalrep,zwalrep,xwalrad,ywalrad,zwalrad, &
              xfa,yfa,zfa,xmyo,ymyo,zmyo,xlk,ylk,zlk,xmemb,ymemb,zmemb,xwall,ywall,zwall)

!stop

12 format(f10.1,14x,i8,59x)

   open(junit2,file='measures.dat')

   if(runtime<1.0)then

      write(junit2,'(a91)')'      time    wallrate      natp  ringthick   ringwid  ringrad   G-actin  actconc  nxlinker'
      write(junit2,'(a91)')'       (s)      (nm/s)                 (nm)      (nm)     (nm)               (uM)          '

      natp=0

   else

      length=2

      read(junit2,'(a91)')a91
      read(junit2,'(a91)')a91

      do n=1,10000000

         read(junit2,12)timecheck,natp

         if(timecheck>runtime*0.000001)then

            exit
         end if

         length=length+1

      end do

      close(junit2)
      open(junit2,file='measures.dat')

      do n=1,length
         read(junit2,'(a91)')a91
      end do



   end if


   wallrate=0.0d0


   jprint=10000

   printtime=100000.0d0

   nstep=5000000000


   jsignal=0


   oldtime=runtime

!  initial setup solid angle indices:

   call solidset(nxsol,nphisol,nmemb,nfa,nmyo,nlk,jmbsol,jfasol,jmysol,jlksol, &
               pi,delta,dxsol,dphisol,xmemb,ymemb,zmemb,xfa,yfa,zfa,xmyo,ymyo,zmyo, &
                xlk,ylk,zlk,xwall,ywall,zwall,xnorwall,ynorwall,znorwall, &
                 xsurf,ysurf,zsurf,xnorsurf,ynorsurf,znorsurf,jsursol)

!=============================================================


!  Dynamics of the ring:

   print*,'start dynamics'

!  start time step:
   dt=0.0d0

!  ring boundary:

   xrmin=-rwid/2
   xrmax=-xrmin

!  for connecting two boundaries:

   shift=xboundmax-xboundmin+l_mem



   do jstep=1,nstep

      if(runtime>printtime)then
         p1_hydr=p1_hydr0
      end if

      if(mod(jstep-1,10)==0)then
         jforce1=1
      else
         jforce1=0
      end if

      if(mod(jstep-1,50)==0)then
         jforce2=1
      else
         jforce2=0
      end if

      if(mod(jstep-1,100)==0)then
         jforce3=1
      else
         jforce3=0
      end if


!     set random force


      if(jrforce==nrforce)then
         call rforceset(jrforce,nrforce,nmemb,nact,nmyo,nlk,pi,rxmemb,rymemb,rzmemb, &
              rxfa,ryfa,rzfa,rxmyo,rymyo,rzmyo,rxlk,rylk,rzlk,k_scale,k_scale_lk)
      end if


!     reset random numbers:

      if(jrand>=nrand)then
         call random_number(rands(:))
         jrand=0
      end if




      if(mod(jstep-1,nturnover)==0)then


!        add new F-actin:

         nmono=nfamax1-nfa1

         if(nmono>falen)then

            jdir=1

            call newactin(nmemb,myonum,crlknum,falen,nmono,jdir,nxsol,nphisol,ntether, &
               fanum,fanum1,fanum2,nfa,nfa1,nfa2,jsursol,apar,fa1stbound,astart,alen, &
                fadist,apos,filid,a2mem,jfasol,fa2myo,fa2lk,ltether,l_a,l_mem, &
                 pi,p_tether,xfa,yfa,zfa,xmemb,ymemb,zmemb, &
                  xnorsurf,ynorsurf,znorsurf,xsurf,ysurf,zsurf)





         end if

         nmono=nfamax2-nfa2

         if(nmono>falen)then

            jdir=-1

            call newactin(nmemb,myonum,crlknum,falen,nmono,jdir,nxsol,nphisol,ntether, &
               fanum,fanum1,fanum2,nfa,nfa1,nfa2,jsursol,apar,fa1stbound,astart,alen, &
                fadist,apos,filid,a2mem,jfasol,fa2myo,fa2lk,ltether,l_a,l_mem, &
                 pi,p_tether,xfa,yfa,zfa,xmemb,ymemb,zmemb, &
                  xnorsurf,ynorsurf,znorsurf,xsurf,ysurf,zsurf)





         end if


!        depolymerization of F-actin:

         if(j_dep==1)then
            call depoly(nfa,nfa1,nfa2,fanum,fanum1,fanum2,apar,apos,filid,a2mem,fadist,jfasol, &
                  astart,alen,fa1stbound,xfa,yfa,zfa, &
                  myonum,mytyp,fa2myo,crlknum,fa2lk,p_dep,dt,tension1,tension2)
         end if

!        myosin turnover:

         if(jmyoturn==1)then
            call myoturnover(nmyoturn,nfa,nmemb,myonum,mybodylen,nxsol,nphisol,mytyp,jmysol,mybody,myhead, &
               jsursol,jfasol,mydist,l_mb,l_mem,pi,pmyturn,dt,dphisol,xfa,yfa,zfa, &
                xmyo,ymyo,zmyo,xmemb,ymemb,zmemb,xnorsurf,ynorsurf,znorsurf,xsurf,ysurf,zsurf)
         end if

!        crosslinker turnover:

         call xlturnover(crlknum,crlknummax,jxldel,jxlturnover,nfa,crlknumactive,fa2lk,jfasol, &
               lkstart,lktyp,lkdist,jlksol,plk_remove,plk_turn,dt,dphisol, &
                xnorsurf,ynorsurf,znorsurf,xsurf,ysurf,zsurf,xfa,yfa,zfa,xlk,ylk,zlk)




!        setup pairs for long run:

         call allpairs(nmemb,xmemb,ymemb,zmemb,neinum5,mynei5,lknei5, &
               nfa,xfa,yfa,zfa,nmyo,xmyo,ymyo,zmyo,nlk,xlk,ylk,zlk,fanum,astart,alen, &
                myonum,myhead,mybody,mybodylen,crlknum,lklen,lkstart, &
                 npair5_myac,pair5_myac,npair5_lkac,pair5_lkac,npair5_mylk,pair5_mylk, &
                  npair5_ac2,pair5_ac2,npair5_my2,pair5_my2,npair5_lk2,pair5_lk2, &
                   npair5_mb2,pair5_mb2,5*r_off,cos_t_2)


      end if

      if(mod(jstep-1,10000)==0)then

!        cofilin breaking F-actin:

         if(cofilin==1)then
            call breakfa(fanum,fanum1,fanum2,astart,alen,apar,fa1stbound,filid,apos,fa2lk,xfa,yfa,zfa,l_a)
         end if



!        count neighbors:

         call neighbors(myonum,neinum,crlknum,lklen,neinum5,mybodylen,filid, &
              astart,alen,lkstart,myhead,mybody,mynei,lknei,mynei5,lknei5, &
              xfa,yfa,zfa,xmyo,ymyo,zmyo,xlk,ylk,zlk)



         call solidupdate(nxsol,nphisol,nmemb,nfa,nmyo,nlk,jmbsol,jfasol,jmysol,jlksol, &
               pi,delta,dxsol,dphisol,xmemb,ymemb,zmemb,xfa,yfa,zfa,xmyo,ymyo,zmyo, &
                xlk,ylk,zlk,xwall,ywall,zwall,xnorwall,ynorwall,znorwall, &
                 xsurf,ysurf,zsurf,xnorsurf,ynorsurf,znorsurf,jsursol)



!        pair setup for exclusion effect:

         call setpair(nmemb,xmemb,ymemb,zmemb,xfa,yfa,zfa,xmyo,ymyo,zmyo,xlk,ylk,zlk, &
               npair5_myac,pair5_myac,npair5_lkac,pair5_lkac,npair5_mylk,pair5_mylk, &
                npair_myac,pair_myac,npair_lkac,pair_lkac,npair_mylk,pair_mylk, &
                 npair5_ac2,pair5_ac2,npair5_my2,pair5_my2,npair5_lk2,pair5_lk2, &
                  npair_ac2,pair_ac2,npair_my2,pair_my2,npair_lk2,pair_lk2, &
                   npair5_mb2,pair5_mb2,npair_mb2,pair_mb2,r_off,l_pair,thet2by2, &
                    pairpart,boundtyp,l_mem,xboundmin,xboundmax,shift)



      end if

      if(mod(jstep-1,1000)==0)then


         call surfremod(nmemb,nthreads,nphisol,nxsol,nwall,jsursol,jmbsol,nsurf,wthick, &
               gap,dphisol,wallrate,xmemb,ymemb,zmemb,xsurf,ysurf,zsurf,xnorsurf, &
                ynorsurf,znorsurf,xwall,ywall,zwall,xnorwall,ynorwall,znorwall)


         jupdate=1

      else
         jupdate=0

      end if



!----------------------------------
!     changing myosin head status:

      if(mod(jstep-1,jskip)==0)then

         call myocycle(natp,myonum,neinum,tension2,nrand,jrand,mybodylen,apar,apos, &
              filid,fa1stbound,mytyp,fa2myo,myhead,mybody,mynei,p1_hydr,p2_bind, &
              p3_pi_rele,p4_adp_rele,p5_ubind,dt,invmylen2,lbind,rands,xmyo,ymyo,zmyo,xfa,yfa,zfa)


!---------------------------

!     update binding status of crosslinkers to actin

         call crlkcycle(jrand,nrand,crlknum,crlknum1,fanum,lklen,neinum,fa2lk,apar, &
              filid,astart,apos,alen,lkstart,lktyp,fa1stbound,lknei,plk_ubind1,plk_ubind2, &
              dt,plk_bind,lbind,rands,xlk,ylk,zlk,xfa,yfa,zfa)

      end if

!     crosslink release due to misorientation:

      if(crlkorient==1.and.mod(jstep-1,10000)==0)then

         do nl=1,crlknum

            if(fa2lk(1,nl)==0.or.fa2lk(2,nl)==0)then
               cycle
            end if

            ja1=fa2lk(1,nl)
            ja2=fa2lk(2,nl)

            j1=ja1-1
            j2=ja2-1

            dx1=xfa(ja1)-xfa(j1)
            dy1=yfa(ja1)-yfa(j1)
            dz1=zfa(ja1)-zfa(j1)

            dx2=xfa(ja2)-xfa(j2)
            dy2=yfa(ja2)-yfa(j2)
            dz2=zfa(ja2)-zfa(j2)

            prob=0.5d0-abs(dx1*dx2+dy1*dy2+dz1*dz2)*invl_a*invl_a

            if(prob>0.0d0)then

               call random_number(r)

               if(prob>r)then

                  fa2lk(1:2,nl)=0

                  apar(1:2,ja1)=0

                  apar(1:2,ja2)=0

                  jf1=filid(ja1)
                  jf2=filid(ja2)

                  jstart=astart(jf1)

                  if(apos(ja1)==fa1stbound(jf1))then

                     fa1stbound(jf1)=alen(jf1)

                     do j=1,alen(jf1)

                        if(apar(2,jstart+j-1)<0)then

                           nltem=-apar(2,jstart+j-1)

                           if(fa2lk(1,nltem)>0.and.fa2lk(2,nltem)>0)then
                              fa1stbound(jf1)=j
                              exit
                           end if

                        end if

                     end do

                  end if

                  jstart=astart(jf2)

                  if(apos(ja2)==fa1stbound(jf2))then

                     fa1stbound(jf2)=alen(jf2)

                     do j=1,alen(jf2)

                        if(apar(2,jstart+j-1)<0)then

                           nltem=-apar(2,jstart+j-1)

                           if(fa2lk(1,nltem)>0.and.fa2lk(2,nltem)>0)then
                              fa1stbound(jf2)=j
                              exit
                           end if

                        end if

                     end do

                  end if

               end if

            end if

         end do

      end if
!----------------------------------------------------------------

!     calculate forces:


!     rigidity of myosin:

      call myoforce(jforce2,fxmyo,fymyo,fzmyo,fxanglmyo,fyanglmyo,fzanglmyo, &
                    myonum,myhead,mytyp,mybody,mybodylen,xmyo,ymyo,zmyo,k_a,l_mh,l_mb, &
                    kmh_thet,kmb_thet,thet_mh1,thet_mh2,thet_mb,delta,invdelta,beta)


!     rigidity of crosslinkers:

      call lkforce(jforce3,fxlk,fylk,fzlk,fxangllk,fyangllk,fzangllk,crlknum1,crlknum2, &
           lkstart,lklen,xlk,ylk,zlk,k_lk,l_lk1,l_lk2,klk_thet,thet_lk,delta,invdelta,beta)


      if(jforce1==1)then

         jrforce=jrforce+1

         fxmemb(1:nmemb)=rxmemb(1:nmemb,jrforce)
         fymemb(1:nmemb)=rymemb(1:nmemb,jrforce)
         fzmemb(1:nmemb)=rzmemb(1:nmemb,jrforce)

         fxfarep(1:nfa)=rxfa(1:nfa,jrforce)
         fyfarep(1:nfa)=ryfa(1:nfa,jrforce)
         fzfarep(1:nfa)=rzfa(1:nfa,jrforce)

         fxmyorep(1:nmyo)=fxanglmyo(1:nmyo)+rxmyo(1:nmyo,jrforce)
         fymyorep(1:nmyo)=fyanglmyo(1:nmyo)+rymyo(1:nmyo,jrforce)
         fzmyorep(1:nmyo)=fzanglmyo(1:nmyo)+rzmyo(1:nmyo,jrforce)

         fxlkrep(1:nlk)=fxangllk(1:nlk)+rxlk(1:nlk,jrforce)
         fylkrep(1:nlk)=fyangllk(1:nlk)+rylk(1:nlk,jrforce)
         fzlkrep(1:nlk)=fzangllk(1:nlk)+rzlk(1:nlk,jrforce)




!        constraint forces: turgor, and actin-membrane tether, wall blocking, boundaries


         call constraints(nfa,nmemb,nmyo,nlk,jupdate,nxsol,nphisol,nthreads,a2mem, &
               fadist,mydist,lkdist,jmbsol,jsursol,jfasol,jmysol,jlksol,nsurf,ktether, &
                ltether,kwall,lsqueez,wthick,k_mem,l_mem,xwall,ywall,zwall,xnorwall, &
                 ynorwall,znorwall,xmemb,ymemb,zmemb,xboundmin,xboundmax,xrmin,xrmax,xfa,yfa,zfa, &
                  xmyo,ymyo,zmyo,xlk,ylk,zlk,xsurf,ysurf,zsurf,xnorsurf,ynorsurf, &
                   znorsurf,fxmemb,fymemb,fzmemb,fxmyorep,fymyorep,fzmyorep, &
                    fxlkrep,fylkrep,fzlkrep,fxfarep,fyfarep,fzfarep)



!        non-bond interaction: exclusion effect and coulomb

         call nonbond(npair_myac,pair_myac,npair_lkac,pair_lkac,npair_mylk,pair_mylk, &
            npair_ac2,pair_ac2,npair_my2,pair_my2,npair_lk2,pair_lk2, &
             npair_mb2,pair_mb2,pairpart,boundtyp, &
              xfa,yfa,zfa,fxfarep,fyfarep,fzfarep, &
               xmyo,ymyo,zmyo,fxmyorep,fymyorep,fzmyorep, &
                xlk,ylk,zlk,fxlkrep,fylkrep,fzlkrep, &
                 xmemb,ymemb,zmemb,fxmemb,fymemb,fzmemb, &
                  rvdw,kvdw,r_off,r_on2,r_off2,fvdwmax, &
                   kpair,l_pair,l_mem,kmemb,shift)


      end if

!     rigidity of F-actin:

      call faforce(jforce1,fxfa,fyfa,fzfa,fxfarep,fyfarep,fzfarep,fanum, &
           astart,alen,xfa,yfa,zfa,k_a,l_a,ka_thet,thet_fa,delta,invdelta,beta)


!     binding force between myosin and actin:

      call fabindmyo(fxfa,fyfa,fzfa,xfa,yfa,zfa,fxmyo,fymyo,fzmyo,xmyo,ymyo,zmyo, &
                     myhead,fa2myo,myonum,kbind,lbind,invl_a)

!     binding force between crosslinker and actin:

      call fabindlk(fxfa,fyfa,fzfa,xfa,yfa,zfa,fxlk,fylk,fzlk,xlk,ylk,zlk, &
                     lkstart,lklen,fa2lk,crlknum,kbind,lbind)




!--------------------------------------------------------


!     combining forces:


      fxfa(1:nfa)=fxfa(1:nfa)+fxfarep(1:nfa)
      fyfa(1:nfa)=fyfa(1:nfa)+fyfarep(1:nfa)
      fzfa(1:nfa)=fzfa(1:nfa)+fzfarep(1:nfa)

      fxmyo(1:nmyo)=fxmyo(1:nmyo)+fxmyorep(1:nmyo)
      fymyo(1:nmyo)=fymyo(1:nmyo)+fymyorep(1:nmyo)
      fzmyo(1:nmyo)=fzmyo(1:nmyo)+fzmyorep(1:nmyo)

      fxlk(1:nlk)=fxlk(1:nlk)+fxlkrep(1:nlk)
      fylk(1:nlk)=fylk(1:nlk)+fylkrep(1:nlk)
      fzlk(1:nlk)=fzlk(1:nlk)+fzlkrep(1:nlk)




!----------------------------------------------------

!     update coordinates:

      if(jforce1==1)then

         fmaxmb2=maxval(fxmemb*fxmemb+fymemb*fymemb+fzmemb*fzmemb)

         fmaxfa2=maxval(fxfa(1:nfa)*fxfa(1:nfa)+fyfa(1:nfa)*fyfa(1:nfa)+fzfa(1:nfa)*fzfa(1:nfa))

         fmaxmyo2=maxval(fxmyo*fxmyo+fymyo*fymyo+fzmyo*fzmyo)

         fmaxlk2=maxval(fxlk*fxlk+fylk*fylk+fzlk*fzlk)

         fmax2=max(fmaxmb2,fmaxfa2,fmaxmyo2,fmaxlk2)

         gam=0.01d0/sqrt(fmax2)

         dt=600*gam*jskip

         runtime=runtime+dt

      end if


      xmemb=xmemb+gam*fxmemb
      ymemb=ymemb+gam*fymemb
      zmemb=zmemb+gam*fzmemb


      xfa(1:nfa)=xfa(1:nfa)+gam*fxfa(1:nfa)
      yfa(1:nfa)=yfa(1:nfa)+gam*fyfa(1:nfa)
      zfa(1:nfa)=zfa(1:nfa)+gam*fzfa(1:nfa)

      xmyo=xmyo+gam*fxmyo
      ymyo=ymyo+gam*fymyo
      zmyo=zmyo+gam*fzmyo

!     treat free crosslinkers differently

!$omp parallel &
!$omp default(none) &
!$omp private(nl,jl) &
!$omp shared(crlknum,lkstart,lktyp,xlk,ylk,zlk,fxlk,fylk,fzlk,gam)
!$omp do

      do nl=1,crlknum

         jl=lkstart(nl)

         xlk(jl:jl+2)=xlk(jl:jl+2)+gam*fxlk(jl:jl+2)*lktyp(nl)
         ylk(jl:jl+2)=ylk(jl:jl+2)+gam*fylk(jl:jl+2)*lktyp(nl)
         zlk(jl:jl+2)=zlk(jl:jl+2)+gam*fzlk(jl:jl+2)*lktyp(nl)

      end do

!$omp end do
!$omp end parallel


      if(runtime-oldtime>printtime)then

         oldtime=runtime

         write(*,15)'step=',jstep,'run time=',runtime*0.000001,'sec','frame',nframe+1

         myorate=1.0*myonum/nmyoturn*runtime*0.000001

         write(*,*)'myo turnover rate',myorate

         call writedcd(junit,nframe,natom,nxsol,nphisol,nwalsurf,nwalrad,nmemb,fanum1,fanum2, &
              natom_a1,natom_a2,nfa,natom_a2m,nmyo,nlk,alen,a2mem,astart,waltyp,wposid, &
              xwalsurf,ywalsurf,zwalsurf,xwalrep,ywalrep,zwalrep,xwalrad,ywalrad,zwalrad, &
              xfa,yfa,zfa,xmyo,ymyo,zmyo,xlk,ylk,zlk,xmemb,ymemb,zmemb,xwall,ywall,zwall)




         if(jstep>10000000.and.jsignal==0)then
            jsignal=1
         elseif(jsignal==1)then
            jsignal=2
         end if

         call measures(junit2,nfa,crlknumactive,natp,jsignal,nmemb,nfamax1,nfamax2,crlknummax, &
              wallrate,aconc1,aconc2,lkconc,thet2by2,cos_t_2,printtime,runtime, &
              pi,l_mem,xfa,yfa,zfa,ymemb,zmemb)



         if(nframe>=200)then
            close(junit)

            call ringout(nstart,natom,natom_a1,natom_a2,natom_a2m,nwall,nxsol,nphisol,nwalsurf, &
               nwalrad,nmemb,nact1,nact2,nfamax1,nfamax2,nfa1,nfa2,fanum1,fanum2,ntether,nfa,nmyoturn, &
                nmyo,myonum,mybodylen,nlk,crlknum1,crlknum2,crlknummax,crlknumactive,astart,alen,filid,a2mem,apos, &
                 fa1stbound,waltyp,lkstart,lktyp,mybody,myhead,mytyp,fa2myo,fa2lk,apar,wposid, &
                  dxsol,dphisol,xboundmin,xboundmax,xsurf,ysurf,zsurf,xnorsurf,ynorsurf,znorsurf, &
                   xfa,yfa,zfa,xmyo,ymyo,zmyo,xlk,ylk,zlk,xmemb,ymemb,zmemb,xwall,ywall,zwall, &
                    xnorwall,ynorwall,znorwall,xwalsurf,ywalsurf,zwalsurf, &
                     xwalrad,ywalrad,zwalrad,xwalrep,ywalrep,zwalrep)




            if(runtime*0.000001>140.0) exit

            call system_clock(timerun)

            OPEN(junit1,FILE='restart.inp')
            WRITE(junit1,*)NSTART,JFILE
            WRITE(junit1,*)TIMERUN-TIMESTART+time0,runtime
            CLOSE(junit1)

            if((nstep-jstep)/jprint>0)then
               jfile=jfile+1
               call dcdheader(junit,jfile,natom)
               nframe=0
            end if
         end if


      end if
   end do

15 format(a5,2x,i10,2x,a9,2x,f10.6,1x,a3,2x,a5,2x,i3)
!-----------------------------------------------------------



   call system_clock(timerun)

   TIMERUN=(TIMERUN-TIMESTART)/rate+time0

   days=timerun/86400
   timerun=timerun-86400*days

   hours=timerun/3600
   timerun=timerun-3600*hours

   mins=timerun/60
   secs=timerun-60*mins

   write(*,1)'RUNNING TIME =',DAYS,'days : ',HOURS,'hours : ',MINS,'mins : ',secs,'secs'

1  FORMAT(A14,2X,I2,1X,A7,I2,1X,A8,1X,I2,1X,A7,1x,I2,1x,A4)


2  FORMAT(A4,2X,I9,4X,A15,2X,F6.3,2x,F6.3)

end
