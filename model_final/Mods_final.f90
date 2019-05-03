module cdk

contains


!---------------------------------------------------------

   SUBROUTINE WHATTIME(TIMERUN)

   IMPLICIT NONE
   INTEGER TIMES(10),YEARS,MONS,TIMERUN,DAYS,HOURS,MINS,SECS,LASTTIME,TIMESTART,NM

   call date_and_time(values=times)

   mins=times(6); hours=times(5); days=times(3)-1; mons=times(2); years=times(1)
   timerun=mins+60*hours+60*24*days

   if(mons>1)then
      do nm=1,mons-1
         if((mod(nm,2)==1.and.nm<=7).or.(mod(nm,2)==0.and.nm>=8))then
            timerun=timerun+60*24*31
         elseif(nm==2.and.mod(years,4)==0)then
            timerun=timerun+60*24*29
         elseif(nm==2.and.mod(years,4)/=0)then
            timerun=timerun+60*24*28
         else
            timerun=timerun+60*24*30
         end if
      end do
   end if

   if(mod(years-1,4)==0)then
      timerun=timerun+60*24*366*(years-1)
   else
      timerun=timerun+60*24*365*(years-1)
   end if

   END SUBROUTINE

!==================================================

   subroutine setjob(nthreads,jobid,mbrad,mwid,rrad,rwid,rthick,falen,fanum,myp2num,crlknum1,crlknum2, &
              myo2num,kmemb,tension1,tension2,crlkorient,j_dep,cofilin,jmyoturn,jxldel,jxlturnover,p_scale,p_rate,halftime)

   implicit none

   integer nthreads,jobid,fanum,falen,myp2num,crlknum1,crlknum2,myo2num
   integer tension1,tension2,crlkorient,j_dep,cofilin,jmyoturn,jxldel,jxlturnover

   character (len=100) chara
   double precision mbrad,mwid,rrad,rwid,rthick,kmemb,p_scale,p_rate,halftime

2  format(a50)
   open(1,file='inputjob.inp')

13  read(1,2)chara
   if(chara(1:4)/='CPUs')then
      goto 13
   end if
   read(1,*)nthreads

3  read(1,2)chara
   if(chara(1:4)/='JOBI')then
      goto 3
   end if
   read(1,*)jobid

4  read(1,2)chara
   if(chara(1:4)/='RRAD')then
      goto 4
   end if
   read(1,*)mbrad

!  Matt's tomograms show 30 nm gap:
   rrad=mbrad-30.0d0

5  read(1,2)chara
   if(chara(1:4)/='RWID')then
      goto 5
   end if
   read(1,*)rwid

35 read(1,2)chara
   if(chara(1:4)/='MWID')then
      goto 35
   end if
   read(1,*)mwid

!   mwid=2*rwid

15 read(1,2)chara
   if(chara(1:4)/='RTHI')then
      goto 15
   end if
   read(1,*)rthick

7  read(1,2)chara
   if(chara(1:4)/='FALE')then
      goto 7
   end if
   read(1,*)falen

9  read(1,2)chara
   if(chara(1:4)/='FANU')then
      goto 9
   end if
   read(1,*)fanum

19 read(1,2)chara
   if(chara(1:4)/='MYO2')then
      goto 19
   end if
   read(1,*)myo2num

10 read(1,2)chara
   if(chara(1:4)/='MYP2')then
      goto 10
   end if
   read(1,*)myp2num

11 read(1,2)chara
   if(chara(1:4)/='CRLK')then
      goto 11
   end if
   read(1,*)crlknum1
   read(1,*)crlknum2


23 read(1,*)chara

   if(chara(1:4)/='BEND')then
      goto 23
   end if

   read(1,*)kmemb


24 read(1,*)chara

   if(chara(1:4)/='TEDE')then
      goto 24
   end if

   read(1,*)tension1

25 read(1,*)chara

   if(chara(1:4)/='TEMY')then
      goto 25
   end if

   read(1,*)tension2

26 read(1,*)chara

   if(chara(1:4)/='CROR')then
      goto 26
   end if

   read(1,*)crlkorient

27 read(1,*)chara

   if(chara(1:4)/='DEPO')then
      goto 27
   end if

   read(1,*)j_dep

17 read(1,*)chara

   if(chara(1:4)/='MYOT')then
      goto 17
   end if

   read(1,*)jmyoturn

18 read(1,*)chara

   if(chara(1:4)/='XLDE')then
      goto 18
   end if

   read(1,*)jxldel

38 read(1,*)chara

   if(chara(1:4)/='XLTU')then
      goto 38
   end if

   read(1,*)jxlturnover

28 read(1,*)chara

   if(chara(1:4)/='COFI')then
      goto 28
   end if

   read(1,*)cofilin

34 read(1,*)chara

   if(chara(1:4)/='RATE')then
      goto 34
   end if

   read(1,*)p_rate

36 read(1,*)chara

   if(chara(1:4)/='SCAL')then
      goto 36
   end if

   read(1,*)p_scale


48 read(1,*)chara

   if(chara(1:4)/='HALF')then
      goto 48
   end if

   read(1,*)halftime

   close(1)

   end subroutine

!=========================================================

   subroutine makering(nwall,nxsol,nphisol,nwalsurf,nwalrad,nmemb,fanum1,fanum2,nfa1,nfa2,myp2num,myp2len, &
               nmyp2,nbondmyp2,myo2num,myo2len,nmyo2,nbondmyo2,crlknum1,crlknum2,lklen,nlk,nbondlk, &
                astart,alen,filid,apos,lkstart,my2mem,myp2body,myp2head,myo2body,myo2head,bondmyp2,bondmyo2, &
                 bondlk,wposid,xboundmin,xboundmax,dxsol,dphisol,xwall,ywall,zwall,xnorwall,ynorwall,znorwall, &
                  xwalsurf,ywalsurf,zwalsurf,xwalrad,ywalrad,zwalrad,xwalrep,ywalrep,zwalrep, &
                   xmemb,ymemb,zmemb,xsurf,ysurf,zsurf,xfa,yfa,zfa,xmyp2,ymyp2,zmyp2,xmyo2,ymyo2,zmyo2,xlk,ylk,zlk)




   implicit none

   integer,value::nwall,nxsol,nphisol,nwalsurf,nwalrad,nmemb,fanum1,fanum2,nfa1,nfa2
   integer,value::myp2num,myp2len,nmyp2,nbondmyp2,myo2num,myo2len,nmyo2,nbondmyo2
   integer,value::crlknum1,crlknum2,lklen,nlk,nbondlk

   integer natom,natom_a1,natom_a2,nact1,nact2,nfamax1,nfamax2,nmyoturn

   integer fanum,nfa,nbond_a1,nbond_a2,nbond_a2m,nbond0,jmb
   integer jstart,na,nwalltemp
   integer jm,nm,j,crlknum,crlknummax,crlknumactive
   integer izero,nbond,ires,n,jatom
   integer j1,j2,j3,j4,j5,j6,j7,j8,nline,i,jx,jp

   integer,allocatable,dimension(:),intent(in)::astart,alen,filid,apos,lkstart,my2mem
   integer,allocatable,dimension(:)::fa1stbound,waltyp,lktyp
   integer,allocatable,dimension(:,:),intent(in)::myp2body,myp2head,myo2body,myo2head,bondmyp2,bondmyo2,bondlk,wposid
   integer,allocatable,dimension(:,:)::bond,myp2typ,myo2typ,apar
   integer,allocatable,dimension(:,:)::fa2myp2,fa2myo2,fa2lk

   double precision,value::dxsol,dphisol,xboundmin,xboundmax
   double precision charg,mass,w1,w2

   double precision,allocatable,dimension(:,:),intent(in)::xwall,ywall,zwall
   double precision,allocatable,dimension(:,:),intent(in)::xnorwall,ynorwall,znorwall
   real,allocatable,dimension(:),intent(in)::xwalsurf,ywalsurf,zwalsurf
   real,allocatable,dimension(:),intent(in)::xwalrad,ywalrad,zwalrad,xwalrep,ywalrep,zwalrep
   double precision,allocatable,dimension(:),intent(in)::xmemb,ymemb,zmemb
   double precision,allocatable,dimension(:,:),intent(in)::xsurf,ysurf,zsurf
   double precision,allocatable,dimension(:),intent(in)::xfa,yfa,zfa,xmyp2,ymyp2,zmyp2,xmyo2,ymyo2,zmyo2,xlk,ylk,zlk

   real zero,x0,y0,z0
   real,allocatable,dimension(:,:)::xtem,ytem,ztem

   character(8) tex,typ,res,segid,resno


   open(1,file='ring.psf')

20  FORMAT(I8,1X,A4,1X,A4,1X,A3,2X,A3,2X,A4,2X,F9.6,6X,F8.4,11X,I1)
21  FORMAT(I8,1X,A)

   open(2,file='startring.pdb')

42 FORMAT(A6,I5,2X,A3,1X,A3,1X,I5,4X,3F8.3,2F6.2,6X,A4)
52 FORMAT(A50)


!--------------------------------------------

!--------------------------------------------

!  cell wall and membrane:

   nwalltemp=4*nwall-2*nphisol

   natom=nwalsurf+nwalltemp+nwalrad+nmemb

   nbond=0

   nfa=nfa1+nfa2

!  actin clockwise polarity:

   nact1=nfa1


   natom=natom+2*nact1

   natom_a1=natom

   nfamax1=nact1

   nbond=nbond+nact1

   nbond_a1=nbond

!  actin counter-clockwise polarity:

   nact2=nfa2

   natom=natom+2*nact2

   natom_a2=natom

   nfamax2=nact2

   nbond=nbond+nact2

   nbond_a2=nbond

!  tethering of actin to membrane:

!   natom=natom+2*ntether

!   natom_a2m=natom

!   nbond=nbond+ntether

!   nbond_a2m=nbond


!  Myp2 myosins:

   natom=natom+nmyp2

   nbond=nbond+nbondmyp2

!  Myo2 myosins:

   natom=natom+nmyo2

   nbond=nbond+nbondmyo2

!  to visualize tethering of Myo2 myosin to membrane:

   natom=natom+2*myo2num
   nbond=nbond+myo2num


!  crosslinkers:

   natom=natom+nlk

   nbond=nbond+nbondlk

   allocate(bond(2,nbond))



   write(1,21)natom,'!NATOM'

   charg=0.0d0
   mass=1.0d0
   izero=0
   tex='ATOM'
   w1=1.0d0
   w2=0.0d0
   zero=0.0

!  for cell wall

   res='WAL'
   typ='WAL'
   segid='WAL'

   ires=1
   write(resno,'(i1)')ires

   jatom=0

!  start with the surface:

   do n=1,nwalsurf

      jatom=jatom+1

      write(1,20)jatom,segid,resno,res,typ,tex,charg,mass,izero

      write(2,42)tex,jatom,typ,res,ires,xwalsurf(n),ywalsurf(n),zwalsurf(n),w1,w2,segid

   end do

   allocate(xtem(nphisol,nxsol),ytem(nphisol,nxsol),ztem(nphisol,nxsol))

   do jx=1,nxsol

      x0=xwall(nphisol,jx)
      y0=ywall(nphisol,jx)
      z0=zwall(nphisol,jx)

      do jp=1,nphisol

         jatom=jatom+1

         write(1,20)jatom,segid,resno,res,typ,tex,charg,mass,izero

         write(2,42)tex,jatom,typ,res,ires,xwall(jp,jx)/10,ywall(jp,jx)/10,zwall(jp,jx)/10,w1,w2,segid

         xtem(jp,jx)=0.05*(x0+xwall(jp,jx))
         ytem(jp,jx)=0.05*(y0+ywall(jp,jx))
         ztem(jp,jx)=0.05*(z0+zwall(jp,jx))

         jatom=jatom+1

         write(1,20)jatom,segid,resno,res,typ,tex,charg,mass,izero

         write(2,42)tex,jatom,typ,res,ires,xtem(jp,jx),ytem(jp,jx),ztem(jp,jx),w1,w2,segid

         x0=xwall(jp,jx)
         y0=ywall(jp,jx)
         z0=zwall(jp,jx)

      end do

   end do

   do jx=1,nxsol-1

      do jp=1,nphisol

         jatom=jatom+1

         write(1,20)jatom,segid,resno,res,typ,tex,charg,mass,izero

         x0=0.05*(xwall(jp,jx)+xwall(jp,jx+1))
         y0=0.05*(ywall(jp,jx)+ywall(jp,jx+1))
         z0=0.05*(zwall(jp,jx)+zwall(jp,jx+1))

         write(2,42)tex,jatom,typ,res,ires,x0,y0,z0,w1,w2,segid

         jatom=jatom+1

         write(1,20)jatom,segid,resno,res,typ,tex,charg,mass,izero

         x0=0.5*(xtem(jp,jx)+xtem(jp,jx+1))
         y0=0.5*(ytem(jp,jx)+ytem(jp,jx+1))
         z0=0.5*(ztem(jp,jx)+ztem(jp,jx+1))

         write(2,42)tex,jatom,typ,res,ires,x0,y0,z0,w1,w2,segid

      end do

   end do

   deallocate(xtem,ytem,ztem)

!  radial wall:

   do n=1,nwalrad

      jatom=jatom+1

      write(1,20)jatom,segid,resno,res,typ,tex,charg,mass,izero

      write(2,42)tex,jatom,typ,res,ires,xwalrad(n)/10,ywalrad(n)/10,zwalrad(n)/10,w1,w2,segid

   end do



!  for membrane:

   res='MBR'
   typ='MBR'
   segid='MBR'

   ires=2
   write(resno,'(i1)')ires

   do n=1,nmemb

      jatom=jatom+1

      write(1,20)jatom,segid,resno,res,typ,tex,charg,mass,izero

      write(2,42)tex,jatom,typ,res,ires,xmemb(n)/10,ymemb(n)/10,zmemb(n)/10,w1,w2,segid


   end do

   nbond=0

!  for F-actin

   res='FAC'
   typ='FAC'

   ires=3
   write(resno,'(i1)')ires


!  clockwise:

   segid='FA1'

   do n=1,fanum1

      jstart=astart(n)

      do j=1,alen(n)-1

         jatom=jatom+1

         na=jstart+j-1

         write(1,20)jatom,segid,resno,res,typ,tex,charg,mass,izero

         write(2,42)tex,jatom,typ,res,ires,xfa(na)/10,yfa(na)/10,zfa(na)/10,w1,w2,segid

         jatom=jatom+1

         na=jstart+j

         write(1,20)jatom,segid,resno,res,typ,tex,charg,mass,izero

         write(2,42)tex,jatom,typ,res,ires,xfa(na)/10,yfa(na)/10,zfa(na)/10,w1,w2,segid

         nbond=nbond+1

         bond(1,nbond)=jatom-1

         bond(2,nbond)=jatom

      end do

   end do


   if(nbond<nbond_a1)then

      nbond0=nbond

      do n=nbond0+1,nbond_a1

         jatom=jatom+1

         write(1,20)jatom,segid,resno,res,typ,tex,charg,mass,izero

         write(2,42)tex,jatom,typ,res,ires,zero,zero,zero,w1,w2,segid

         jatom=jatom+1

         write(1,20)jatom,segid,resno,res,typ,tex,charg,mass,izero

         write(2,42)tex,jatom,typ,res,ires,zero,zero,zero,w1,w2,segid

         nbond=nbond+1

         bond(1,nbond)=jatom-1

         bond(2,nbond)=jatom

      end do

   end if

!print*,'b1',jatom

!  counter-clockwise:


   segid='FA2'

!   fanum=fanum1+fanum2

   do n=1,fanum2


      jstart=astart(n+fanum1)

      do j=1,alen(n+fanum1)-1

         jatom=jatom+1

         na=jstart+j-1

         write(1,20)jatom,segid,resno,res,typ,tex,charg,mass,izero

         write(2,42)tex,jatom,typ,res,ires,xfa(na)/10,yfa(na)/10,zfa(na)/10,w1,w2,segid

         jatom=jatom+1

         na=jstart+j

         write(1,20)jatom,segid,resno,res,typ,tex,charg,mass,izero

         write(2,42)tex,jatom,typ,res,ires,xfa(na)/10,yfa(na)/10,zfa(na)/10,w1,w2,segid

         nbond=nbond+1

         bond(1,nbond)=jatom-1

         bond(2,nbond)=jatom

      end do

   end do


   if(nbond<nbond_a2)then

      nbond0=nbond

      do n=nbond0+1,nbond_a2

         jatom=jatom+1

         write(1,20)jatom,segid,resno,res,typ,tex,charg,mass,izero

         write(2,42)tex,jatom,typ,res,ires,zero,zero,zero,w1,w2,segid

         jatom=jatom+1

         write(1,20)jatom,segid,resno,res,typ,tex,charg,mass,izero

         write(2,42)tex,jatom,typ,res,ires,zero,zero,zero,w1,w2,segid

         nbond=nbond+1

         bond(1,nbond)=jatom-1

         bond(2,nbond)=jatom

      end do

   end if

!print*,'b2',jatom


!  to visualize membrane-tethered F-actin:

!   nfa=nfa1+nfa2


!   res='A2M'
!   typ='A2M'
!   segid='A2M'

!   ires=4
!   write(resno,'(i1)')ires


!   do n=1,nfa

!      if(a2mem(n)==0)then
!         cycle
!      end if

!      jatom=jatom+1

!      write(1,20)jatom,segid,resno,res,typ,tex,charg,mass,izero

!      write(2,42)tex,jatom,typ,res,ires,xfa(n)/10,yfa(n)/10,zfa(n)/10,w1,w2,segid

!      jatom=jatom+1

!      write(1,20)jatom,segid,resno,res,typ,tex,charg,mass,izero

!      jm=a2mem(n)

!      write(2,42)tex,jatom,typ,res,ires,xmemb(jm)/10,ymemb(jm)/10,zmemb(jm)/10,w1,w2,segid

!      nbond=nbond+1

!      bond(1,nbond)=jatom-1

!      bond(2,nbond)=jatom

!   end do


!   if(nbond<nbond_a2m)then

!      nbond0=nbond

!      do n=nbond0+1,nbond_a2m

!         jatom=jatom+1

!         write(1,20)jatom,segid,resno,res,typ,tex,charg,mass,izero

!         write(2,42)tex,jatom,typ,res,ires,zero,zero,zero,w1,w2,segid

!         jatom=jatom+1

!         write(1,20)jatom,segid,resno,res,typ,tex,charg,mass,izero

!         write(2,42)tex,jatom,typ,res,ires,zero,zero,zero,w1,w2,segid

!         nbond=nbond+1

!         bond(1,nbond)=jatom-1

!         bond(2,nbond)=jatom

!      end do

!   end if

!print*,'a2m',jatom

!  for Myp2 myosin

   res='MYP'
   typ='MYP'

   ires=4
   write(resno,'(i1)')ires


   n=0

   do nm=1,myp2num

      segid='MYB'

      do j=1,myp2len

         n=n+1

         jatom=jatom+1

         write(1,20)jatom,segid,resno,res,typ,tex,charg,mass,izero

         write(2,42)tex,jatom,typ,res,ires,xmyp2(n)/10,ymyp2(n)/10,zmyp2(n)/10,w1,w2,segid

      end do

      segid='MYH'

      do j=1,4

         jatom=jatom+1

         n=n+1

         write(1,20)jatom,segid,resno,res,typ,tex,charg,mass,izero

         write(2,42)tex,jatom,typ,res,ires,xmyp2(n)/10,ymyp2(n)/10,zmyp2(n)/10,w1,w2,segid

      end do

   end do

   do n=1,nbondmyp2

      nbond=nbond+1

      bond(1,nbond)=bondmyp2(1,n)+jatom-nmyp2
      bond(2,nbond)=bondmyp2(2,n)+jatom-nmyp2

   end do


!  for Myo2 myosin

   res='MYO'
   typ='MYO'

   ires=5
   write(resno,'(i1)')ires


   n=0

   do nm=1,myo2num

      segid='MYB'

      do j=1,myo2len

         n=n+1

         jatom=jatom+1

         write(1,20)jatom,segid,resno,res,typ,tex,charg,mass,izero

         write(2,42)tex,jatom,typ,res,ires,xmyo2(n)/10,ymyo2(n)/10,zmyo2(n)/10,w1,w2,segid

      end do

      segid='MYH'

      do j=1,2

         jatom=jatom+1

         n=n+1

         write(1,20)jatom,segid,resno,res,typ,tex,charg,mass,izero

         write(2,42)tex,jatom,typ,res,ires,xmyo2(n)/10,ymyo2(n)/10,zmyo2(n)/10,w1,w2,segid

      end do

   end do

   do n=1,nbondmyo2

      nbond=nbond+1

      bond(1,nbond)=bondmyo2(1,n)+jatom-nmyo2
      bond(2,nbond)=bondmyo2(2,n)+jatom-nmyo2

   end do


!  binding of myosin to membrane:

   segid='M2M'
   typ='M2M'
   res='M2M'
   ires=6

   write(resno,'(i1)')ires

   do nm=1,myo2num

      jatom=jatom+1

      write(1,20)jatom,segid,resno,res,typ,tex,charg,mass,izero

      jm=myo2body(1,nm)

      write(2,42)tex,jatom,typ,res,ires,xmyo2(jm)/10,ymyo2(jm)/10,zmyo2(jm)/10,w1,w2,segid

      jatom=jatom+1

      write(1,20)jatom,segid,resno,res,typ,tex,charg,mass,izero

      jmb=my2mem(nm)

      write(2,42)tex,jatom,typ,res,ires,xmemb(jmb)/10,ymemb(jmb)/10,zmemb(jmb)/10,w1,w2,segid

      nbond=nbond+1

      bond(1,nbond)=jatom-1

      bond(2,nbond)=jatom

   end do




!  for crosslinkers:

   res='CLK'
   typ='CLK'
   segid='CLK'

   ires=7
   write(resno,'(i1)')ires

   do n=1,nlk

      jatom=jatom+1

      write(1,20)jatom,segid,resno,res,typ,tex,charg,mass,izero

      write(2,42)tex,jatom,typ,res,ires,xlk(n)/10,ylk(n)/10,zlk(n)/10,w1,w2,segid

   end do

   do n=1,nbondlk

      nbond=nbond+1

      bond(1,nbond)=bondlk(1,n)+jatom-nlk
      bond(2,nbond)=bondlk(2,n)+jatom-nlk

   end do

   tex='END'

   write(2,'(a3)')'END'

   close(2)

!-------------------------------------------
! -- WRITE LIST OF BONDS:

   NLINE=NBOND/4

   WRITE(1,*)
   WRITE(1,21)NBOND,'!NBOND: bonds'

   DO N=1,NLINE
      I=1+4*(N-1)
      J1=BOND(1,I); J2=BOND(2,I); J3=BOND(1,I+1); J4=BOND(2,I+1)
      J5=BOND(1,I+2); J6=BOND(2,I+2); J7=BOND(1,I+3); J8=BOND(2,I+3)

      WRITE(1,'(8I8)')J1,J2,J3,J4,J5,J6,J7,J8
   END DO

   IF(MOD(NBOND,4)==1)THEN
      J1=BOND(1,NBOND); J2=BOND(2,NBOND)
      WRITE(1,'(2I8)')J1,J2

   ELSE IF(MOD(NBOND,4)==2)THEN
      J1=BOND(1,NBOND-1); J2=BOND(2,NBOND-1)
      J3=BOND(1,NBOND); J4=BOND(2,NBOND)
      WRITE(1,'(4I8)')J1,J2,J3,J4

   ELSE IF(MOD(NBOND,4)==3)THEN
      J1=BOND(1,NBOND-2); J2=BOND(2,NBOND-2)
      J3=BOND(1,NBOND-1); J4=BOND(2,NBOND-1)
      J5=BOND(1,NBOND); J6=BOND(2,NBOND)
      WRITE(1,'(6I8)')J1,J2,J3,J4,J5,J6

   END IF

   CLOSE(1)

!---------------------------------------------

!  write coordinate

   open(1,file='rcoor000.inp',form='unformatted')

   write(1)dxsol,dphisol

   write(1)xwall(1:nphisol,1:nxsol)
   write(1)ywall(1:nphisol,1:nxsol)
   write(1)zwall(1:nphisol,1:nxsol)
   write(1)xnorwall(1:nphisol,1:nxsol)
   write(1)ynorwall(1:nphisol,1:nxsol)
   write(1)znorwall(1:nphisol,1:nxsol)



   write(1)xwalsurf(1:nwalsurf)
   write(1)ywalsurf(1:nwalsurf)
   write(1)zwalsurf(1:nwalsurf)

   write(1)xwalrad(1:nwalrad)
   write(1)ywalrad(1:nwalrad)
   write(1)zwalrad(1:nwalrad)

   write(1)xwalrep(1:nwalrad)
   write(1)ywalrep(1:nwalrad)
   write(1)zwalrep(1:nwalrad)



   write(1)nmemb
   write(1)xmemb(1:nmemb)
   write(1)ymemb(1:nmemb)
   write(1)zmemb(1:nmemb)
   write(1)xboundmin,xboundmax
   write(1)xsurf(1:nphisol,1:nxsol)
   write(1)ysurf(1:nphisol,1:nxsol)
   write(1)zsurf(1:nphisol,1:nxsol)
   write(1)xnorwall(1:nphisol,1:nxsol)
   write(1)ynorwall(1:nphisol,1:nxsol)
   write(1)znorwall(1:nphisol,1:nxsol)


   write(1)nfa
   write(1)xfa(1:nfa)
   write(1)yfa(1:nfa)
   write(1)zfa(1:nfa)

   write(1)nmyp2
   write(1)xmyp2(1:nmyp2)
   write(1)ymyp2(1:nmyp2)
   write(1)zmyp2(1:nmyp2)

   write(1)nmyo2
   write(1)xmyo2(1:nmyo2)
   write(1)ymyo2(1:nmyo2)
   write(1)zmyo2(1:nmyo2)


   write(1)nlk
   write(1)xlk(1:nlk)
   write(1)ylk(1:nlk)
   write(1)zlk(1:nlk)

   close(1)

!---------------------------------------------

!  write configuration

   open(1,file='rconf000.inp')

   write(1,*)'VISUAL'

   write(1,*)natom,natom_a1,natom_a2

   write(1,*)'CELLWALL'
   write(1,*)nwall,nxsol,nphisol,nwalsurf,nwalrad

   allocate(waltyp(nwalrad))
   waltyp=0

   write(1,*)waltyp(1:nwalrad)

   deallocate(waltyp)

   write(1,*)wposid(1,1:nwalrad)
   write(1,*)wposid(2,1:nwalrad)

   write(1,*)

   write(1,*)'MEMBRANE'
   write(1,*)nmemb
   write(1,*)


!print*,'nfa=',nfa
   write(1,*)'FACTIN'
   write(1,*)nact1,nact2,nfamax1,nfamax2,nfa1,nfa2,fanum1,fanum2
   write(1,*)filid(1:nfa)
!   write(1,*)a2mem(1:nfa)

   allocate(apar(2,nfa))
   apar=0

   fanum=fanum1+fanum2

   apar(1,astart(1:fanum))=-1

   write(1,*)apar(1:2,1:nfa)
   deallocate(apar)

   write(1,*)apos(1:nfa)

   allocate(fa1stbound(fanum))
   fa1stbound(1:fanum)=alen(1:fanum)

   write(1,*)fa1stbound(1:fanum)
   deallocate(fa1stbound)


   write(1,*)alen(1:fanum)

   write(1,*)astart(1:fanum)

   write(1,*)


   write(1,*)'MYP2 Myosin'
   write(1,*)nmyp2,myp2num

   allocate(myp2typ(4,myp2num),fa2myp2(4,myp2num))
   myp2typ=1; fa2myp2=0

   do n=1,myp2num
      write(1,*)myp2head(1:4,n),myp2body(1:myp2len,n)
      write(1,*)myp2typ(1:4,n),fa2myp2(1:4,n)
   end do
   write(1,*)

   deallocate(myp2typ,fa2myp2)

   nmyoturn=0

   write(1,*)'MYO2 Myosin'
   write(1,*)nmyo2,myo2num,nmyoturn

   allocate(myo2typ(2,myo2num),fa2myo2(2,myo2num))
   myo2typ=1; fa2myo2=0

   do n=1,myo2num
      write(1,*)my2mem(n),myo2head(1:2,n),myo2body(1:myo2len,n)
      write(1,*)myo2typ(1:2,n),fa2myo2(1:2,n)
   end do
   write(1,*)

   deallocate(myo2typ,fa2myo2)



   crlknum=crlknum1+crlknum2

   crlknummax=crlknum

   crlknumactive=crlknum

   write(1,*)'CROSSLINKER'
   write(1,*)nlk,crlknum1,crlknum2,crlknummax,crlknumactive

   allocate(fa2lk(2,crlknum),lktyp(crlknum))
   fa2lk=0
   lktyp=1

   write(1,*)lkstart(1:crlknum)
   write(1,*)lktyp(1:crlknum)
   write(1,*)fa2lk(1,1:crlknum)
   write(1,*)fa2lk(2,1:crlknum)



   deallocate(fa2lk,lktyp)


   close(1)

   open(1,file='restart.inp')

   write(1,*)0,0
   write(1,*)0,0.0
   close(1)

   end subroutine

!=========================================================
   subroutine writedcd(junit,nframe,natom,nxsol,nphisol,nwalsurf,nwalrad,nmemb,fanum1,fanum2, &
              natom_a1,natom_a2,nfa,nmyp2,nmyo2,myo2num,nlk,alen,astart,my2mem,waltyp,wposid,myo2body, &
              xwalsurf,ywalsurf,zwalsurf,xwalrep,ywalrep,zwalrep,xwalrad,ywalrad,zwalrad, &
              xfa,yfa,zfa,xmyp2,ymyp2,zmyp2,xmyo2,ymyo2,zmyo2,xlk,ylk,zlk,xmemb,ymemb,zmemb,xwall,ywall,zwall)



   implicit none

   integer,value::junit,nxsol,natom,nphisol,nwalsurf,nwalrad
   integer,value::nmemb,fanum1,fanum2,natom_a1,natom_a2,nfa,nmyp2,nmyo2,myo2num,nlk
   integer nframe
   integer n,j,jw,jm,jp,jx,jstart,jmemb

   integer,allocatable,intent(in),dimension(:)::alen,astart,my2mem
   integer,allocatable,dimension(:)::waltyp
   integer,allocatable,intent(in),dimension(:,:)::wposid,myo2body

!   double precision,value::wallbead
!   double precision dx,dy,dz,dist,invdist

   real rad2,x0,y0,z0,x1,y1,z1

   real,allocatable,dimension(:)::xw,yw,zw

   real,allocatable,dimension(:),intent(in)::xwalsurf,ywalsurf,zwalsurf
   real,allocatable,dimension(:),intent(in)::xwalrep,ywalrep,zwalrep
   real,allocatable,dimension(:)::xwalrad,ywalrad,zwalrad
   double precision,allocatable,dimension(:),intent(in)::xfa,yfa,zfa

   double precision,allocatable,intent(in),dimension(:)::xmyp2,ymyp2,zmyp2,xmyo2,ymyo2,zmyo2
   double precision,allocatable,intent(in),dimension(:)::xlk,ylk,zlk,xmemb,ymemb,zmemb
   double precision,allocatable,intent(in),dimension(:,:)::xwall,ywall,zwall

   real,allocatable,dimension(:,:)::xtem,ytem,ztem

   nframe=nframe+1

   allocate(xw(natom),yw(natom),zw(natom))

!  surface wall beads:


   xw(1:nwalsurf)=xwalsurf(1:nwalsurf)
   yw(1:nwalsurf)=ywalsurf(1:nwalsurf)
   zw(1:nwalsurf)=zwalsurf(1:nwalsurf)

   jw=nwalsurf


   allocate(xtem(nphisol,nxsol),ytem(nphisol,nxsol),ztem(nphisol,nxsol))



   do jx=1,nxsol

      x0=xwall(nphisol,jx)
      y0=ywall(nphisol,jx)
      z0=zwall(nphisol,jx)

      do jp=1,nphisol

         jw=jw+1

         x1=xwall(jp,jx)
         y1=ywall(jp,jx)
         z1=zwall(jp,jx)

         xw(jw)=0.1*x1
         yw(jw)=0.1*y1
         zw(jw)=0.1*z1

         xtem(jp,jx)=0.05*(x0+x1)
         ytem(jp,jx)=0.05*(y0+y1)
         ztem(jp,jx)=0.05*(z0+z1)

         x0=x1
         y0=y1
         z0=z1

         jw=jw+1

         xw(jw)=xtem(jp,jx)
         yw(jw)=ytem(jp,jx)
         zw(jw)=ztem(jp,jx)

      end do

   end do

   do jx=1,nxsol-1

      do jp=1,nphisol

         jw=jw+1

         xw(jw)=0.05*(xwall(jp,jx)+xwall(jp,jx+1))
         yw(jw)=0.05*(ywall(jp,jx)+ywall(jp,jx+1))
         zw(jw)=0.05*(zwall(jp,jx)+zwall(jp,jx+1))

         jw=jw+1

         xw(jw)=0.5*(xtem(jp,jx)+xtem(jp,jx+1))
         yw(jw)=0.5*(ytem(jp,jx)+ytem(jp,jx+1))
         zw(jw)=0.5*(ztem(jp,jx)+ztem(jp,jx+1))

      end do

   end do

   deallocate(xtem,ytem,ztem)

!  radial wall beads:

   do n=1,nwalrad

      if(waltyp(n)==1)then
         cycle
      end if

      jp=wposid(1,n)

      jx=wposid(2,n)

      rad2=ywalrep(n)*ywalrep(n)+zwalrep(n)*zwalrep(n)

      if(rad2>ywall(jp,jx)*ywall(jp,jx)+zwall(jp,jx)*zwall(jp,jx))then
         waltyp(n)=1
         xwalrad(n)=xwalrep(n)
         ywalrad(n)=ywalrep(n)
         zwalrad(n)=zwalrep(n)
      end if

   end do

   xw(jw+1:jw+nwalrad)=0.1*xwalrad(1:nwalrad)
   yw(jw+1:jw+nwalrad)=0.1*ywalrad(1:nwalrad)
   zw(jw+1:jw+nwalrad)=0.1*zwalrad(1:nwalrad)

   jw=jw+nwalrad

!  membrane:

   xw(jw+1:jw+nmemb)=0.1*xmemb(1:nmemb)
   yw(jw+1:jw+nmemb)=0.1*ymemb(1:nmemb)
   zw(jw+1:jw+nmemb)=0.1*zmemb(1:nmemb)

   jw=jw+nmemb



!  visualize F-actin:


   do n=1,fanum1

      jstart=astart(n)

      do j=1,alen(n)-1

         jw=jw+1

         xw(jw)=0.1*xfa(jstart+j-1)
         yw(jw)=0.1*yfa(jstart+j-1)
         zw(jw)=0.1*zfa(jstart+j-1)

         jw=jw+1

         xw(jw)=0.1*xfa(jstart+j)
         yw(jw)=0.1*yfa(jstart+j)
         zw(jw)=0.1*zfa(jstart+j)


      end do

   end do

   if(jw<natom_a1)then

!print*,'a1',jw,natom_a1

      xw(jw+1:natom_a1)=0.0
      yw(jw+1:natom_a1)=0.0
      zw(jw+1:natom_a1)=0.0

      jw=natom_a1


   end if


   do n=1+fanum1,fanum1+fanum2

      jstart=astart(n)

      do j=1,alen(n)-1
   
         jw=jw+1

         xw(jw)=0.1*xfa(jstart+j-1)
         yw(jw)=0.1*yfa(jstart+j-1)
         zw(jw)=0.1*zfa(jstart+j-1)

         jw=jw+1

         xw(jw)=0.1*xfa(jstart+j)
         yw(jw)=0.1*yfa(jstart+j)
         zw(jw)=0.1*zfa(jstart+j)



      end do

   end do

   if(jw<natom_a2)then

!print*,'a2',jw,natom_a2

      xw(jw+1:natom_a2)=0.0
      yw(jw+1:natom_a2)=0.0
      zw(jw+1:natom_a2)=0.0

      jw=natom_a2

   end if

!print*,'natom_a2',natom_a2,jw



!   do n=1,nfa

!      if(a2mem(n)>0)then

!         jw=jw+1
!         xw(jw)=0.1*xfa(n)
!         yw(jw)=0.1*yfa(n)
!         zw(jw)=0.1*zfa(n)

!         jm=a2mem(n)

!         jw=jw+1
!         xw(jw)=0.1*xmemb(jm)
!         yw(jw)=0.1*ymemb(jm)
!         zw(jw)=0.1*zmemb(jm)


!      end if

!   end do

!   if(jw<natom_a2m)then

!      xw(jw+1:natom_a2m)=0.0
!      yw(jw+1:natom_a2m)=0.0
!      zw(jw+1:natom_a2m)=0.0

!      jw=natom_a2m

!   end if

!print*,'a2m',jw,natom_a2m

!   jw=jw+2*ntether

!  Myp2 myosin:

   xw(jw+1:jw+nmyp2)=0.1*xmyp2(1:nmyp2)
   yw(jw+1:jw+nmyp2)=0.1*ymyp2(1:nmyp2)
   zw(jw+1:jw+nmyp2)=0.1*zmyp2(1:nmyp2)

   jw=jw+nmyp2


!  visualize Myo2 myosin:


   xw(jw+1:jw+nmyo2)=0.1*xmyo2(1:nmyo2)
   yw(jw+1:jw+nmyo2)=0.1*ymyo2(1:nmyo2)
   zw(jw+1:jw+nmyo2)=0.1*zmyo2(1:nmyo2)

   jw=jw+nmyo2

!  tethering Myo2 myosin to membrane:

   do n=1,myo2num

      jmemb=my2mem(n)

      jm=myo2body(1,n)

      jw=jw+1

      xw(jw)=0.1*xmemb(jmemb)
      yw(jw)=0.1*ymemb(jmemb)
      zw(jw)=0.1*zmemb(jmemb)

      jw=jw+1

      xw(jw)=0.1*xmyo2(jm)
      yw(jw)=0.1*ymyo2(jm)
      zw(jw)=0.1*zmyo2(jm)

   end do

!  crosslinkers:

   xw(jw+1:jw+nlk)=0.1*xlk(1:nlk)
   yw(jw+1:jw+nlk)=0.1*ylk(1:nlk)
   zw(jw+1:jw+nlk)=0.1*zlk(1:nlk)

   jw=jw+nlk

   if(jw<natom)then

      xw(jw+1:natom)=0.0
      yw(jw+1:natom)=0.0
      zw(jw+1:natom)=0.0

   end if

   write(junit)xw(1:natom)
   write(junit)yw(1:natom)
   write(junit)zw(1:natom)

   deallocate(xw,yw,zw)

   end subroutine

!=========================================================

   subroutine getinfo(nstart,jfile,time0,runtime,natom,natom_a1,natom_a2,nwall,nxsol,nphisol, &
               nwalsurf,nwalrad,nmemb,nact1,nact2,nfamax1,nfamax2,nfa1,nfa2,fanum1,fanum2, &
                nmyp2,myp2num,nmyo2,myo2num,nlk,crlknum1,crlknum2,crlknummax,crlknumactive,nmyoturn)


   implicit none

   integer(kind=8)::time0

   integer nstart,jfile
   double precision runtime

   integer natom,natom_a1,natom_a2,nwall,nxsol,nphisol,nwalsurf,nwalrad
   integer nmemb,nact1,nact2,nfamax1,nfamax2,nfa1,nfa2,fanum1,fanum2
   integer nmyp2,myp2num,nmyo2,myo2num,nlk,crlknum1,crlknum2,crlknummax,crlknumactive,nmyoturn

   character zero*1,charid1*1,charid2*2,charid3*3
   character (len=64) fileconf,chara

   open(1,file='restart.inp')
   read(1,*)nstart,jfile
   read(1,*)time0,runtime
   close(1)

   write(zero,'(i1)')0

   if(nstart<10)then
      write(charid1,'(i1)')nstart
      fileconf='rconf'//zero//zero//charid1//'.inp'
   elseif(nstart<100)then
      write(charid2,'(i2)')nstart
      fileconf='rconf'//zero//charid2//'.inp'
   else
      write(charid3,'(i3)')nstart
      fileconf='rconf'//charid3//'.inp'
   end if

   open(1,file=fileconf)

11 read(1,*)chara

   if(chara(1:4)/='VISU')then
      goto 11
   end if

   read(1,*)natom,natom_a1,natom_a2

2  read(1,*)chara

   if(chara(1:4)/='CELL')then
      goto 2
   end if

   read(1,*)nwall,nxsol,nphisol,nwalsurf,nwalrad

12 read(1,*)chara

   if(chara(1:4)/='MEMB')then
      goto 12
   end if

   read(1,*)nmemb

3  read(1,*)chara

   if(chara(1:4)/='FACT')then
      goto 3
   end if

   read(1,*)nact1,nact2,nfamax1,nfamax2,nfa1,nfa2,fanum1,fanum2!,ntether

!   nfa=nfa1+nfa2

!   fanum=fanum1+fanum2


5  read(1,*)chara

   if(chara(1:4)/='MYP2')then
      goto 5
   end if

   read(1,*)nmyp2,myp2num

15 read(1,*)chara

   if(chara(1:4)/='MYO2')then
      goto 15
   end if

   read(1,*)nmyo2,myo2num,nmyoturn

8  read(1,*)chara

   if(chara(1:4)/='CROS')then
      goto 8
   end if

   read(1,*)nlk,crlknum1,crlknum2,crlknummax,crlknumactive

!   crlknum=crlknum1+crlknum2

   close(1)

   end subroutine

!=========================================================

   subroutine ringin(nstart,nxsol,nphisol,nmemb,fanum,nfa,myp2num,nmyp2,myo2num,nmyo2,nwalrad,nwalsurf,nwall, &
               myp2len,myo2len,crlknum,nlk,astart,alen,filid,apos,fa1stbound,my2mem,waltyp,lkstart,lktyp, &
                 myp2body,myp2head,myp2typ,myo2body,myo2head,myo2typ,apar,wposid,fa2myp2,fa2myo2,fa2lk, &
                  xboundmin,xboundmax,dxsol,dphisol,xwall,ywall,zwall,xnorwall,ynorwall,znorwall,xwalsurf, &
                   ywalsurf,zwalsurf,xwalrad,ywalrad,zwalrad,xwalrep,ywalrep,zwalrep,xfa,yfa,zfa,xmyp2,ymyp2,zmyp2, &
                    xmyo2,ymyo2,zmyo2,xlk,ylk,zlk,xmemb,ymemb,zmemb,xsurf,ysurf,zsurf,xnorsurf,ynorsurf,znorsurf)



   implicit none

   integer,value::nstart,nxsol,nphisol,nmemb,fanum,nfa,myp2num,nmyp2,myo2num,nmyo2,nwalrad,nwalsurf,nwall
   integer,value::myp2len,myo2len,crlknum,nlk

   integer n,jread(1:20),nmyoturn

   integer,allocatable,dimension(:)::astart,alen,filid,apos,fa1stbound,my2mem,waltyp,lkstart,lktyp
   integer,allocatable,dimension(:,:)::myp2body,myp2head,myp2typ,myo2body,myo2head,myo2typ,apar,wposid
   integer,allocatable,dimension(:,:)::fa2myp2,fa2myo2,fa2lk


   double precision dxsol,dphisol,xboundmin,xboundmax

   double precision,allocatable,dimension(:,:)::xwall,ywall,zwall,xnorwall,ynorwall,znorwall
   real,allocatable,dimension(:)::xwalsurf,ywalsurf,zwalsurf,xwalrad,ywalrad,zwalrad
   real,allocatable,dimension(:)::xwalrep,ywalrep,zwalrep

   double precision,allocatable,dimension(:)::xfa,yfa,zfa,xmyp2,ymyp2,zmyp2,xmyo2,ymyo2,zmyo2
   double precision,allocatable,dimension(:)::xlk,ylk,zlk,xmemb,ymemb,zmemb
   double precision,allocatable,dimension(:,:)::xsurf,ysurf,zsurf,xnorsurf,ynorsurf,znorsurf

   character zero*1,charid1*1,charid2*2,charid3*3
   character (len=64) filecoor,fileconf,chara

   write(zero,'(i1)')0

   if(nstart<10)then
      write(charid1,'(i1)')nstart
      fileconf='rconf'//zero//zero//charid1//'.inp'
      filecoor='rcoor'//zero//zero//charid1//'.inp'
   elseif(nstart<100)then
      write(charid2,'(i2)')nstart
      fileconf='rconf'//zero//charid2//'.inp'
      filecoor='rcoor'//zero//charid2//'.inp'
   else
      write(charid3,'(i3)')nstart
      fileconf='rconf'//charid3//'.inp'
      filecoor='rcoor'//charid3//'.inp'
   end if

!  read coordinates

   open(1,file=filecoor,form='unformatted')

   read(1)dxsol,dphisol

   read(1)xwall(1:nphisol,1:nxsol)
   read(1)ywall(1:nphisol,1:nxsol)
   read(1)zwall(1:nphisol,1:nxsol)
   read(1)xnorwall(1:nphisol,1:nxsol)
   read(1)ynorwall(1:nphisol,1:nxsol)
   read(1)znorwall(1:nphisol,1:nxsol)

   read(1)xwalsurf(1:nwalsurf)
   read(1)ywalsurf(1:nwalsurf)
   read(1)zwalsurf(1:nwalsurf)

   read(1)xwalrad(1:nwalrad)
   read(1)ywalrad(1:nwalrad)
   read(1)zwalrad(1:nwalrad)

   read(1)xwalrep(1:nwalrad)
   read(1)ywalrep(1:nwalrad)
   read(1)zwalrep(1:nwalrad)


!---------------------------      

   read(1)jread(1)!nmemb
   read(1)xmemb(1:nmemb)
   read(1)ymemb(1:nmemb)
   read(1)zmemb(1:nmemb)
   read(1)xboundmin,xboundmax
   read(1)xsurf(1:nphisol,1:nxsol)
   read(1)ysurf(1:nphisol,1:nxsol)
   read(1)zsurf(1:nphisol,1:nxsol)
   read(1)xnorsurf(1:nphisol,1:nxsol)
   read(1)ynorsurf(1:nphisol,1:nxsol)
   read(1)znorsurf(1:nphisol,1:nxsol)


   read(1)jread(1)!nfa
   read(1)xfa(1:nfa)
   read(1)yfa(1:nfa)
   read(1)zfa(1:nfa)

   read(1)jread(1)!nmyp2
   read(1)xmyp2(1:nmyp2)
   read(1)ymyp2(1:nmyp2)
   read(1)zmyp2(1:nmyp2)

   read(1)jread(1)!nmyo2
   read(1)xmyo2(1:nmyo2)
   read(1)ymyo2(1:nmyo2)
   read(1)zmyo2(1:nmyo2)

   read(1)jread(1)!nlk
   read(1)xlk(1:nlk)
   read(1)ylk(1:nlk)
   read(1)zlk(1:nlk)

   close(1)

!  read configuration

   open(1,file=fileconf)

9  read(1,*)chara

   if(chara(1:4)/='CELL')then
      goto 9
   end if

   read(1,*)jread(1:5)!nwall,nxsol,nphisol,nwalsurf,nwalrad

   read(1,*)waltyp(1:nwalrad)

   read(1,*)wposid(1,1:nwalrad)
   read(1,*)wposid(2,1:nwalrad)



3  read(1,*)chara

   if(chara(1:4)/='FACT')then
      goto 3
   end if

   read(1,*)jread(1:8)!nact1,nact2,nfamax1,nfamax2,nfa1,nfa2,fanum1,fanum2

   read(1,*)filid(1:nfa)
!   read(1,*)a2mem(1:nfa)

   read(1,*)apar(1:2,1:nfa)
   read(1,*)apos(1:nfa)
   read(1,*)fa1stbound(1:fanum)
   read(1,*)alen(1:fanum)
   read(1,*)astart(1:fanum)


!print*,'myo'

5  read(1,*)chara

   if(chara(1:4)/='MYP2')then
      goto 5
   end if

   read(1,*)jread(1:2)!!nmyp2,myp2num
   do n=1,myp2num
      read(1,*)myp2head(1:4,n),myp2body(1:myp2len,n)
      read(1,*)myp2typ(1:4,n),fa2myp2(1:4,n)
   end do


15 read(1,*)chara

   if(chara(1:4)/='MYO2')then
      goto 15
   end if

   read(1,*)jread(1:3)!nmyo2,myo2num
   do n=1,myo2num
      read(1,*)my2mem(n),myo2head(1:2,n),myo2body(1:myo2len,n)
      read(1,*)myo2typ(1:2,n),fa2myo2(1:2,n)
   end do

!print*,'crlk'

8  read(1,*)chara

   if(chara(1:4)/='CROS')then
      goto 8
   end if

   read(1,*)jread(1:5)!nlk,crlknum1,crlknum2,crlknummax,crlknumactive


   read(1,*)lkstart(1:crlknum)
   read(1,*)lktyp(1:crlknum)
   read(1,*)fa2lk(1,1:crlknum)
   read(1,*)fa2lk(2,1:crlknum)



   close(1)

   end subroutine

!=========================================================

   subroutine measures(junit2,nfa,crlknumactive,natp,jsignal,nmemb,nfamax1,nfamax2,crlknummax, &
              wallrate,aconc1,aconc2,lkconc,thet2by2,cos_t_2,printtime,runtime, &
              pi,l_mem,xfa,yfa,zfa,ymemb,zmemb)

   integer,value::junit2,nfa,crlknumactive,natp,jsignal,nmemb
   integer nfamax1,nfamax2,crlknummax

   integer ncount(4),n,nxlinker,j
!   integer,allocatable,dimension(:,:),intent(in)::jfasol
!   integer,allocatable,dimension(:),intent(in)::lktyp

   real(kind=8)::wallrate,ringthick,aconc1,aconc2,lkconc,thet2by2,cos_t_2
   real(kind=8),value::printtime,runtime,pi,l_mem
   real(kind=8)::e_tk(4),e_tk2(4),e_wid(4),e_wid2(4),ringrad,rad,rad2,ringwid,a_conc,ringvol,thet,y0,z0
   real(kind=8),allocatable,dimension(:),intent(in)::xfa,yfa,zfa,ymemb,zmemb
!   real(kind=8),allocatable,dimension(:,:),intent(in)::xsurf,ysurf,zsurf,xnorsurf,ynorsurf,znorsurf

!  wall remodeling rate in second:
   wallrate=wallrate/printtime*1000000

   ncount(1:4)=0

   e_tk(1:4)=0.0d0

   e_tk2(1:4)=0.0d0

   e_wid(1:4)=0.0d0

   e_wid2(1:4)=0.0d0

   ringrad=0.0d0


!  find the center:

!   y0=sum(yfa(1:nfa))/nfa
!   z0=sum(zfa(1:nfa))/nfa


   do n=1,nfa,4

!      jx=jfasol(1,n)

!      jp=jfasol(2,n)

!      dx=xfa(n)-xsurf(jp,jx)
!      dy=yfa(n)-ysurf(jp,jx)
!      dz=zfa(n)-zsurf(jp,jx)


!      thick=dx*xnorsurf(jp,jx)+dy*ynorsurf(jp,jx)+dz*znorsurf(jp,jx)

!      rad=sqrt(yfa(n)*yfa(n)+zfa(n)*zfa(n))

!      thick=rmemb(jp,jx)-rad

      if(yfa(n)>0.0d0.and.zfa(n)>0.0d0)then

         j=1

      elseif(yfa(n)<0.0d0.and.zfa(n)>0.0d0)then

         j=2

      elseif(yfa(n)<0.0d0.and.zfa(n)<0.0d0)then

         j=3

      else

         j=4

      end if

!      rad2=(yfa(n)-y0)*(yfa(n)-y0)+(zfa(n)-z0)*(zfa(n)-z0)

      rad2=yfa(n)*yfa(n)+zfa(n)*zfa(n)

      rad=sqrt(rad2)

      ncount(j)=ncount(j)+1

      e_tk(j)=e_tk(j)+rad

      e_tk2(j)=e_tk2(j)+rad2

      e_wid(j)=e_wid(j)+xfa(n)

      e_wid2(j)=e_wid2(j)+xfa(n)*xfa(n)


!      ncount=ncount+1

!      e_tk=e_tk+thick

!      e_tk2=e_tk2+thick*thick

!      e_wid=e_wid+xfa(n)

!      e_wid2=e_wid2+xfa(n)*xfa(n)

      ringrad=ringrad+rad

   end do

!  use formular std^2 = E(x^2) - (E(x))^2

   ringthick=sum(sqrt(e_tk2(1:4)*ncount(1:4)-e_tk(1:4)*e_tk(1:4))/ncount(1:4))

   ringwid=sum(sqrt(e_wid2(1:4)*ncount(1:4)-e_wid(1:4)*e_wid(1:4))/ncount(1:4))

   ringrad=ringrad/sum(ncount(1:4))

   ringvol=2*pi*ringrad*ringthick*ringwid

!  actin conc in uM:
   a_conc=(nfa*2/6.02/ringvol)*10000000


!   nxlinker=0

!   do n=1,crlknum

!      if(lktyp(n)==1)then
!         nxlinker=nxlinker+1
!      end if

!   end do

   write(junit2,2)runtime*0.000001,wallrate,natp,ringthick,ringwid,ringrad,2*nfa,a_conc,crlknumactive

   wallrate=0.0d0

2  format(f10.1,2x,f10.6,2x,i8,5x,f6.1,4x,f6.1,3x,f6.1,4x,i6,2x,f7.1,5x,i5)


   if(jsignal==1)then

      aconc1=nfamax1/ringrad

      aconc2=nfamax2/ringrad

      lkconc=crlknummax/ringrad

   elseif(jsignal==2)then

      n=aconc1*ringrad
      nfamax1=min(nfamax1,n)

      n=aconc2*ringrad
      nfamax2=min(nfamax2,n)

      n=lkconc*ringrad
      crlknummax=min(crlknummax,n)

   end if

!  update membrane radius:

   rad=0.0d0

!$omp parallel &
!$omp default(none) &
!$omp private(n) &
!$omp shared(nmemb,ymemb,zmemb) &
!$omp reduction(+:rad)
!$omp do schedule (guided,64)

   do n=1,nmemb

      rad=rad+sqrt(ymemb(n)*ymemb(n)+zmemb(n)*zmemb(n))

   end do

!$omp end do nowait
!$omp end parallel

   rad=rad/nmemb

   thet=2.5d0*l_mem/rad

   cos_t_2=(1.0d0-(2*thet)**2/2)**2

   thet2by2=thet*thet/2

   end subroutine

!=========================================================

   subroutine ringout(nstart,natom,natom_a1,natom_a2,nwall,nxsol,nphisol,nwalsurf,nwalrad, &
               nmemb,nact1,nact2,nfamax1,nfamax2,nfa1,nfa2,fanum1,fanum2,nfa,nmyp2,myp2num,nmyoturn, &
                myp2len,nmyo2,myo2num,myo2len,nlk,crlknum1,crlknum2,crlknummax,crlknumactive, &
                 astart,alen,filid,apos,fa1stbound,my2mem,waltyp,lkstart,lktyp,myp2body,myp2head, &
                  myp2typ,fa2myp2,myo2body,myo2head,myo2typ,fa2myo2,fa2lk,apar,wposid, &
                   dxsol,dphisol,xboundmin,xboundmax,xsurf,ysurf,zsurf,xnorsurf,ynorsurf,znorsurf, &
                    xfa,yfa,zfa,xmyp2,ymyp2,zmyp2,xmyo2,ymyo2,zmyo2,xlk,ylk,zlk,xmemb,ymemb,zmemb, &
                     xwall,ywall,zwall,xnorwall,ynorwall,znorwall,xwalsurf,ywalsurf,zwalsurf, &
                      xwalrad,ywalrad,zwalrad,xwalrep,ywalrep,zwalrep)



   implicit none

   integer nstart

   integer,value::natom,natom_a1,natom_a2
   integer,value::nwall,nxsol,nphisol,nwalsurf,nwalrad
   integer,value::nmemb,nact1,nact2,nfamax1,nfamax2,nfa1,nfa2,fanum1,fanum2!,ntether
   integer,value::nfa,nmyp2,myp2num,myp2len,nmyo2,myo2num,myo2len,nmyoturn
   integer,value::nlk,crlknum1,crlknum2,crlknummax,crlknumactive

   integer lentemp,n,fanum,crlknum

   integer,allocatable,intent(in),dimension(:)::astart,alen,filid,apos,fa1stbound,my2mem,waltyp,lkstart,lktyp
   integer,allocatable,intent(in),dimension(:,:)::myp2body,myp2head,myp2typ,fa2myp2
   integer,allocatable,intent(in),dimension(:,:)::myo2body,myo2head,myo2typ,fa2myo2,fa2lk,apar,wposid

   double precision,value::dxsol,dphisol,xboundmin,xboundmax
   double precision,allocatable,intent(in),dimension(:)::xfa,yfa,zfa,xmyp2,ymyp2,zmyp2,xmyo2,ymyo2,zmyo2
   double precision,allocatable,intent(in),dimension(:)::xlk,ylk,zlk,xmemb,ymemb,zmemb
   double precision,allocatable,intent(in),dimension(:,:)::xsurf,ysurf,zsurf,xnorsurf,ynorsurf,znorsurf
   double precision,allocatable,intent(in),dimension(:,:)::xwall,ywall,zwall,xnorwall,ynorwall,znorwall
   real,allocatable,intent(in),dimension(:)::xwalsurf,ywalsurf,zwalsurf
   real,allocatable,dimension(:),intent(in)::xwalrad,ywalrad,zwalrad,xwalrep,ywalrep,zwalrep


   character zero*1,charid1*1,charid2*2,charid3*3
   character (len=64) filecoor,fileconf

   nstart=nstart+1

   write(zero,'(i1)')0

   if(nstart<10)then
      write(charid1,'(i1)')nstart
      fileconf='rconf'//zero//zero//charid1//'.inp'
      filecoor='rcoor'//zero//zero//charid1//'.inp'
   elseif(nstart<100)then
      write(charid2,'(i2)')nstart
      fileconf='rconf'//zero//charid2//'.inp'
      filecoor='rcoor'//zero//charid2//'.inp'
   else
      write(charid3,'(i3)')nstart
      fileconf='rconf'//charid3//'.inp'
      filecoor='rcoor'//charid3//'.inp'
   end if

!  coordinates

   open(1,file=filecoor,form='unformatted')

   write(1)dxsol,dphisol

   write(1)xwall(1:nphisol,1:nxsol)
   write(1)ywall(1:nphisol,1:nxsol)
   write(1)zwall(1:nphisol,1:nxsol)
   write(1)xnorwall(1:nphisol,1:nxsol)
   write(1)ynorwall(1:nphisol,1:nxsol)
   write(1)znorwall(1:nphisol,1:nxsol)
!   write(1)rwall(1:nphisol,1:nxsol)

   write(1)xwalsurf(1:nwalsurf)
   write(1)ywalsurf(1:nwalsurf)
   write(1)zwalsurf(1:nwalsurf)

   write(1)xwalrad(1:nwalrad)
   write(1)ywalrad(1:nwalrad)
   write(1)zwalrad(1:nwalrad)

   write(1)xwalrep(1:nwalrad)
   write(1)ywalrep(1:nwalrad)
   write(1)zwalrep(1:nwalrad)

   write(1)nmemb
   write(1)xmemb(1:nmemb)
   write(1)ymemb(1:nmemb)
   write(1)zmemb(1:nmemb)
   write(1)xboundmin,xboundmax
   write(1)xsurf(1:nphisol,1:nxsol)
   write(1)ysurf(1:nphisol,1:nxsol)
   write(1)zsurf(1:nphisol,1:nxsol)
   write(1)xnorsurf(1:nphisol,1:nxsol)
   write(1)ynorsurf(1:nphisol,1:nxsol)
   write(1)znorsurf(1:nphisol,1:nxsol)


   write(1)nfa
   write(1)xfa(1:nfa)
   write(1)yfa(1:nfa)
   write(1)zfa(1:nfa)

   write(1)nmyp2
   write(1)xmyp2(1:nmyp2)
   write(1)ymyp2(1:nmyp2)
   write(1)zmyp2(1:nmyp2)

   write(1)nmyo2
   write(1)xmyo2(1:nmyo2)
   write(1)ymyo2(1:nmyo2)
   write(1)zmyo2(1:nmyo2)

   write(1)nlk
   write(1)xlk(1:nlk)
   write(1)ylk(1:nlk)
   write(1)zlk(1:nlk)

   close(1)

!  write configuration

   open(1,file=fileconf)

   write(1,*)'VISUAL'

   write(1,*)natom,natom_a1,natom_a2!,natom_a2m

   write(1,*)'CELLWALL'

   write(1,*)nwall,nxsol,nphisol,nwalsurf,nwalrad

   write(1,*)waltyp(1:nwalrad)

   write(1,*)wposid(1,1:nwalrad)
   write(1,*)wposid(2,1:nwalrad)

   write(1,*)


   write(1,*)'MEMBRANE'
   write(1,*)nmemb
!   write(1,*)xboundmin,xboundmax
   write(1,*)



   fanum=fanum1+fanum2

   write(1,*)'FACTIN'
   write(1,*)nact1,nact2,nfamax1,nfamax2,nfa1,nfa2,fanum1,fanum2!,ntether
   write(1,*)filid(1:nfa)
!   write(1,*)a2mem(1:nfa)

   write(1,*)apar(1:2,1:nfa)

   write(1,*)apos(1:nfa)

   write(1,*)fa1stbound(1:fanum)

   write(1,*)alen(1:fanum)

   write(1,*)astart(1:fanum)


   write(1,*)'MYP2 Myosin'
   write(1,*)nmyp2,myp2num
   do n=1,myp2num
      write(1,*)myp2head(1:4,n),myp2body(1:myp2len,n)
      write(1,*)myp2typ(1:4,n),fa2myp2(1:4,n)
   end do
   write(1,*)


   write(1,*)'MYO2 myosin'
   write(1,*)nmyo2,myo2num,nmyoturn
   do n=1,myo2num
      write(1,*)my2mem(n),myo2head(1:2,n),myo2body(1:myo2len,n)
      write(1,*)myo2typ(1:2,n),fa2myo2(1:2,n)
   end do
   write(1,*)

   write(1,*)'CROSSLINKER'
   write(1,*)nlk,crlknum1,crlknum2,crlknummax,crlknumactive

   crlknum=crlknum1+crlknum2

   write(1,*)lkstart(1:crlknum)
   write(1,*)lktyp(1:crlknum)
   write(1,*)fa2lk(1,1:crlknum)
   write(1,*)fa2lk(2,1:crlknum)



   close(1)

   end subroutine




!=========================================================

          SUBROUTINE init_random_seed()
            INTEGER :: i, n, clock
            INTEGER, DIMENSION(:), ALLOCATABLE :: seed

            CALL RANDOM_SEED(size = n)
            ALLOCATE(seed(n))

            CALL SYSTEM_CLOCK(COUNT=clock)

            seed = clock + 37 * (/ (i - 1, i = 1, n) /)
            CALL RANDOM_SEED(PUT = seed)

            DEALLOCATE(seed)
          END SUBROUTINE


!=========================================================

function r8_uniform_01 ( seed )

!*****************************************************************************80
!
!! R8_UNIFORM_01 returns a unit pseudorandom R8.
!
!  Discussion:
!
!    An R8 is a real ( kind = 8 ) value.
!
!    For now, the input quantity SEED is an integer variable.
!
!    This routine implements the recursion
!
!      seed = ( 16807 * seed ) mod ( 2^31 - 1 )
!      r8_uniform_01 = seed / ( 2^31 - 1 )
!
!    The integer arithmetic never requires more than 32 bits,
!    including a sign bit.
!
!    If the initial seed is 12345, then the first three computations are
!
!      Input     Output      R8_UNIFORM_01
!      SEED      SEED
!
!         12345   207482415  0.096616
!     207482415  1790989824  0.833995
!    1790989824  2035175616  0.947702
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Second Edition,
!    Springer, 1987,
!    ISBN: 0387964673,
!    LC: QA76.9.C65.B73.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, December 1986, pages 362-376.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley, 1998,
!    ISBN: 0471134031,
!    LC: T57.62.H37.
!
!    Peter Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, Number 2, 1969, pages 136-143.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which should
!    NOT be 0. On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R8_UNIFORM_01, a new pseudorandom variate,
!    strictly between 0 and 1.
!
  implicit none

  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) k
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop 1
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + i4_huge
  end if

  r8_uniform_01 = real ( seed, kind = 8 ) * 4.656612875D-10

  return
end function


!=========================================================

subroutine r4vec_uniform_01 ( n, seed, r )

!*****************************************************************************80
!
!! R4VEC_UNIFORM_01 returns a unit pseudorandom R4VEC.
!
!  Discussion:
!
!    An R4VEC is an array of real ( kind = 4 ) values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Second Edition,
!    Springer, 1987,
!    ISBN: 0387964673,
!    LC: QA76.9.C65.B73.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, December 1986, pages 362-376.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley, 1998,
!    ISBN: 0471134031,
!    LC: T57.62.H37.
!
!    Peter Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, Number 2, 1969, pages 136-143.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value,
!    which should NOT be 0.
!    On output, SEED has been updated.
!
!    Output, real ( kind = 4 ) R(N), the vector of pseudorandom values.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed
  real ( kind = 8 ) r(n)

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4VEC_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop 1
  end if

  do i = 1, n

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + i4_huge
    end if

    r(i) = real ( seed, kind = 8 ) * 4.656612875E-10

  end do

  return
end subroutine

!=========================================================

   subroutine neighbors(myp2num,myo2num,neinum,crlknum,lklen,neinum5,myp2len,myo2len,filid, &
              astart,alen,lkstart,myp2head,myp2body,myo2head,myo2body,myp2nei,myo2nei,lknei, &
              myp2nei5,myo2nei5,lknei5,xfa,yfa,zfa,xmyp2,ymyp2,zmyp2,xmyo2,ymyo2,zmyo2,xlk,ylk,zlk)


   implicit none

   integer,value:: myp2num,myo2num,neinum,crlknum,lklen,neinum5,myp2len,myo2len
   integer nm,jmh,jmyo,jnei,ja,nl,jlk,jn,j1,j2,nf

   integer,allocatable,intent(in),dimension(:)::filid,astart,alen,lkstart
   integer,allocatable,intent(in),dimension(:,:)::myp2head,myp2body,myo2head,myo2body
   integer,allocatable,dimension(:,:,:)::myp2nei,myo2nei,lknei
   integer,allocatable,dimension(:,:,:),intent(in)::myp2nei5,myo2nei5,lknei5
   double precision dx,dy,dz,d2,d2max,dx1,dy1,dz1,dx2,dy2,dz2

   double precision,allocatable,intent(in),dimension(:)::xfa,yfa,zfa,xlk,ylk,zlk
   double precision,allocatable,intent(in),dimension(:)::xmyp2,ymyp2,zmyp2,xmyo2,ymyo2,zmyo2

   d2max=15.0d0*15.0d0
   myp2nei=0
   myo2nei=0
   lknei=0

!$omp parallel &
!$omp default(none) &
!$omp private(nm,j1,j2,nf,dx1,dy1,dz1,dx2,dy2,dz2,jmh,jmyo,jnei,jn,ja,dx,dy,dz,d2) &
!$omp shared(myp2num,myp2head,myp2body,myp2len,filid,astart,alen,xfa,yfa,zfa,xmyp2,ymyp2,zmyp2) &
!$omp shared(myp2nei,neinum,myp2nei5,neinum5,d2max) &
!$omp shared(myo2num,myo2head,myo2body,myo2len,myo2nei,myo2nei5,xmyo2,ymyo2,zmyo2) &
!$omp private(nl,jlk) &
!$omp shared(crlknum,lkstart,lklen,xlk,ylk,zlk,lknei,lknei5)


!  for Myp2 myosin:

!$omp do schedule(guided,32)

   do nm=1,myp2num

      j1=myp2body(1,nm)

      j2=myp2body(myp2len,nm)

      dx1=xmyp2(j1)-xmyp2(j2)
      dy1=ymyp2(j1)-ymyp2(j2)
      dz1=zmyp2(j1)-zmyp2(j2)

      do jmh=1,2

         jmyo=myp2head(jmh,nm)

         jnei=0

         do jn=1,neinum5

            if(myp2nei5(jn,jmh,nm)==0)then
               exit
            end if

            ja=myp2nei5(jn,jmh,nm)

            nf=filid(ja)

            if(ja==astart(nf)+alen(nf)-1)then

               dx2=xfa(ja-1)-xfa(ja)
               dy2=yfa(ja-1)-yfa(ja)
               dz2=zfa(ja-1)-zfa(ja)

            else

               dx2=xfa(ja)-xfa(ja+1)
               dy2=yfa(ja)-yfa(ja+1)
               dz2=zfa(ja)-zfa(ja+1)

            end if

            if(dx1*dx2+dy1*dy2+dz1*dz2<0.0d0)then
               cycle
            end if





!            if(adir(ja)/=mydir(nm))then
!               cycle
!            end if

            dx=xfa(ja)-xmyp2(jmyo)
            dy=yfa(ja)-ymyp2(jmyo)
            dz=zfa(ja)-zmyp2(jmyo)

            d2=dx*dx+dy*dy+dz*dz

            if(d2<d2max)then
               jnei=jnei+1
               myp2nei(jnei,jmh,nm)=ja
               if(jnei==neinum)then
                  print*,'too crowded at head',jmh,nm
                  exit
               end if
            end if

         end do

      end do

      do jmh=3,4

         jmyo=myp2head(jmh,nm)

         jnei=0

         do jn=1,neinum5

            if(myp2nei5(jn,jmh,nm)==0)then
               exit
            end if

            ja=myp2nei5(jn,jmh,nm)

            nf=filid(ja)


            if(ja==astart(nf)+alen(nf)-1)then
               dx2=xfa(ja-1)-xfa(ja)
               dy2=yfa(ja-1)-yfa(ja)
               dz2=zfa(ja-1)-zfa(ja)

            else
               dx2=xfa(ja)-xfa(ja+1)
               dy2=yfa(ja)-yfa(ja+1)
               dz2=zfa(ja)-zfa(ja+1)
            end if

            if(dx1*dx2+dy1*dy2+dz1*dz2>0.0d0)then
               cycle
            end if

!            if(adir(ja)==mydir(nm))then
!               cycle
!            end if

            dx=xfa(ja)-xmyp2(jmyo)
            dy=yfa(ja)-ymyp2(jmyo)
            dz=zfa(ja)-zmyp2(jmyo)

            d2=dx*dx+dy*dy+dz*dz

            if(d2<d2max)then
               jnei=jnei+1
               myp2nei(jnei,jmh,nm)=ja
               if(jnei==neinum)then
                  print*,'too crowded at head',jmh,nm
                  exit
               end if

            end if

         end do

      end do

   end do

!$omp end do nowait

!-----------------------------


!  for Myo2 myosin:

!$omp do schedule(guided,32)

   do nm=1,myo2num

      j1=myo2body(1,nm)

      j2=myo2body(myo2len,nm)

      dx1=xmyo2(j1)-xmyo2(j2)
      dy1=ymyo2(j1)-ymyo2(j2)
      dz1=zmyo2(j1)-zmyo2(j2)

      do jmh=1,2

         jmyo=myo2head(jmh,nm)

         jnei=0

         do jn=1,neinum5

            if(myo2nei5(jn,jmh,nm)==0)then
               exit
            end if

            ja=myo2nei5(jn,jmh,nm)

            nf=filid(ja)

            if(ja==astart(nf)+alen(nf)-1)then

               dx2=xfa(ja-1)-xfa(ja)
               dy2=yfa(ja-1)-yfa(ja)
               dz2=zfa(ja-1)-zfa(ja)

            else

               dx2=xfa(ja)-xfa(ja+1)
               dy2=yfa(ja)-yfa(ja+1)
               dz2=zfa(ja)-zfa(ja+1)

            end if

            if(dx1*dx2+dy1*dy2+dz1*dz2>0.0d0)then
               cycle
            end if


            dx=xfa(ja)-xmyo2(jmyo)
            dy=yfa(ja)-ymyo2(jmyo)
            dz=zfa(ja)-zmyo2(jmyo)

            d2=dx*dx+dy*dy+dz*dz

            if(d2<d2max)then
               jnei=jnei+1
               myo2nei(jnei,jmh,nm)=ja
               if(jnei==neinum)then
                  print*,'too crowded at head',jmh,nm
                  exit
               end if
            end if

         end do

      end do

   end do

!$omp end do nowait

!-----------------------------

!  for crosslinkers:

!$omp do schedule(guided,32)

   do nl=1,crlknum

!      jlk=crlk(1,nl)

      jlk=lkstart(nl)

      jnei=0

      do jn=1,neinum5

         if(lknei5(jn,1,nl)==0)then
            exit
         end if

         ja=lknei5(jn,1,nl)

         dx=xfa(ja)-xlk(jlk)
         dy=yfa(ja)-ylk(jlk)
         dz=zfa(ja)-zlk(jlk)

         d2=dx*dx+dy*dy+dz*dz

         if(d2<d2max.and.jnei<neinum)then
            jnei=jnei+1
            lknei(jnei,1,nl)=ja
         end if

      end do

      jlk=jlk+lklen-1  !crlk(lklen,nl)

      jnei=0

      do jn=1,neinum5

         if(lknei5(jn,2,nl)==0)then
            exit
         end if
      
         ja=lknei5(jn,2,nl)

         dx=xfa(ja)-xlk(jlk)
         dy=yfa(ja)-ylk(jlk)
         dz=zfa(ja)-zlk(jlk)

         d2=dx*dx+dy*dy+dz*dz

         if(d2<d2max.and.jnei<neinum)then
            jnei=jnei+1
            lknei(jnei,2,nl)=ja
         end if

      end do

   end do

!$omp end do nowait
!$omp end parallel

   end subroutine

!=========================================================

   subroutine dcdheader(junit,jfile,natom)

   implicit none

   character coor*4,zero*1,charid1*1,charid2*2,charid3*3
   character (len=80) string1,string2,filedcd
   integer ifirst,nframe,nfreq,zeros5(5),peroff,zeros7(7),two,twentyfour,ntot
   integer,value:: junit,jfile,natom
   integer*8 jdelta

   write(zero,'(i1)')0

   if(jfile<10)then
      write(charid1,'(i1)')jfile
      filedcd='traj'//zero//zero//charid1//'.dcd'
   elseif(jfile<100)then
      write(charid2,'(i2)')jfile
      filedcd='traj'//zero//charid2//'.dcd'
   elseif(jfile<1000)then
      write(charid3,'(i3)')jfile
      filedcd='traj'//charid3//'.dcd'
   else
      print*,'too many dcd file already, stopping now...'
      stop
   end if

   open(junit,file=filedcd,form='unformatted')

   COOR='CORD'; NFRAME=10000; IFIRST=0; NFREQ=1; NTOT=100
   ZEROS5=0; JDELTA=1; PEROFF=0; ZEROS7=0; TWENTYFOUR=24; TWO=2
   STRING1='HELLOOOOOOO'; STRING2='WHAT THE HELL!'

   WRITE(JUNIT)COOR,NFRAME,IFIRST,NFREQ,NTOT,ZEROS5,JDELTA,PEROFF,ZEROS7,TWENTYFOUR
   WRITE(JUNIT)TWO,STRING1,STRING2
   WRITE(JUNIT)natom

   end subroutine

!=========================================================

   subroutine myocycle(natp,myp2num,myo2num,neinum,tension2,nrand,jrand,myp2len,myo2len,apar,apos, &
              filid,fa1stbound,myp2typ,fa2myp2,myo2typ,fa2myo2,myp2head,myp2body,myo2head,myo2body,myp2nei, &
              myo2nei,p1_hydr,p2_bind,p3_pi_rele,p4_adp_rele,p5_ubind,p_scale,dt,invmyp2len2,invmyo2len2, &
              lbind,rands,xmyp2,ymyp2,zmyp2,xmyo2,ymyo2,zmyo2,xfa,yfa,zfa)

   implicit none

   integer,value::myp2num,myo2num,neinum,tension2,nrand,myp2len,myo2len
   integer natp,jrand

   integer nm,jh,jtem,jmyo,ncan,jn,ja,jf,j,jexit,j1,j2,omp_get_thread_num,seed,iseed
   integer,allocatable,dimension(:),intent(in)::apos,filid,fa1stbound
   integer,allocatable,dimension(:,:)::myp2typ,fa2myp2,myo2typ,fa2myo2,apar
   integer,allocatable,dimension(:,:),intent(in)::myp2head,myp2body,myo2head,myo2body
   integer,allocatable,dimension(:,:,:),intent(in)::myp2nei,myo2nei
   integer acan(neinum)!,jtem0(nthreads)

   real(kind=8)::d2max,dx,dy,dz,d2,psum,rtem,rat0,normal,pact,ratio,lbind2
   real(kind=8)::pcan(neinum)
   real(kind=8),value::p1_hydr,p2_bind,p3_pi_rele,p4_adp_rele,p5_ubind,p_scale,dt,invmyp2len2,invmyo2len2,lbind
   real(kind=8),allocatable,dimension(:),intent(in)::rands,xmyp2,ymyp2,zmyp2,xmyo2,ymyo2,zmyo2,xfa,yfa,zfa


!     changing Myp2 myosin head status:

!   do tid=1,nthreads

!      jtem0(tid)=jrand+4*myp2num+tid-1

!   end do

   rat0=0.8

   normal=1.0d0/(1.0d0-rat0)

   normal=normal*normal


   lbind2=lbind*lbind


   CALL SYSTEM_CLOCK(COUNT=seed)



!$omp parallel &
!$omp default(none) &
!$omp private(nm,jh,jtem,jmyo,ncan,jn,ja,jf,j,jexit,acan,pcan,j1,j2,iseed) &
!$omp private(d2max,dx,dy,dz,d2,psum,rtem,pact,ratio) &
!$omp shared(myp2num,neinum,tension2,nrand,rat0,normal,seed) &
!$omp shared(jrand,apos,filid,fa1stbound,myp2typ,fa2myp2,apar,myp2head,myp2nei) &
!$omp shared(p1_hydr,p2_bind,p3_pi_rele,p4_adp_rele,p5_ubind,dt,lbind2) &
!$omp shared(rands,xmyp2,ymyp2,zmyp2,xfa,yfa,zfa,invmyp2len2,myp2body,myp2len) &
!$omp reduction(+:natp)

   iseed=seed+omp_get_thread_num( )


!$omp do schedule (guided,32)


   do nm=1,myp2num

      do jh=1,4

         jtem=(nm-1)*4+jh+jrand

         if(jtem<=nrand)then
            rtem=rands(jtem)
         else

            rtem=r8_uniform_01(iseed)

!            rtem=0.0d0
         end if

         if(myp2typ(jh,nm)==1)then

            j1=myp2body(1,nm)

            j2=myp2body(myp2len,nm)


            dx=xmyp2(j1)-xmyp2(j2)
            dy=ymyp2(j1)-ymyp2(j2)
            dz=zmyp2(j1)-zmyp2(j2)

            d2=dx*dx+dy*dy+dz*dz

            ratio=d2*invmyp2len2

            if(ratio>rat0)then
               pact=normal*(ratio-rat0)**2
            else
               pact=0.0d0
            endif


            if(pact*p1_hydr*dt>rtem)then
               myp2typ(jh,nm)=2
!               natp=natp+1
            end if

         else if(myp2typ(jh,nm)==2)then

            if(p2_bind*dt>rtem)then

!              max binding distance square:

               d2max=15.0d0*15.0d0

               jmyo=myp2head(jh,nm)

               ncan=0

               do jn=1,neinum

                  if(myp2nei(jn,jh,nm)==0)then
                     exit
                  end if

                  ja=myp2nei(jn,jh,nm)

                  if(apar(1,ja)/=0)then
                     cycle
                  end if

                  jf=filid(ja)

                  if(apos(ja)<=fa1stbound(jf).and.tension2==1)then
                     cycle
                  end if

                  dx=xfa(ja)-xmyp2(jmyo)
                  dy=yfa(ja)-ymyp2(jmyo)
                  dz=zfa(ja)-zmyp2(jmyo)

                  d2=dx*dx+dy*dy+dz*dz

                  if(d2<d2max)then

                     ncan=ncan+1

                     pcan(ncan)=lbind2*(d2max-d2)/d2/(d2max-lbind2)!1.0d0-d2/d2max

                     pcan(ncan)=min(1.0d0,pcan(ncan))

!                     pcan(ncan)=1.0d0-d2/d2max

                     acan(ncan)=ja

                  end if

               end do

               if(ncan>0)then


!                  tid=omp_get_thread_num()+1
!                  jtem0(tid)=jtem0(tid)+nthreads
!                  jtem=jtem0(tid)


!                  if(jtem<=nrand)then
!                     rtem=rands(jtem)
!                  else
!                     rtem=0.5d0
!                  end if

!                  call random_number(rtem)

                  rtem=r8_uniform_01(iseed)

                  psum=sum(pcan(1:ncan))

                  if(psum>1.0d0)then
                     pcan(1:ncan)=pcan(1:ncan)/psum
                  end if

                  psum=0.0d0

                  do j=1,ncan

                     psum=psum+pcan(j)

                     if(psum>rtem)then

                        jexit=0

                        ja=acan(j)

                        !$omp critical

                        if(apar(1,ja)==0)then

                           myp2typ(jh,nm)=3
                           fa2myp2(jh,nm)=ja
                           apar(1,ja)=jh
                           apar(2,ja)=nm
                           jexit=1
                        end if

                        !$omp end critical

                        if(jexit==1)then
                           exit
                        end if

                     end if

                  end do

               end if


            end if

         else if(myp2typ(jh,nm)==3)then

            if(p3_pi_rele*dt>rtem)then
               myp2typ(jh,nm)=4
               natp=natp+1
            end if

         else if(myp2typ(jh,nm)==4)then

            if(p4_adp_rele*dt>rtem)then
               myp2typ(jh,nm)=5
            end if

         else 

            if(p5_ubind*dt>rtem)then
               myp2typ(jh,nm)=1
               ja=fa2myp2(jh,nm)
               apar(1:2,ja)=0
               fa2myp2(jh,nm)=0

            end if

         end if

      end do

   end do

!$omp end do nowait
!$omp end parallel

!   jrand=maxval(jtem0(1:nthreads))!jrand+2*myp2num

   jrand=jrand+4*myp2num

!---------------------------------------------

!     changing Myo2 myosin head status:

!   do tid=1,nthreads
!      jtem0(tid)=jrand+2*myo2num+tid-1
!   end do

   CALL SYSTEM_CLOCK(COUNT=seed)


!  scaling down the rates:
   dt=dt*p_scale

!$omp parallel &
!$omp default(none) &
!$omp private(nm,jh,jtem,jmyo,ncan,jn,ja,jf,j,jexit,acan,pcan,j1,j2,iseed) &
!$omp private(d2max,dx,dy,dz,d2,psum,rtem,ratio,pact) &
!$omp shared(myp2num,myo2num,neinum,rat0,normal,nrand,seed) &
!$omp shared(jrand,apos,filid,fa1stbound,myo2typ,fa2myo2,apar,myo2head,myo2nei) &
!$omp shared(p1_hydr,p2_bind,p3_pi_rele,p4_adp_rele,p5_ubind,dt,lbind2) &
!$omp shared(rands,xmyo2,ymyo2,zmyo2,xfa,yfa,zfa,invmyo2len2,myo2body,myo2len,tension2) &
!$omp reduction(+:natp)


   iseed=seed+omp_get_thread_num( )


!$omp do schedule (guided,32)

   do nm=1,myo2num

      do jh=1,2

         jtem=(nm-1)*2+jh+jrand

         if(jtem<=nrand)then

            rtem=rands(jtem)
         else

            rtem=r8_uniform_01(iseed)

!            rtem=0.0d0
         end if


         if(myo2typ(jh,nm)==1)then

            j1=myo2body(1,nm)

            j2=myo2body(myo2len,nm)

            dx=xmyo2(j1)-xmyo2(j2)
            dy=ymyo2(j1)-ymyo2(j2)
            dz=zmyo2(j1)-zmyo2(j2)

            d2=dx*dx+dy*dy+dz*dz

            ratio=d2*invmyo2len2

            if(ratio>rat0)then
               pact=normal*(ratio-rat0)**2
            else
               pact=0.0d0
            endif


            if(pact*p1_hydr*dt>rtem)then
               myo2typ(jh,nm)=2
!               natp=natp+1

            end if

         else if(myo2typ(jh,nm)==2)then

            if(p2_bind*dt>rtem)then

!              max binding distance square:

               d2max=15.0d0*15.0d0

               jmyo=myo2head(jh,nm)

               ncan=0

               do jn=1,neinum

                  if(myo2nei(jn,jh,nm)==0)then
                     exit
                  end if

                  ja=myo2nei(jn,jh,nm)

                  if(apar(1,ja)/=0)then
                     cycle
                  end if

                  jf=filid(ja)

                  if(apos(ja)<=fa1stbound(jf).and.tension2==1)then
                     cycle
                  end if

                  dx=xfa(ja)-xmyo2(jmyo)
                  dy=yfa(ja)-ymyo2(jmyo)
                  dz=zfa(ja)-zmyo2(jmyo)

                  d2=dx*dx+dy*dy+dz*dz

                  if(d2<d2max)then

                     ncan=ncan+1

                     pcan(ncan)=lbind2*(d2max-d2)/d2/(d2max-lbind2)!1.0d0-d2/d2max

                     pcan(ncan)=min(1.0d0,pcan(ncan))

!                     pcan(ncan)=1.0d0-d2/d2max

                     acan(ncan)=ja

                  end if

               end do

               if(ncan>0)then


!                  tid=omp_get_thread_num()+1
!                  jtem0(tid)=jtem0(tid)+nthreads
!                  jtem=jtem0(tid)!jrand+4*myo2num+1


!                  if(jtem<=nrand)then
!                     rtem=rands(jtem)
!                  else
!                     rtem=0.5d0
!                  end if

                  rtem=r8_uniform_01(iseed)

                  psum=sum(pcan(1:ncan))

                  if(psum>1.0d0)then
                     pcan(1:ncan)=pcan(1:ncan)/psum
                  end if

                  psum=0.0d0

                  do j=1,ncan

                     psum=psum+pcan(j)

                     if(psum>rtem)then

                        jexit=0

                        ja=acan(j)

                        !$omp critical

                        if(apar(1,ja)==0)then

                           myo2typ(jh,nm)=3
                           fa2myo2(jh,nm)=ja
                           apar(1,ja)=jh
                           apar(2,ja)=nm+myp2num
                           jexit=1
                        end if

                        !$omp end critical

                        if(jexit==1)then
                           exit
                        end if

                     end if

                  end do

               end if


            end if

         else if(myo2typ(jh,nm)==3)then

            if(p3_pi_rele*dt>rtem)then
               myo2typ(jh,nm)=4
               natp=natp+1
            end if

         else if(myo2typ(jh,nm)==4)then

            if(p4_adp_rele*dt>rtem)then
               myo2typ(jh,nm)=5
            end if

         else

            if(p5_ubind*dt>rtem)then
               myo2typ(jh,nm)=1
               ja=fa2myo2(jh,nm)
               apar(1:2,ja)=0
               fa2myo2(jh,nm)=0
            end if

         end if

      end do

   end do

!$omp end do nowait
!$omp end parallel

   jrand=jrand+2*myo2num



   end subroutine
!=========================================================

   subroutine crlkcycle(jrand,nrand,crlknum,crlknum1,fanum,lklen,neinum,fa2lk,apar, &
              filid,astart,apos,alen,lkstart,lktyp,fa1stbound,lknei,plk_ubind1,plk_ubind2, &
              dt,plk_bind,lbind,rands,xlk,ylk,zlk,xfa,yfa,zfa)

   implicit none

   integer nl,jl,jtem,ja,jf,jstart,jchange,jlk,filid0,ncan,jn,j,ja1,jexit,omp_get_thread_num,seed,iseed
   integer jrand
   integer,value::nrand,crlknum,crlknum1,fanum,lklen,neinum
   integer,allocatable,dimension(:,:)::fa2lk,apar
   integer,allocatable,dimension(:),intent(in)::filid,astart,apos,alen,lkstart!,lktyp
   integer,allocatable,dimension(:)::fa1stbound,lktyp
   integer,allocatable,dimension(:,:,:),intent(in)::lknei
   integer acan(neinum)!,jtem0(nthreads)

   real(kind=8)::rtem,plk_off,d2max,dx,dy,dz,d2,psum,lbind2
   real(kind=8),value::plk_ubind1,plk_ubind2,dt,plk_bind,lbind
   real(kind=8),allocatable,dimension(:),intent(in)::rands,xlk,ylk,zlk,xfa,yfa,zfa
   real(kind=8) pcan(neinum)

!     update binding status of crosslinkers to actin

   jchange=0

!   do tid=1,nthreads
!      jtem0(tid)=jrand+2*crlknum+tid-1
!   end do

!   jtem0=jrand+2*crlknum

   lbind2=lbind*lbind


   CALL SYSTEM_CLOCK(COUNT=seed)


!$omp parallel &
!$omp default(none) &
!$omp private(nl,jl,jtem,ja,jf,jstart,jlk,filid0,ncan,jn,j,ja1,jexit,iseed) &
!$omp private(acan,pcan,rtem,plk_off,d2max,dx,dy,dz,d2,psum) &
!$omp shared(jrand,nrand,crlknum,crlknum1,lklen,neinum,seed) &
!$omp shared(fa2lk,apar,lkstart,lktyp,filid,astart,apos,fa1stbound,lknei) &
!$omp shared(plk_ubind1,plk_ubind2,dt,plk_bind,lbind2,rands,xlk,ylk,zlk,xfa,yfa,zfa) &
!$omp reduction(+:jchange)

   iseed=seed+omp_get_thread_num( )

!$omp do schedule (guided, 32)


   do nl=1,crlknum

      if(lktyp(nl)==0)then
         cycle
      end if

      do jl=1,2

         jtem=(nl-1)*2+jl+jrand

         if(jtem<nrand)then
            rtem=rands(jtem)
         else

            rtem=r8_uniform_01(iseed)

!            rtem=0.0d0
         end if

         if(fa2lk(jl,nl)>0)then

            if(nl<=crlknum1)then
               plk_off=plk_ubind1*dt
            else
               plk_off=plk_ubind2*dt
            end if

            if(plk_off>rtem)then

               ja=fa2lk(jl,nl)

               apar(1:2,ja)=0

!              update 1st crosslinked bead on a filament:

               jf=filid(ja)

               jstart=astart(jf)


               if(apos(ja)==fa1stbound(jf))then

                  jchange=jchange+1

                  !$omp critical

                  fa1stbound(jf)=-1

                  !$omp end critical

               end if


               fa2lk(jl,nl)=0

!              update first crosslinked bead on the former crosslinking partner filament:

               if(jl==1)then
                  ja=fa2lk(2,nl)
               else
                  ja=fa2lk(1,nl)
               end if

               if(ja>0)then

                  jf=filid(ja)


                  if(apos(ja)==fa1stbound(jf))then

                     jchange=jchange+1

                     !$omp critical

                     fa1stbound(jf)=-1

                     !$omp end critical


                  end if

               end if


            end if

!        crosslink binds to filament:

         else

            if(plk_bind*dt>rtem)then

!              max binding distance square:

               d2max=15.0d0*15.0d0

               if(jl==1)then

                  jlk=lkstart(nl)   !crlk(1,nl)

                  if(fa2lk(2,nl)==0)then
                     filid0=0
                  else
                     filid0=filid(fa2lk(2,nl))
                  end if

               else

                  jlk=lkstart(nl)+lklen-1   !crlk(lklen,nl)

                  if(fa2lk(1,nl)==0)then
                     filid0=0
                  else
                     filid0=filid(fa2lk(1,nl))
                  end if

               end if

               ncan=0

               do jn=1,neinum

                  if(lknei(jn,jl,nl)==0)then
                     exit
                  end if

                  ja=lknei(jn,jl,nl)

                  if(apar(1,ja)/=0)then
                     cycle
                  end if

                  if(filid(ja)==filid0)then
                     cycle
                  end if

                  dx=xfa(ja)-xlk(jlk)
                  dy=yfa(ja)-ylk(jlk)
                  dz=zfa(ja)-zlk(jlk)

                  d2=dx*dx+dy*dy+dz*dz

                  if(d2<d2max)then

                     ncan=ncan+1

                     pcan(ncan)=lbind2*(d2max-d2)/d2/(d2max-lbind2)!1.0d0-d2/d2max

                     pcan(ncan)=min(1.0d0,pcan(ncan))

!                     pcan(ncan)=1.0d0-d2/d2max

                     acan(ncan)=ja

                  end if

               end do

               if(ncan==0)then
                  cycle
               end if


!               tid=omp_get_thread_num()+1
!               jtem0(tid)=jtem0(tid)+nthreads!jrand+2*crlknum+1
!               jtem=jtem0(tid)

!               if(jtem<=nrand)then
!                  rtem=rands(jtem)
!               else
!                  rtem=0.5d0
!               end if

!               call random_number(rtem)


               rtem=r8_uniform_01(iseed)

               psum=sum(pcan(1:ncan))

               if(psum>1.0d0)then
                  pcan(1:ncan)=pcan(1:ncan)/psum
               end if

               psum=0.0d0

               do j=1,ncan

                  psum=psum+pcan(j)

                  if(psum>rtem)then

                     jexit=0

                     ja=acan(j)

                     !$omp critical

                     if(apar(1,ja)==0)then

                        fa2lk(jl,nl)=ja
                        apar(1,ja)=jl
                        apar(2,ja)=-nl

                        if(jl==1)then
                           ja1=fa2lk(2,nl)
                        else
                           ja1=fa2lk(1,nl)
                        end if


                        if(ja1>0)then

                           jf=filid(ja)

                           if(apos(ja)<fa1stbound(jf))then
                              fa1stbound(jf)=apos(ja)
                           end if

                           jf=filid(ja1)

                           if(apos(ja1)<fa1stbound(jf))then
                              fa1stbound(jf)=apos(ja1)
                           end if


                        end if

                        jexit=1

                     end if

                     !$omp end critical

                     if(jexit==1)then
                        exit
                     end if

                  end if

               end do


            end if

         end if

      end do

      if(fa2lk(1,nl)+fa2lk(2,nl)==0)then
         lktyp(nl)=10
      else
         lktyp(nl)=1
      end if

   end do

!$omp end do nowait
!$omp end parallel 

!   jrand=maxval(jtem0(1:nthreads))!jrand+2*crlknum

   jrand=jrand+2*crlknum

   if(jchange==0)then
      return
   end if


!$omp parallel &
!$omp default(none) &
!$omp private(jf,jstart,j,nl) &
!$omp shared(fanum,fa1stbound,astart,alen,apar,fa2lk)
!$omp do schedule (guided, 32)

   do jf=1,fanum

      if(fa1stbound(jf)>0)then
         cycle
      end if

      fa1stbound(jf)=alen(jf)

      jstart=astart(jf)


      do j=1,alen(jf)

         if(apar(2,jstart+j-1)<0)then

            nl=-apar(2,jstart+j-1)

            if(fa2lk(1,nl)>0.and.fa2lk(2,nl)>0)then
               fa1stbound(jf)=j
               exit
            end if

         end if

      end do

   end do

!$omp end do nowait
!$omp end parallel 


   end subroutine

!=========================================================

   subroutine faforce(jforce,fxbond,fybond,fzbond,fxangl,fyangl,fzangl,fanum, &
                 astart,alen,xfa,yfa,zfa,k_a,l_a,kthet,thet0,delta,invdelta,beta)

   implicit none

   integer,value:: jforce,fanum
   integer nf,jf,n1,n2,n3,jstart

   integer,allocatable,intent(in),dimension(:)::astart,alen

   double precision,value::k_a,l_a,kthet,thet0,delta,invdelta,beta
   double precision dx,dy,dz,dx1,dy1,dz1,dx3,dy3,dz3,invdist1,invdist3,invdist
   double precision f0,cos_t0,cos_t,thet,drep,dfx1,dfy1,dfz1,dfx3,dfy3,dfz3,f
   double precision x1,y1,z1,x2,y2,z2,x3,y3,z3,dist
   double precision,allocatable,intent(in),dimension(:)::xfa,yfa,zfa
   double precision,allocatable,dimension(:)::fxbond,fybond,fzbond,fxangl,fyangl,fzangl






!$omp parallel &
!$omp default(none) &
!$omp private(nf,jf,n1,n2,jstart,dx,dy,dz,dist,f,x1,y1,z1,x2,y2,z2) &
!$omp shared(fanum,alen,astart) &
!$omp shared(k_a,l_a) &
!$omp shared(xfa,yfa,zfa,fxbond,fybond,fzbond)

!$omp do schedule(guided,32)

   do nf=1,fanum

      jstart=astart(nf)

      n2=jstart  ! used to be afil(1,nf)

      fxbond(n2)=0.0d0
      fybond(n2)=0.0d0
      fzbond(n2)=0.0d0

      x2=xfa(n2)
      y2=yfa(n2)
      z2=zfa(n2)

      do jf=2,alen(nf)

         n1=n2

         x1=x2
         y1=y2
         z1=z2


         n2=jstart+jf-1  ! used to be afil(jf,nf)

!         if(atyp(n2)==0)then
!            exit
!         end if

         x2=xfa(n2)
         y2=yfa(n2)
         z2=zfa(n2)


         dx=x1-x2
         dy=y1-y2
         dz=z1-z2

         dist=sqrt(dx*dx+dy*dy+dz*dz)

!         invdist=1.0d0/dist

         f=k_a*(dist-l_a)/dist

         fxbond(n2)=f*dx!*invdist
         fybond(n2)=f*dy!*invdist
         fzbond(n2)=f*dz!*invdist

         fxbond(n1)=fxbond(n1)-fxbond(n2)
         fybond(n1)=fybond(n1)-fybond(n2)
         fzbond(n1)=fzbond(n1)-fzbond(n2)

!         fxbond(n2)=fx
!         fybond(n2)=fy
!         fzbond(n2)=fz


      end do

   end do


!$omp end do nowait
!$omp end parallel


!-------------------------------------

   if(jforce==0)then
      return
   end if


!$omp parallel &
!$omp default(none) &
!$omp private(nf,jf,n1,n2,n3,dx1,dy1,dz1,dx3,dy3,dz3,invdist1,invdist3,x2,y2,z2,x3,y3,z3) &
!$omp private(jstart,dx,dy,dz,f0,cos_t0,cos_t,thet,drep,dfx1,dfy1,dfz1,dfx3,dfy3,dfz3) &
!$omp shared(fanum,alen,astart) &
!$omp shared(kthet,thet0,delta,invdelta,beta) &
!$omp shared(xfa,yfa,zfa,fxangl,fyangl,fzangl)

!$omp do schedule(guided,32)

   do nf=1,fanum

      jstart=astart(nf)

      n2=jstart  ! used to be afil(1,nf)

!      fxangl(n2)=0.0d0
!      fyangl(n2)=0.0d0
!      fzangl(n2)=0.0d0

      n3=jstart+1  ! used to be afil(2,nf)

!      fxangl(n3)=0.0d0
!      fyangl(n3)=0.0d0
!      fzangl(n3)=0.0d0

      x3=xfa(n3)
      y3=yfa(n3)
      z3=zfa(n3)

      dx3=x3-xfa(n2)
      dy3=y3-yfa(n2)
      dz3=z3-zfa(n2)

      invdist3=1.0d0/sqrt(dx3*dx3+dy3*dy3+dz3*dz3)

      do jf=2,alen(nf)-1


         dx1=-dx3
         dy1=-dy3
         dz1=-dz3

         invdist1=invdist3

         x2=x3
         y2=y3
         z2=z3

         n1=n2

         n2=n3


         n3=jstart+jf  ! used to be afil(jf+1,nf)

!         if(atyp(n3)==0)then
!            exit
!         end if

         x3=xfa(n3)
         y3=yfa(n3)
         z3=zfa(n3)

         dx3=x3-x2
         dy3=y3-y2
         dz3=z3-z2

         invdist3=1.0d0/sqrt(dx3*dx3+dy3*dy3+dz3*dz3)

         cos_t0=(dx1*dx3+dy1*dy3+dz1*dz3)*invdist1*invdist3


         thet=acos((1.0d0-beta)*cos_t0)

         f0=kthet*(thet-thet0)/sin(thet)*invdelta

!if(n2==na)then
!print*,'angle'
!print*,cos_t0,thet,f0
!endif
!        force on n1 along x:

         dx=dx1+delta

         drep=sqrt(dx*dx+dy1*dy1+dz1*dz1)

         cos_t=(dx*dx3+dy1*dy3+dz1*dz3)/drep*invdist3


         dfx1=f0*(cos_t-cos_t0)


         fxangl(n1)=fxangl(n1)+dfx1

!        force on n1 along y:

         dy=dy1+delta
   
         drep=sqrt(dx1*dx1+dy*dy+dz1*dz1)

         cos_t=(dx1*dx3+dy*dy3+dz1*dz3)/drep*invdist3
   

         dfy1=f0*(cos_t-cos_t0)

         fyangl(n1)=fyangl(n1)+dfy1

!if(n2==na)then
!print*,'force1'
!print*,dfy1,cos_t
!endif

!        force on n1 along z:

         dz=dz1+delta

         drep=sqrt(dx1*dx1+dy1*dy1+dz*dz)

         cos_t=(dx1*dx3+dy1*dy3+dz*dz3)/drep*invdist3


         dfz1=f0*(cos_t-cos_t0)

         fzangl(n1)=fzangl(n1)+dfz1

!        force on n3 along x:

         dx=dx3+delta

         drep=sqrt(dx*dx+dy3*dy3+dz3*dz3)

         cos_t=(dx1*dx+dy1*dy3+dz1*dz3)/drep*invdist1


         dfx3=f0*(cos_t-cos_t0)

         fxangl(n3)=fxangl(n3)+dfx3

!        force on n3 along y:
   
         dy=dy3+delta

         drep=sqrt(dx3*dx3+dy*dy+dz3*dz3)

         cos_t=(dx1*dx3+dy1*dy+dz1*dz3)/drep*invdist1


         dfy3=f0*(cos_t-cos_t0)

         fyangl(n3)=fyangl(n3)+dfy3


!        force on n3 along z:

         dz=dz3+delta

         drep=sqrt(dx3*dx3+dy3*dy3+dz*dz)

         cos_t=(dx1*dx3+dy1*dy3+dz1*dz)/drep*invdist1


         dfz3=f0*(cos_t-cos_t0)

         fzangl(n3)=fzangl(n3)+dfz3

!        forces on n2:

         fxangl(n2)=fxangl(n2)-dfx1-dfx3
         fyangl(n2)=fyangl(n2)-dfy1-dfy3
         fzangl(n2)=fzangl(n2)-dfz1-dfz3


      end do

   end do

!$omp end do nowait
!$omp end parallel

   end subroutine

!=========================================================

   subroutine myoforce(jforce,fxmyp2bond,fymyp2bond,fzmyp2bond,fxmyp2angl,fymyp2angl,fzmyp2angl, &
                       fxmyo2bond,fymyo2bond,fzmyo2bond,fxmyo2angl,fymyo2angl,fzmyo2angl, &
                       myp2num,myp2head,myp2typ,myp2body,myp2len,xmyp2,ymyp2,zmyp2, &
                       myo2num,myo2head,myo2typ,myo2body,myo2len,xmyo2,ymyo2,zmyo2, &
                       k_a,l_mh,l_mb,kmh_thet,kmb_thet,thet_mh1,thet_mh2,thet_mb,delta,invdelta,beta)

   implicit none

   integer,value:: jforce,myp2num,myp2len,myo2num,myo2len
   integer nm,jm,n1,n2,n3

   integer,allocatable,intent(in),dimension(:,:)::myp2head,myp2body,myp2typ,myo2head,myo2body,myo2typ

   double precision,value::k_a,l_mh,l_mb,kmh_thet,kmb_thet,thet_mh1,thet_mh2,thet_mb,delta,invdelta,beta
   double precision dx,dy,dz,dx1,dy1,dz1,dx3,dy3,dz3,invdist1,invdist3,thet0
   double precision f0,cos_t,cos_t0,thet,drep,dfx1,dfy1,dfz1,dfx3,dfy3,dfz3,f,fx,fy,fz
   double precision x1,y1,z1,x2,y2,z2,x3,y3,z3,dist

   real(kind=8),allocatable,intent(in),dimension(:)::xmyp2,ymyp2,zmyp2,xmyo2,ymyo2,zmyo2
   real(kind=8),allocatable,dimension(:)::fxmyp2bond,fymyp2bond,fzmyp2bond,fxmyp2angl,fymyp2angl,fzmyp2angl
   real(kind=8),allocatable,dimension(:)::fxmyo2bond,fymyo2bond,fzmyo2bond,fxmyo2angl,fymyo2angl,fzmyo2angl


!$omp parallel &
!$omp default(none) &
!$omp private(nm,jm,n1,n2,dx,dy,dz,dist,f,x1,y1,z1,x2,y2,z2) &
!$omp shared(myp2num,myp2len,myp2head,myp2body) &
!$omp shared(myo2num,myo2len,myo2head,myo2body) &
!$omp shared(k_a,l_mh,l_mb) &
!$omp shared(xmyp2,ymyp2,zmyp2,fxmyp2bond,fymyp2bond,fzmyp2bond) &
!$omp shared(xmyo2,ymyo2,zmyo2,fxmyo2bond,fymyo2bond,fzmyo2bond)

!$omp do schedule(guided,32)

!  linear stiffness of Myp2:

   do nm=1,myp2num

!     force on the body:

      n2=myp2body(1,nm)

      fxmyp2bond(n2)=0.0d0
      fymyp2bond(n2)=0.0d0
      fzmyp2bond(n2)=0.0d0

      x2=xmyp2(n2)
      y2=ymyp2(n2)
      z2=zmyp2(n2)

      do jm=2,myp2len

         x1=x2
         y1=y2
         z1=z2

         n1=n2


         n2=myp2body(jm,nm)

         x2=xmyp2(n2)
         y2=ymyp2(n2)
         z2=zmyp2(n2)


         dx=x1-x2
         dy=y1-y2
         dz=z1-z2

         dist=sqrt(dx*dx+dy*dy+dz*dz)

!         invdist=1.0d0/dist

         f=k_a*(dist-l_mb)/dist

         fxmyp2bond(n2)=f*dx!*invdist
         fymyp2bond(n2)=f*dy!*invdist
         fzmyp2bond(n2)=f*dz!*invdist

         fxmyp2bond(n1)=fxmyp2bond(n1)-fxmyp2bond(n2)
         fymyp2bond(n1)=fymyp2bond(n1)-fymyp2bond(n2)
         fzmyp2bond(n1)=fzmyp2bond(n1)-fzmyp2bond(n2)

!         fxbond(n2)=fxbond(n2)+fx
!         fybond(n2)=fybond(n2)+fy
!         fzbond(n2)=fzbond(n2)+fz

      end do

!     force on the head:

      do jm=1,2

         n1=myp2head(jm,nm)

         n2=myp2body(1,nm)


         dx=xmyp2(n1)-xmyp2(n2)
         dy=ymyp2(n1)-ymyp2(n2)
         dz=zmyp2(n1)-zmyp2(n2)

         dist=sqrt(dx*dx+dy*dy+dz*dz)

!         invdist=1.0d0/dist

         f=-k_a*(dist-l_mh)/dist

         fxmyp2bond(n1)=f*dx!*invdist
         fymyp2bond(n1)=f*dy!*invdist
         fzmyp2bond(n1)=f*dz!*invdist

         fxmyp2bond(n2)=fxmyp2bond(n2)-fxmyp2bond(n1)
         fymyp2bond(n2)=fymyp2bond(n2)-fymyp2bond(n1)
         fzmyp2bond(n2)=fzmyp2bond(n2)-fzmyp2bond(n1)

      end do


      do jm=3,4

         n1=myp2head(jm,nm)

         n2=myp2body(myp2len,nm)


         dx=xmyp2(n1)-xmyp2(n2)
         dy=ymyp2(n1)-ymyp2(n2)
         dz=zmyp2(n1)-zmyp2(n2)

         dist=sqrt(dx*dx+dy*dy+dz*dz)

!         invdist=1.0d0/dist

         f=-k_a*(dist-l_mh)/dist

         fxmyp2bond(n1)=f*dx!*invdist
         fymyp2bond(n1)=f*dy!*invdist
         fzmyp2bond(n1)=f*dz!*invdist

         fxmyp2bond(n2)=fxmyp2bond(n2)-fxmyp2bond(n1)
         fymyp2bond(n2)=fymyp2bond(n2)-fymyp2bond(n1)
         fzmyp2bond(n2)=fzmyp2bond(n2)-fzmyp2bond(n1)

      end do

   end do

!$omp end do nowait


!  linear stiffness of Myo2:


!$omp do schedule(guided,32)

   do nm=1,myo2num

!     force on the body:

      n2=myo2body(1,nm)

      fxmyo2bond(n2)=0.0d0
      fymyo2bond(n2)=0.0d0
      fzmyo2bond(n2)=0.0d0

      x2=xmyo2(n2)
      y2=ymyo2(n2)
      z2=zmyo2(n2)

      do jm=2,myo2len

         x1=x2
         y1=y2
         z1=z2

         n1=n2


         n2=n2+1!myo2body(jm,nm)

         x2=xmyo2(n2)
         y2=ymyo2(n2)
         z2=zmyo2(n2)


         dx=x1-x2
         dy=y1-y2
         dz=z1-z2

         dist=sqrt(dx*dx+dy*dy+dz*dz)


         f=k_a*(dist-l_mb)/dist

         fxmyo2bond(n2)=f*dx!*invdist
         fymyo2bond(n2)=f*dy!*invdist
         fzmyo2bond(n2)=f*dz!*invdist

         fxmyo2bond(n1)=fxmyo2bond(n1)-fxmyo2bond(n2)
         fymyo2bond(n1)=fymyo2bond(n1)-fymyo2bond(n2)
         fzmyo2bond(n1)=fzmyo2bond(n1)-fzmyo2bond(n2)


      end do

!     force on the head:

      n2=myo2body(myo2len,nm)

      do jm=1,2

         n1=n2+jm!myo2head(jm,nm)

         dx=xmyo2(n1)-xmyo2(n2)
         dy=ymyo2(n1)-ymyo2(n2)
         dz=zmyo2(n1)-zmyo2(n2)

         dist=sqrt(dx*dx+dy*dy+dz*dz)

         f=-k_a*(dist-l_mh)/dist

         fxmyo2bond(n1)=f*dx!*invdist
         fymyo2bond(n1)=f*dy!*invdist
         fzmyo2bond(n1)=f*dz!*invdist

         fxmyo2bond(n2)=fxmyo2bond(n2)-fxmyo2bond(n1)
         fymyo2bond(n2)=fymyo2bond(n2)-fymyo2bond(n1)
         fzmyo2bond(n2)=fzmyo2bond(n2)-fzmyo2bond(n1)

      end do

   end do

!$omp end do nowait

!$omp end parallel


!--------------------------------------------------
!  calculating angle term:

   if(jforce==0)then
      return
   end if


!$omp parallel &
!$omp default(none) &
!$omp private(nm,jm,n1,n2,n3,dx,dy,dz,dx1,dy1,dz1,dx3,dy3,dz3,invdist1,invdist3,thet0) &
!$omp private(f0,cos_t,cos_t0,thet,drep,dfx1,dfy1,dfz1,dfx3,dfy3,dfz3,x2,y2,z2,x3,y3,z3) &
!$omp shared(myp2num,myp2len,myp2head,myp2body,myp2typ) &
!$omp shared(myo2num,myo2len,myo2head,myo2body,myo2typ) &
!$omp shared(kmh_thet,kmb_thet,thet_mh1,thet_mh2,thet_mb,delta,invdelta,beta) &
!$omp shared(xmyp2,ymyp2,zmyp2,fxmyp2angl,fymyp2angl,fzmyp2angl) &
!$omp shared(xmyo2,ymyo2,zmyo2,fxmyo2angl,fymyo2angl,fzmyo2angl)


!  bending stiffness of Myp2

!$omp do schedule(guided,32)


   do nm=1,myp2num

      n2=myp2body(1,nm)

      fxmyp2angl(n2)=0.0d0
      fymyp2angl(n2)=0.0d0
      fzmyp2angl(n2)=0.0d0

      n3=myp2body(2,nm)

      fxmyp2angl(n3)=0.0d0
      fymyp2angl(n3)=0.0d0
      fzmyp2angl(n3)=0.0d0

      x3=xmyp2(n3)
      y3=ymyp2(n3)
      z3=zmyp2(n3)

      dx3=x3-xmyp2(n2)
      dy3=y3-ymyp2(n2)
      dz3=z3-zmyp2(n2)

      invdist3=1.0d0/sqrt(dx3*dx3+dy3*dy3+dz3*dz3)

      do jm=2,myp2len-1


         dx1=-dx3
         dy1=-dy3
         dz1=-dz3

         invdist1=invdist3

         x2=x3
         y2=y3
         z2=z3

         n1=n2

         n2=n3


         n3=myp2body(jm+1,nm)

         x3=xmyp2(n3)
         y3=ymyp2(n3)
         z3=zmyp2(n3)

         dx3=x3-x2
         dy3=y3-y2
         dz3=z3-z2


         invdist3=1.0d0/sqrt(dx3*dx3+dy3*dy3+dz3*dz3)

         cos_t0=(dx1*dx3+dy1*dy3+dz1*dz3)*invdist1*invdist3

         thet=acos((1.0d0-beta)*cos_t0)

         f0=kmb_thet*(thet-thet_mb)/sin(thet)*invdelta


!        force on n1 along x:

         dx=dx1+delta

         drep=sqrt(dx*dx+dy1*dy1+dz1*dz1)

         cos_t=(dx*dx3+dy1*dy3+dz1*dz3)/drep*invdist3


         dfx1=f0*(cos_t-cos_t0)

         fxmyp2angl(n1)=fxmyp2angl(n1)+dfx1

!        force on n1 along y:

         dy=dy1+delta

         drep=sqrt(dx1*dx1+dy*dy+dz1*dz1)

         cos_t=(dx1*dx3+dy*dy3+dz1*dz3)/drep*invdist3


         dfy1=f0*(cos_t-cos_t0)

         fymyp2angl(n1)=fymyp2angl(n1)+dfy1

!        force on n1 along z:

         dz=dz1+delta

         drep=sqrt(dx1*dx1+dy1*dy1+dz*dz)

         cos_t=(dx1*dx3+dy1*dy3+dz*dz3)/drep*invdist3

         dfz1=f0*(cos_t-cos_t0)

         fzmyp2angl(n1)=fzmyp2angl(n1)+dfz1

!        force on n3 along x:

         dx=dx3+delta

         drep=sqrt(dx*dx+dy3*dy3+dz3*dz3)

         cos_t=(dx1*dx+dy1*dy3+dz1*dz3)/drep*invdist1

         dfx3=f0*(cos_t-cos_t0)

         fxmyp2angl(n3)=dfx3

!        force on n3 along y:

         dy=dy3+delta

         drep=sqrt(dx3*dx3+dy*dy+dz3*dz3)

         cos_t=(dx1*dx3+dy1*dy+dz1*dz3)/drep*invdist1

         dfy3=f0*(cos_t-cos_t0)

         fymyp2angl(n3)=dfy3

!        force on n3 along z:

         dz=dz3+delta

         drep=sqrt(dx3*dx3+dy3*dy3+dz*dz)

         cos_t=(dx1*dx3+dy1*dy3+dz1*dz)/drep*invdist1

         dfz3=f0*(cos_t-cos_t0)

         fzmyp2angl(n3)=dfz3

!        forces on n2:

         fxmyp2angl(n2)=fxmyp2angl(n2)-dfx1-dfx3
         fymyp2angl(n2)=fymyp2angl(n2)-dfy1-dfy3
         fzmyp2angl(n2)=fzmyp2angl(n2)-dfz1-dfz3

      end do

!     force on the head:

      do jm=1,2

         n1=myp2head(jm,nm)


         n2=myp2body(1,nm)

         n3=myp2body(2,nm)


         dx1=xmyp2(n1)-xmyp2(n2)
         dy1=ymyp2(n1)-ymyp2(n2)
         dz1=zmyp2(n1)-zmyp2(n2)

         invdist1=1.0d0/sqrt(dx1*dx1+dy1*dy1+dz1*dz1)


         if(myp2typ(jm,nm)==2.or.myp2typ(jm,nm)==3)then
            thet0=thet_mh2
         else
            thet0=thet_mh1
         end if

         dx3=xmyp2(n3)-xmyp2(n2)
         dy3=ymyp2(n3)-ymyp2(n2)
         dz3=zmyp2(n3)-zmyp2(n2)

         invdist3=1.0d0/sqrt(dx3*dx3+dy3*dy3+dz3*dz3)

         cos_t0=(dx1*dx3+dy1*dy3+dz1*dz3)*invdist1*invdist3

         thet=acos((1.0d0-beta)*cos_t0)

         f0=kmh_thet*(thet-thet0)/sin(thet)*invdelta

!        force on n1 along x:

         dx=dx1+delta

         drep=sqrt(dx*dx+dy1*dy1+dz1*dz1)

         cos_t=(dx*dx3+dy1*dy3+dz1*dz3)/drep*invdist3

         dfx1=f0*(cos_t-cos_t0)

         fxmyp2angl(n1)=dfx1

!        force on n1 along y:

         dy=dy1+delta

         drep=sqrt(dx1*dx1+dy*dy+dz1*dz1)

         cos_t=(dx1*dx3+dy*dy3+dz1*dz3)/drep*invdist3

         dfy1=f0*(cos_t-cos_t0)

         fymyp2angl(n1)=dfy1

!        force on n1 along z:

         dz=dz1+delta

         drep=sqrt(dx1*dx1+dy1*dy1+dz*dz)

         cos_t=(dx1*dx3+dy1*dy3+dz*dz3)/drep*invdist3

         dfz1=f0*(cos_t-cos_t0)

         fzmyp2angl(n1)=dfz1

!        force on n3 along x:

         dx=dx3+delta

         drep=sqrt(dx*dx+dy3*dy3+dz3*dz3)

         cos_t=(dx1*dx+dy1*dy3+dz1*dz3)/drep*invdist1

         dfx3=f0*(cos_t-cos_t0)

         fxmyp2angl(n3)=fxmyp2angl(n3)+dfx3

!        force on n3 along y:

         dy=dy3+delta

         drep=sqrt(dx3*dx3+dy*dy+dz3*dz3)

         cos_t=(dx1*dx3+dy1*dy+dz1*dz3)/drep*invdist1

         dfy3=f0*(cos_t-cos_t0)

         fymyp2angl(n3)=fymyp2angl(n3)+dfy3

!        force on n3 along z:

         dz=dz3+delta

         drep=sqrt(dx3*dx3+dy3*dy3+dz*dz)

         cos_t=(dx1*dx3+dy1*dy3+dz1*dz)/drep*invdist1

         dfz3=f0*(cos_t-cos_t0)

         fzmyp2angl(n3)=fzmyp2angl(n3)+dfz3

!        forces on n2:

         fxmyp2angl(n2)=fxmyp2angl(n2)-dfx1-dfx3
         fymyp2angl(n2)=fymyp2angl(n2)-dfy1-dfy3
         fzmyp2angl(n2)=fzmyp2angl(n2)-dfz1-dfz3

      end do



      do jm=3,4

         n1=myp2head(jm,nm)

         n2=myp2body(myp2len,nm)

         n3=myp2body(myp2len-1,nm)


         dx1=xmyp2(n1)-xmyp2(n2)
         dy1=ymyp2(n1)-ymyp2(n2)
         dz1=zmyp2(n1)-zmyp2(n2)

         invdist1=1.0d0/sqrt(dx1*dx1+dy1*dy1+dz1*dz1)


         if(myp2typ(jm,nm)==2.or.myp2typ(jm,nm)==3)then
            thet0=thet_mh2
         else
            thet0=thet_mh1
         end if

         dx3=xmyp2(n3)-xmyp2(n2)
         dy3=ymyp2(n3)-ymyp2(n2)
         dz3=zmyp2(n3)-zmyp2(n2)

         invdist3=1.0d0/sqrt(dx3*dx3+dy3*dy3+dz3*dz3)

         cos_t0=(dx1*dx3+dy1*dy3+dz1*dz3)*invdist1*invdist3

         thet=acos((1.0d0-beta)*cos_t0)

         f0=kmh_thet*(thet-thet0)/sin(thet)*invdelta

!        force on n1 along x:

         dx=dx1+delta

         drep=sqrt(dx*dx+dy1*dy1+dz1*dz1)

         cos_t=(dx*dx3+dy1*dy3+dz1*dz3)/drep*invdist3

         dfx1=f0*(cos_t-cos_t0)

         fxmyp2angl(n1)=dfx1

!        force on n1 along y:

         dy=dy1+delta

         drep=sqrt(dx1*dx1+dy*dy+dz1*dz1)

         cos_t=(dx1*dx3+dy*dy3+dz1*dz3)/drep*invdist3

         dfy1=f0*(cos_t-cos_t0)

         fymyp2angl(n1)=dfy1

!        force on n1 along z:

         dz=dz1+delta

         drep=sqrt(dx1*dx1+dy1*dy1+dz*dz)

         cos_t=(dx1*dx3+dy1*dy3+dz*dz3)/drep*invdist3

         dfz1=f0*(cos_t-cos_t0)

         fzmyp2angl(n1)=dfz1

!        force on n3 along x:

         dx=dx3+delta

         drep=sqrt(dx*dx+dy3*dy3+dz3*dz3)

         cos_t=(dx1*dx+dy1*dy3+dz1*dz3)/drep*invdist1

         dfx3=f0*(cos_t-cos_t0)

         fxmyp2angl(n3)=fxmyp2angl(n3)+dfx3

!        force on n3 along y:

         dy=dy3+delta

         drep=sqrt(dx3*dx3+dy*dy+dz3*dz3)

         cos_t=(dx1*dx3+dy1*dy+dz1*dz3)/drep*invdist1

         dfy3=f0*(cos_t-cos_t0)

         fymyp2angl(n3)=fymyp2angl(n3)+dfy3

!        force on n3 along z:

         dz=dz3+delta

         drep=sqrt(dx3*dx3+dy3*dy3+dz*dz)

         cos_t=(dx1*dx3+dy1*dy3+dz1*dz)/drep*invdist1

         dfz3=f0*(cos_t-cos_t0)

         fzmyp2angl(n3)=fzmyp2angl(n3)+dfz3

!        forces on n2:

         fxmyp2angl(n2)=fxmyp2angl(n2)-dfx1-dfx3
         fymyp2angl(n2)=fymyp2angl(n2)-dfy1-dfy3
         fzmyp2angl(n2)=fzmyp2angl(n2)-dfz1-dfz3

      end do



   end do

!$omp end do nowait


!  bending stiffness of Myp2


!$omp do schedule(guided,32)


   do nm=1,myo2num

      n2=myo2body(1,nm)

      fxmyo2angl(n2)=0.0d0
      fymyo2angl(n2)=0.0d0
      fzmyo2angl(n2)=0.0d0

      n3=myo2body(2,nm)

      fxmyo2angl(n3)=0.0d0
      fymyo2angl(n3)=0.0d0
      fzmyo2angl(n3)=0.0d0

      x3=xmyo2(n3)
      y3=ymyo2(n3)
      z3=zmyo2(n3)

      dx3=x3-xmyo2(n2)
      dy3=y3-ymyo2(n2)
      dz3=z3-zmyo2(n2)

      invdist3=1.0d0/sqrt(dx3*dx3+dy3*dy3+dz3*dz3)

      do jm=2,myo2len-1


         dx1=-dx3
         dy1=-dy3
         dz1=-dz3

         invdist1=invdist3

         x2=x3
         y2=y3
         z2=z3

         n1=n2

         n2=n3


         n3=n3+1!myo2body(jm+1,nm)

         x3=xmyo2(n3)
         y3=ymyo2(n3)
         z3=zmyo2(n3)

         dx3=x3-x2
         dy3=y3-y2
         dz3=z3-z2

         invdist3=1.0d0/sqrt(dx3*dx3+dy3*dy3+dz3*dz3)

         cos_t0=(dx1*dx3+dy1*dy3+dz1*dz3)*invdist1*invdist3

         thet=acos((1.0d0-beta)*cos_t0)

         f0=kmb_thet*(thet-thet_mb)/sin(thet)*invdelta


!        force on n1 along x:

         dx=dx1+delta

         drep=sqrt(dx*dx+dy1*dy1+dz1*dz1)

         cos_t=(dx*dx3+dy1*dy3+dz1*dz3)/drep*invdist3


         dfx1=f0*(cos_t-cos_t0)

         fxmyo2angl(n1)=fxmyo2angl(n1)+dfx1

!        force on n1 along y:

         dy=dy1+delta

         drep=sqrt(dx1*dx1+dy*dy+dz1*dz1)

         cos_t=(dx1*dx3+dy*dy3+dz1*dz3)/drep*invdist3


         dfy1=f0*(cos_t-cos_t0)

         fymyo2angl(n1)=fymyo2angl(n1)+dfy1

!        force on n1 along z:

         dz=dz1+delta

         drep=sqrt(dx1*dx1+dy1*dy1+dz*dz)

         cos_t=(dx1*dx3+dy1*dy3+dz*dz3)/drep*invdist3

         dfz1=f0*(cos_t-cos_t0)

         fzmyo2angl(n1)=fzmyo2angl(n1)+dfz1

!        force on n3 along x:

         dx=dx3+delta

         drep=sqrt(dx*dx+dy3*dy3+dz3*dz3)

         cos_t=(dx1*dx+dy1*dy3+dz1*dz3)/drep*invdist1

         dfx3=f0*(cos_t-cos_t0)

         fxmyo2angl(n3)=dfx3

!        force on n3 along y:

         dy=dy3+delta

         drep=sqrt(dx3*dx3+dy*dy+dz3*dz3)

         cos_t=(dx1*dx3+dy1*dy+dz1*dz3)/drep*invdist1

         dfy3=f0*(cos_t-cos_t0)

         fymyo2angl(n3)=dfy3

!        force on n3 along z:

         dz=dz3+delta

         drep=sqrt(dx3*dx3+dy3*dy3+dz*dz)

         cos_t=(dx1*dx3+dy1*dy3+dz1*dz)/drep*invdist1

         dfz3=f0*(cos_t-cos_t0)

         fzmyo2angl(n3)=dfz3

!        forces on n2:

         fxmyo2angl(n2)=fxmyo2angl(n2)-dfx1-dfx3
         fymyo2angl(n2)=fymyo2angl(n2)-dfy1-dfy3
         fzmyo2angl(n2)=fzmyo2angl(n2)-dfz1-dfz3

      end do

!     force on the head:

      n2=myo2body(myo2len,nm)
      n3=n2-1

      do jm=1,2

         n1=n2+jm!myo2head(jm,nm)


         dx1=xmyo2(n1)-xmyo2(n2)
         dy1=ymyo2(n1)-ymyo2(n2)
         dz1=zmyo2(n1)-zmyo2(n2)

         invdist1=1.0d0/sqrt(dx1*dx1+dy1*dy1+dz1*dz1)


         if(myo2typ(jm,nm)==2.or.myo2typ(jm,nm)==3)then
            thet0=thet_mh2
         else
            thet0=thet_mh1
         end if

         dx3=xmyo2(n3)-xmyo2(n2)
         dy3=ymyo2(n3)-ymyo2(n2)
         dz3=zmyo2(n3)-zmyo2(n2)

         invdist3=1.0d0/sqrt(dx3*dx3+dy3*dy3+dz3*dz3)

         cos_t0=(dx1*dx3+dy1*dy3+dz1*dz3)*invdist1*invdist3

         thet=acos((1.0d0-beta)*cos_t0)

         f0=kmh_thet*(thet-thet0)/sin(thet)*invdelta

!        force on n1 along x:

         dx=dx1+delta

         drep=sqrt(dx*dx+dy1*dy1+dz1*dz1)

         cos_t=(dx*dx3+dy1*dy3+dz1*dz3)/drep*invdist3

         dfx1=f0*(cos_t-cos_t0)

         fxmyo2angl(n1)=dfx1

!        force on n1 along y:

         dy=dy1+delta

         drep=sqrt(dx1*dx1+dy*dy+dz1*dz1)

         cos_t=(dx1*dx3+dy*dy3+dz1*dz3)/drep*invdist3

         dfy1=f0*(cos_t-cos_t0)

         fymyo2angl(n1)=dfy1

!        force on n1 along z:

         dz=dz1+delta

         drep=sqrt(dx1*dx1+dy1*dy1+dz*dz)

         cos_t=(dx1*dx3+dy1*dy3+dz*dz3)/drep*invdist3

         dfz1=f0*(cos_t-cos_t0)

         fzmyo2angl(n1)=dfz1

!        force on n3 along x:

         dx=dx3+delta

         drep=sqrt(dx*dx+dy3*dy3+dz3*dz3)

         cos_t=(dx1*dx+dy1*dy3+dz1*dz3)/drep*invdist1

         dfx3=f0*(cos_t-cos_t0)

         fxmyo2angl(n3)=fxmyo2angl(n3)+dfx3

!        force on n3 along y:

         dy=dy3+delta

         drep=sqrt(dx3*dx3+dy*dy+dz3*dz3)

         cos_t=(dx1*dx3+dy1*dy+dz1*dz3)/drep*invdist1

         dfy3=f0*(cos_t-cos_t0)

         fymyo2angl(n3)=fymyo2angl(n3)+dfy3

!        force on n3 along z:

         dz=dz3+delta

         drep=sqrt(dx3*dx3+dy3*dy3+dz*dz)

         cos_t=(dx1*dx3+dy1*dy3+dz1*dz)/drep*invdist1

         dfz3=f0*(cos_t-cos_t0)

         fzmyo2angl(n3)=fzmyo2angl(n3)+dfz3

!        forces on n2:

         fxmyo2angl(n2)=fxmyo2angl(n2)-dfx1-dfx3
         fymyo2angl(n2)=fymyo2angl(n2)-dfy1-dfy3
         fzmyo2angl(n2)=fzmyo2angl(n2)-dfz1-dfz3

      end do

   end do

!$omp end do nowait

!$omp end parallel

   end subroutine

!=========================================================

   subroutine lkforce(jforce,fxbond,fybond,fzbond,fxangl,fyangl,fzangl,crlknum1,crlknum2, &
           lkstart,lklen,xlk,ylk,zlk,k_lk,l_lk1,l_lk2,kthet,thet0,delta,invdelta,beta)

   implicit none

   integer,value:: jforce,crlknum1,crlknum2,lklen
   integer nl,jl,n1,n2,n3

   integer,allocatable,intent(in),dimension(:)::lkstart

   double precision,value::k_lk,l_lk1,l_lk2,kthet,thet0,delta,invdelta,beta
   double precision dx,dy,dz,dx1,dy1,dz1,dx3,dy3,dz3,invdist1,invdist3,l_lk,dist
   double precision f0,cos_t,cos_t0,thet,drep,dfx1,dfy1,dfz1,dfx3,dfy3,dfz3,f

   double precision,allocatable,intent(in),dimension(:)::xlk,ylk,zlk
   double precision,allocatable,dimension(:)::fxbond,fybond,fzbond,fxangl,fyangl,fzangl


!$omp parallel &
!$omp default(none) &
!$omp private(nl,n1,n2,n3,dx,dy,dz,dist,l_lk,f) &
!$omp shared(crlknum1,crlknum2,lkstart,k_lk,l_lk1,l_lk2) &
!$omp shared(xlk,ylk,zlk,fxbond,fybond,fzbond)

!$omp do schedule(guided,32)

   do nl=1,crlknum1+crlknum2

      if(nl<=crlknum1)then
         l_lk=l_lk1
      else
         l_lk=l_lk2
      end if


      n1=lkstart(nl)   !crlk(1,nl)
      n2=n1+1          !crlk(2,nl)
      n3=n2+1          !crlk(3,nl)

      dx=xlk(n1)-xlk(n2)
      dy=ylk(n1)-ylk(n2)
      dz=zlk(n1)-zlk(n2)

      dist=sqrt(dx*dx+dy*dy+dz*dz)

!      invdist=1.0d0/dist

      f=k_lk*(dist-l_lk)/dist

      fxbond(n2)=f*dx!*invdist
      fybond(n2)=f*dy!*invdist
      fzbond(n2)=f*dz!*invdist

      fxbond(n1)=-fxbond(n2)
      fybond(n1)=-fybond(n2)
      fzbond(n1)=-fzbond(n2)



      dx=xlk(n3)-xlk(n2)
      dy=ylk(n3)-ylk(n2)
      dz=zlk(n3)-zlk(n2)

      dist=sqrt(dx*dx+dy*dy+dz*dz)

!      invdist=1.0d0/dist

      f=-k_lk*(dist-l_lk)/dist

      fxbond(n3)=f*dx!*invdist
      fybond(n3)=f*dy!*invdist
      fzbond(n3)=f*dz!*invdist

      fxbond(n2)=fxbond(n2)-fxbond(n3)
      fybond(n2)=fybond(n2)-fybond(n3)
      fzbond(n2)=fzbond(n2)-fzbond(n3)


   end do

!$omp end do nowait
!$omp end parallel

!--------------------------------------------------
!  calculating angle term:

   if(jforce==0)then
      return
   end if

!$omp parallel &
!$omp default(none) &
!$omp private(nl,n1,n2,n3,dx,dy,dz,dx1,dy1,dz1,dx3,dy3,dz3,invdist1,invdist3) &
!$omp private(f0,cos_t,cos_t0,thet,drep,dfx1,dfy1,dfz1,dfx3,dfy3,dfz3) &
!$omp shared(crlknum1,crlknum2,lkstart,kthet,thet0,delta,invdelta,beta) &
!$omp shared(xlk,ylk,zlk,fxangl,fyangl,fzangl)

!$omp do schedule(guided,32)

   do nl=1,crlknum1+crlknum2

      n1=lkstart(nl)   !crlk(1,nl)
      n2=n1+1          !crlk(2,nl)
      n3=n2+1          !crlk(3,nl)


      dx1=xlk(n1)-xlk(n2)
      dy1=ylk(n1)-ylk(n2)
      dz1=zlk(n1)-zlk(n2)

      invdist1=1.0d0/sqrt(dx1*dx1+dy1*dy1+dz1*dz1)

      dx3=xlk(n3)-xlk(n2)
      dy3=ylk(n3)-ylk(n2)
      dz3=zlk(n3)-zlk(n2)

      invdist3=1.0d0/sqrt(dx3*dx3+dy3*dy3+dz3*dz3)

      cos_t0=(dx1*dx3+dy1*dy3+dz1*dz3)*invdist1*invdist3

      thet=acos((1.0d0-beta)*cos_t0)

      f0=kthet*(thet-thet0)/sin(thet)*invdelta

!        force on n1 along x:

      dx=dx1+delta

      drep=sqrt(dx*dx+dy1*dy1+dz1*dz1)

      cos_t=(dx*dx3+dy1*dy3+dz1*dz3)/drep*invdist3

      fxangl(n1)=f0*(cos_t-cos_t0)

!         fxangl(n1)=fxangl(n1)+dfx1

!        force on n1 along y:

      dy=dy1+delta

      drep=sqrt(dx1*dx1+dy*dy+dz1*dz1)

      cos_t=(dx1*dx3+dy*dy3+dz1*dz3)/drep*invdist3

      fyangl(n1)=f0*(cos_t-cos_t0)

!         fyangl(n1)=fyangl(n1)+dfy1

!        force on n1 along z:

      dz=dz1+delta

      drep=sqrt(dx1*dx1+dy1*dy1+dz*dz)

      cos_t=(dx1*dx3+dy1*dy3+dz*dz3)/drep*invdist3

      fzangl(n1)=f0*(cos_t-cos_t0)

!         fzangl(n1)=fzangl(n1)+dfz1

!        force on n3 along x:

      dx=dx3+delta

      drep=sqrt(dx*dx+dy3*dy3+dz3*dz3)

      cos_t=(dx1*dx+dy1*dy3+dz1*dz3)/drep*invdist1

      fxangl(n3)=f0*(cos_t-cos_t0)

!         fxangl(n3)=fxangl(n3)+dfx3

!        force on n3 along y:

      dy=dy3+delta

      drep=sqrt(dx3*dx3+dy*dy+dz3*dz3)

      cos_t=(dx1*dx3+dy1*dy+dz1*dz3)/drep*invdist1

      fyangl(n3)=f0*(cos_t-cos_t0)

!         fyangl(n3)=fyangl(n3)+dfy3

!        force on n3 along z:

      dz=dz3+delta

      drep=sqrt(dx3*dx3+dy3*dy3+dz*dz)

      cos_t=(dx1*dx3+dy1*dy3+dz1*dz)/drep*invdist1

      fzangl(n3)=f0*(cos_t-cos_t0)

!         fzangl(n3)=fzangl(n3)+dfz3

!        forces on n2:

      fxangl(n2)=-fxangl(n1)-fxangl(n3)
      fyangl(n2)=-fyangl(n1)-fyangl(n3)
      fzangl(n2)=-fzangl(n1)-fzangl(n3)


   end do

!$omp end do nowait
!$omp end parallel

   end subroutine

!=========================================================


   subroutine fabindmyo(myp2num,myo2num,myp2head,fa2myp2,myo2head,fa2myo2,kbind,lbind,invl_a, &
               xfa,yfa,zfa,xmyp2,ymyp2,zmyp2,xmyo2,ymyo2,zmyo2, &
                fxfa,fyfa,fzfa,fxmyp2,fymyp2,fzmyp2,fxmyo2,fymyo2,fzmyo2)




   implicit none

   integer,value::myp2num,myo2num
   integer nm,jm,ja,jmyo

   integer,allocatable,intent(in),dimension(:,:)::myp2head,fa2myp2,myo2head,fa2myo2

   double precision,value::kbind,lbind,invl_a
   double precision dx,dy,dz,dist,f,fx,fy,fz,dx1,dy1,dz1

   double precision,allocatable,intent(in),dimension(:)::xfa,yfa,zfa,xmyp2,ymyp2,zmyp2,xmyo2,ymyo2,zmyo2
   double precision,allocatable,dimension(:)::fxfa,fyfa,fzfa,fxmyp2,fymyp2,fzmyp2,fxmyo2,fymyo2,fzmyo2

!integer na

!na=643

!$omp parallel &
!$omp default(none) &
!$omp private(nm,jm,ja,jmyo,dx,dy,dz,dx1,dy1,dz1,dist,f,fx,fy,fz) &
!$omp shared(myp2num,myp2head,myo2num,myo2head,fa2myp2,kbind,lbind,invl_a) &
!$omp shared(xfa,yfa,zfa,xmyp2,ymyp2,zmyp2,fxfa,fyfa,fzfa,fxmyp2,fymyp2,fzmyp2) &
!$omp shared(fa2myo2,xmyo2,ymyo2,zmyo2,fxmyo2,fymyo2,fzmyo2)


!  binding of Myp2 to actin:


!$omp do schedule(guided,32)


   do nm=1,myp2num

      do jm=1,4

         ja=fa2myp2(jm,nm)

         if(ja==0)then
            cycle
         end if

         jmyo=myp2head(jm,nm)

         dx=xfa(ja)-xmyp2(jmyo)
         dy=yfa(ja)-ymyp2(jmyo)
         dz=zfa(ja)-zmyp2(jmyo)

         dist=sqrt(dx*dx+dy*dy+dz*dz)

         f=kbind*(dist-lbind)/dist

!         invdist=1.0d0/dist

         fx=f*dx!*invdist
         fy=f*dy!*invdist
         fz=f*dz!*invdist

         fxfa(ja)=fxfa(ja)-fx
         fyfa(ja)=fyfa(ja)-fy
         fzfa(ja)=fzfa(ja)-fz

         fxmyp2(jmyo)=fxmyp2(jmyo)+fx
         fymyp2(jmyo)=fymyp2(jmyo)+fy
         fzmyp2(jmyo)=fzmyp2(jmyo)+fz


!        constraining to right angle:

         dx1=(xfa(ja-1)-xfa(ja))*invl_a
         dy1=(yfa(ja-1)-yfa(ja))*invl_a
         dz1=(zfa(ja-1)-zfa(ja))*invl_a

         f=kbind*(dx*dx1+dy*dy1+dz*dz1)

         fx=f*dx1
         fy=f*dy1
         fz=f*dz1

         fxfa(ja)=fxfa(ja)-fx
         fyfa(ja)=fyfa(ja)-fy
         fzfa(ja)=fzfa(ja)-fz

         fxmyp2(jmyo)=fxmyp2(jmyo)+fx
         fymyp2(jmyo)=fymyp2(jmyo)+fy
         fzmyp2(jmyo)=fzmyp2(jmyo)+fz



      end do

   end do

!$omp end do nowait


!--------------------------

!  binding of Myo2 to actin:

!$omp do schedule(guided,32)

   do nm=1,myo2num

      do jm=1,2

         ja=fa2myo2(jm,nm)

         if(ja==0)then
            cycle
         end if

         jmyo=myo2head(jm,nm)

         dx=xfa(ja)-xmyo2(jmyo)
         dy=yfa(ja)-ymyo2(jmyo)
         dz=zfa(ja)-zmyo2(jmyo)

         dist=sqrt(dx*dx+dy*dy+dz*dz)

         f=kbind*(dist-lbind)/dist

         fx=f*dx!*invdist
         fy=f*dy!*invdist
         fz=f*dz!*invdist

         fxfa(ja)=fxfa(ja)-fx
         fyfa(ja)=fyfa(ja)-fy
         fzfa(ja)=fzfa(ja)-fz

         fxmyo2(jmyo)=fxmyo2(jmyo)+fx
         fymyo2(jmyo)=fymyo2(jmyo)+fy
         fzmyo2(jmyo)=fzmyo2(jmyo)+fz

!        constraining to right angle:

         dx1=(xfa(ja-1)-xfa(ja))*invl_a
         dy1=(yfa(ja-1)-yfa(ja))*invl_a
         dz1=(zfa(ja-1)-zfa(ja))*invl_a

         f=kbind*(dx*dx1+dy*dy1+dz*dz1)

         fx=f*dx1
         fy=f*dy1
         fz=f*dz1

         fxfa(ja)=fxfa(ja)-fx
         fyfa(ja)=fyfa(ja)-fy
         fzfa(ja)=fzfa(ja)-fz

         fxmyo2(jmyo)=fxmyo2(jmyo)+fx
         fymyo2(jmyo)=fymyo2(jmyo)+fy
         fzmyo2(jmyo)=fzmyo2(jmyo)+fz

      end do

   end do

!$omp end do nowait

!$omp end parallel

   end subroutine

!=========================================================

   subroutine fabindlk(fxfa,fyfa,fzfa,xfa,yfa,zfa,fxlk,fylk,fzlk,xlk,ylk,zlk, &
                     lkstart,lklen,fa2lk,crlknum,kbind,lbind)

   implicit none

   integer,value::crlknum,lklen
   integer nl,jl,ja,jlk

   integer,allocatable,intent(in),dimension(:)::lkstart
   integer,allocatable,intent(in),dimension(:,:)::fa2lk

   double precision,value::kbind,lbind
   double precision dx,dy,dz,dist,f,fx,fy,fz

   double precision,allocatable,intent(in),dimension(:)::xfa,yfa,zfa,xlk,ylk,zlk
   double precision,allocatable,dimension(:)::fxfa,fyfa,fzfa,fxlk,fylk,fzlk


!$omp parallel &
!$omp default(none) &
!$omp private(nl,jl,ja,jlk,dx,dy,dz,dist,f,fx,fy,fz) &
!$omp shared(crlknum,lklen,lkstart,fa2lk,kbind,lbind) &
!$omp shared(xfa,yfa,zfa,xlk,ylk,zlk,fxfa,fyfa,fzfa,fxlk,fylk,fzlk)
!$omp do schedule(guided,32)

   do nl=1,crlknum

      do jl=1,2

         ja=fa2lk(jl,nl)

         if(ja==0)then
            cycle
         end if

         if(jl==1)then
            jlk=lkstart(nl)   !crlk(1,nl)
         else
            jlk=lkstart(nl)+lklen-1   !crlk(lklen,nl)
         end if

         dx=xfa(ja)-xlk(jlk)
         dy=yfa(ja)-ylk(jlk)
         dz=zfa(ja)-zlk(jlk)

         dist=sqrt(dx*dx+dy*dy+dz*dz)

         f=kbind*(dist-lbind)/dist

!         invdist=1.0d0/dist

         fx=f*dx!*invdist
         fy=f*dy!*invdist
         fz=f*dz!*invdist

         fxfa(ja)=fxfa(ja)-fx
         fyfa(ja)=fyfa(ja)-fy
         fzfa(ja)=fzfa(ja)-fz

         fxlk(jlk)=fxlk(jlk)+fx
         fylk(jlk)=fylk(jlk)+fy
         fzlk(jlk)=fzlk(jlk)+fz

      end do

   end do

!$omp end do nowait
!$omp end parallel

   end subroutine

!=========================================================

   subroutine rforceset(jrforce,nrforce,nmemb,nfa,nmyp2,nmyo2,nlk,pi,rxmemb,rymemb,rzmemb, &
              rxfa,ryfa,rzfa,rxmyp2,rymyp2,rzmyp2,rxmyo2,rymyo2,rzmyo2,rxlk,rylk,rzlk,k_scale,k_scale_lk)

   implicit none

   integer,value::nrforce,nmemb,nfa,nmyp2,nmyo2,nlk
   integer jrforce,j,nmax,n,nh1,nh2

   double precision,value::pi,k_scale,k_scale_lk
   double precision,allocatable,dimension(:,:)::rxfa,ryfa,rzfa,rxmyp2,rymyp2,rzmyp2,rxlk,rylk,rzlk
   double precision,allocatable,dimension(:,:)::rxmemb,rymemb,rzmemb,rxmyo2,rymyo2,rzmyo2
   double precision,allocatable,dimension(:)::rh1,rh2,rx,ry,rz
!   double precision k,kmb
   integer seed,iseed,omp_get_thread_num


   nmax=nmemb+nfa+nmyp2+nmyo2+nlk
   nh1=(nmax+1)/2
   nh2=nmax-nh1

   jrforce=0

!   k=0.1d0
!   kmb=0.1d0


   allocate(rh1(nh1),rh2(nh1),rx(nmax),ry(nmax),rz(nmax))

   CALL SYSTEM_CLOCK(COUNT=seed)


!$omp parallel &
!$omp default(none) &
!$omp private(j,rh1,rh2,rx,ry,rz,n,iseed) &
!$omp shared(nrforce,nmax,nh1,nh2,nmemb,nfa,nmyp2,nmyo2,nlk,pi,k_scale,k_scale_lk) &
!$omp shared(rxfa,ryfa,rzfa,rxmyp2,rymyp2,rzmyp2,rxmyo2,rymyo2,rzmyo2,rxlk,rylk,rzlk,rxmemb,rymemb,rzmemb,seed) 

   iseed=seed+omp_get_thread_num( )

!$omp do schedule(guided,32)


   do j=1,nrforce

!      iseed=seed/127773+j*539


      call r4vec_uniform_01 ( nh1, iseed, rh1 )
      call r4vec_uniform_01 ( nh1, iseed, rh2 )

      rx(1:nh1)=sqrt(-2*log(rh1(1:nh1)))*sin(2*pi*rh2(1:nh1))
      rx(nh1+1:nmax)=sqrt(-2*log(rh1(1:nh2)))*cos(2*pi*rh2(1:nh2))


      call r4vec_uniform_01 ( nh1, iseed, rh1 )
      call r4vec_uniform_01 ( nh1, iseed, rh2 )

      ry(1:nh1)=sqrt(-2*log(rh1(1:nh1)))*sin(2*pi*rh2(1:nh1))
      ry(nh1+1:nmax)=sqrt(-2*log(rh1(1:nh2)))*cos(2*pi*rh2(1:nh2))


      call r4vec_uniform_01 ( nh1, iseed, rh1 )
      call r4vec_uniform_01 ( nh1, iseed, rh2 )

      rz(1:nh1)=sqrt(-2*log(rh1(1:nh1)))*sin(2*pi*rh2(1:nh1))
      rz(nh1+1:nmax)=sqrt(-2*log(rh1(1:nh2)))*cos(2*pi*rh2(1:nh2))

      rxmemb(1:nmemb,j)=rx(1:nmemb)*k_scale
      rymemb(1:nmemb,j)=ry(1:nmemb)*k_scale
      rzmemb(1:nmemb,j)=rz(1:nmemb)*k_scale

      n=nmemb

      rxfa(1:nfa,j)=k_scale_lk*rx(n+1:n+nfa)
      ryfa(1:nfa,j)=k_scale_lk*ry(n+1:n+nfa)
      rzfa(1:nfa,j)=k_scale_lk*rz(n+1:n+nfa)

      n=n+nfa

      rxmyp2(1:nmyp2,j)=k_scale_lk*rx(n+1:n+nmyp2)
      rymyp2(1:nmyp2,j)=k_scale_lk*ry(n+1:n+nmyp2)
      rzmyp2(1:nmyp2,j)=k_scale_lk*rz(n+1:n+nmyp2)

      n=n+nmyp2

      rxmyo2(1:nmyo2,j)=k_scale_lk*rx(n+1:n+nmyo2)
      rymyo2(1:nmyo2,j)=k_scale_lk*ry(n+1:n+nmyo2)
      rzmyo2(1:nmyo2,j)=k_scale_lk*rz(n+1:n+nmyo2)

      n=n+nmyo2

      rxlk(1:nlk,j)=k_scale_lk*rx(n+1:n+nlk)
      rylk(1:nlk,j)=k_scale_lk*ry(n+1:n+nlk)
      rzlk(1:nlk,j)=k_scale_lk*rz(n+1:n+nlk)


   end do


!$omp enddo nowait
!$omp end parallel



   deallocate(rh1,rh2)


   end subroutine

!=========================================================

   subroutine randforce(nmemb,fxmemb,fymemb,fzmemb, &
               nfa,fxfa,fyfa,fzfa,nmyp2,fxmyp2,fymyp2,fzmyp2,nlk,fxlk,fylk,fzlk, &
                jrforce,rxmemb,rymemb,rzmemb, &
                 rxfa,ryfa,rzfa,rxmyp2,rymyp2,rzmyp2,rxlk,rylk,rzlk)

   implicit none

   integer,value::nmemb,nfa,nmyp2,nlk
   integer jrforce

   double precision,allocatable,intent(in),dimension(:,:)::rxfa,ryfa,rzfa,rxmyp2,rymyp2,rzmyp2,rxlk,rylk,rzlk
   double precision,allocatable,intent(in),dimension(:,:)::rxmemb,rymemb,rzmemb

   double precision,allocatable,dimension(:)::fxfa,fyfa,fzfa,fxmyp2,fymyp2,fzmyp2,fxlk,fylk,fzlk
   double precision,allocatable,dimension(:)::fxmemb,fymemb,fzmemb

!   double precision,allocatable,dimension(:)::rh1,rh2,rx,ry,rz



   jrforce=jrforce+1

   fxmemb(1:nmemb)=rxmemb(1:nmemb,jrforce)
   fymemb(1:nmemb)=rymemb(1:nmemb,jrforce)
   fzmemb(1:nmemb)=rzmemb(1:nmemb,jrforce)


   fxfa(1:nfa)=rxfa(1:nfa,jrforce)
   fyfa(1:nfa)=ryfa(1:nfa,jrforce)
   fzfa(1:nfa)=rzfa(1:nfa,jrforce)

   fxmyp2(1:nmyp2)=rxmyp2(1:nmyp2,jrforce)
   fymyp2(1:nmyp2)=rymyp2(1:nmyp2,jrforce)
   fzmyp2(1:nmyp2)=rzmyp2(1:nmyp2,jrforce)

   fxlk(1:nlk)=rxlk(1:nlk,jrforce)
   fylk(1:nlk)=rylk(1:nlk,jrforce)
   fzlk(1:nlk)=rzlk(1:nlk,jrforce)

   end subroutine

!=========================================================

   subroutine newactin(myp2num,myo2num,crlknum,falen,nmono,jdir,nxsol,nphisol, &
               fanum,fanum1,fanum2,nfa,nfa1,nfa2,jsursol,apar,fa1stbound,astart,alen, &
                fadist,apos,filid,jfasol,fa2myp2,fa2myo2,fa2lk,l_a,l_mem,pi, &
                 xfa,yfa,zfa,xmemb,ymemb,zmemb,xnorsurf,ynorsurf,znorsurf,xsurf,ysurf,zsurf)

   implicit none

   integer,value::myp2num,myo2num,crlknum,falen,nmono,jdir,nxsol,nphisol!,ntether
   integer fanum,fanum1,fanum2,nfa,nfa1,nfa2

   integer length,na,nnew,nm,jm,nl,jl,jx0,jp0,j,jx,jp,jxget,jpget,j1,j2,n

   integer ncan!,acan(500)

   integer,allocatable,dimension(:,:),intent(in)::jsursol

   integer,allocatable,dimension(:)::apar(:,:)
   integer,allocatable,dimension(:)::fa1stbound,astart,alen,fadist,apos,filid!,a2mem
   integer,allocatable,dimension(:,:)::jfasol,fa2myp2,fa2myo2,fa2lk

!   integer,allocatable,dimension(:,:)::apartem
!   integer,allocatable,dimension(:)::mark
!   integer,allocatable,dimension(:,:)::jfasoltem

   real(kind=8),value::l_a,l_mem,pi!,p_tether

!   real(kind=8)::distmin,distmax!,pcan(100),dsmall,dp,psum

   real(kind=8)::d2max,dshift,r,xn0,yn0,zn0,distance,x0,y0,z0,phi,rad,dphi,xmb,ymb,zmb
   real(kind=8)::dist2,d2,dx,dy,dz,xn,yn,zn,proj,proj0

   real(kind=8),allocatable,dimension(:)::xfa,yfa,zfa
   real(kind=8),allocatable,dimension(:),intent(in)::xmemb,ymemb,zmemb
   real(kind=8),allocatable,dimension(:,:),intent(in)::xnorsurf,ynorsurf,znorsurf,xsurf,ysurf,zsurf
!   real(kind=8),allocatable,dimension(:)::xfatem,yfatem,zfatem


   d2max=1e10!16*ltether*ltether

!  pick a solid angle

   call random_number(r)

   jp0=nphisol*r+1


20 call random_number(r)

   jx0=nxsol*r+1

   x0=0.0d0; y0=0.0d0;z0=0.0d0

   ncan=0

   do na=1,nfa

      if(jfasol(1,na)==jx0.and.jfasol(2,na)==jp0)then

         ncan=ncan+1

         x0=x0+xfa(na)
         y0=y0+yfa(na)
         z0=z0+zfa(na)

      end if

   end do

   if(ncan<5) goto 20

   x0=x0/ncan
   y0=y0/ncan
   z0=z0/ncan


!   ncan=0

!   distance=0.0d0

!   do na=1,nfa

!      if(jfasol(2,na)==jp0)then
!         ncan=ncan+1
!         acan(ncan)=na

!         jx=jfasol(1,na)

!         jp=jfasol(2,na)

!         dx=xfa(na)-xsurf(jp,jx)
!         dy=yfa(na)-ysurf(jp,jx)
!         dz=zfa(na)-zsurf(jp,jx)

!         xn=xnorsurf(jp,jx)
!         yn=ynorsurf(jp,jx)
!         zn=znorsurf(jp,jx)

!         distance=distance+dx*xn+dy*yn+dz*zn

!      end if

!      if(ncan==500) exit

!   end do

!   distance=distance/ncan

!  pick an actin bead as a reference:


!   call random_number(r)

!   j=ncan*r+1

!   na=acan(j)


!  then the reference membrane bead:

!   jx0=jfasol(1,na)

!   jp0=jfasol(2,na)


   xn0=xnorsurf(jp0,jx0)
   yn0=ynorsurf(jp0,jx0)
   zn0=znorsurf(jp0,jx0)

   xmb=xsurf(jp0,jx0)
   ymb=ysurf(jp0,jx0)
   zmb=zsurf(jp0,jx0)

!  distance away from membrane:

!   distance=(xfa(na)-xmb)*xn0+(yfa(na)-ymb)*yn0+(zfa(na)-zmb)*zn0


!  position of the first bead:

!   call random_number(r)

!   x0=xfa(na)+2*l_a*(r-0.5d0)
!   y0=yfa(na)
!   z0=zfa(na)

!   x0=xmb+distance*xn0
!   y0=ymb+distance*yn0
!   z0=zmb+distance*zn0


!  pick a length for filament:

   call random_number(r)

   length=falen+(r-0.5)*(falen/2)

   if(length>nmono)then
      length=nmono
   end if


!  total tethering of actin to membrane:

!   nsum=0

!   do j=1,nfa

!      if(a2mem(j)>0) nsum=nsum+1

!   end do


   if(jdir==1)then

      do n=nfa,nfa1+1,-1

         filid(n+length)=filid(n)+1

         apar(1:2,n+length)=apar(1:2,n)

         jfasol(1:2,n+length)=jfasol(1:2,n)

         fadist(n+length)=fadist(n)

         apos(n+length)=apos(n)

!         a2mem(n+length)=a2mem(n)

         xfa(n+length)=xfa(n)
         yfa(n+length)=yfa(n)
         zfa(n+length)=zfa(n)

      end do

      do n=fanum,fanum1+1,-1

         astart(n+1)=astart(n)+length

         alen(n+1)=alen(n)

         fa1stbound(n+1)=fa1stbound(n)
      end do


      na=nfa1

      nnew=fanum1+1

      do nm=1,myp2num

         do jm=1,4

            if(fa2myp2(jm,nm)>nfa1)then

               fa2myp2(jm,nm)=fa2myp2(jm,nm)+length

            end if

         end do

      end do

      do nm=1,myo2num

         do jm=1,2

            if(fa2myo2(jm,nm)>nfa1)then

               fa2myo2(jm,nm)=fa2myo2(jm,nm)+length

            end if

         end do

      end do

      do nl=1,crlknum

         do jl=1,2

            if(fa2lk(jl,nl)>nfa1)then

               fa2lk(jl,nl)=fa2lk(jl,nl)+length

            end if

         end do

      end do


   else

      na=nfa

      nnew=fanum+1

   end if



!--------------------------------


   alen(nnew)=length

   astart(nnew)=na+1

   filid(na+1:na+length)=nnew

   apar(1:2,na+1:na+length)=0
   apar(1,na+1)=-1

   fa1stbound(nnew)=length

   fadist(na+1:na+length)=1

!   a2mem(na+1:na+length)=0



!   if(distance<ltether.and.nsum<ntether)then
!      jtether=1

!      allocate(mark(nmemb))
!      mark=0

!      do j=1,nfa
!         if(a2mem(j)>0) mark(a2mem(j))=1
!      end do

!   else
!      jtether=0
!   end if




   phi=atan(z0/y0)

   if(y0<0.0d0)then
      phi=phi+pi
   end if

   rad=sqrt(y0*y0+z0*z0)

   dphi=l_a/rad*jdir




!  now assign coordinates:

   do j=1,length

      na=na+1

      xfa(na)=x0
      yfa(na)=y0
      zfa(na)=z0

      jfasol(1,na)=jx0

      jfasol(2,na)=jp0

      apos(na)=j

!     assign tethering on the beads:

!      if(jtether==1.and.nsum<ntether)then

!         call random_number(r)


!         if(p_tether>r)then

!            dist2=d2max

!            npick=0

!            do jm=1,nmemb

!               if(mark(jm)==1)then
!                  cycle
!               end if

!               jx=jsursol(1,jm)

!               if(jx<jx0-1.or.jx>jx0+1)then
!                  cycle
!               end if

!               jp=jsursol(2,jm)

!               if(jp<jp0-1.or.jp>jp0+1)then
!                  cycle
!               end if

!               dx=x0-xmemb(jm)

!               dy=y0-ymemb(jm)

!               dz=z0-zmemb(jm)

!               d2=dx*dx+dy*dy+dz*dz

!               if(d2<dist2)then

!                  dist2=d2

!                  npick=jm

!                  jxget=jx

!                  jpget=jp

!               end if

!            end do

!            if(npick==0)then
!               print*,'could not find npick'

!            else

!               a2mem(na)=npick

!               mark(npick)=1

!               nsum=nsum+1

!            end if

!         end if

!      end if

!     update coordinates for the next bead:

      rad=sqrt(y0*y0+z0*z0)

      phi=phi+dphi

      y0=rad*cos(phi)

      z0=rad*sin(phi)

      jxget=0

      dist2=d2max

      do j1=1,3

         jx=jx0+j1-2

         if(jx<1)then
            cycle
         end if

         if(jx>nxsol)then
            cycle
         end if

         do j2=1,3

            jp=jp0+j2-2

            if(jp<1)then
               jp=jp+nphisol
            end if

            if(jp>nphisol)then
               jp=jp-nphisol
            end if

            dx=x0-xsurf(jp,jx)
            dy=y0-ysurf(jp,jx)
            dz=z0-zsurf(jp,jx)

            xn=xnorsurf(jp,jx)
            yn=ynorsurf(jp,jx)
            zn=znorsurf(jp,jx)

            proj=dx*xn+dy*yn+dz*zn

            d2=dx*dx+dy*dy+dz*dz-proj*proj

            if(dist2>d2)then

               dist2=d2

               jxget=jx

               jpget=jp

               proj0=proj

               xn0=xn
               yn0=yn
               zn0=zn

            end if

         end do

      end do

      if(jxget==0)then
         print*,'error in jxget',jx0,jp0
         stop
      end if

      jx0=jxget

      jp0=jpget

!      if(distance-proj0>l_memby2)then

!         x0=x0+xn0*dshift
!         y0=y0+yn0*dshift
!         z0=z0+zn0*dshift

!      elseif(proj0-distance>l_memby2)then

!         x0=x0-xn0*dshift
!         y0=y0-yn0*dshift
!         z0=z0-zn0*dshift

!      end if

!     make sure distance from new bead to membrane not too close:

      if(proj0<l_mem)then

         dshift=l_mem-proj0

         x0=x0+xn0*dshift
         y0=y0+yn0*dshift
         z0=z0+zn0*dshift

      end if


   end do


!--------------------------


!   if(jtether==1)then
!      deallocate(mark)
!   end if

   if(jdir==1)then

      fanum1=nnew

      fanum=fanum+1

      nfa1=na

      nfa=nfa1+nfa2



   else

      fanum=nnew

      fanum2=fanum2+1

      nfa2=nfa2+length

      nfa=na

   end if




   end subroutine

!=========================================================

   subroutine depoly(myp2num,myo2num,crlknum,tension1,tension2,fanum,fanum1,fanum2,nfa,nfa1,nfa2, &
               astart,alen,apos,filid,fa1stbound,fadist,fa2myp2,fa2myo2,fa2lk,myp2typ,myo2typ, &
                jfasol,apar,p_dep,dt,xfa,yfa,zfa)



   implicit none

   integer,value::myp2num,myo2num,crlknum,tension1,tension2
   integer fanum,fanum1,fanum2,nfa,nfa1,nfa2

   integer nfil,newfanum1,nnew,check
   integer j,n,ja,jstop,length,nm,jh,nl,jap,jf,jl,jstart,jstart_jf

   integer,allocatable,dimension(:)::astart,alen,apos,filid,fa1stbound,fadist,map
   integer,allocatable,dimension(:,:)::fa2myp2,fa2myo2,fa2lk,myp2typ,myo2typ,jfasol,apar

   double precision,value::p_dep,dt
   double precision r

   double precision,allocatable,dimension(:)::xfa,yfa,zfa

   check=0

   do n=1,fanum


      length=alen(n)


      jstart=astart(n)

      jstop=jstart+length-1     ! used to be afil(j,n)

      if(tension1==1)then

         if(apar(1,jstop)>0.or.apar(1,jstop-1)>0)then
            cycle
         end if

      end if



      call random_number(r)

      if(p_dep*dt>r)then

         if(length>5)then

            jstart=jstop

            alen(n)=length-1

            if(fa1stbound(n)==length)then
               fa1stbound(n)=length-1
            end if

         else
            alen(n)=0
         end if

         do ja=jstart,jstop

            if(apar(1,ja)>0)then


!              check myosin binding:

               if(apar(2,ja)>0)then

                  jh=apar(1,ja)

                  nm=apar(2,ja)

                  if(nm<=myp2num)then

                     fa2myp2(jh,nm)=0

                     myp2typ(jh,nm)=1

                  else

                     nm=nm-myp2num

                     fa2myo2(jh,nm)=0

                     myo2typ(jh,nm)=1

                  end if

!              check crosslinker binding:

               else

                  jl=apar(1,ja)

                  nl=-apar(2,ja)

                  fa2lk(jl,nl)=0

                  if(jl==1)then

                     jap=fa2lk(2,nl)

                  else

                     jap=fa2lk(1,nl)

                  end if

!                 update 1st bound bead on the partner F-actin:


                  if(jap>0.and.tension2==1)then


                     jf=filid(jap)

                     jstart_jf=astart(jf)

                     if(apos(jap)==fa1stbound(jf))then

                        fa1stbound(jf)=alen(jf)

                        do jl=1,alen(jf)

                           if(apar(2,jstart_jf+jl-1)<0)then

                              nl=-apar(2,jstart_jf+jl-1)

                              if(fa2lk(1,nl)>0.and.fa2lk(2,nl)>0)then

                                 fa1stbound(jf)=jl

                                 exit
                              end if


                           end if

                        end do

                     end if


                  end if


               end if

            end if


!            a2mem(ja)=0

            check=check+1

         end do



      end if


   end do

   if(check==0)then
      return
   end if

!  cleaning up the system:

   allocate(map(nfa))

   nnew=0

   nfil=0

   newfanum1=0

   do n=1,fanum

      if(alen(n)==0)then
         cycle
      end if

      if(n>fanum1.and.newfanum1==0)then
         newfanum1=nfil
      end if

      nfil=nfil+1

      jstart=astart(n)

      astart(nfil)=nnew+1

      alen(nfil)=alen(n)

      fa1stbound(nfil)=fa1stbound(n)

      do j=1,alen(n)

         ja=jstart+j-1

         nnew=nnew+1

         map(ja)=nnew

         apar(1:2,nnew)=apar(1:2,ja)

         apos(nnew)=apos(ja)

!         adir(nnew)=adir(ja)

         filid(nnew)=nfil!filid(ja)

!         a2mem(nnew)=a2mem(ja)

         fadist(nnew)=fadist(ja)

         jfasol(1:2,nnew)=jfasol(1:2,ja)

         xfa(nnew)=xfa(ja)

         yfa(nnew)=yfa(ja)

         zfa(nnew)=zfa(ja)


      end do

   end do

   nfa=nnew

   fanum=nfil

   fanum1=newfanum1

   fanum2=fanum-fanum1


   nfa1=astart(fanum1+1)-1

   nfa2=nfa-nfa1

   do n=1,myp2num

      do j=1,4

         if(fa2myp2(j,n)>0)then
            fa2myp2(j,n)=map(fa2myp2(j,n))
         end if

      end do

   end do

   do n=1,myo2num

      do j=1,2

         if(fa2myo2(j,n)>0)then
            fa2myo2(j,n)=map(fa2myo2(j,n))
         end if

      end do

   end do


   do n=1,crlknum

      do j=1,2

         if(fa2lk(j,n)>0)then
            fa2lk(j,n)=map(fa2lk(j,n))
         end if

      end do

   end do

   deallocate(map)



   end subroutine

!=========================================================

   subroutine myoturnover(nmyoturn,nfa,nmemb,myp2num,myp2len,myo2num,myo2len,nxsol,nphisol,myp2typ,jmyp2sol, &
               myp2body,myp2head,myo2typ,jmyo2sol,myo2body,myo2head,jsursol,jfasol,myp2dist,myo2dist, &
                my2mem,l_mb,l_mem,pi,pmyturn,dt,dphisol,xrmin,xrmax,lmy2mem,xfa,yfa,zfa, &
                 xmyp2,ymyp2,zmyp2,xmyo2,ymyo2,zmyo2,xmemb,ymemb,zmemb, &
                  xnorsurf,ynorsurf,znorsurf,xsurf,ysurf,zsurf)

   implicit none

   integer,value::nfa,nmemb,myp2num,myp2len,myo2num,myo2len,nxsol,nphisol

   integer nmyoturn

   integer na,jdir,jb,jh,nm,jx0,jp0,jm,jx,jp,j,jxget,jpget,j1,j2

   integer,allocatable,dimension(:,:)::myp2typ,jmyp2sol,myo2typ,jmyo2sol
   integer,allocatable,dimension(:,:),intent(in)::myp2body,myp2head,myo2body,myo2head,jsursol,jfasol
   integer,allocatable,dimension(:)::myp2dist,myo2dist,my2mem

   integer,allocatable,dimension(:)::mark

   real(kind=8),value::l_mb,l_mem,pi,pmyturn,dt,dphisol,xrmin,xrmax,lmy2mem
!   real(kind=8),value::xrmin,xrmax,distmin,distmax

   real(kind=8)::dshift,r,xn0,yn0,zn0,distance,x0,y0,z0,phi,rad,dphi
   real(kind=8)::dist2,d2,dx,dy,dz,xn,yn,zn,proj,proj0,xmb,ymb,zmb
   real(kind=8)::dx1,dy1,dz1,dx2,dy2,dz2,dx3,dy3,dz3,dx4,dy4,dz4,cosp,sinp
   real(kind=8)::xdir,ydir,zdir,invdist

   real(kind=8),allocatable,dimension(:)::xmyp2,ymyp2,zmyp2,xmyo2,ymyo2,zmyo2,xfa,yfa,zfa
   real(kind=8),allocatable,dimension(:),intent(in)::xmemb,ymemb,zmemb
   real(kind=8),allocatable,dimension(:,:),intent(in)::xnorsurf,ynorsurf,znorsurf,xsurf,ysurf,zsurf


!   l_memby2=0.5d0*l_mem

!   dshift=0.1d0*l_mb

   allocate(mark(nmemb))
   mark=0
   mark(my2mem(1:myo2num))=1


!$omp parallel &
!$omp default(none) &
!$omp private(na,jdir,jb,jh,nm,jx0,jp0,jm,jx,jp,j,jxget,jpget,j1,j2) &
!$omp private(dshift,r,xn0,yn0,zn0,distance,x0,y0,z0,phi,rad,dphi) &
!$omp private(dist2,d2,dx,dy,dz,xn,yn,zn,proj,proj0,xmb,ymb,zmb) &
!$omp private(dx1,dy1,dz1,dx2,dy2,dz2,dx3,dy3,dz3,dx4,dy4,dz4,cosp,sinp) &
!$omp private(xdir,ydir,zdir,invdist) &
!$omp shared(nfa,nmemb,myp2num,myp2len,myo2num,myo2len,nxsol,nphisol) &
!$omp shared(myp2typ,jmyp2sol,myo2typ,jmyo2sol) &
!$omp shared(myp2body,myp2head,myo2body,myo2head,jsursol,jfasol) &
!$omp shared(myp2dist,myo2dist,my2mem,mark) &
!$omp shared(l_mb,l_mem,pi,pmyturn,dt,dphisol,xrmin,xrmax,lmy2mem) &
!$omp shared(xmyp2,ymyp2,zmyp2,xmyo2,ymyo2,zmyo2,xfa,yfa,zfa) &
!$omp shared(xmemb,ymemb,zmemb,xnorsurf,ynorsurf,znorsurf,xsurf,ysurf,zsurf) &

!$omp reduction(+:nmyoturn)

!  Myp2 turnover:

!$omp do

   do nm=1,myp2num

      if(myp2typ(1,nm)/=2.or.myp2typ(2,nm)/=2.or.myp2typ(3,nm)/=2.or.myp2typ(4,nm)/=2) cycle

      call random_number(r)

      if(pmyturn*dt<r) cycle

      nmyoturn=nmyoturn+1

      myp2typ(1:4,nm)=1

!     preserve heads configuration:

      jb=myp2body(1,nm)

      jh=myp2head(1,nm)

      dx1=xmyp2(jh)-xmyp2(jb)
      dy1=ymyp2(jh)-ymyp2(jb)
      dz1=zmyp2(jh)-zmyp2(jb)

      jh=myp2head(2,nm)

      dx2=xmyp2(jh)-xmyp2(jb)
      dy2=ymyp2(jh)-ymyp2(jb)
      dz2=zmyp2(jh)-zmyp2(jb)

      jb=myp2body(myp2len,nm)

      jh=myp2head(3,nm)

      dx3=xmyp2(jh)-xmyp2(jb)
      dy3=ymyp2(jh)-ymyp2(jb)
      dz3=zmyp2(jh)-zmyp2(jb)

      jh=myp2head(4,nm)

      dx4=xmyp2(jh)-xmyp2(jb)
      dy4=ymyp2(jh)-ymyp2(jb)
      dz4=zmyp2(jh)-zmyp2(jb)


!     pick an actin bead as a reference:

      call random_number(r)

      na=nfa*r+1


!     membrane reference:

!12    call random_number(r)

!      jm=nmemb*r+1

!      if(xmemb(jm)<xrmin.or.xmemb(jm)>xrmax)then
!         goto 12
!      end if


!      jx0=jsursol(1,jm)

!      jp0=jsursol(2,jm)

      jx0=jfasol(1,na)

      jp0=jfasol(2,na)

      xn0=xnorsurf(jp0,jx0)
      yn0=ynorsurf(jp0,jx0)
      zn0=znorsurf(jp0,jx0)

      xmb=xsurf(jp0,jx0)
      ymb=ysurf(jp0,jx0)
      zmb=zsurf(jp0,jx0)

!     a distance away from the membrane:

!      call random_number(r)

!      distance=(xfa(na)-xmb)*xn0+(yfa(na)-ymb)*yn0+(zfa(na)-zmb)*zn0+r

!     angular shift of myosin:

      jm=myp2body(1,nm)

      jp=jmyp2sol(2,jm)

      phi=(jp0-jp)*dphisol

      cosp=cos(phi)

      sinp=sin(phi)

!     update distance signal:

      myp2dist(jm:jm+myp2len+3)=1

!     pick a distance away from the membrane:

!      call random_number(r)

!      distance=(distmax-distmin)*r+distmin

!     position of the first bead:

      call random_number(r)

      x0=xfa(na)+10*(r-0.5d0)
      y0=yfa(na)
      z0=zfa(na)

!      x0=xmb+distance*xn0
!      y0=ymb+distance*yn0
!      z0=zmb+distance*zn0

!     pick "direction"

      call random_number(r)

      if(r>0.5)then
         jdir=1
      else
         jdir=-1
      end if

      phi=atan(z0/y0)

      if(y0<0.0d0)then
         phi=phi+pi
      end if

      rad=sqrt(y0*y0+z0*z0)

      dphi=l_mb/rad*jdir

!     now assign coordinates for the body:

      do j=1,myp2len

         jm=myp2body(j,nm)

         xmyp2(jm)=x0
         ymyp2(jm)=y0
         zmyp2(jm)=z0

         jmyp2sol(1,jm)=jx0

         jmyp2sol(2,jm)=jp0

         if(j==1)then

            jmyp2sol(1,myp2head(1:2,nm))=jx0

            jmyp2sol(2,myp2head(1:2,nm))=jp0

         end if

         if(j==myp2len)then

            jmyp2sol(1,myp2head(3:4,nm))=jx0

            jmyp2sol(2,myp2head(3:4,nm))=jp0

         end if


!        update coordinates for the next bead:

         phi=phi+dphi

         y0=rad*cos(phi)

         z0=rad*sin(phi)

         jxget=0

         dist2=1e10

         do j1=1,3

            jx=jx0+j1-2

            if(jx<1)then
               cycle
            end if

            if(jx>nxsol)then
               cycle
            end if

            do j2=1,3

               jp=jp0+j2-2

               if(jp<1)then
                  jp=jp+nphisol
               end if

               if(jp>nphisol)then
                  jp=jp-nphisol
               end if

               dx=x0-xsurf(jp,jx)
               dy=y0-ysurf(jp,jx)
               dz=z0-zsurf(jp,jx)

               xn=xnorsurf(jp,jx)
               yn=ynorsurf(jp,jx)
               zn=znorsurf(jp,jx)

               proj=dx*xn+dy*yn+dz*zn

               d2=dx*dx+dy*dy+dz*dz-proj*proj


               if(dist2>d2)then

                  dist2=d2

                  jxget=jx

                  jpget=jp

                  proj0=proj

                  xn0=xn
                  yn0=yn
                  zn0=zn

               end if

            end do

         end do

         if(jxget==0)then
            print*,'error in jxget for Myp2',jx0,jp0
            stop
         end if

         jx0=jxget

         jp0=jpget

!         if(distance-proj0>l_memby2)then

!            x0=x0+xn0*dshift
!            y0=y0+yn0*dshift
!            z0=z0+zn0*dshift

!         elseif(proj0-distance>l_memby2)then

!            x0=x0-xn0*dshift
!            y0=y0-yn0*dshift
!            z0=z0-zn0*dshift

!         end if


!        make sure distance from new bead to membrane not too close:

         if(proj0<l_mem)then

            dshift=l_mem-proj0

            x0=x0+xn0*dshift
            y0=y0+yn0*dshift
            z0=z0+zn0*dshift

         end if


      end do

!     coordinates for the heads:

      jb=myp2body(1,nm)

      x0=xmyp2(jb)
      y0=ymyp2(jb)
      z0=zmyp2(jb)

      jh=myp2head(1,nm)

      xmyp2(jh)=x0+dx1

      ymyp2(jh)=y0+dy1*cosp+dz1*sinp

      zmyp2(jh)=z0+dz1*cosp-dy1*sinp

      jh=myp2head(2,nm)

      xmyp2(jh)=x0+dx2

      ymyp2(jh)=y0+dy2*cosp+dz2*sinp

      zmyp2(jh)=z0+dz2*cosp-dy2*sinp

      jb=myp2body(myp2len,nm)

      x0=xmyp2(jb)
      y0=ymyp2(jb)
      z0=zmyp2(jb)

      jh=myp2head(3,nm)

      xmyp2(jh)=x0+dx3

      ymyp2(jh)=y0+dy3*cosp+dz3*sinp

      zmyp2(jh)=z0+dz3*cosp-dy3*sinp

      jh=myp2head(4,nm)

      xmyp2(jh)=x0+dx4

      ymyp2(jh)=y0+dy4*cosp+dz4*sinp

      zmyp2(jh)=z0+dz4*cosp-dy4*sinp

   end do

!$omp enddo nowait


!  Myo2 turnover:

!   allocate(mark(nmemb))
!   mark=0
!   mark(my2mem(1:myo2num))=1

!$omp do schedule(guided,64)

   do nm=1,myo2num

      if(myo2typ(1,nm)/=2.or.myo2typ(2,nm)/=2) cycle

      call random_number(r)

      if(pmyturn*dt<r) cycle

      nmyoturn=nmyoturn+1

      myo2typ(1:2,nm)=1

!     preserve heads configuration:

      jb=myo2body(myo2len,nm)

      jh=myo2head(1,nm)

      dx1=xmyo2(jh)-xmyo2(jb)
      dy1=ymyo2(jh)-ymyo2(jb)
      dz1=zmyo2(jh)-zmyo2(jb)

      jh=myo2head(2,nm)

      dx2=xmyo2(jh)-xmyo2(jb)
      dy2=ymyo2(jh)-ymyo2(jb)
      dz2=zmyo2(jh)-zmyo2(jb)

!     pick a tethering membrane bead:

12    call random_number(r)

      jm=nmemb*r+1

      if(xmemb(jm)<xrmin.or.xmemb(jm)>xrmax) goto 12

!$omp atomic
      mark(jm)=mark(jm)+1

      if(mark(jm)>1) goto 12

      my2mem(nm)=jm


      jx0=jsursol(1,jm)

      jp0=jsursol(2,jm)

      xn0=xnorsurf(jp0,jx0)
      yn0=ynorsurf(jp0,jx0)
      zn0=znorsurf(jp0,jx0)

!      xmb=xsurf(jp0,jx0)
!      ymb=ysurf(jp0,jx0)
!      zmb=zsurf(jp0,jx0)

      xmb=xmemb(jm)
      ymb=ymemb(jm)
      zmb=zmemb(jm)


!     angular shift of myosin:

      jm=myo2body(1,nm)

      jp=jmyo2sol(2,jm)

      phi=(jp0-jp)*dphisol

      cosp=cos(phi)

      sinp=sin(phi)

!     update distance signal:

      myo2dist(jm:jm+myo2len+1)=1


!     position of the first bead:

      call random_number(r)

      distance=lmy2mem*(1.0d0+r)*0.5d0

      x0=xmb+distance*xn0
      y0=ymb+distance*yn0
      z0=zmb+distance*zn0

!     pick "direction"


      call random_number(r)

      xdir=r-0.5d0

      call random_number(r)

      ydir=r-0.5d0

      call random_number(r)

      zdir=r-0.5d0

      invdist=1.0d0/sqrt(ydir*ydir+zdir*zdir)

      rad=0.8d0*sqrt(y0*y0+z0*z0)

      ydir=ydir*rad*invdist
      zdir=zdir*rad*invdist

      ydir=ydir-y0
      zdir=zdir-z0

      invdist=1.0d0/sqrt(ydir*ydir+zdir*zdir)

      ydir=ydir*invdist
      zdir=zdir*invdist

      invdist=l_mb/sqrt(xdir*xdir+ydir*ydir+zdir*zdir)

      xdir=xdir*invdist
      ydir=ydir*invdist
      zdir=zdir*invdist


!     now assign coordinates for the body:

      do j=1,myo2len

         jm=myo2body(j,nm)

         xmyo2(jm)=x0
         ymyo2(jm)=y0
         zmyo2(jm)=z0

         jmyo2sol(1,jm)=jx0

         jmyo2sol(2,jm)=jp0

         if(j==myo2len)then

            jmyo2sol(1,myo2head(1:2,nm))=jx0

            jmyo2sol(2,myo2head(1:2,nm))=jp0

         end if


!        update coordinates for the next bead:

         x0=x0+xdir
         y0=y0+ydir
         z0=z0+zdir

         jxget=0

         dist2=1e10

         do j1=1,3

            jx=jx0+j1-2

            if(jx<1)then
               cycle
            end if

            if(jx>nxsol)then
               cycle
            end if

            do j2=1,3

               jp=jp0+j2-2

               if(jp<1)then
                  jp=jp+nphisol
               end if

               if(jp>nphisol)then
                  jp=jp-nphisol
               end if

               dx=x0-xsurf(jp,jx)
               dy=y0-ysurf(jp,jx)
               dz=z0-zsurf(jp,jx)

               xn=xnorsurf(jp,jx)
               yn=ynorsurf(jp,jx)
               zn=znorsurf(jp,jx)

               proj=dx*xn+dy*yn+dz*zn

               d2=dx*dx+dy*dy+dz*dz-proj*proj

               if(dist2>d2)then

                  dist2=d2

                  jxget=jx

                  jpget=jp

                  proj0=proj

                  xn0=xn
                  yn0=yn
                  zn0=zn

               end if

            end do

         end do

         if(jxget==0)then
            print*,'error in jxget for Myo2',jx0,jp0
            stop
         end if

         jx0=jxget

         jp0=jpget

!        make sure distance from new bead to membrane not too close:

         if(proj0<l_mem)then

            dshift=l_mem-proj0

            x0=x0+xn0*dshift
            y0=y0+yn0*dshift
            z0=z0+zn0*dshift

         end if

      end do

!     coordinates for the heads:

      jb=myo2body(myo2len,nm)

      x0=xmyo2(jb)
      y0=ymyo2(jb)
      z0=zmyo2(jb)

      jh=myo2head(1,nm)

      xmyo2(jh)=x0+dx1

      ymyo2(jh)=y0+dy1*cosp+dz1*sinp

      zmyo2(jh)=z0+dz1*cosp-dy1*sinp

      jh=myo2head(2,nm)

      xmyo2(jh)=x0+dx2

      ymyo2(jh)=y0+dy2*cosp+dz2*sinp

      zmyo2(jh)=z0+dz2*cosp-dy2*sinp

   end do

!$omp enddo nowait

!$omp end parallel


   deallocate(mark)

   end subroutine

!=========================================================

   subroutine xlturnover(crlknum,crlknummax,jxldel,jxlturnover,nfa,crlknumactive,fa2lk,jfasol, &
               lkstart,lktyp,lkdist,jlksol,plk_remove,plk_turn,dt,dphisol, &
                xnorsurf,ynorsurf,znorsurf,xsurf,ysurf,zsurf,xfa,yfa,zfa,xlk,ylk,zlk)


   implicit none
   integer,value::crlknum,crlknummax,jxldel,jxlturnover,nfa
   integer crlknumactive

   integer nl,n1,n2,n3,na,jx0,jp0,jp

   integer,allocatable,dimension(:,:),intent(in)::fa2lk,jfasol
   integer,allocatable,dimension(:),intent(in)::lkstart
   integer,allocatable,dimension(:)::lktyp,lkdist
   integer,allocatable,dimension(:,:)::jlksol

   real(kind=8),value::plk_remove,plk_turn,dt,dphisol
   real(kind=8)::r,prob1,prob2,xori,yori,zori,dy,dz,xn0,yn0,zn0,xmb,ymb,zmb,distance,phi,cosp,sinp

   real(kind=8),allocatable,dimension(:,:),intent(in)::xnorsurf,ynorsurf,znorsurf,xsurf,ysurf,zsurf
   real(kind=8),allocatable,dimension(:),intent(in)::xfa,yfa,zfa
   real(kind=8),allocatable,dimension(:)::xlk,ylk,zlk

   if(crlknumactive<=crlknummax) jxldel=0

   prob1=plk_remove*dt

   prob2=plk_turn*dt

!   inv2l_lk=0.5d0/l_lk

!$omp parallel &
!$omp default(none) &
!$omp private(nl,n1,n2,n3,na,jx0,jp0,jp,r,xori,yori,zori,dy,dz) &
!$omp private(xn0,yn0,zn0,xmb,ymb,zmb,distance,phi,cosp,sinp) &
!$omp shared(crlknum,jxldel,jxlturnover,nfa,fa2lk,jfasol,lkstart,lktyp,lkdist,jlksol) &
!$omp shared(prob1,prob2,dphisol,xnorsurf,ynorsurf,znorsurf,xsurf,ysurf,zsurf) &
!$omp shared(xfa,yfa,zfa,xlk,ylk,zlk) &
!$omp reduction(+:crlknumactive)
!$omp do

   do nl=1,crlknum

      if(jxldel==1)then

         if(fa2lk(1,nl)==0.and.fa2lk(2,nl)==0.and.lktyp(nl)==1)then

            call random_number(r)

            if(prob1>r)then

               lktyp(nl)=0

               crlknumactive=crlknumactive-1

            end if

         end if

      end if

      if(jxlturnover==1)then

         if(fa2lk(1,nl)==0.and.fa2lk(2,nl)==0.and.lktyp(nl)==1)then

            call random_number(r)

            if(prob2<r) cycle


!           orientation of crosslinker:

            n1=lkstart(nl)

            n2=n1+1

            n3=n2+1

            xori=xlk(2)-xlk(1)
            dy=ylk(2)-ylk(1)
            dz=zlk(2)-zlk(1)

!           pick an actin bead as a reference:

            call random_number(r)

            na=nfa*r+1

            jx0=jfasol(1,na)

            jp0=jfasol(2,na)

            xn0=xnorsurf(jp0,jx0)
            yn0=ynorsurf(jp0,jx0)
            zn0=znorsurf(jp0,jx0)

            xmb=xsurf(jp0,jx0)
            ymb=ysurf(jp0,jx0)
            zmb=zsurf(jp0,jx0)

!           a distance away from the membrane:

            call random_number(r)

            distance=(xfa(na)-xmb)*xn0+(yfa(na)-ymb)*yn0+(zfa(na)-zmb)*zn0+r

!           rotation of crosslinker:

            jp=jlksol(2,n2)

            phi=(jp0-jp)*dphisol

            cosp=cos(phi)

            sinp=sin(phi)

            yori=dy*cosp+dz*sinp

            zori=dz*cosp-dy*sinp

!           update distance signal:

            lkdist(n1:n3)=1

!           position of the middle bead:

            xlk(n2)=xmb+distance*xn0
            ylk(n2)=ymb+distance*yn0
            zlk(n2)=zmb+distance*zn0

!           position of the first bead:

            xlk(n1)=xlk(n2)-xori
            ylk(n1)=ylk(n2)-yori
            zlk(n1)=zlk(n2)-zori

!           position of the third bead:

            xlk(n3)=xlk(n2)+xori
            ylk(n3)=ylk(n2)+yori
            zlk(n3)=zlk(n2)+zori

            jlksol(1,n1:n3)=jx0

            jlksol(2,n1:n3)=jp0

         end if

      end if


   end do

!$omp end do nowait
!$omp end parallel


   end subroutine

!=========================================================

   subroutine breakfa(fanum,fanum1,fanum2,astart,alen,apar,fa1stbound,filid,apos,fa2lk,xfa,yfa,zfa,l_a)

   implicit none

   integer fanum,fanum1,fanum2,fanumtem,fanum1tem,lengthtem
   integer n,length,jstart,ja,j,jcheck,jbreak,nl
   integer,allocatable,dimension(:)::astart,alen,j_cof,fa1stbound,filid,apos
   integer,allocatable,dimension(:)::astarttem,alentem,fa1stboundtem
   integer,allocatable,dimension(:,:)::apar
   integer,allocatable,dimension(:,:),intent(in)::fa2lk
   double precision, value::l_a
   double precision half,dist0,dist,dx,dy,dz,dcheck,dchange,prob,r
   double precision,allocatable,dimension(:)::xfa,yfa,zfa


   allocate(j_cof(fanum))

   j_cof=0

   half=0.5d0*l_a

!$omp parallel &
!$omp default(none) &
!$omp private(n,length,jstart,ja,j,jcheck,dx,dy,dz,dist0,dist,dcheck,dchange,prob,r) &
!$omp shared(fanum,astart,alen,apar,j_cof,l_a,half,xfa,yfa,zfa)
!$omp do schedule(guided,32)

   do n=1,fanum

      length=alen(n)-20

      if(length<25)then
         cycle
      end if

      jstart=astart(n)

      ja=jstart+19

      dx=xfa(ja)-xfa(jstart)
      dy=yfa(ja)-yfa(jstart)
      dz=zfa(ja)-zfa(jstart)

      dist0=sqrt(dx*dx+dy*dy+dz*dz)

      dcheck=l_a

      do j=20,length

         ja=jstart+j

         dx=xfa(ja)-xfa(jstart)
         dy=yfa(ja)-yfa(jstart)
         dz=zfa(ja)-zfa(jstart)

         dist=sqrt(dx*dx+dy*dy+dz*dz)

         dchange=dist-dist0

         dist0=dist


         if(dchange<dcheck.and.apar(1,ja)==0.and.apar(2,ja)==0)then
            dcheck=dchange
            jcheck=j
         end if


      end do

      if(dcheck>half)then
         cycle
      end if

      prob=1.0d0-dcheck/l_a

      call random_number(r)

      if(prob>r)then
         j_cof(n)=jcheck
      end if

   end do

!$omp end do nowait
!$omp end parallel

   jbreak=sum(j_cof(1:fanum))

   if(jbreak==0)then
     return
   end if

!-------------------------

!  updating up the system:

   allocate(fa1stboundtem(fanum+jbreak),alentem(fanum+jbreak),astarttem(fanum+jbreak))

   fanumtem=0

   do n=1,fanum

      fanumtem=fanumtem+1

      jstart=astart(n)

      length=alen(n)

      astarttem(fanumtem)=jstart

      filid(jstart:jstart+length-1)=fanumtem

      if(j_cof(n)==0)then

         fa1stboundtem(fanumtem)=fa1stbound(n)

         alentem(fanumtem)=alen(n)

      else

         lengthtem=j_cof(n)

         alentem(fanumtem)=lengthtem

         if(fa1stbound(n)>lengthtem)then
            fa1stboundtem(fanumtem)=lengthtem
         else
            fa1stboundtem(fanumtem)=fa1stbound(n)
         end if

         fanumtem=fanumtem+1

         astarttem(fanumtem)=jstart+lengthtem

         apar(1,jstart+lengthtem)=-1

         apar(2,jstart+lengthtem)=0

         filid(jstart+lengthtem:jstart+length-1)=fanumtem

         apos(jstart+lengthtem:jstart+length-1)=apos(jstart+lengthtem:jstart+length-1)-lengthtem

         alentem(fanumtem)=length-lengthtem

         if(fa1stbound(n)>lengthtem)then
            fa1stboundtem(fanumtem)=fa1stbound(n)-lengthtem
         else

            fa1stboundtem(fanumtem)=length-lengthtem

            do j=1,length-lengthtem

               if(apar(2,jstart+lengthtem+j-1)<0)then

                  nl=-apar(2,jstart+lengthtem+j-1)

                  if(fa2lk(1,nl)>0.and.fa2lk(2,nl)>0)then
                     fa1stboundtem(fanumtem)=j
                     exit
                  end if
 
               end if

            end do

         end if

      end if

      if(n==fanum1)then
         fanum1tem=fanumtem
      end if

   end do

   fanum=fanumtem

   fanum1=fanum1tem

   fanum2=fanum-fanum1


   fa1stbound(1:fanum)=fa1stboundtem(1:fanum)

   alen(1:fanum)=alentem(1:fanum)

   astart(1:fanum)=astarttem(1:fanum)

   deallocate(j_cof,fa1stboundtem,alentem,astarttem)

   end subroutine

!=========================================================
   subroutine memb_crowd(rcl_off,nmemb,xmemb,ymemb,zmemb,ncrowd)
   implicit none

   integer, value::nmemb
   integer ncrowd,jmb,jo
   integer,allocatable,dimension(:)::part
   double precision, value::rcl_off
   double precision d2max_cl,d2,dx,dy,dz
   double precision,allocatable,dimension(:)::xmemb,ymemb,zmemb

   allocate(part(nmemb))

   part=0


   d2max_cl=2*rcl_off*rcl_off

   do jmb=1,nmemb-1

      do jo=jmb+1,nmemb

         dx=xmemb(jmb)-xmemb(jo)
         dy=ymemb(jmb)-ymemb(jo)
         dz=zmemb(jmb)-zmemb(jo)

         d2=dx*dx+dy*dy+dz*dz

         if(d2<d2max_cl)then
            part(jmb)=part(jmb)+1
            part(jo)=part(jo)+1
         end if

      end do

   end do

   ncrowd=maxval(part(1:nmemb))

   deallocate(part)

   end subroutine

!=========================================================

   subroutine allpairs(nmemb,xmemb,ymemb,zmemb,neinum,myp2nei,myo2nei,lknei,nfa,xfa,yfa,zfa, &
               nmyp2,xmyp2,ymyp2,zmyp2,nmyo2,xmyo2,ymyo2,zmyo2,nlk,xlk,ylk,zlk,fanum,astart,alen, &
                myp2num,myp2head,myp2body,myp2len,myo2num,myo2head,myo2body,myo2len,crlknum,lklen,lkstart, &
                 npair_myp2ac,pair_myp2ac,npair_lkac,pair_lkac,npair_myp2lk,pair_myp2lk,npair_ac2,pair_ac2, &
                  npair_myp2,pair_myp2,npair_lk2,pair_lk2,npair_mb2,pair_mb2,r_off,cos_t_2, &
                    npair_myo2ac,pair_myo2ac,npair_myo2lk,pair_myo2lk,npair_myo2,pair_myo2,npair_mymy,pair_mymy)

   implicit none

   integer,value::nmemb,nfa,nmyp2,nmyo2,nlk,fanum,myp2num,myo2num,crlknum,lklen,myp2len,myo2len,neinum
   integer npair_myp2ac,npair_lkac,npair_myp2lk,npair_ac2,npair_myp2,npair_lk2
   integer npair_mb2,npair_myo2ac,npair_myo2lk,npair_myo2,npair_mymy

   integer ja,ja1,ja2,jm,jm1,jm2,jl,jl1,jl2,nf1,nf2,jf1,jf2,nm1,nm2,jb1,jb2,jnei,jp
   integer jh1,jh2,nl1,nl2,jc1,jc2,jmb,jo,ji,j,n,jstart1,jstart2

   integer,allocatable,intent(in),dimension(:)::alen,astart,lkstart
   integer,allocatable,intent(in),dimension(:,:)::myp2body,myp2head,myo2body,myo2head
   integer,allocatable,dimension(:,:)::pair_myp2ac,pair_lkac,pair_myp2lk,pair_ac2,pair_myp2,pair_lk2
   integer,allocatable,dimension(:,:)::pair_mb2,pair_myo2ac,pair_myo2lk,pair_myo2,pair_mymy
   integer,allocatable,dimension(:)::pairtyp
   integer,allocatable,dimension(:,:,:)::myp2nei,myo2nei,lknei

   double precision,value::r_off,cos_t_2
   double precision dx,dx2,rad1_2,rad2_2,dy,dz,d2,d2max,dmax,dist,cos_2

   double precision,allocatable,intent(in),dimension(:)::xfa,yfa,zfa,xmyp2,ymyp2,zmyp2,xlk,ylk,zlk
   double precision,allocatable,intent(in),dimension(:)::xmemb,ymemb,zmemb,xmyo2,ymyo2,zmyo2
   double precision,allocatable,dimension(:)::pairlen
   integer,allocatable,dimension(:)::pairnum
   integer,allocatable,dimension(:,:)::pairtem
   integer nmax,nnei





   d2max=2*r_off

   nmax=max(nfa,nmyp2,nmyo2,nlk,nmemb)

   nnei=500

   allocate(pairtem(nnei,nmax),pairnum(nmax))


!  actin-Myp2 pairs:

   myp2nei=0

!$omp parallel &
!$omp default(none) &
!$omp private(ja,n,jm,nm1,jb1,jh1,jnei,dx,dy,dz,d2,jl) &
!$omp shared(nfa,myp2num,myp2len,myp2body,myp2head,xfa,yfa,zfa,d2max,pairtem,pairnum) &
!$omp shared(xmyp2,ymyp2,zmyp2,neinum,myp2nei,nnei)
!$omp do schedule(guided,64)



   do nm1=1,myp2num

      do jb1=1,myp2len

         jm=myp2body(jb1,nm1)

         n=0

         do ja=1,nfa

            dx=xfa(ja)-xmyp2(jm)
            dy=yfa(ja)-ymyp2(jm)
            dz=zfa(ja)-zmyp2(jm)

            d2=abs(dx)+abs(dy)+abs(dz)!dx*dx+dy*dy+dz*dz

            if(d2<d2max)then

               if(n==nnei)then
                  exit
               end if

               n=n+1
               pairtem(n,jm)=ja
            end if

         end do

         pairnum(jm)=n

      end do

      do jh1=1,4

         jm=myp2head(jh1,nm1)

         n=0

         jnei=0

         do ja=1,nfa

            dx=xfa(ja)-xmyp2(jm)
            dy=yfa(ja)-ymyp2(jm)
            dz=zfa(ja)-zmyp2(jm)

            d2=abs(dx)+abs(dy)+abs(dz)!dx*dx+dy*dy+dz*dz

            if(d2<d2max)then

               if(n==nnei)then
                  exit
               end if

               n=n+1
               pairtem(n,jm)=ja


               jnei=jnei+1

               if(jnei>neinum)then
                  print*,'too crowded at head',jm,nm1
               else
                  myp2nei(jnei,jh1,nm1)=ja
               end if

            end if

         end do

         pairnum(jm)=n

      end do

   end do

!$omp end do nowait
!$omp end parallel

   npair_myp2ac=0

   do jm=1,nmyp2

      if(pairnum(jm)==0)then
         cycle
      end if

      do j=1,pairnum(jm)

         npair_myp2ac=npair_myp2ac+1

         pair_myp2ac(1,npair_myp2ac)=jm
         pair_myp2ac(2,npair_myp2ac)=pairtem(j,jm)

      end do

   end do

!print*,'actin-myosin',npair_myp2ac

!-------------------------------------


!  actin-Myo2 pairs:

   myo2nei=0

!$omp parallel &
!$omp default(none) &
!$omp private(ja,n,jm,nm1,jb1,jh1,jnei,dx,dy,dz,d2,jl) &
!$omp shared(nfa,myo2num,myo2len,myo2body,myo2head,xfa,yfa,zfa,d2max,pairtem,pairnum) &
!$omp shared(xmyo2,ymyo2,zmyo2,neinum,myo2nei,nnei)
!$omp do schedule(guided,64)



   do nm1=1,myo2num

      do jb1=1,myo2len

         jm=myo2body(jb1,nm1)

         n=0

         do ja=1,nfa

            dx=xfa(ja)-xmyo2(jm)
            dy=yfa(ja)-ymyo2(jm)
            dz=zfa(ja)-zmyo2(jm)

            d2=abs(dx)+abs(dy)+abs(dz)!dx*dx+dy*dy+dz*dz

            if(d2<d2max)then

               if(n==nnei)then
                  exit
               end if

               n=n+1
               pairtem(n,jm)=ja
            end if

         end do

         pairnum(jm)=n

      end do

      do jh1=1,2

         jm=myo2head(jh1,nm1)

         n=0

         jnei=0

         do ja=1,nfa

            dx=xfa(ja)-xmyo2(jm)
            dy=yfa(ja)-ymyo2(jm)
            dz=zfa(ja)-zmyo2(jm)

            d2=abs(dx)+abs(dy)+abs(dz)!dx*dx+dy*dy+dz*dz

            if(d2<d2max)then

               if(n==nnei)then
                  exit
               end if

               n=n+1
               pairtem(n,jm)=ja


               jnei=jnei+1

               if(jnei>neinum)then
                  print*,'too crowded at head',jm,nm1
               else
                  myo2nei(jnei,jh1,nm1)=ja
               end if

            end if

         end do

         pairnum(jm)=n

      end do

   end do

!$omp end do nowait
!$omp end parallel

   npair_myo2ac=0

   do jm=1,nmyo2

      if(pairnum(jm)==0)then
         cycle
      end if

      do j=1,pairnum(jm)

         npair_myo2ac=npair_myo2ac+1

         pair_myo2ac(1,npair_myo2ac)=jm
         pair_myo2ac(2,npair_myo2ac)=pairtem(j,jm)

      end do

   end do

!-------------------------------------


!     Myp2-crosslinker pairs:



!$omp parallel &
!$omp default(none) &
!$omp private(n,jm,dx,dy,dz,d2,jl) &
!$omp shared(nmyp2,d2max,pairtem,nlk,pairnum) &
!$omp shared(xmyp2,ymyp2,zmyp2,xlk,ylk,zlk,nnei)
!$omp do schedule(guided,64)

   do jm=1,nmyp2



      n=0

      do jl=1,nlk

         dx=xmyp2(jm)-xlk(jl)
         dy=ymyp2(jm)-ylk(jl)
         dz=zmyp2(jm)-zlk(jl)

!         d2=dx*dx+dy*dy+dz*dz
         d2=abs(dx)+abs(dy)+abs(dz)!dx*dx+dy*dy+dz*dz

         if(d2<d2max)then

            if(n==nnei)then
               exit
            end if

            n=n+1
            pairtem(n,jm)=jl
         end if

      end do

      pairnum(jm)=n


   end do

!$omp end do nowait
!$omp end parallel


   npair_myp2lk=0


   do jm=1,nmyp2

      if(pairnum(jm)==0)then
         cycle
      end if

      do j=1,pairnum(jm)!nnei


         npair_myp2lk=npair_myp2lk+1

         pair_myp2lk(1,npair_myp2lk)=jm
         pair_myp2lk(2,npair_myp2lk)=pairtem(j,jm)

      end do

   end do


!print*,'myosin-crlk',npair_myp2lk

!-------------------------------------------------------

!   Myo2-crosslinker pairs:

!$omp parallel &
!$omp default(none) &
!$omp private(n,jm,dx,dy,dz,d2,jl) &
!$omp shared(nmyo2,d2max,pairtem,nlk,pairnum) &
!$omp shared(xmyo2,ymyo2,zmyo2,xlk,ylk,zlk,nnei)
!$omp do schedule(guided,64)

   do jm=1,nmyo2



      n=0

      do jl=1,nlk

         dx=xmyo2(jm)-xlk(jl)
         dy=ymyo2(jm)-ylk(jl)
         dz=zmyo2(jm)-zlk(jl)

         d2=abs(dx)+abs(dy)+abs(dz)!dx*dx+dy*dy+dz*dz

         if(d2<d2max)then

            if(n==nnei)then
               exit
            end if

            n=n+1
            pairtem(n,jm)=jl
         end if

      end do

      pairnum(jm)=n


   end do

!$omp end do nowait
!$omp end parallel

   npair_myo2lk=0


   do jm=1,nmyo2

      if(pairnum(jm)==0)then
         cycle
      end if

      do j=1,pairnum(jm)!nnei


         npair_myo2lk=npair_myo2lk+1

         pair_myo2lk(1,npair_myo2lk)=jm
         pair_myo2lk(2,npair_myo2lk)=pairtem(j,jm)

      end do

   end do

!-------------------------------------------------------


!   actin-crosslinker pairs:

   lknei=0

!$omp parallel &
!$omp default(none) &
!$omp private(n,ja,dx,dy,dz,d2,jl,nl1,jc1,jnei) &
!$omp shared(nfa,d2max,pairtem,lkstart,crlknum,lklen,pairnum) &
!$omp shared(xfa,yfa,zfa,xlk,ylk,zlk,lknei,neinum,nnei)
!$omp do schedule(guided,64)

   do nl1=1,crlknum

      do jc1=1,lklen

         jl=lkstart(nl1)+jc1-1   !crlk(jc1,nl1)


         n=0

         jnei=0

         do ja=1,nfa

            dx=xfa(ja)-xlk(jl)
            dy=yfa(ja)-ylk(jl)
            dz=zfa(ja)-zlk(jl)

!         d2=dx*dx+dy*dy+dz*dz
            d2=abs(dx)+abs(dy)+abs(dz)!dx*dx+dy*dy+dz*dz

            if(d2<d2max)then

               if(n==nnei)then
                  exit
               end if

               n=n+1
               pairtem(n,jl)=ja

               if(jc1==1)then

                  jnei=jnei+1

                  if(jnei>neinum)then
                     print*,'too crowded at crlk',jl
                  else
                     lknei(jnei,1,nl1)=ja
                  end if

               end if

               if(jc1==lklen)then

                  jnei=jnei+1

                  if(jnei>neinum)then
                     print*,'too crowded at crlk',jl
                  else
                     lknei(jnei,2,nl1)=ja
                  end if

               end if

            end if

         end do

         pairnum(jl)=n

      end do


   end do

!$omp end do nowait
!$omp end parallel

   npair_lkac=0


   do jl=1,nlk

      if(pairnum(jl)==0)then
         cycle
      end if

      do j=1,pairnum(jl)!nnei

         npair_lkac=npair_lkac+1

         pair_lkac(1,npair_lkac)=jl

         pair_lkac(2,npair_lkac)=pairtem(j,jl)

      end do

   end do


!print*,'actin-crlk',npair_lkac

!write(*,*)

!----------------------------------------------
!  actin-actin pairs:



!$omp parallel &
!$omp default(none) &
!$omp private(nf1,jf1,ja1,n,nf2,ja2,dx,dy,dz,d2,jstart1,jstart2) &
!$omp shared(fanum,d2max,pairtem,alen,astart,pairnum) &
!$omp shared(xfa,yfa,zfa,nnei)
!$omp do schedule(guided,32)

   do nf1=1,fanum-1

      jstart1=astart(nf1)

      do jf1=1,alen(nf1)

         ja1=jstart1+jf1-1     ! used to be afil(jf1,nf1)

         n=0


         do nf2=nf1+1,fanum

            jstart2=astart(nf2)

            do jf2=1,alen(nf2)

               ja2=jstart2+jf2-1     ! used to be afil(jf2,nf2)

               dx=xfa(ja1)-xfa(ja2)
               dy=yfa(ja1)-yfa(ja2)
               dz=zfa(ja1)-zfa(ja2)

!               d2=dx*dx+dy*dy+dz*dz
               d2=abs(dx)+abs(dy)+abs(dz)!dx*dx+dy*dy+dz*dz

               if(d2<d2max)then

                  if(n==nnei)then
                     exit
                  end if

                  n=n+1
                  pairtem(n,ja1)=ja2

               end if

            end do


         end do

         pairnum(ja1)=n

      end do

   end do

!$omp end do nowait
!$omp end parallel

   npair_ac2=0


   do ja=1,nfa-alen(fanum) ! the last filament is not in pairs

      if(pairnum(ja)==0)then
         cycle
      end if

      do j=1,pairnum(ja)!nnei

         npair_ac2=npair_ac2+1

         pair_ac2(1,npair_ac2)=ja

         pair_ac2(2,npair_ac2)=pairtem(j,ja)

      end do

   end do

!omp end do nowait
!omp end parallel

!print*,'actin-actin',npair_ac2

!write(*,*)

!--------------------------------------------------
!  Myp2-Myp2 pairs:



!$omp parallel &
!$omp default(none) &
!$omp private(nm1,jb1,jm1,n,nm2,jb2,jm2,dx,dy,dz,d2,jh1,jh2) &
!$omp shared(myp2num,d2max,pairtem,myp2body,myp2head,myp2len,pairnum) &
!$omp shared(xmyp2,ymyp2,zmyp2,nnei)
!$omp do schedule(guided,32)

   do nm1=1,myp2num-1

      do jb1=1,myp2len

         jm1=myp2body(jb1,nm1)

         n=0


         do nm2=nm1+1,myp2num

!           body-body:

            do jb2=1,myp2len

               jm2=myp2body(jb2,nm2)

               dx=xmyp2(jm1)-xmyp2(jm2)
               dy=ymyp2(jm1)-ymyp2(jm2)
               dz=zmyp2(jm1)-zmyp2(jm2)

!               d2=dx*dx+dy*dy+dz*dz
               d2=abs(dx)+abs(dy)+abs(dz)!dx*dx+dy*dy+dz*dz

               if(d2<d2max)then

                  if(n==nnei)then
                     exit
                  end if

                  n=n+1
                  pairtem(n,jm1)=jm2

               end if

            end do

!           body-head:

            do jh2=1,4

               jm2=myp2head(jh2,nm2)

               dx=xmyp2(jm1)-xmyp2(jm2)
               dy=ymyp2(jm1)-ymyp2(jm2)
               dz=zmyp2(jm1)-zmyp2(jm2)

!               d2=dx*dx+dy*dy+dz*dz
               d2=abs(dx)+abs(dy)+abs(dz)!dx*dx+dy*dy+dz*dz

               if(d2<d2max)then

                  if(n==nnei)then
                     exit
                  end if

                  n=n+1
                  pairtem(n,jm1)=jm2

               end if

            end do

         end do

         pairnum(jm1)=n

      end do

      do jh1=1,4

         jm1=myp2head(jh1,nm1)


         n=0

         do nm2=nm1+1,myp2num

!           head-body:

            do jb2=1,myp2len

               jm2=myp2body(jb2,nm2)

               dx=xmyp2(jm1)-xmyp2(jm2)
               dy=ymyp2(jm1)-ymyp2(jm2)
               dz=zmyp2(jm1)-zmyp2(jm2)

!               d2=dx*dx+dy*dy+dz*dz
               d2=abs(dx)+abs(dy)+abs(dz)!dx*dx+dy*dy+dz*dz

               if(d2<d2max)then

                  if(n==nnei)then
                     exit
                  end if

                  n=n+1
                  pairtem(n,jm1)=jm2

               end if

            end do

!           head-head:

            do jh2=1,4

               jm2=myp2head(jh2,nm2)

               dx=xmyp2(jm1)-xmyp2(jm2)
               dy=ymyp2(jm1)-ymyp2(jm2)
               dz=zmyp2(jm1)-zmyp2(jm2)

!               d2=dx*dx+dy*dy+dz*dz
               d2=abs(dx)+abs(dy)+abs(dz)!dx*dx+dy*dy+dz*dz

               if(d2<d2max)then

                  if(n==nnei)then
                     exit
                  end if

                  n=n+1
                  pairtem(n,jm1)=jm2

               end if

            end do

         end do

         pairnum(jm1)=n

      end do

   end do

!$omp end do nowait
!$omp end parallel


   npair_myp2=0


   do jm=1,nmyp2-4-myp2len ! the last myosin beads don't have pairs

      if(pairnum(jm)==0)then
         cycle
      end if


      do j=1,pairnum(jm)!nnei

         npair_myp2=npair_myp2+1

         pair_myp2(1,npair_myp2)=jm

         pair_myp2(2,npair_myp2)=pairtem(j,jm)


      end do

   end do


!print*,'myo-myo',npair_myp2

!--------------------------------------------------


!  Myo2-Myo2 pairs:


!$omp parallel &
!$omp default(none) &
!$omp private(nm1,jb1,jm1,n,nm2,jb2,jm2,dx,dy,dz,d2,jh1,jh2) &
!$omp shared(myo2num,d2max,pairtem,myo2body,myo2head,myo2len,pairnum) &
!$omp shared(xmyo2,ymyo2,zmyo2,nnei)
!$omp do schedule(guided,32)

   do nm1=1,myo2num-1

      do jb1=1,myo2len

         jm1=myo2body(jb1,nm1)

         n=0


         do nm2=nm1+1,myo2num

!           body-body:

            do jb2=1,myo2len

               jm2=myo2body(jb2,nm2)

               dx=xmyo2(jm1)-xmyo2(jm2)
               dy=ymyo2(jm1)-ymyo2(jm2)
               dz=zmyo2(jm1)-zmyo2(jm2)

               d2=abs(dx)+abs(dy)+abs(dz)!dx*dx+dy*dy+dz*dz

               if(d2<d2max)then

                  if(n==nnei)then
                     exit
                  end if

                  n=n+1
                  pairtem(n,jm1)=jm2

               end if

            end do

!           body-head:

            do jh2=1,2

               jm2=myo2head(jh2,nm2)

               dx=xmyo2(jm1)-xmyo2(jm2)
               dy=ymyo2(jm1)-ymyo2(jm2)
               dz=zmyo2(jm1)-zmyo2(jm2)

               d2=abs(dx)+abs(dy)+abs(dz)!dx*dx+dy*dy+dz*dz

               if(d2<d2max)then

                  if(n==nnei)then
                     exit
                  end if

                  n=n+1
                  pairtem(n,jm1)=jm2

               end if

            end do

         end do

         pairnum(jm1)=n

      end do

      do jh1=1,2

         jm1=myo2head(jh1,nm1)


         n=0

         do nm2=nm1+1,myo2num

!           head-body:

            do jb2=1,myo2len

               jm2=myo2body(jb2,nm2)

               dx=xmyo2(jm1)-xmyo2(jm2)
               dy=ymyo2(jm1)-ymyo2(jm2)
               dz=zmyo2(jm1)-zmyo2(jm2)

               d2=abs(dx)+abs(dy)+abs(dz)!dx*dx+dy*dy+dz*dz

               if(d2<d2max)then

                  if(n==nnei)then
                     exit
                  end if

                  n=n+1
                  pairtem(n,jm1)=jm2

               end if

            end do

!           head-head:

            do jh2=1,2

               jm2=myo2head(jh2,nm2)

               dx=xmyo2(jm1)-xmyo2(jm2)
               dy=ymyo2(jm1)-ymyo2(jm2)
               dz=zmyo2(jm1)-zmyo2(jm2)

               d2=abs(dx)+abs(dy)+abs(dz)!dx*dx+dy*dy+dz*dz

               if(d2<d2max)then

                  if(n==nnei)then
                     exit
                  end if

                  n=n+1
                  pairtem(n,jm1)=jm2

               end if

            end do

         end do

         pairnum(jm1)=n

      end do

   end do

!$omp end do nowait
!$omp end parallel

   npair_myo2=0


   do jm=1,nmyo2-2-myo2len ! the last myosin beads don't have pairs

      if(pairnum(jm)==0)then
         cycle
      end if


      do j=1,pairnum(jm)!nnei

         npair_myo2=npair_myo2+1

         pair_myo2(1,npair_myo2)=jm

         pair_myo2(2,npair_myo2)=pairtem(j,jm)


      end do

   end do

!--------------------------------------------------


!  Myp2-Myo2 pairs:


!$omp parallel &
!$omp default(none) &
!$omp private(n,jm,dx,dy,dz,d2,jp) &
!$omp shared(nmyo2,nmyp2,d2max,pairtem,pairnum) &
!$omp shared(xmyo2,ymyo2,zmyo2,xmyp2,ymyp2,zmyp2,nnei)
!$omp do schedule(guided,64)


   do jm=1,nmyo2

      n=0

      do jp=1,nmyp2

         dx=xmyo2(jm)-xmyp2(jp)
         dy=ymyo2(jm)-ymyp2(jp)
         dz=zmyo2(jm)-zmyp2(jp)

         d2=abs(dx)+abs(dy)+abs(dz)

         if(d2<d2max)then

            if(n==nnei) exit

            n=n+1
            pairtem(n,jm)=jp

         end if

      end do

      pairnum(jm)=n

   end do

!$omp end do nowait
!$omp end parallel


   npair_mymy=0

   do jm=1,nmyo2

      if(pairnum(jm)==0)then
         cycle
      end if

      do j=1,pairnum(jm)!nnei

         npair_mymy=npair_mymy+1

         pair_mymy(1,npair_mymy)=jm

         pair_mymy(2,npair_mymy)=pairtem(j,jm)

      end do

   end do

!-------------------------------------------------------


!  crosslinker-crosslinker pairs:


!$omp parallel &
!$omp default(none) &
!$omp private(nl1,jc1,jl1,n,nl2,jc2,jl2,dx,dy,dz,d2) &
!$omp shared(crlknum,d2max,pairtem,lkstart,lklen,pairnum) &
!$omp shared(xlk,ylk,zlk,nnei)
!$omp do schedule(guided,32)

   do nl1=1,crlknum-1

      do jc1=1,lklen

         jl1=lkstart(nl1)+jc1-1   !crlk(jc1,nl1)

         n=0

         do nl2=nl1+1,crlknum

            do jc2=1,lklen

               jl2=lkstart(nl2)+jc2-1   !crlk(jc2,nl2)

               dx=xlk(jl1)-xlk(jl2)
               dy=ylk(jl1)-ylk(jl2)
               dz=zlk(jl1)-zlk(jl2)

!               d2=dx*dx+dy*dy+dz*dz
               d2=abs(dx)+abs(dy)+abs(dz)!dx*dx+dy*dy+dz*dz

               if(d2<d2max)then

                  if(n==nnei)then
                     exit
                  end if

                  n=n+1
                  pairtem(n,jl1)=jl2

               end if

            end do

         end do

         pairnum(jl1)=n

      end do

   end do

!$omp end do nowait
!$omp end parallel


   npair_lk2=0


   do jl=1,nlk-lklen ! the last crosslinker is not in pairs

      if(pairnum(jl)==0)then
         cycle
      end if

      do j=1,pairnum(jl)!nnei


         npair_lk2=npair_lk2+1

         pair_lk2(1,npair_lk2)=jl

         pair_lk2(2,npair_lk2)=pairtem(j,jl)

      end do

   end do


!print*,'crlk-crlk',npair_lk2

!-------------------------------------------------------

!  membrane-membrane pairs:


   d2max=100.0d0!5*d2max

!$omp parallel &
!$omp default(none) &
!$omp private(jmb,n,jo,dx,dx2,d2,rad1_2,rad2_2,cos_2) &
!$omp shared(nmemb,xmemb,ymemb,zmemb,cos_t_2,pairtem,pairnum,d2max,nnei)

!$omp do schedule(guided,64)





   do jmb=1,nmemb-1

      n=0

      do jo=jmb+1,nmemb

         dx=0.5d0*(xmemb(jmb)-xmemb(jo))

         d2=abs(2*dx)+abs(ymemb(jmb)-ymemb(jo))+abs(zmemb(jmb)-zmemb(jo))

         if(d2>d2max)then
            cycle
         end if

         dx2=dx*dx

         rad1_2=dx2+ymemb(jmb)*ymemb(jmb)+zmemb(jmb)*zmemb(jmb)

         rad2_2=dx2+ymemb(jo)*ymemb(jo)+zmemb(jo)*zmemb(jo)

         cos_2=(dx2-ymemb(jmb)*ymemb(jo)-zmemb(jmb)*zmemb(jo))**2/rad1_2/rad2_2

         if(cos_2>cos_t_2)then

            if(n==nnei)then
               exit
            end if

            n=n+1
            pairtem(n,jmb)=jo


         end if

      end do

      pairnum(jmb)=n

   end do

!$omp end do nowait
!$omp end parallel

   npair_mb2=0


   do jmb=1,nmemb-1

      if(pairnum(jmb)==0)then
         cycle
      end if

      do j=1,pairnum(jmb)!nnei

         npair_mb2=npair_mb2+1

         pair_mb2(1,npair_mb2)=jmb

         pair_mb2(2,npair_mb2)=pairtem(j,jmb)


      end do

   end do


!print*,'membrane',npair_mb2

   deallocate(pairtem,pairnum)


   end subroutine

!============================================================================

   subroutine setpair(nmemb,npair5_myp2ac,npair5_lkac,npair5_myp2lk,npair5_ac2,npair5_myp2,npair5_lk2, &
               npair5_myo2ac,npair5_myo2lk,npair5_myo2,npair5_mymy,npair5_mb2,npair_myp2ac,npair_lkac, &
                npair_myp2lk,npair_ac2,npair_myp2,npair_lk2,npair_myo2ac,npair_myo2lk,npair_myo2, &
                 npair_mymy,npair_mb2,pair5_myp2ac,pair5_lkac,pair5_myp2lk,pair5_ac2,pair5_myp2, &
                  pair5_lk2,pair5_myo2ac,pair5_myo2lk,pair5_myo2,pair5_mymy,pair5_mb2,pair_myp2ac, &
                   pair_lkac,pair_myp2lk,pair_ac2,pair_myp2,pair_lk2,pair_myo2ac,pair_myo2lk,pair_myo2, &
                    pair_mymy,pair_mb2,pairpart,boundtyp,r_off,l_pair,thet2by2,l_mem,xboundmin,xboundmax, &
                     shift,xfa,yfa,zfa,xmyp2,ymyp2,zmyp2,xmyo2,ymyo2,zmyo2,xlk,ylk,zlk,xmemb,ymemb,zmemb)



   implicit none

   integer,value::nmemb,npair5_myp2ac,npair5_lkac,npair5_myp2lk,npair5_ac2,npair5_myp2,npair5_lk2
   integer,value::npair5_myo2ac,npair5_myo2lk,npair5_myo2,npair5_mymy,npair5_mb2

   integer npair_myp2ac,npair_lkac,npair_myp2lk,npair_ac2,npair_myp2,npair_lk2
   integer npair_myo2ac,npair_myo2lk,npair_myo2,npair_mymy,npair_mb2

   integer ja,ja1,ja2,jm,jm1,jm2,jl,jl1,jl2
   integer jmb,jo,n
   integer np1,np2,jb,jc,jd,nnew,jcycle
   integer nmax

   integer n1,n2,n3,n4,jexit,j1,j2,j3,j4,m1,m2,jp,np

   integer,allocatable,intent(in),dimension(:,:)::pair5_myp2ac,pair5_lkac,pair5_myp2lk,pair5_ac2,pair5_myp2
   integer,allocatable,intent(in),dimension(:,:)::pair5_lk2,pair5_myo2ac,pair5_myo2lk,pair5_myo2,pair5_mymy,pair5_mb2
   integer,allocatable,dimension(:,:)::pair_myp2ac,pair_lkac,pair_myp2lk,pair_ac2,pair_myp2,pair_lk2
   integer,allocatable,dimension(:,:)::pair_myo2ac,pair_myo2lk,pair_myo2,pair_mymy,pair_mb2,pairpart

   integer,allocatable,dimension(:)::boundtyp

   integer,allocatable,dimension(:)::pairtyp,mark,partlen

   integer,allocatable,dimension(:,:)::partner

   double precision,value::r_off,l_pair,thet2by2,l_mem,xboundmin,xboundmax,shift
   double precision dx,dx2,arg,rad1,rad2,dy,dz,d2,d2max,dist,xc,yc,zc,xsum,ysum,zsum,xmin,xmax
   double precision xl1,yl1,zl1,xl2,yl2,zl2,xp1,yp1,zp1,xp2,yp2,zp2,xp3,yp3,zp3

   double precision,allocatable,intent(in),dimension(:)::xfa,yfa,zfa,xmyp2,ymyp2,zmyp2,xmyo2,ymyo2,zmyo2
   double precision,allocatable,intent(in),dimension(:)::xlk,ylk,zlk,xmemb,ymemb,zmemb
   double precision,allocatable,dimension(:)::pairlen,lentem


   DOUBLE PRECISION EX1,EY1,EZ1,EX2,EY2,EZ2,TX,TY,TZ,PX,PY,PZ,QX,QY,QZ
   DOUBLE PRECISION DET,INVDET,T,U,V,delta,low,up,up1
!--------------------------



   d2max=2*r_off*r_off


   nmax=max(npair5_myp2ac,npair5_lkac,npair5_myp2lk,npair5_ac2,npair5_myp2,npair5_lk2,npair5_mb2)
   nmax=max(nmax,npair5_myo2ac,npair5_myo2lk,npair5_myo2,npair5_mymy)


   allocate(mark(nmax))


!  Myp2-actin pairs:

   mark(1:npair5_myp2ac)=0

!$omp parallel &
!$omp default(none) &
!$omp private(n,ja,jm,dx,dy,dz,d2) &
!$omp shared(npair5_myp2ac,pair5_myp2ac,xfa,yfa,zfa,xmyp2,ymyp2,zmyp2,d2max,mark)
!$omp do

   do n=1,npair5_myp2ac

      jm=pair5_myp2ac(1,n)

      ja=pair5_myp2ac(2,n)


      dx=xfa(ja)-xmyp2(jm)
      dy=yfa(ja)-ymyp2(jm)
      dz=zfa(ja)-zmyp2(jm)

      d2=dx*dx+dy*dy+dz*dz

      if(d2<d2max)then
         mark(n)=1
      end if

   end do

!$omp end do nowait
!$omp end parallel


   npair_myp2ac=0


   do n=1,npair5_myp2ac

      if(mark(n)==1)then

         npair_myp2ac=npair_myp2ac+1

         pair_myp2ac(1:2,npair_myp2ac)=pair5_myp2ac(1:2,n)

      end if

   end do

!-------------------------------------

!  Myo2-actin pairs:

   mark(1:npair5_myo2ac)=0

!$omp parallel &
!$omp default(none) &
!$omp private(n,ja,jm,dx,dy,dz,d2) &
!$omp shared(npair5_myo2ac,pair5_myo2ac,xfa,yfa,zfa,xmyo2,ymyo2,zmyo2,d2max,mark)
!$omp do

   do n=1,npair5_myo2ac

      jm=pair5_myo2ac(1,n)

      ja=pair5_myo2ac(2,n)


      dx=xfa(ja)-xmyo2(jm)
      dy=yfa(ja)-ymyo2(jm)
      dz=zfa(ja)-zmyo2(jm)

      d2=dx*dx+dy*dy+dz*dz

      if(d2<d2max)then
         mark(n)=1
      end if

   end do

!$omp end do nowait
!$omp end parallel

   npair_myo2ac=0


   do n=1,npair5_myo2ac

      if(mark(n)==1)then

         npair_myo2ac=npair_myo2ac+1

         pair_myo2ac(1:2,npair_myo2ac)=pair5_myo2ac(1:2,n)

      end if

   end do

!-------------------------------------

!  crosslinker-actin pairs:

   mark(1:npair5_lkac)=0


!$omp parallel &
!$omp default(none) &
!$omp private(n,ja,jl,dx,dy,dz,d2) &
!$omp shared(npair5_lkac,pair5_lkac,xfa,yfa,zfa,xlk,ylk,zlk,d2max,mark)
!$omp do

   do n=1,npair5_lkac

      jl=pair5_lkac(1,n)

      ja=pair5_lkac(2,n)


      dx=xfa(ja)-xlk(jl)
      dy=yfa(ja)-ylk(jl)
      dz=zfa(ja)-zlk(jl)

      d2=dx*dx+dy*dy+dz*dz

      if(d2<d2max)then
         mark(n)=1
      end if

   end do

!$omp end do nowait
!$omp end parallel


   npair_lkac=0


   do n=1,npair5_lkac

      if(mark(n)==1)then

         npair_lkac=npair_lkac+1

         pair_lkac(1:2,npair_lkac)=pair5_lkac(1:2,n)

      end if

   end do

!------------------------------------

!  Myp2-crosslinker pairs:


   mark(1:npair5_myp2lk)=0


!$omp parallel &
!$omp default(none) &
!$omp private(n,jm,jl,dx,dy,dz,d2) &
!$omp shared(npair5_myp2lk,pair5_myp2lk,xmyp2,ymyp2,zmyp2,xlk,ylk,zlk,d2max,mark)
!$omp do

   do n=1,npair5_myp2lk

      jm=pair5_myp2lk(1,n)

      jl=pair5_myp2lk(2,n)


      dx=xlk(jl)-xmyp2(jm)
      dy=ylk(jl)-ymyp2(jm)
      dz=zlk(jl)-zmyp2(jm)

      d2=dx*dx+dy*dy+dz*dz

      if(d2<d2max)then
         mark(n)=1
      end if

   end do


!$omp end do nowait
!$omp end parallel

   npair_myp2lk=0


   do n=1,npair5_myp2lk

      if(mark(n)==1)then

         npair_myp2lk=npair_myp2lk+1

         pair_myp2lk(1:2,npair_myp2lk)=pair5_myp2lk(1:2,n)

      end if

   end do


!--------------------------------------

!  Myo2-crosslinker pairs:

   mark(1:npair5_myo2lk)=0


!$omp parallel &
!$omp default(none) &
!$omp private(n,jm,jl,dx,dy,dz,d2) &
!$omp shared(npair5_myo2lk,pair5_myo2lk,xmyo2,ymyo2,zmyo2,xlk,ylk,zlk,d2max,mark)
!$omp do

   do n=1,npair5_myo2lk

      jm=pair5_myo2lk(1,n)

      jl=pair5_myo2lk(2,n)


      dx=xlk(jl)-xmyo2(jm)
      dy=ylk(jl)-ymyo2(jm)
      dz=zlk(jl)-zmyo2(jm)

      d2=dx*dx+dy*dy+dz*dz

      if(d2<d2max)then
         mark(n)=1
      end if

   end do


!$omp end do nowait
!$omp end parallel

   npair_myo2lk=0


   do n=1,npair5_myo2lk

      if(mark(n)==1)then

         npair_myo2lk=npair_myo2lk+1

         pair_myo2lk(1:2,npair_myo2lk)=pair5_myo2lk(1:2,n)

      end if

   end do


!--------------------------------------

!  actin-actin pairs:

   mark(1:npair5_ac2)=0


!$omp parallel &
!$omp default(none) &
!$omp private(n,ja1,ja2,dx,dy,dz,d2) &
!$omp shared(npair5_ac2,pair5_ac2,xfa,yfa,zfa,d2max,mark)
!$omp do


   do n=1,npair5_ac2

      ja1=pair5_ac2(1,n)

      ja2=pair5_ac2(2,n)


      dx=xfa(ja1)-xfa(ja2)
      dy=yfa(ja1)-yfa(ja2)
      dz=zfa(ja1)-zfa(ja2)

      d2=dx*dx+dy*dy+dz*dz

      if(d2<d2max)then
         mark(n)=1
      end if

   end do

!$omp end do nowait
!$omp end parallel


   npair_ac2=0


   do n=1,npair5_ac2

      if(mark(n)==1)then

         npair_ac2=npair_ac2+1

         pair_ac2(1:2,npair_ac2)=pair5_ac2(1:2,n)

      end if

   end do

!--------------------------------------

!  Myp2-Myp2 pairs:

   mark(1:npair5_myp2)=0

!$omp parallel &
!$omp default(none) &
!$omp private(n,jm1,jm2,dx,dy,dz,d2) &
!$omp shared(npair5_myp2,pair5_myp2,xmyp2,ymyp2,zmyp2,d2max,mark)
!$omp do

   do n=1,npair5_myp2

      jm1=pair5_myp2(1,n)

      jm2=pair5_myp2(2,n)

      dx=xmyp2(jm1)-xmyp2(jm2)
      dy=ymyp2(jm1)-ymyp2(jm2)
      dz=zmyp2(jm1)-zmyp2(jm2)

      d2=dx*dx+dy*dy+dz*dz

      if(d2<d2max)then
         mark(n)=1
      end if

   end do

!$omp end do nowait
!$omp end parallel

   npair_myp2=0


   do n=1,npair5_myp2

      if(mark(n)==1)then

         npair_myp2=npair_myp2+1

         pair_myp2(1:2,npair_myp2)=pair5_myp2(1:2,n)

      end if

   end do


!--------------------------------------

!  Myo2-Myo2 pairs:

   mark(1:npair5_myo2)=0

!$omp parallel &
!$omp default(none) &
!$omp private(n,jm1,jm2,dx,dy,dz,d2) &
!$omp shared(npair5_myo2,pair5_myo2,xmyo2,ymyo2,zmyo2,d2max,mark)
!$omp do

   do n=1,npair5_myo2

      jm1=pair5_myo2(1,n)

      jm2=pair5_myo2(2,n)

      dx=xmyo2(jm1)-xmyo2(jm2)
      dy=ymyo2(jm1)-ymyo2(jm2)
      dz=zmyo2(jm1)-zmyo2(jm2)

      d2=dx*dx+dy*dy+dz*dz

      if(d2<d2max)then
         mark(n)=1
      end if

   end do

!$omp end do nowait
!$omp end parallel

   npair_myo2=0

   do n=1,npair5_myo2

      if(mark(n)==1)then

         npair_myo2=npair_myo2+1

         pair_myo2(1:2,npair_myo2)=pair5_myo2(1:2,n)

      end if

   end do


!--------------------------------------

!  crosslinker-crosslinker pairs:

   mark(1:npair5_lk2)=0


!$omp parallel &
!$omp default(none) &
!$omp private(n,jl1,jl2,dx,dy,dz,d2) &
!$omp shared(npair5_lk2,pair5_lk2,xlk,ylk,zlk,d2max,mark)
!$omp do

   do n=1,npair5_lk2

      jl1=pair5_lk2(1,n)

      jl2=pair5_lk2(2,n)


      dx=xlk(jl1)-xlk(jl2)
      dy=ylk(jl1)-ylk(jl2)
      dz=zlk(jl1)-zlk(jl2)

      d2=dx*dx+dy*dy+dz*dz

      if(d2<d2max)then
         mark(n)=1
      end if

   end do

!$omp end do nowait
!$omp end parallel

   npair_lk2=0


   do n=1,npair5_lk2

      if(mark(n)==1)then

         npair_lk2=npair_lk2+1

         pair_lk2(1:2,npair_lk2)=pair5_lk2(1:2,n)

      end if

   end do

!-----------------------------------------

!  Myp2-Myo2 pairs:

   mark(1:npair5_mymy)=0

!$omp parallel &
!$omp default(none) &
!$omp private(n,jm,jp,dx,dy,dz,d2) &
!$omp shared(npair5_mymy,pair5_mymy,xmyo2,ymyo2,zmyo2,xmyp2,ymyp2,zmyp2,d2max,mark)
!$omp do

   do n=1,npair5_mymy

      jm=pair5_mymy(1,n)

      jp=pair5_mymy(2,n)

      dx=xmyp2(jp)-xmyo2(jm)
      dy=ymyp2(jp)-ymyo2(jm)
      dz=zmyp2(jp)-zmyo2(jm)

      d2=dx*dx+dy*dy+dz*dz

      if(d2<d2max)then
         mark(n)=1
      end if

   end do

!$omp end do nowait
!$omp end parallel

   npair_mymy=0

   do n=1,npair5_mymy

      if(mark(n)==1)then

         npair_mymy=npair_mymy+1

         pair_mymy(1:2,npair_mymy)=pair5_mymy(1:2,n)

      end if

   end do



!-----------------------------------------

   mark(1:npair5_mb2)=0

   allocate(lentem(npair5_mb2),pairlen(npair5_mb2))

!$omp parallel &
!$omp default(none) &
!$omp private(n,jmb,jo,dx,dx2,rad1,rad2,arg) &
!$omp shared(npair5_mb2,pair5_mb2,xmemb,ymemb,zmemb,thet2by2,mark,lentem)
!$omp do


   do n=1,npair5_mb2

      jmb=pair5_mb2(1,n)

      jo=pair5_mb2(2,n)

      dx=0.5d0*(xmemb(jmb)-xmemb(jo))

      dx2=dx*dx

      rad1=sqrt(dx2+ymemb(jmb)*ymemb(jmb)+zmemb(jmb)*zmemb(jmb))

      rad2=sqrt(dx2+ymemb(jo)*ymemb(jo)+zmemb(jo)*zmemb(jo))

      arg=1.0d0+(dx2-ymemb(jmb)*ymemb(jo)-zmemb(jmb)*zmemb(jo))/rad1/rad2

      if(arg<thet2by2)then

         mark(n)=1

         lentem(n)=arg

      end if

   end do

!$omp end do nowait
!$omp end parallel

   npair_mb2=0

   do n=1,npair5_mb2

      if(mark(n)==1)then

         npair_mb2=npair_mb2+1

         pair_mb2(1:2,npair_mb2)=pair5_mb2(1:2,n)

         pairlen(npair_mb2)=lentem(n)

      end if

   end do

   deallocate(lentem,mark)

!  boundary problem:

   boundtyp=0


   xmin=xboundmin+l_mem

   xmax=xboundmax-l_mem

   do n=1,nmemb


      if(xmemb(n)<xmin)then

         boundtyp(n)=1


      end if

      if(xmemb(n)>xmax)then

         boundtyp(n)=2

      end if

   end do

   d2max=100.0d0

   do n1=1,nmemb

      if(boundtyp(n1)/=1) cycle

      do n2=1,nmemb

         if(boundtyp(n2)/=2) cycle

         dx=0.5d0*(xmemb(n1)+shift-xmemb(n2))

         d2=abs(2*dx)+abs(ymemb(n1)-ymemb(n2))+abs(zmemb(n1)-zmemb(n2))

         if(d2>d2max) cycle

         dx2=dx*dx

         rad1=sqrt(dx2+ymemb(n1)*ymemb(n1)+zmemb(n1)*zmemb(n1))

         rad2=sqrt(dx2+ymemb(n2)*ymemb(n2)+zmemb(n2)*zmemb(n2))

         arg=1.0d0+(dx2-ymemb(n1)*ymemb(n2)-zmemb(n1)*zmemb(n2))/rad1/rad2

         if(arg<thet2by2)then

            npair_mb2=npair_mb2+1

            pair_mb2(1,npair_mb2)=n1

            pair_mb2(2,npair_mb2)=n2

            pairlen(npair_mb2)=arg


         end if

      end do

   end do



!=============================================


   call sorting(pairlen(:npair_mb2),pair_mb2(:,:npair_mb2))


!  remove cross-over pairs:


   allocate(pairtyp(npair_mb2))
   pairtyp=1

   yp3=0.0d0; zp3=0.0d0


   dist=4*l_pair

   delta=0.000000000001d0

   low=-delta
   up=1.0d0+delta
   up1=1.02d0

   do np1=1,npair_mb2-1

      if(pairtyp(np1)==0)then
         cycle
      end if

      ja=pair_mb2(1,np1)

      jb=pair_mb2(2,np1)


      xp1=xmemb(ja)
      yp1=ymemb(ja)
      zp1=zmemb(ja)

      if(boundtyp(ja)==1.and.boundtyp(jb)==2) xp1=xp1+shift



      xp2=xmemb(jb)
      yp2=ymemb(jb)
      zp2=zmemb(jb)

      xc=0.5d0*(xp1+xp2)
      yc=0.5d0*(yp1+yp2)
      zc=0.5d0*(zp1+zp2)

      xp1=xc+(xp1-xc)*up1
      yp1=yc+(yp1-yc)*up1
      zp1=zc+(zp1-zc)*up1

      xp2=xc+(xp2-xc)*up1
      yp2=yc+(yp2-yc)*up1
      zp2=zc+(zp2-zc)*up1

      xp3=0.5d0*(xp1+xp2)

!     original sum:
      xsum=xp1+xp2
      ysum=yp1+yp2
      zsum=zp1+zp2

!     then scaled:
      xp1=xp1+xp1-xp3
      yp1=yp1+yp1
      zp1=zp1+zp1

      xp2=xp2+xp2-xp3
      yp2=yp2+yp2
      zp2=zp2+zp2


!$omp parallel &
!$omp default(none) &
!$omp private(np2,jc,jd) &
!$omp private(dx,dy,dz,xl1,yl1,zl1,xl2,yl2,zl2,xc,yc,zc) &
!$omp private(EX1,EY1,EZ1,EX2,EY2,EZ2,TX,TY,TZ,PX,PY,PZ,QX,QY,QZ,DET,INVDET,T,U,V)&

!$omp shared(np1,npair_mb2,pairtyp,ja,jb,xp1,yp1,zp1,xp2,yp2,zp2,xp3,yp3,zp3) &
!$omp shared(xmemb,ymemb,zmemb,dist,pairlen,pair_mb2,delta,low,up,up1,l_mem,xsum,ysum,zsum) &
!$omp shared(shift,boundtyp)

!$omp do schedule(guided,64)



      do np2=np1+1,npair_mb2

         if(pairtyp(np2)==0)then
            cycle
         end if


         jc=pair_mb2(1,np2)

         jd=pair_mb2(2,np2)

         if(jc==ja.or.jc==jb.or.jd==ja.or.jd==jb)then
            cycle
         end if

         if(boundtyp(ja)==1.and.boundtyp(jb)==2.and.boundtyp(jc)==0.and.boundtyp(jd)==0) cycle

         if(boundtyp(ja)==0.and.boundtyp(jb)==0.and.boundtyp(jc)==1.and.boundtyp(jd)==2) cycle


         xl1=xmemb(jc)
         yl1=ymemb(jc)
         zl1=zmemb(jc)

         xl2=xmemb(jd)
         yl2=ymemb(jd)
         zl2=zmemb(jd)

         if(boundtyp(ja)==1.and.boundtyp(jb)==2)then

            if(boundtyp(jc)==1) xl1=xl1+shift

            if(boundtyp(jd)==1) xl2=xl2+shift

         else

            if(boundtyp(ja)==1.and.boundtyp(jc)==1.and.boundtyp(jd)==2) xl2=xl2-shift

            if(boundtyp(ja)==2.and.boundtyp(jc)==1.and.boundtyp(jd)==2) xl1=xl1+shift

         end if


         dy=abs(ysum-yl1-yl2)
!         dy=abs(ymemb(ja)+ymemb(jb)-ymemb(jc)-ymemb(jd))

         if(dy>dist)then
            cycle
         end if

         dz=abs(zsum-zl1-zl2)

!         dz=abs(zmemb(ja)+zmemb(jb)-zmemb(jc)-zmemb(jd))

         if(dz>dist)then
            cycle
         end if

         dx=abs(xsum-xl1-xl2)

!         dx=abs(xmemb(ja)+xmemb(jb)-xmemb(jc)-xmemb(jd))

         if(dx>dist)then
            cycle
         end if


!         xl1=xmemb(jc)
!         yl1=ymemb(jc)
!         zl1=zmemb(jc)

!         xl2=xmemb(jd)
!         yl2=ymemb(jd)
!         zl2=zmemb(jd)

         xc=0.5d0*(xl1+xl2)
         yc=0.5d0*(yl1+yl2)
         zc=0.5d0*(zl1+zl2)

         xl1=xc+(xl1-xc)*up1
         yl1=yc+(yl1-yc)*up1
         zl1=zc+(zl1-zc)*up1

         xl2=xc+(xl2-xc)*up1
         yl2=yc+(yl2-yc)*up1
         zl2=zc+(zl2-zc)*up1

!------------------------------------------


         DX=XL2-XL1
         DY=YL2-YL1
         DZ=ZL2-ZL1

         EX1=XP2-XP1
         EY1=YP2-YP1
         EZ1=ZP2-ZP1

         EX2=XP3-XP1
         EY2=YP3-YP1
         EZ2=ZP3-ZP1

         TX=XL1-XP1
         TY=YL1-YP1
         TZ=ZL1-ZP1

         PX=DY*EZ2-EY2*DZ
         PY=DZ*EX2-EZ2*DX
         PZ=DX*EY2-EX2*DY

         QX=TY*EZ1-EY1*TZ
         QY=TZ*EX1-EZ1*TX
         QZ=TX*EY1-EX1*TY

         DET=PX*EX1+PY*EY1+PZ*EZ1

         if(abs(det)<delta)then
            cycle
         endif

         INVDET=1.0D0/DET

         T=(QX*EX2+QY*EY2+QZ*EZ2)*INVDET


         IF(T<low.OR.T>up)THEN
            cycle
         END IF

         U=(PX*TX+PY*TY+PZ*TZ)*INVDET

         IF(U<low.OR.U>up)THEN
            cycle
         END IF

         V=(QX*DX+QY*DY+QZ*DZ)*INVDET

         if(u+v<low.or.u+v>up)then
            cycle
         end if



         pairtyp(np2)=0


      end do

!$omp end do nowait
!$omp end parallel

   end do

   nnew=0

   do n=1,npair_mb2

      if(pairtyp(n)==0)then
         cycle
      end if

      n1=pair_mb2(1,n)
      n2=pair_mb2(2,n)

      if(n1==n2)then
         cycle
      end if

      nnew=nnew+1

      pair_mb2(1,nnew)=n1!pair_mb2(1:2,n)

      pair_mb2(2,nnew)=n2

   end do

   npair_mb2=nnew

!print*,'pairs',nnew

   deallocate(pairlen,pairtyp)


!-------------------------------

!  list of tetrahedrons

!  partners of the pairs in tetrahedrons:

!   pairtyp(1:npair_mb2)=0

   allocate(partner(20,nmemb),partlen(nmemb))

   partlen=0

   do n=1,npair_mb2

      n1=pair_mb2(1,n)

      n2=pair_mb2(2,n)

      partlen(n1)=partlen(n1)+1

      partner(partlen(n1),n1)=n2

      partlen(n2)=partlen(n2)+1

      partner(partlen(n2),n2)=n1

!      pairidlen(n1)=pairidlen(n1)+1

!      pairid(pairidlen(n1),n1)=n

!      pairidlen(n2)=pairidlen(n2)+1

!      pairid(pairidlen(n2),n2)=n

   end do

!   ntetrahed=0

   pairpart=0

   do n=1,npair_mb2

!      if(pairtyp(n)==1)then
!         cycle
!      end if

      n1=pair_mb2(1,n)

      n2=pair_mb2(2,n)

      n3=0

      n4=0

      jexit=0

      do j1=1,partlen(n1)

         m1=partner(j1,n1)

         do j2=1,partlen(n2)

            m2=partner(j2,n2)

            if(m1==m2)then

               if(n3==0)then

                  n3=m1

               else

                  n4=m1

                  jexit=1

                  exit

               end if

            end if

         end do

         if(jexit==1)then

            exit

         end if

      end do

      if(n4>0)then

         pairpart(1,n)=n3
         pairpart(2,n)=n4



      end if

   end do


   deallocate(partner,partlen)



   end subroutine

!==========================================

recursive   subroutine sorting(pairlen,pair)

   implicit none
   integer iq
   integer, intent(in out),dimension(:,:)::pair
   double precision,intent(in out),dimension(:)::pairlen

   if(size(pairlen)>1)then

      call partition(pairlen,pair,iq)

      call sorting(pairlen(:iq-1),pair(:,:iq-1))

      call sorting(pairlen(iq:),pair(:,iq:))

   end if

   end subroutine sorting

   subroutine partition(pairlen,pair,marker)

   integer marker,i,j,n1,n2
   integer, intent(in out),dimension(:,:)::pair
   double precision,intent(in out),dimension(:)::pairlen
   double precision x,temp

   x=pairlen(1)

   i=0

   j=size(pairlen)+1

   do

      j=j-1

      do

         if(pairlen(j)<=x) exit
         j=j-1

      end do

      i=i+1

      do

        if(pairlen(i)>=x) exit
        i=i+1

      end do

      if(i<j)then

         temp=pairlen(i)

         n1=pair(1,i)
         n2=pair(2,i)

         pairlen(i)=pairlen(j)

         pair(1:2,i)=pair(1:2,j)

         pairlen(j)=temp

         pair(1,j)=n1

         pair(2,j)=n2

      elseif(i==j)then
         marker=i+1
         return
      else
         marker=i
         return
      end if

   end do

   end subroutine partition

!==========================================

   SUBROUTINE INTERSEC(CHECK,XL1,YL1,ZL1,XL2,YL2,ZL2,XP1,YP1,ZP1,XP2,YP2,ZP2,XP3,YP3,ZP3)

   IMPLICIT NONE

   INTEGER CHECK

   DOUBLE PRECISION,VALUE:: XL1,YL1,ZL1,XL2,YL2,ZL2,XP1,YP1,ZP1,XP2,YP2,ZP2,XP3,YP3,ZP3
   DOUBLE PRECISION M11,M12,M13,M21,M22,M23,M31,M32,M33
   DOUBLE PRECISION I11,I12,I13,I21,I22,I23,I31,I32,I33
   DOUBLE PRECISION DET,INVDET,T,U,V,XB,YB,ZB,amin,amax

! -- FORM THE MATRIX FROM THE LINE THROUGH L1 L2 AND PLANE THROUGH P1 P2 P3:


   M11=XL1-XL2
   M21=YL1-YL2
   M31=ZL1-ZL2

   M12=XP2-XP1
   M22=YP2-YP1
   M32=ZP2-ZP1

   M13=XP3-XP1
   M23=YP3-YP1
   M33=ZP3-ZP1

   DET=M11*M22*M33+M12*M23*M31+M13*M21*M32-M11*M23*M32-M12*M21*M33-M13*M22*M31


   if(abs(det)<0.000000001d0)then
      check=0
      return
   end if

   INVDET=1.0D0/DET

!--- MATRIX INVERSE:

   I11=(M22*M33-M23*M32)*INVDET
   I12=(M13*M32-M12*M33)*INVDET
   I13=(M12*M23-M13*M22)*INVDET

   I21=(M23*M31-M21*M33)*INVDET
   I22=(M11*M33-M13*M31)*INVDET
   I23=(M13*M21-M11*M23)*INVDET

   I31=(M21*M32-M22*M31)*INVDET
   I32=(M12*M31-M11*M32)*INVDET
   I33=(M11*M22-M12*M21)*INVDET

! --- BASE VECTOR:
   XB=XL1-XP1
   YB=YL1-YP1
   ZB=ZL1-ZP1

! -- INTERSECTION IS REPRESENTED BY (T U V):

   T=I11*XB+I12*YB+I13*ZB

   U=I21*XB+I22*YB+I23*ZB

   V=I31*XB+I32*YB+I33*ZB


! --- INTERSECTION OCCURS IF T =(0,1); U,V = (0,1); U+V = (0,1)
! --- FIRST ASSUME CHECK = 1:

   CHECK = 1

!   IF(T<0.0D0.OR.T>1.0D0)THEN
!      CHECK=0
!   END IF

!   IF(U<0.0D0.OR.U>1.0D0)THEN
!      CHECK=0
!   END IF

!   IF(V<0.0D0.OR.V>1.0D0)THEN
!      CHECK=0
!   END IF

!   IF(U+V<0.0D0.OR.U+V>1.0D0)THEN
!      CHECK=0
!   END IF

   amin=min(T,U,V,U+V)
   amax=max(T,U,V,U+V)

!   IF(T<0.0D0.OR.T>1.0D0.OR.U<0.0D0.OR.U>1.0D0.OR.V<0.0D0.OR.V>1.0D0.OR.U+V<0.0D0.OR.U+V>1.0D0)THEN
!      CHECK=0
!   END IF

   if(amin<0.0d0.or.amax>1.0d0)then
      check=0
   end if

   END SUBROUTINE

!=========================================================

   SUBROUTINE INTERSEC1(CHECK,XL1,YL1,ZL1,XL2,YL2,ZL2,XP1,YP1,ZP1,XP2,YP2,ZP2,XP3,YP3,ZP3)

   IMPLICIT NONE

   INTEGER CHECK

   DOUBLE PRECISION,VALUE:: XL1,YL1,ZL1,XL2,YL2,ZL2,XP1,YP1,ZP1,XP2,YP2,ZP2,XP3,YP3,ZP3
   DOUBLE PRECISION DX,DY,DZ,EX1,EY1,EZ1,EX2,EY2,EZ2,TX,TY,TZ,PX,PY,PZ,QX,QY,QZ
   DOUBLE PRECISION DET,INVDET,T,U,V

   DX=XL2-XL1
   DY=YL2-YL1
   DZ=ZL2-ZL1

   EX1=XP2-XP1
   EY1=YP2-YP1
   EZ1=ZP2-ZP1

   EX2=XP3-XP1
   EY2=YP3-YP1
   EZ2=ZP3-ZP1

   TX=XL1-XP1
   TY=YL1-YP1
   TZ=ZL1-ZP1

   PX=DY*EZ2-EY2*DZ
   PY=DZ*EX2-EZ2*DX
   PZ=DX*EY2-EX2*DY

   QX=TY*EZ1-EY1*TZ
   QY=TZ*EX1-EZ1*TX
   QZ=TX*EY1-EX1*TY

   DET=PX*EX1+PY*EY1+PZ*EZ1

   INVDET=1.0D0/DET

   T=(QX*EX2+QY*EY2+QZ*EZ2)*INVDET
   U=(PX*TX+PY*TY+PZ*TZ)*INVDET
   V=(QX*DX+QY*DY+QZ*DZ)*INVDET

   CHECK = 1

   IF(T<0.0D0.OR.T>1.0D0)THEN
      CHECK=0
      RETURN
   END IF

   IF(U<0.0D0.OR.U>1.0D0)THEN
      CHECK=0
      RETURN
   END IF

   IF(V<0.0D0.OR.U+V>1.0D0)THEN
      CHECK=0
      RETURN
   END IF

   END SUBROUTINE


!=========================================================

   subroutine solidset(nxsol,nphisol,nmemb,nfa,nmyp2,nmyo2,nlk,jmbsol,jfasol,jmyp2sol,jmyo2sol,jlksol, &
               pi,delta,dxsol,dphisol,xmemb,ymemb,zmemb,xfa,yfa,zfa,xmyp2,ymyp2,zmyp2, &
                xmyo2,ymyo2,zmyo2,xlk,ylk,zlk,xwall,ywall,zwall,xnorwall,ynorwall,znorwall, &
                 xsurf,ysurf,zsurf,xnorsurf,ynorsurf,znorsurf,jsursol)



   implicit none

   integer,value:: nxsol,nphisol,nmemb,nfa,nmyp2,nmyo2,nlk
   integer n,jx,jp,jx0,jp0,j1,j2,jxget,jpget

   integer,allocatable,dimension(:,:)::jmbsol,jfasol,jmyp2sol,jmyo2sol,jlksol,jsursol

   double precision,value:: pi,delta,dxsol,dphisol
   double precision xmin,dxsolby2,piby2,dphisolby2,twopi,phi,arg
   double precision dx,dy,dz,d2,dist2

   double precision,allocatable,intent(in),dimension(:)::xmemb,ymemb,zmemb,xfa,yfa,zfa
   double precision,allocatable,intent(in),dimension(:)::xmyp2,ymyp2,zmyp2,xmyo2,ymyo2,zmyo2,xlk,ylk,zlk
   double precision,allocatable,intent(in),dimension(:,:)::xwall,ywall,zwall,xnorwall,ynorwall,znorwall
   double precision,allocatable,intent(in),dimension(:,:)::xsurf,ysurf,zsurf,xnorsurf,ynorsurf,znorsurf

   xmin=-(nxsol-1)/2*dxsol

   dxsolby2=0.5d0*dxsol

   piby2=0.5d0*pi

   dphisolby2=0.5d0*dphisol

   twopi=2*pi

!$omp parallel &
!$omp default(none) &
!$omp private(n,jx,phi,arg,jp,jx0,jp0,j1,j2,jxget,jpget,dx,dy,dz,d2,dist2) &
!$omp shared(nmemb,xmemb,xmin,dxsol,dxsolby2,nxsol,jmbsol,ymemb,delta,zmemb,piby2) &
!$omp shared(dphisolby2,nphisol,pi,twopi,dphisol) &
!$omp shared(xwall,ywall,zwall,xnorwall,ynorwall,znorwall,jsursol) &

!$omp shared(nfa,xfa,jfasol,yfa,zfa) &
!$omp shared(xsurf,ysurf,zsurf,xnorsurf,ynorsurf,znorsurf) &
!$omp shared(nmyp2,xmyp2,jmyp2sol,ymyp2,zmyp2) &
!$omp shared(nmyo2,xmyo2,jmyo2sol,ymyo2,zmyo2) &

!$omp shared(nlk,xlk,jlksol,ylk,zlk)

!$omp do schedule(guided,64)

   do n=1,nmemb

      jx=1+(xmemb(n)-xmin)/dxsol

      if(xmemb(n)-xmin-(jx-1)*dxsol>dxsolby2)then
         jx=jx+1
      end if

      if(jx<1)then
         jx=1
      end if

      if(jx>nxsol)then
         jx=nxsol
      end if

      jsursol(1,n)=jx

      jx0=jx

      if(abs(ymemb(n))<delta)then
         if(zmemb(n)>0.0d0)then
            phi=piby2
         else
            phi=-piby2
         end if

      else

         arg=zmemb(n)/ymemb(n)

         phi=atan(arg)

         if(ymemb(n)<0.0d0)then
            phi=phi+pi
         end if

         if(arg<0.0d0.and.ymemb(n)>0.0d0)then
            phi=phi+twopi
         end if


      end if

      if(phi<dphisolby2)then
         jp=nphisol
      else

         jp=phi/dphisol

         if(phi-jp*dphisol>dphisolby2)then
            jp=jp+1
         end if

      end if

      jsursol(2,n)=jp


      jp0=jp


      dist2=1000000.0d0

      jxget=0

      do j1=1,20

         jx=jx0+j1-10

         if(jx<1) cycle

         if(jx>nxsol) cycle

         do j2=1,20

            jp=jp0+j2-10

            if(jp<1) jp=jp+nphisol

            if(jp>nphisol) jp=jp-nphisol

            if(ymemb(n)*ywall(jp,jx)+zmemb(n)*zwall(jp,jx)<0.0) cycle

            dx=xmemb(n)-xwall(jp,jx)
            dy=ymemb(n)-ywall(jp,jx)
            dz=zmemb(n)-zwall(jp,jx)

            arg=dx*xnorwall(jp,jx)+dy*ynorwall(jp,jx)+dz*znorwall(jp,jx)

            if(arg<0.0d0) cycle

            dx=dx-arg*xnorwall(jp,jx)
            dy=dy-arg*ynorwall(jp,jx)
            dz=dz-arg*znorwall(jp,jx)

            d2=dx*dx+dy*dy+dz*dz

            if(dist2>d2)then

               dist2=d2

               jxget=jx

               jpget=jp

            end if

         end do

      end do

      if(jxget==0)then
         print*,'error: could not get indices',n,jx0,jp0
         stop
      end if

      jmbsol(1,n)=jxget

      jmbsol(2,n)=jpget


   end do

!$omp end do nowait



!$omp do schedule(guided,64)

   do n=1,nfa

      jx=1+(xfa(n)-xmin)/dxsol

      if(xfa(n)-xmin-(jx-1)*dxsol>dxsolby2)then
         jx=jx+1
      end if

      if(jx<1)then
         jx=1
      end if

      if(jx>nxsol)then
         jx=nxsol
      end if

!      jfasol(1,n)=jx

      jx0=jx

      if(abs(yfa(n))<delta)then
         if(zfa(n)>0.0d0)then
            phi=piby2
         else
            phi=-piby2
         end if

      else

         arg=zfa(n)/yfa(n)

         phi=atan(arg)

         if(yfa(n)<0.0d0)then
            phi=phi+pi
         end if

         if(arg<0.0d0.and.yfa(n)>0.0d0)then
            phi=phi+twopi
         end if

      end if

      if(phi<dphisolby2)then
         jp=nphisol
      else

         jp=phi/dphisol

         if(phi-jp*dphisol>dphisolby2)then
            jp=jp+1
         end if

      end if

!      jfasol(2,n)=jp

      jp0=jp


      dist2=1000000.0d0

      jxget=0

      do j1=1,20

         jx=jx0+j1-10

         if(jx<1) cycle

         if(jx>nxsol) cycle

         do j2=1,20

            jp=jp0+j2-10

            if(jp<1) jp=jp+nphisol

            if(jp>nphisol) jp=jp-nphisol

            if(yfa(n)*ysurf(jp,jx)+zfa(n)*zsurf(jp,jx)<0.0) cycle

            dx=xfa(n)-xsurf(jp,jx)
            dy=yfa(n)-ysurf(jp,jx)
            dz=zfa(n)-zsurf(jp,jx)

            arg=dx*xnorsurf(jp,jx)+dy*ynorsurf(jp,jx)+dz*znorsurf(jp,jx)

            if(arg<0.0d0) cycle

            dx=dx-arg*xnorsurf(jp,jx)
            dy=dy-arg*ynorsurf(jp,jx)
            dz=dz-arg*znorsurf(jp,jx)

            d2=dx*dx+dy*dy+dz*dz

            if(dist2>d2)then

               dist2=d2

               jxget=jx

               jpget=jp

            end if

         end do

      end do

      if(jxget==0)then
         print*,'error: could not get indices',n,jx0,jp0
         stop
      end if

      jfasol(1,n)=jxget

      jfasol(2,n)=jpget

   end do

!$omp end do nowait


!$omp do schedule(guided,64)

   do n=1,nmyp2

      jx=1+(xmyp2(n)-xmin)/dxsol

      if(xmyp2(n)-xmin-(jx-1)*dxsol>dxsolby2)then
         jx=jx+1
      end if

      if(jx<1)then
         jx=1
      end if

      if(jx>nxsol)then
         jx=nxsol
      end if

!      jmyp2sol(1,n)=jx

      jx0=jx

      if(abs(ymyp2(n))<delta)then
         if(zmyp2(n)>0.0d0)then
            phi=piby2
         else
            phi=-piby2
         end if

      else

         arg=zmyp2(n)/ymyp2(n)

         phi=atan(arg)

         if(ymyp2(n)<0.0d0)then
            phi=phi+pi
         end if

         if(arg<0.0d0.and.ymyp2(n)>0.0d0)then
            phi=phi+twopi
         end if

      end if

      if(phi<dphisolby2)then
         jp=nphisol
      else

         jp=phi/dphisol

         if(phi-jp*dphisol>dphisolby2)then
            jp=jp+1
         end if

      end if

!      jmyp2sol(2,n)=jp

      jp0=jp

      dist2=1000000.0d0

      jxget=0

      do j1=1,20

         jx=jx0+j1-10

         if(jx<1) cycle

         if(jx>nxsol) cycle

         do j2=1,20

            jp=jp0+j2-10

            if(jp<1) jp=jp+nphisol

            if(jp>nphisol) jp=jp-nphisol

            if(ymyp2(n)*ysurf(jp,jx)+zmyp2(n)*zsurf(jp,jx)<0.0) cycle

            dx=xmyp2(n)-xsurf(jp,jx)
            dy=ymyp2(n)-ysurf(jp,jx)
            dz=zmyp2(n)-zsurf(jp,jx)

            arg=dx*xnorsurf(jp,jx)+dy*ynorsurf(jp,jx)+dz*znorsurf(jp,jx)

            if(arg<0.0d0) cycle

            dx=dx-arg*xnorsurf(jp,jx)
            dy=dy-arg*ynorsurf(jp,jx)
            dz=dz-arg*znorsurf(jp,jx)

            d2=dx*dx+dy*dy+dz*dz

            if(dist2>d2)then

               dist2=d2

               jxget=jx

               jpget=jp

            end if

         end do

      end do

      if(jxget==0)then
         print*,'error: could not get indices',n,jx0,jp0
         stop
      end if

      jmyp2sol(1,n)=jxget

      jmyp2sol(2,n)=jpget

   end do

!$omp end do nowait


!$omp do schedule(guided,64)

   do n=1,nmyo2

      jx=1+(xmyo2(n)-xmin)/dxsol

      if(xmyo2(n)-xmin-(jx-1)*dxsol>dxsolby2)then
         jx=jx+1
      end if

      if(jx<1)then
         jx=1
      end if

      if(jx>nxsol)then
         jx=nxsol
      end if

      jx0=jx

      if(abs(ymyo2(n))<delta)then
         if(zmyo2(n)>0.0d0)then
            phi=piby2
         else
            phi=-piby2
         end if

      else

         arg=zmyo2(n)/ymyo2(n)

         phi=atan(arg)

         if(ymyo2(n)<0.0d0)then
            phi=phi+pi
         end if

         if(arg<0.0d0.and.ymyo2(n)>0.0d0)then
            phi=phi+twopi
         end if

      end if

      if(phi<dphisolby2)then
         jp=nphisol
      else

         jp=phi/dphisol

         if(phi-jp*dphisol>dphisolby2)then
            jp=jp+1
         end if

      end if

!      jmyo2sol(2,n)=jp

      jp0=jp

      dist2=1000000.0d0

      jxget=0

      do j1=1,20

         jx=jx0+j1-10

         if(jx<1) cycle

         if(jx>nxsol) cycle

         do j2=1,20

            jp=jp0+j2-10

            if(jp<1) jp=jp+nphisol

            if(jp>nphisol) jp=jp-nphisol

            if(ymyo2(n)*ysurf(jp,jx)+zmyo2(n)*zsurf(jp,jx)<0.0) cycle

            dx=xmyo2(n)-xsurf(jp,jx)
            dy=ymyo2(n)-ysurf(jp,jx)
            dz=zmyo2(n)-zsurf(jp,jx)

            arg=dx*xnorsurf(jp,jx)+dy*ynorsurf(jp,jx)+dz*znorsurf(jp,jx)

            if(arg<0.0d0) cycle

            dx=dx-arg*xnorsurf(jp,jx)
            dy=dy-arg*ynorsurf(jp,jx)
            dz=dz-arg*znorsurf(jp,jx)

            d2=dx*dx+dy*dy+dz*dz

            if(dist2>d2)then

               dist2=d2

               jxget=jx

               jpget=jp

            end if

         end do

      end do

      if(jxget==0)then
         print*,'error: could not get indices',n,jx0,jp0
         stop
      end if

      jmyo2sol(1,n)=jxget

      jmyo2sol(2,n)=jpget

   end do

!$omp end do nowait



!$omp do schedule(guided,32)

   do n=1,nlk

      jx=1+(xlk(n)-xmin)/dxsol

      if(xlk(n)-xmin-(jx-1)*dxsol>dxsolby2)then
         jx=jx+1
      end if

      if(jx<1)then
         jx=1
      end if

      if(jx>nxsol)then
         jx=nxsol
      end if

!      jlksol(1,n)=jx

      jx0=jx

      if(abs(ylk(n))<delta)then
         if(zlk(n)>0.0d0)then
            phi=piby2
         else
            phi=-piby2
         end if

      else

         arg=zlk(n)/ylk(n)

         phi=atan(arg)

         if(ylk(n)<0.0d0)then
            phi=phi+pi
         end if

         if(arg<0.0d0.and.ylk(n)>0.0d0)then
            phi=phi+twopi
         end if

      end if

      if(phi<dphisolby2)then
         jp=nphisol
      else

         jp=phi/dphisol

         if(phi-jp*dphisol>dphisolby2)then
            jp=jp+1
         end if

      end if

!      jlksol(2,n)=jp

      jp0=jp


      dist2=1000000.0d0

      jxget=0

      do j1=1,20

         jx=jx0+j1-10

         if(jx<1) cycle

         if(jx>nxsol) cycle

         do j2=1,20

            jp=jp0+j2-10

            if(jp<1) jp=jp+nphisol

            if(jp>nphisol) jp=jp-nphisol

            if(ylk(n)*ysurf(jp,jx)+zlk(n)*zsurf(jp,jx)<0.0) cycle

            dx=xlk(n)-xsurf(jp,jx)
            dy=ylk(n)-ysurf(jp,jx)
            dz=zlk(n)-zsurf(jp,jx)

            arg=dx*xnorsurf(jp,jx)+dy*ynorsurf(jp,jx)+dz*znorsurf(jp,jx)

            if(arg<0.0d0) cycle

            dx=dx-arg*xnorsurf(jp,jx)
            dy=dy-arg*ynorsurf(jp,jx)
            dz=dz-arg*znorsurf(jp,jx)

            d2=dx*dx+dy*dy+dz*dz

            if(dist2>d2)then

               dist2=d2

               jxget=jx

               jpget=jp

            end if

         end do

      end do

      if(jxget==0)then
         print*,'error: could not get indices',n,jx0,jp0
         stop
      end if

      jlksol(1,n)=jxget

      jlksol(2,n)=jpget


   end do

!$omp end do nowait
!$omp end parallel

   end subroutine


!=========================================================

   subroutine solidupdate(nxsol,nphisol,nmemb,nfa,nmyp2,nmyo2,nlk,jmbsol,jfasol,jmyp2sol,jmyo2sol,jlksol, &
               pi,delta,dxsol,dphisol,xmemb,ymemb,zmemb,xfa,yfa,zfa,xmyp2,ymyp2,zmyp2,xmyo2,ymyo2,zmyo2, &
                xlk,ylk,zlk,xwall,ywall,zwall,xnorwall,ynorwall,znorwall, &
                 xsurf,ysurf,zsurf,xnorsurf,ynorsurf,znorsurf,jsursol)



   implicit none

   integer,value:: nxsol,nphisol,nmemb,nfa,nmyp2,nmyo2,nlk
   integer n,jx,jp,jx0,jp0,j1,j2,jxget,jpget

   integer,allocatable,dimension(:,:)::jmbsol,jfasol,jmyp2sol,jmyo2sol,jlksol,jsursol

   double precision,value:: pi,delta,dxsol,dphisol
   double precision xmin,dxsolby2,piby2,dphisolby2,twopi,phi,arg
   double precision dx,dy,dz,d2,dist2

   double precision,allocatable,intent(in),dimension(:)::xmemb,ymemb,zmemb,xfa,yfa,zfa
   double precision,allocatable,intent(in),dimension(:)::xmyp2,ymyp2,zmyp2,xmyo2,ymyo2,zmyo2,xlk,ylk,zlk
   double precision,allocatable,intent(in),dimension(:,:)::xwall,ywall,zwall,xnorwall,ynorwall,znorwall
   double precision,allocatable,intent(in),dimension(:,:)::xsurf,ysurf,zsurf,xnorsurf,ynorsurf,znorsurf


   xmin=-(nxsol-1)/2*dxsol

   dxsolby2=0.5d0*dxsol

   piby2=0.5d0*pi

   dphisolby2=0.5d0*dphisol

   twopi=2*pi

!$omp parallel &
!$omp default(none) &
!$omp private(n,jx,phi,arg,jp,jx0,jp0,j1,j2,jxget,jpget,dx,dy,dz,d2,dist2) &
!$omp shared(nmemb,xmemb,xmin,dxsol,dxsolby2,nxsol,jsursol,ymemb,delta,zmemb,piby2) &
!$omp shared(dphisolby2,nphisol,pi,twopi,dphisol) &

!$omp shared(xwall,ywall,zwall,xnorwall,ynorwall,znorwall,jmbsol) &

!$omp shared(nfa,xfa,jfasol,yfa,zfa) &
!$omp shared(xsurf,ysurf,zsurf,xnorsurf,ynorsurf,znorsurf) &
!$omp shared(nmyp2,xmyp2,jmyp2sol,ymyp2,zmyp2) &
!$omp shared(nmyo2,xmyo2,jmyo2sol,ymyo2,zmyo2) &

!$omp shared(nlk,xlk,jlksol,ylk,zlk)


!  for membrane:

!$omp do schedule(guided,64)


   do n=1,nmemb

!     for membrane to interact with the ring:

      jx=1+(xmemb(n)-xmin)/dxsol

      if(xmemb(n)-xmin-(jx-1)*dxsol>dxsolby2)then
         jx=jx+1
      end if

      if(jx<1)then
         jx=1
      end if

      if(jx>nxsol)then
         jx=nxsol
      end if

      jsursol(1,n)=jx

      if(abs(ymemb(n))<delta)then
         if(zmemb(n)>0.0d0)then
            phi=piby2
         else
            phi=-piby2
         end if

      else

         arg=zmemb(n)/ymemb(n)

         phi=atan(arg)

         if(ymemb(n)<0.0d0)then
            phi=phi+pi
         end if

         if(arg<0.0d0.and.ymemb(n)>0.0d0)then
            phi=phi+twopi
         end if


      end if

      if(phi<dphisolby2)then
         jp=nphisol
      else

         jp=phi/dphisol

         if(phi-jp*dphisol>dphisolby2)then
            jp=jp+1
         end if

      end if

      jsursol(2,n)=jp

!------------------
!     for membrane to interact with the wall:

      jx0=jmbsol(1,n)

      jp0=jmbsol(2,n)

      jxget=jx0
      jpget=jp0

      dist2=10000000.0d0

      do j1=1,3

         jx=jx0-2+j1

         if(jx<1) cycle

         if(jx>nxsol) cycle

         do j2=1,3

            jp=jp0-2+j2

            if(jp<1) jp=jp+nphisol

            if(jp>nphisol) jp=jp-nphisol


            dx=xmemb(n)-xwall(jp,jx)
            dy=ymemb(n)-ywall(jp,jx)
            dz=zmemb(n)-zwall(jp,jx)

            arg=dx*xnorwall(jp,jx)+dy*ynorwall(jp,jx)+dz*znorwall(jp,jx)

            if(arg<0.0d0) cycle

            dx=dx-arg*xnorwall(jp,jx)
            dy=dy-arg*ynorwall(jp,jx)
            dz=dz-arg*znorwall(jp,jx)

            d2=dx*dx+dy*dy+dz*dz

            if(dist2>d2)then

               dist2=d2

               jxget=jx

               jpget=jp

            end if

         end do

      end do

      jmbsol(1,n)=jxget

      jmbsol(2,n)=jpget

   end do

!$omp enddo nowait


!----------------------------

!  for actin:


!$omp do schedule(guided,64)

   do n=1,nfa

      jx0=jfasol(1,n)

      jp0=jfasol(2,n)


      jxget=jx0
      jpget=jp0

      dist2=10000000.0d0

      do j1=1,3

         jx=jx0-2+j1

         if(jx<1) cycle

         if(jx>nxsol) cycle

         do j2=1,3

            jp=jp0-2+j2

            if(jp<1) jp=jp+nphisol

            if(jp>nphisol) jp=jp-nphisol

            dx=xfa(n)-xsurf(jp,jx)
            dy=yfa(n)-ysurf(jp,jx)
            dz=zfa(n)-zsurf(jp,jx)

            arg=dx*xnorsurf(jp,jx)+dy*ynorsurf(jp,jx)+dz*znorsurf(jp,jx)

            if(arg<0.0d0) cycle

            dx=dx-arg*xnorwall(jp,jx)
            dy=dy-arg*ynorwall(jp,jx)
            dz=dz-arg*znorwall(jp,jx)

            d2=dx*dx+dy*dy+dz*dz

            if(dist2>d2)then

               dist2=d2

               jxget=jx

               jpget=jp

            end if

         end do

      end do

      jfasol(1,n)=jxget

      jfasol(2,n)=jpget

   end do

!$omp enddo nowait

!----------------------------

!  for Myp2 myosin:


!$omp do schedule(guided,64)

   do n=1,nmyp2

      jx0=jmyp2sol(1,n)

      jp0=jmyp2sol(2,n)

      jxget=jx0
      jpget=jp0

      dist2=10000000.0d0

      do j1=1,3

         jx=jx0-2+j1

         if(jx<1) cycle

         if(jx>nxsol) cycle

         do j2=1,3

            jp=jp0-2+j2

            if(jp<1) jp=jp+nphisol

            if(jp>nphisol) jp=jp-nphisol

            dx=xmyp2(n)-xsurf(jp,jx)
            dy=ymyp2(n)-ysurf(jp,jx)
            dz=zmyp2(n)-zsurf(jp,jx)

            arg=dx*xnorsurf(jp,jx)+dy*ynorsurf(jp,jx)+dz*znorsurf(jp,jx)

            if(arg<0.0d0) cycle

            dx=dx-arg*xnorwall(jp,jx)
            dy=dy-arg*ynorwall(jp,jx)
            dz=dz-arg*znorwall(jp,jx)

            d2=dx*dx+dy*dy+dz*dz

            if(dist2>d2)then

               dist2=d2

               jxget=jx

               jpget=jp

            end if

         end do

      end do

      jmyp2sol(1,n)=jxget

      jmyp2sol(2,n)=jpget

   end do

!$omp enddo nowait

!----------------------------


!  for Myo2 myosin:

!$omp do schedule(guided,64)

   do n=1,nmyo2

      jx0=jmyo2sol(1,n)

      jp0=jmyo2sol(2,n)

      jxget=jx0
      jpget=jp0

      dist2=10000000.0d0

      do j1=1,3

         jx=jx0-2+j1

         if(jx<1) cycle

         if(jx>nxsol) cycle

         do j2=1,3

            jp=jp0-2+j2

            if(jp<1) jp=jp+nphisol

            if(jp>nphisol) jp=jp-nphisol

            dx=xmyo2(n)-xsurf(jp,jx)
            dy=ymyo2(n)-ysurf(jp,jx)
            dz=zmyo2(n)-zsurf(jp,jx)

            arg=dx*xnorsurf(jp,jx)+dy*ynorsurf(jp,jx)+dz*znorsurf(jp,jx)

            if(arg<0.0d0) cycle

            dx=dx-arg*xnorwall(jp,jx)
            dy=dy-arg*ynorwall(jp,jx)
            dz=dz-arg*znorwall(jp,jx)

            d2=dx*dx+dy*dy+dz*dz

            if(dist2>d2)then

               dist2=d2

               jxget=jx

               jpget=jp

            end if

         end do

      end do

      jmyo2sol(1,n)=jxget

      jmyo2sol(2,n)=jpget

   end do

!$omp enddo nowait

!----------------------------

!  for crosslinkers:


!$omp do schedule(guided,64)

   do n=1,nlk

      jx0=jlksol(1,n)

      jp0=jlksol(2,n)

      jxget=jx0
      jpget=jp0

      dist2=10000000.0d0

      do j1=1,3

         jx=jx0-2+j1

         if(jx<1) cycle

         if(jx>nxsol) cycle

         do j2=1,3

            jp=jp0-2+j2

            if(jp<1) jp=jp+nphisol

            if(jp>nphisol) jp=jp-nphisol

            dx=xlk(n)-xsurf(jp,jx)
            dy=ylk(n)-ysurf(jp,jx)
            dz=zlk(n)-zsurf(jp,jx)

            arg=dx*xnorsurf(jp,jx)+dy*ynorsurf(jp,jx)+dz*znorsurf(jp,jx)

            if(arg<0.0d0) cycle

            dx=dx-arg*xnorwall(jp,jx)
            dy=dy-arg*ynorwall(jp,jx)
            dz=dz-arg*znorwall(jp,jx)

            d2=dx*dx+dy*dy+dz*dz

            if(dist2>d2)then

               dist2=d2

               jxget=jx

               jpget=jp

            end if

         end do

      end do

      jlksol(1,n)=jxget

      jlksol(2,n)=jpget

   end do

!$omp enddo nowait
!$omp end parallel

   end subroutine

!=========================================================

   subroutine surfremod(nmemb,nthreads,nphisol,nxsol,nwall,jsursol,jmbsol,nsurf,wthick, &
               gap,dphisol,wallrate,xmemb,ymemb,zmemb,xsurf,ysurf,zsurf,xnorsurf, &
                ynorsurf,znorsurf,xwall,ywall,zwall,xnorwall,ynorwall,znorwall)


   implicit none

   integer,value::nmemb,nthreads,nphisol,nxsol,nwall
   integer n,jx,jp,jx1,jx2,jp1,jp2,tid,omp_get_thread_num,nsum

   integer,allocatable,intent(in),dimension(:,:)::jsursol,jmbsol
   integer,allocatable,dimension(:,:)::nsurf

   integer,allocatable,dimension(:,:,:)::nsurfcount,nwallcount
   integer,allocatable,dimension(:,:)::surfmark,wallmark

   double precision,value::wthick,gap,dphisol
   double precision wallrate

   double precision dx1,dy1,dz1,dx2,dy2,dz2,xn,yn,zn,invdist
   double precision dx,dy,dz,arg,d2,rad,dwallrate,rad2surf,radsurf,radwall,rad2wall,radmaxmemb

   double precision,allocatable,intent(in),dimension(:)::xmemb,ymemb,zmemb

   double precision,allocatable,dimension(:,:)::xsurf,ysurf,zsurf,xnorsurf,ynorsurf,znorsurf
   double precision,allocatable,dimension(:,:)::xwall,ywall,zwall,xnorwall,ynorwall,znorwall
   double precision,allocatable,dimension(:,:,:)::rsurftemp,dist,radmemb
   double precision,allocatable,dimension(:)::radsurfmax,radwallmax


   allocate(rsurftemp(nthreads,nphisol,nxsol))
   allocate(nsurfcount(nthreads,nphisol,nxsol),nwallcount(nthreads,nphisol,nxsol))
   allocate(dist(nthreads,nphisol,nxsol),radmemb(nthreads,nphisol,nxsol))

   nsurfcount=0

   rsurftemp=0.0d0

   nwallcount=0

   dist=10000.0d0

   radmemb=0.0d0


!$omp parallel &
!$omp default(none) &
!$omp private(n,jx,jp,tid,arg,dx,dy,dz) &
!$omp shared(nmemb,jsursol,nsurfcount,xmemb,ymemb,zmemb,rsurftemp) &
!$omp shared(jmbsol,nwallcount,xwall,ywall,zwall,xnorwall,ynorwall,znorwall,dist,radmemb)
!$omp do schedule(guided,64)

   do n=1,nmemb

      tid=omp_get_thread_num()+1
!tid=1

!     for membrane surface:

      jx=jsursol(1,n)

      jp=jsursol(2,n)

      nsurfcount(tid,jp,jx)=nsurfcount(tid,jp,jx)+1

      radmemb(tid,jp,jx)=sqrt(ymemb(n)*ymemb(n)+zmemb(n)*zmemb(n))

      rsurftemp(tid,jp,jx)=rsurftemp(tid,jp,jx)+radmemb(tid,jp,jx)!sqrt(ymemb(n)*ymemb(n)+zmemb(n)*zmemb(n))

!      xsurftemp(tid,jp,jx)=xsurftemp(tid,jp,jx)+xmemb(n)
!      ysurftemp(tid,jp,jx)=ysurftemp(tid,jp,jx)+ymemb(n)
!      zsurftemp(tid,jp,jx)=zsurftemp(tid,jp,jx)+zmemb(n)


!     for cell wall:

      jx=jmbsol(1,n)

      jp=jmbsol(2,n)

      nwallcount(tid,jp,jx)=nwallcount(tid,jp,jx)+1

      dx=xmemb(n)-xwall(jp,jx)
      dy=ymemb(n)-ywall(jp,jx)
      dz=zmemb(n)-zwall(jp,jx)

      arg=dx*xnorwall(jp,jx)+dy*ynorwall(jp,jx)+dz*znorwall(jp,jx)


      if(arg<dist(tid,jp,jx))then
         dist(tid,jp,jx)=arg
      end if


   end do

!$omp enddo nowait
!$omp end parallel



   allocate(surfmark(nphisol,nxsol),wallmark(nphisol,nxsol))

   surfmark=0

   wallmark=0

   allocate(radsurfmax(nthreads),radwallmax(nthreads))

   radsurfmax=0.0d0

   radwallmax=0.0d0

!$omp parallel &
!$omp default(none) &
!$omp private(jx,jp,nsum,tid,d2,rad,dwallrate,radmaxmemb) &
!$omp shared(nxsol,nphisol,nsurfcount,nthreads,ysurf,zsurf,rsurftemp) &
!$omp shared(nwallcount,dist,radmemb,ywall,zwall,wthick,gap,nwall,dphisol) &
!$omp shared(radsurfmax,surfmark,radwallmax,wallmark,nsurf) &
!$omp reduction(+:wallrate)

!$omp do

   do jx=1,nxsol

      tid=omp_get_thread_num()+1


      do jp=1,nphisol

!        for membrane surface:

         nsum=sum(nsurfcount(1:nthreads,jp,jx))

         nsurf(jp,jx)=nsum

         if(nsum>0)then

            rad=sum(rsurftemp(1:nthreads,jp,jx))/nsum

            ysurf(jp,jx)=rad*cos(jp*dphisol)
            zsurf(jp,jx)=rad*sin(jp*dphisol)

!            xsurf(jp,jx)=sum(xsurftemp(1:nthreads,jp,jx))/nsum
!            ysurf(jp,jx)=sum(ysurftemp(1:nthreads,jp,jx))/nsum
!            zsurf(jp,jx)=sum(zsurftemp(1:nthreads,jp,jx))/nsum

!            rad2=ysurf(jp,jx)*ysurf(jp,jx)+zsurf(jp,jx)*zsurf(jp,jx)

            if(rad>radsurfmax(tid))then

               radsurfmax(tid)=rad

            end if

         else

            surfmark(jp,jx)=1

         end if


!        for cell wall:

         nsum=sum(nwallcount(1:nthreads,jp,jx))

         if(nsum>0)then

            d2=minval(dist(1:nthreads,jp,jx))

            radmaxmemb=maxval(radmemb(1:nthreads,jp,jx))

!            rad2=ywall(jp,jx)*ywall(jp,jx)+zwall(jp,jx)*zwall(jp,jx)

            rad=sqrt(ywall(jp,jx)*ywall(jp,jx)+zwall(jp,jx)*zwall(jp,jx))

            if(d2>wthick+gap.and.rad>radmaxmemb)then

               dwallrate=0.1d0*gap

               wallrate=wallrate+dwallrate/nwall

               rad=rad-dwallrate

               ywall(jp,jx)=rad*cos(jp*dphisol)
               zwall(jp,jx)=rad*sin(jp*dphisol)

            end if

            if(rad>radwallmax(tid))then

               radwallmax(tid)=rad

            end if

         else

            wallmark(jp,jx)=1

         end if


      end do

   end do

!$omp enddo nowait
!$omp end parallel





   radsurf=maxval(radsurfmax(1:nthreads))

   rad2surf=radsurf*radsurf

   radwall=maxval(radwallmax(1:nthreads))

   rad2wall=radwall*radwall


!$omp parallel &
!$omp default(none) &
!$omp private(jx,jp) &
!$omp shared(nxsol,nphisol,surfmark,rad2surf,radsurf,ysurf,zsurf,dphisol) &
!$omp shared(wallmark,rad2wall,ywall,zwall,radwall)

!$omp do


   do jx=1,nxsol

      do jp=1,nphisol

!        for membrane surface:

         if(surfmark(jp,jx)==1)then

            if(rad2surf<ysurf(jp,jx)*ysurf(jp,jx)+zsurf(jp,jx)*zsurf(jp,jx))then

               ysurf(jp,jx)=radsurf*cos(jp*dphisol)
               zsurf(jp,jx)=radsurf*sin(jp*dphisol)

            end if

         end if

!        for cell wall:

         if(wallmark(jp,jx)==1)then

            if(rad2wall<ywall(jp,jx)*ywall(jp,jx)+zwall(jp,jx)*zwall(jp,jx))then

               ywall(jp,jx)=radwall*cos(jp*dphisol)
               zwall(jp,jx)=radwall*sin(jp*dphisol)

            end if

         end if

      end do

   end do

!$omp enddo nowait
!$omp end parallel




!$omp parallel &
!$omp default(none) &
!$omp private(jx,jp,jx1,jx2,jp1,jp2,dx1,dy1,dz1,dx2,dy2,dz2,xn,yn,zn,invdist) &
!$omp shared(nxsol,nphisol,xsurf,ysurf,zsurf,xnorsurf,ynorsurf,znorsurf) &
!$omp shared(xwall,ywall,zwall,xnorwall,ynorwall,znorwall)

!$omp do

   do jx=1,nxsol

      if(jx==1)then
         jx1=1
         jx2=jx+1
      elseif(jx==nxsol)then
         jx2=nxsol
         jx1=jx-1
      else
         jx1=jx-1
         jx2=jx+1
      end if


      do jp=1,nphisol



         if(jp==1)then
            jp1=nphisol
            jp2=jp+1
         elseif(jp==nphisol)then
            jp2=1
            jp1=jp-1
         else
            jp1=jp-1
            jp2=jp+1
         end if


!if(jx==0.or.jx1==0.or.jx2==0)then
!print*,'jx',jx,jx1,jx2
!end if

!if(jp==0.or.jp1==0.or.jp2==0)then
!print*,'jp',jp,jp1,jp2
!end if

!        for membrane surface:

         dx1=xsurf(jp,jx2)-xsurf(jp,jx1)
         dy1=ysurf(jp,jx2)-ysurf(jp,jx1)
         dz1=zsurf(jp,jx2)-zsurf(jp,jx1)

         dx2=xsurf(jp2,jx)-xsurf(jp1,jx)
         dy2=ysurf(jp2,jx)-ysurf(jp1,jx)
         dz2=zsurf(jp2,jx)-zsurf(jp1,jx)

         xn=dy1*dz2-dy2*dz1
         yn=dz1*dx2-dz2*dx1
         zn=dx1*dy2-dx2*dy1

         invdist=1.0d0/sqrt(xn*xn+yn*yn+zn*zn)

         xnorsurf(jp,jx)=xn*invdist
         ynorsurf(jp,jx)=yn*invdist
         znorsurf(jp,jx)=zn*invdist


!        for cell wall:

         dx1=xwall(jp,jx2)-xwall(jp,jx1)
         dy1=ywall(jp,jx2)-ywall(jp,jx1)
         dz1=zwall(jp,jx2)-zwall(jp,jx1)

         dx2=xwall(jp2,jx)-xwall(jp1,jx)
         dy2=ywall(jp2,jx)-ywall(jp1,jx)
         dz2=zwall(jp2,jx)-zwall(jp1,jx)

         xn=dy1*dz2-dy2*dz1
         yn=dz1*dx2-dz2*dx1
         zn=dx1*dy2-dx2*dy1

         invdist=1.0d0/sqrt(xn*xn+yn*yn+zn*zn)

         xnorwall(jp,jx)=xn*invdist
         ynorwall(jp,jx)=yn*invdist
         znorwall(jp,jx)=zn*invdist




      end do

   end do

!$omp enddo nowait
!$omp end parallel

   deallocate(nsurfcount,nwallcount,rsurftemp,dist,radsurfmax,surfmark,radwallmax,wallmark)

   end subroutine



!=========================================================

   subroutine constraints(nfa,nmemb,nmyp2,nmyo2,nlk,jupdate,nxsol,nphisol,nthreads, &
               fadist,myp2dist,myo2dist,lkdist,jmbsol,jsursol,jfasol,jmyp2sol,jmyo2sol,jlksol,nsurf, &
                kwall,lsqueez,wthick,k_mem,l_mem,xwall,ywall,zwall,xnorwall, &
                 ynorwall,znorwall,xmemb,ymemb,zmemb,xboundmin,xboundmax,xrmin,xrmax,xfa,yfa,zfa, &
                  xmyp2,ymyp2,zmyp2,xmyo2,ymyo2,zmyo2,xlk,ylk,zlk,xsurf,ysurf,zsurf,xnorsurf,ynorsurf, &
                   znorsurf,fxmemb,fymemb,fzmemb,fxmyp2,fymyp2,fzmyp2,fxmyo2,fymyo2,fzmyo2, &
                    fxlk,fylk,fzlk,fxfa,fyfa,fzfa)




   implicit none

   integer,value:: nfa,nmemb,nmyp2,nmyo2,nlk,jupdate,nxsol,nphisol,nthreads
   integer n,jx,jp,jm,tid,omp_get_thread_num

!   integer,allocatable,intent(in),dimension(:)::a2mem
   integer,allocatable,dimension(:)::fadist,myp2dist,myo2dist,lkdist
   integer,allocatable,intent(in),dimension(:,:)::jmbsol,jsursol,jfasol,jmyp2sol,jmyo2sol,jlksol,nsurf

   double precision,value::kwall,lsqueez,wthick,k_mem,l_mem,xboundmin,xboundmax,xrmin,xrmax

   double precision dx,dy,dz,dist,f,d2,l2max,rad,l_memby2,dfx,dfy,dfz
   double precision xn,yn,zn

   double precision,allocatable,intent(in),dimension(:,:)::xwall,ywall,zwall,xnorwall,ynorwall,znorwall
   double precision,allocatable,intent(in),dimension(:)::xmemb,ymemb,zmemb
   double precision,allocatable,intent(in),dimension(:)::xfa,yfa,zfa,xmyo2,ymyo2,zmyo2
   double precision,allocatable,intent(in),dimension(:)::xmyp2,ymyp2,zmyp2,xlk,ylk,zlk
   double precision,allocatable,intent(in),dimension(:,:)::xsurf,ysurf,zsurf,xnorsurf,ynorsurf,znorsurf
   double precision,allocatable,dimension(:)::fxmemb,fymemb,fzmemb,fxmyp2,fymyp2,fzmyp2,fxlk,fylk,fzlk
   double precision,allocatable,dimension(:)::fxfa,fyfa,fzfa,fxmyo2,fymyo2,fzmyo2

   double precision,allocatable,dimension(:,:,:)::ftem
   double precision,allocatable,dimension(:,:)::fxsurf,fysurf,fzsurf



   l_memby2=0.5d0*l_mem

!$omp parallel &
!$omp default(none) &
!$omp private(n,jx,jp,xn,yn,zn,dx,dy,dz,dist,f) &
!$omp shared(nmemb,jmbsol,xwall,ywall,zwall,xnorwall,ynorwall,znorwall) &
!$omp shared(kwall,lsqueez,wthick) &
!$omp shared(xmemb,ymemb,zmemb,fxmemb,fymemb,fzmemb) &
!$omp shared(jsursol,xboundmin,xboundmax) 

!$omp do schedule(guided,64)

   do n=1,nmemb

      jx=jmbsol(1,n)

      jp=jmbsol(2,n)

!     normal vector of the wall:

      xn=xnorwall(jp,jx)
      yn=ynorwall(jp,jx)
      zn=znorwall(jp,jx)


      dx=xmemb(n)-xwall(jp,jx)
      dy=ymemb(n)-ywall(jp,jx)
      dz=zmemb(n)-zwall(jp,jx)


      dist=dx*xn+dy*yn+dz*zn-wthick

      if(dist<0.0d0)then


         f=kwall*dist*dist

         fxmemb(n)=fxmemb(n)+f*xn
         fymemb(n)=fymemb(n)+f*yn
         fzmemb(n)=fzmemb(n)+f*zn

      elseif(dist<lsqueez)then

         f=-kwall*dist

         fxmemb(n)=fxmemb(n)+f*xn
         fymemb(n)=fymemb(n)+f*yn
         fzmemb(n)=fzmemb(n)+f*zn




      end if



!     blocking from the x boundaries:

!      jp=jsursol(2,n)

      if(xmemb(n)<xboundmin)then
         fxmemb(n)=fxmemb(n)+xboundmin-xmemb(n)
      else if(xmemb(n)>xboundmax)then
         fxmemb(n)=fxmemb(n)+xboundmax-xmemb(n)
      end if


   end do

!$omp enddo nowait
!$omp end parallel

!----------------------------------------------------

!  constraint on the ring from membrane:

   allocate(ftem(nthreads,nphisol,nxsol))
   ftem=0.0d0

!   l2max=ltether*ltether


!$omp parallel &
!$omp default(none) &
!$omp private(n,jx,jp,jm,dx,dy,dz,xn,yn,zn,dist,f,tid,d2,dfx,dfy,dfz) &
!$omp shared(nfa,jfasol,fadist,jupdate,ftem) &
!$omp shared(k_mem,l_mem,l_memby2) &
!$omp shared(xmemb,ymemb,zmemb,fxmemb,fymemb,fzmemb,xrmin,xrmax) &
!$omp shared(xfa,yfa,zfa,fxfa,fyfa,fzfa) &
!$omp shared(xsurf,ysurf,zsurf,xnorsurf,ynorsurf,znorsurf) &

!$omp shared(nmyp2,jmyp2sol,myp2dist) &
!$omp shared(xmyp2,ymyp2,zmyp2,fxmyp2,fymyp2,fzmyp2) &

!$omp shared(nmyo2,jmyo2sol,myo2dist) &
!$omp shared(xmyo2,ymyo2,zmyo2,fxmyo2,fymyo2,fzmyo2) &

!$omp shared(nlk,jlksol,lkdist) &
!$omp shared(xlk,ylk,zlk,fxlk,fylk,fzlk)


!  constrain on actin:

!$omp do schedule(guided,64)

   do n=1,nfa


      if(fadist(n)==1.or.jupdate==1)then

         jx=jfasol(1,n)

         jp=jfasol(2,n)

!         dist=sqrt(yfa(n)*yfa(n)+zfa(n)*zfa(n))

         dx=xfa(n)-xsurf(jp,jx)
         dy=yfa(n)-ysurf(jp,jx)
         dz=zfa(n)-zsurf(jp,jx)

         xn=xnorsurf(jp,jx)
         yn=ynorsurf(jp,jx)
         zn=znorsurf(jp,jx)

         dist=dx*xn+dy*yn+dz*zn

         if(dist<l_mem)then

            if(dist<0.0d0) dist=0.0d0

            f=k_mem*(l_mem-dist)/(dist+1.0d0)

            fxfa(n)=fxfa(n)+f*xn
            fyfa(n)=fyfa(n)+f*yn
            fzfa(n)=fzfa(n)+f*zn

            fadist(n)=1

            tid=omp_get_thread_num()+1

            ftem(tid,jp,jx)=ftem(tid,jp,jx)-f

         else

            if(dist<l_mem+l_memby2)then
               fadist(n)=1
            else
               fadist(n)=0
            end if

         end if



      end if

!     blocking from the x boundaries:

      if(xfa(n)<xrmin)then
         fxfa(n)=fxfa(n)+xrmin-xfa(n)
      else if(xfa(n)>xrmax)then
         fxfa(n)=fxfa(n)+xrmax-xfa(n)
      end if



!      if(a2mem(n)==0)then
!         cycle
!      end if

!      jm=a2mem(n)

!      dx=xfa(n)-xmemb(jm)
!      dy=yfa(n)-ymemb(jm)
!      dz=zfa(n)-zmemb(jm)

!      d2=dx*dx+dy*dy+dz*dz

!      if(d2>l2max)then
!         dist=sqrt(d2)


!         f=-ktether*(dist-ltether)/dist

!         dfx=f*dx
!         dfy=f*dy
!         dfz=f*dz

!         fxfa(n)=fxfa(n)+dfx!*invdist
!         fyfa(n)=fyfa(n)+dfy!*invdist
!         fzfa(n)=fzfa(n)+dfz!*invdist

!         fxmemb(jm)=fxmemb(jm)-dfx
!         fymemb(jm)=fymemb(jm)-dfy
!         fzmemb(jm)=fzmemb(jm)-dfz


!      end if

   end do

!$omp enddo nowait

!---------------------------------


!  constrain on Myp2:


!$omp do schedule(guided,64)


   do n=1,nmyp2


      jp=jmyp2sol(2,n)

      if(myp2dist(n)==1.or.jupdate==1)then

         jx=jmyp2sol(1,n)

!         jp=jmyp2sol(2,n)


         dx=xmyp2(n)-xsurf(jp,jx)
         dy=ymyp2(n)-ysurf(jp,jx)
         dz=zmyp2(n)-zsurf(jp,jx)

         xn=xnorsurf(jp,jx)
         yn=ynorsurf(jp,jx)
         zn=znorsurf(jp,jx)


         dist=dx*xn+dy*yn+dz*zn

         if(dist<l_mem)then

            if(dist<0.0d0) dist=0.0d0

            f=k_mem*(l_mem-dist)/(dist+1.0d0)

            fxmyp2(n)=fxmyp2(n)+f*xn
            fymyp2(n)=fymyp2(n)+f*yn
            fzmyp2(n)=fzmyp2(n)+f*zn

            myp2dist(n)=1

            tid=omp_get_thread_num()+1

            ftem(tid,jp,jx)=ftem(tid,jp,jx)-f

         else

            if(dist<l_mem+l_memby2)then
               myp2dist(n)=1
            else
               myp2dist(n)=0
            end if

         end if



      end if

!     blocking from the x boundaries:


      if(xmyp2(n)<xrmin)then
         fxmyp2(n)=fxmyp2(n)+xrmin-xmyp2(n)
      else if(xmyp2(n)>xrmax)then
         fxmyp2(n)=fxmyp2(n)+xrmax-xmyp2(n)
      end if

   end do

!$omp enddo nowait


!---------------------------------

!  constrain on Myo2:

!$omp do schedule(guided,64)


   do n=1,nmyo2


      jp=jmyo2sol(2,n)

      if(myo2dist(n)==1.or.jupdate==1)then

         jx=jmyo2sol(1,n)

         dx=xmyo2(n)-xsurf(jp,jx)
         dy=ymyo2(n)-ysurf(jp,jx)
         dz=zmyo2(n)-zsurf(jp,jx)

         xn=xnorsurf(jp,jx)
         yn=ynorsurf(jp,jx)
         zn=znorsurf(jp,jx)


         dist=dx*xn+dy*yn+dz*zn

         if(dist<l_mem)then

            if(dist<0.0d0) dist=0.0d0

            f=k_mem*(l_mem-dist)/(dist+1.0d0)

            fxmyo2(n)=fxmyo2(n)+f*xn
            fymyo2(n)=fymyo2(n)+f*yn
            fzmyo2(n)=fzmyo2(n)+f*zn

            myo2dist(n)=1

            tid=omp_get_thread_num()+1

            ftem(tid,jp,jx)=ftem(tid,jp,jx)-f

         else

            if(dist<l_mem+l_memby2)then
               myo2dist(n)=1
            else
               myo2dist(n)=0
            end if

         end if

      end if

!     blocking from the x boundaries:


      if(xmyo2(n)<xrmin)then
         fxmyo2(n)=fxmyo2(n)+xrmin-xmyo2(n)
      else if(xmyo2(n)>xrmax)then
         fxmyo2(n)=fxmyo2(n)+xrmax-xmyo2(n)
      end if

   end do

!$omp enddo nowait


!---------------------------------

!  constraint on crosslinkers:

!$omp do schedule(guided,32)

   do n=1,nlk

      jp=jlksol(2,n)

      if(lkdist(n)==1.or.jupdate==1)then

         jx=jlksol(1,n)

!         jp=jlksol(2,n)

!         dist=sqrt(ylk(n)*ylk(n)+zlk(n)*zlk(n))

         dx=xlk(n)-xsurf(jp,jx)
         dy=ylk(n)-ysurf(jp,jx)
         dz=zlk(n)-zsurf(jp,jx)

         xn=xnorsurf(jp,jx)
         yn=ynorsurf(jp,jx)
         zn=znorsurf(jp,jx)


         dist=dx*xn+dy*yn+dz*zn

         if(dist<l_mem)then

            if(dist<0.0d0) dist=0.0d0

            f=k_mem*(l_mem-dist)/(dist+1.0d0)

            fxlk(n)=fxlk(n)+f*xn
            fylk(n)=fylk(n)+f*yn
            fzlk(n)=fzlk(n)+f*zn

            lkdist(n)=1

            tid=omp_get_thread_num()+1

            ftem(tid,jp,jx)=ftem(tid,jp,jx)-f

         else

            if(dist<l_mem+l_memby2)then
               lkdist(n)=1
            else
               lkdist(n)=0
            end if

         end if




      end if

!     blocking from the x boundaries:

      if(xlk(n)<xrmin)then
         fxlk(n)=fxlk(n)+xrmin-xlk(n)
      else if(xlk(n)>xrmax)then
         fxlk(n)=fxlk(n)+xrmax-xlk(n)
      end if

   end do

!$omp enddo nowait
!$omp end parallel




   allocate(fxsurf(nphisol,nxsol),fysurf(nphisol,nxsol),fzsurf(nphisol,nxsol))


!$omp parallel &
!$omp default(none) &
!$omp private(jx,jp,f) &
!$omp shared(nxsol,nphisol,ftem,nthreads,nsurf,fxsurf,fysurf,fzsurf,xnorsurf,ynorsurf,znorsurf)

!$omp do

   do jx=1,nxsol

      do jp=1,nphisol

         f=sum(ftem(1:nthreads,jp,jx))/nsurf(jp,jx)

         fxsurf(jp,jx)=f*xnorsurf(jp,jx)
         fysurf(jp,jx)=f*ynorsurf(jp,jx)
         fzsurf(jp,jx)=f*znorsurf(jp,jx)

      end do

   end do

!$omp enddo nowait
!$omp end parallel


!$omp parallel &
!$omp default(none) &
!$omp private(n,jx,jp) &
!$omp shared(nmemb,jsursol,fxmemb,fymemb,fzmemb,fxsurf,fysurf,fzsurf)

!$omp do

   do n=1,nmemb

      jx=jsursol(1,n)
      jp=jsursol(2,n)

      fxmemb(n)=fxmemb(n)+fxsurf(jp,jx)
      fymemb(n)=fymemb(n)+fysurf(jp,jx)
      fzmemb(n)=fzmemb(n)+fzsurf(jp,jx)

   end do

!$omp enddo nowait
!$omp end parallel


   deallocate(fxsurf,fysurf,fzsurf,ftem)


   end subroutine

!=========================================================

!=========================================================

   subroutine nonbond(npair_myp2ac,pair_myp2ac,npair_lkac,pair_lkac,npair_myp2lk,pair_myp2lk, &
               npair_ac2,pair_ac2,npair_myp2,pair_myp2,npair_lk2,pair_lk2,npair_myo2ac,pair_myo2ac, &
                npair_myo2lk,pair_myo2lk,npair_myo2,pair_myo2,npair_mymy,pair_mymy, &
                 npair_mb2,pair_mb2,pairpart,boundtyp,xfa,yfa,zfa,fxfa,fyfa,fzfa, &
                  xmyp2,ymyp2,zmyp2,fxmyp2,fymyp2,fzmyp2,xmyo2,ymyo2,zmyo2,fxmyo2,fymyo2,fzmyo2, &
                   xlk,ylk,zlk,fxlk,fylk,fzlk,xmemb,ymemb,zmemb,fxmemb,fymemb,fzmemb, &
                     rvdw,kvdw,r_off,r_on2,r_off2,fvdwmax,kpair,l_pair,l_mem,kmemb,shift)


   implicit none

   integer,value::npair_myp2ac,npair_lkac,npair_myp2lk,npair_ac2,npair_myp2,npair_lk2
   integer,value::npair_myo2ac,npair_myo2lk,npair_myo2,npair_mymy,npair_mb2

   integer n,ja,ja1,ja2,jm,jm1,jm2,jl,jl1,jl2,n1,n2,n3,n4,jp

   integer,allocatable,intent(in),dimension(:,:)::pair_myp2ac,pair_lkac,pair_myp2lk,pair_ac2,pair_myp2,pair_lk2
   integer,allocatable,intent(in),dimension(:,:)::pair_myo2ac,pair_myo2lk,pair_myo2,pair_mymy,pair_mb2,pairpart

   integer,allocatable,intent(in),dimension(:)::boundtyp

   double precision,value::rvdw,kvdw,r_off,r_on2,r_off2,fvdwmax

   double precision,value::kpair,l_pair,l_mem,kmemb,shift

   double precision dx,dy,dz,d2,dist,invdist,f,dfx,dfy,dfz,ratio
   double precision invd2,dx34,dy34,dz34,proj,xn,yn,zn

   double precision,allocatable,intent(in),dimension(:)::xfa,yfa,zfa,xmyp2,ymyp2,zmyp2,xlk,ylk,zlk
   double precision,allocatable,intent(in),dimension(:)::xmemb,ymemb,zmemb,xmyo2,ymyo2,zmyo2
   double precision,allocatable,dimension(:)::fxmemb,fymemb,fzmemb,fxmyo2,fymyo2,fzmyo2
   double precision,allocatable,dimension(:)::fxfa,fyfa,fzfa,fxmyp2,fymyp2,fzmyp2,fxlk,fylk,fzlk


!   halfbound=0.5d0*shift


!$omp parallel &
!$omp default(none) &

!$omp shared(npair_myp2ac,npair_lkac,npair_myp2lk,npair_ac2,npair_myp2,npair_lk2,npair_mb2) &
!$omp shared(npair_myo2ac,npair_myo2lk,npair_myo2,npair_mymy) &
!$omp private(n,ja,ja1,ja2,jm,jm1,jm2,jl,jl1,jl2,jp) &
!$omp shared(pair_myp2ac,pair_lkac,pair_myp2lk,pair_ac2,pair_myp2,pair_lk2,pair_mb2) &
!$omp shared(pair_myo2ac,pair_myo2lk,pair_myo2,pair_mymy) &
!$omp shared(rvdw,kvdw,r_off,r_on2,r_off2,fvdwmax) &
!$omp shared(kpair,l_pair,l_mem) &
!$omp private(dx,dy,dz,d2,dist,invdist,f,dfx,dfy,dfz,ratio) &
!$omp shared(xfa,yfa,zfa,xmyp2,ymyp2,zmyp2,xmyo2,ymyo2,zmyo2,xlk,ylk,zlk) &
!$omp shared(xmemb,ymemb,zmemb,fxmemb,fymemb,fzmemb) &
!$omp shared(fxfa,fyfa,fzfa,fxmyp2,fymyp2,fzmyp2,fxmyo2,fymyo2,fzmyo2,fxlk,fylk,fzlk) &
!$omp private(n1,n2,n3,n4,invd2,dx34,dy34,dz34,proj,xn,yn,zn) &
!$omp shared(pairpart,kmemb,boundtyp,shift)


!  actin-Myp2 forces:

!$omp do schedule(guided,64)

   do n=1,npair_myp2ac

      jm=pair_myp2ac(1,n)
      ja=pair_myp2ac(2,n)


      dx=xfa(ja)-xmyp2(jm)
      dy=yfa(ja)-ymyp2(jm)
      dz=zfa(ja)-zmyp2(jm)

      d2=dx*dx+dy*dy+dz*dz

      if(d2<r_on2)then

         dist=sqrt(d2)

!         invdist=1/dist
         f=fvdwmax/dist

         dfx=f*dx!*invdist
         dfy=f*dy!*invdist
         dfz=f*dz!*invdist

!$omp atomic
         fxfa(ja)=fxfa(ja)+dfx
!$omp atomic
         fyfa(ja)=fyfa(ja)+dfy
!$omp atomic
         fzfa(ja)=fzfa(ja)+dfz

!$omp atomic
         fxmyp2(jm)=fxmyp2(jm)-dfx
!$omp atomic
         fymyp2(jm)=fymyp2(jm)-dfy
!$omp atomic
         fzmyp2(jm)=fzmyp2(jm)-dfz

      else if(d2<r_off2)then

         dist=sqrt(d2)

!         invdist=1/dist

         ratio=(r_off-dist)/(dist-rvdw)

         f=kvdw*ratio*ratio/dist!(r_off-dist)*(r_off-dist)/(dist-rvdw)/(dist-rvdw)

         dfx=f*dx!*invdist
         dfy=f*dy!*invdist
         dfz=f*dz!*invdist

!$omp atomic
         fxfa(ja)=fxfa(ja)+dfx
!$omp atomic
         fyfa(ja)=fyfa(ja)+dfy
!$omp atomic
         fzfa(ja)=fzfa(ja)+dfz

!$omp atomic
         fxmyp2(jm)=fxmyp2(jm)-dfx
!$omp atomic
         fymyp2(jm)=fymyp2(jm)-dfy
!$omp atomic
         fzmyp2(jm)=fzmyp2(jm)-dfz

      end if

   end do

!$omp enddo nowait


!-------------------------------

!  actin-Myo2 forces:


!$omp do schedule(guided,64)

   do n=1,npair_myo2ac

      jm=pair_myo2ac(1,n)
      ja=pair_myo2ac(2,n)


      dx=xfa(ja)-xmyo2(jm)
      dy=yfa(ja)-ymyo2(jm)
      dz=zfa(ja)-zmyo2(jm)

      d2=dx*dx+dy*dy+dz*dz

      if(d2<r_on2)then

         dist=sqrt(d2)

         f=fvdwmax/dist

         dfx=f*dx!*invdist
         dfy=f*dy!*invdist
         dfz=f*dz!*invdist

!$omp atomic
         fxfa(ja)=fxfa(ja)+dfx
!$omp atomic
         fyfa(ja)=fyfa(ja)+dfy
!$omp atomic
         fzfa(ja)=fzfa(ja)+dfz

!$omp atomic
         fxmyo2(jm)=fxmyo2(jm)-dfx
!$omp atomic
         fymyo2(jm)=fymyo2(jm)-dfy
!$omp atomic
         fzmyo2(jm)=fzmyo2(jm)-dfz

      else if(d2<r_off2)then

         dist=sqrt(d2)

!         invdist=1/dist

         ratio=(r_off-dist)/(dist-rvdw)

         f=kvdw*ratio*ratio/dist!(r_off-dist)*(r_off-dist)/(dist-rvdw)/(dist-rvdw)

         dfx=f*dx!*invdist
         dfy=f*dy!*invdist
         dfz=f*dz!*invdist

!$omp atomic
         fxfa(ja)=fxfa(ja)+dfx
!$omp atomic
         fyfa(ja)=fyfa(ja)+dfy
!$omp atomic
         fzfa(ja)=fzfa(ja)+dfz

!$omp atomic
         fxmyo2(jm)=fxmyo2(jm)-dfx
!$omp atomic
         fymyo2(jm)=fymyo2(jm)-dfy
!$omp atomic
         fzmyo2(jm)=fzmyo2(jm)-dfz

      end if

   end do

!$omp enddo nowait

!-------------------------------


!  actin-crosslinker forces:

!$omp do schedule(guided,64)

   do n=1,npair_lkac

      jl=pair_lkac(1,n)
      ja=pair_lkac(2,n)


      dx=xfa(ja)-xlk(jl)
      dy=yfa(ja)-ylk(jl)
      dz=zfa(ja)-zlk(jl)

      d2=dx*dx+dy*dy+dz*dz

      if(d2<r_on2)then

         dist=sqrt(d2)

!         invdist=1/dist
         f=fvdwmax/dist

         dfx=f*dx!*invdist
         dfy=f*dy!*invdist
         dfz=f*dz!*invdist

!$omp atomic
         fxfa(ja)=fxfa(ja)+dfx
!$omp atomic
         fyfa(ja)=fyfa(ja)+dfy
!$omp atomic
         fzfa(ja)=fzfa(ja)+dfz

!$omp atomic
         fxlk(jl)=fxlk(jl)-dfx
!$omp atomic
         fylk(jl)=fylk(jl)-dfy
!$omp atomic
         fzlk(jl)=fzlk(jl)-dfz

      else if(d2<r_off2)then

         dist=sqrt(d2)

!         invdist=1/dist

         ratio=(r_off-dist)/(dist-rvdw)

         f=kvdw*ratio*ratio/dist!(r_off-dist)*(r_off-dist)/(dist-rvdw)/(dist-rvdw)

         dfx=f*dx!*invdist
         dfy=f*dy!*invdist
         dfz=f*dz!*invdist

!$omp atomic
         fxfa(ja)=fxfa(ja)+dfx
!$omp atomic
         fyfa(ja)=fyfa(ja)+dfy
!$omp atomic
         fzfa(ja)=fzfa(ja)+dfz

!$omp atomic
         fxlk(jl)=fxlk(jl)-dfx
!$omp atomic
         fylk(jl)=fylk(jl)-dfy
!$omp atomic
         fzlk(jl)=fzlk(jl)-dfz

      end if

   end do

!$omp enddo nowait


!-------------------------------


!  Myp2-crosslinker forces:

!$omp do schedule(guided,64)

   do n=1,npair_myp2lk

      jm=pair_myp2lk(1,n)
      jl=pair_myp2lk(2,n)

      dx=xmyp2(jm)-xlk(jl)
      dy=ymyp2(jm)-ylk(jl)
      dz=zmyp2(jm)-zlk(jl)

      d2=dx*dx+dy*dy+dz*dz

      if(d2<r_on2)then

         dist=sqrt(d2)

         f=fvdwmax/dist

         dfx=f*dx!*invdist
         dfy=f*dy!*invdist
         dfz=f*dz!*invdist

!$omp atomic
         fxmyp2(jm)=fxmyp2(jm)+dfx
!$omp atomic
         fymyp2(jm)=fymyp2(jm)+dfy
!$omp atomic
         fzmyp2(jm)=fzmyp2(jm)+dfz

!$omp atomic
         fxlk(jl)=fxlk(jl)-dfx
!$omp atomic
         fylk(jl)=fylk(jl)-dfy
!$omp atomic
         fzlk(jl)=fzlk(jl)-dfz

      else if(d2<r_off2)then

         dist=sqrt(d2)

!         invdist=1/dist

         ratio=(r_off-dist)/(dist-rvdw)

         f=kvdw*ratio*ratio/dist!(r_off-dist)*(r_off-dist)/(dist-rvdw)/(dist-rvdw)

         dfx=f*dx!*invdist
         dfy=f*dy!*invdist
         dfz=f*dz!*invdist

!$omp atomic
         fxmyp2(jm)=fxmyp2(jm)+dfx
!$omp atomic
         fymyp2(jm)=fymyp2(jm)+dfy
!$omp atomic
         fzmyp2(jm)=fzmyp2(jm)+dfz

!$omp atomic
         fxlk(jl)=fxlk(jl)-dfx
!$omp atomic
         fylk(jl)=fylk(jl)-dfy
!$omp atomic
         fzlk(jl)=fzlk(jl)-dfz

      end if

   end do

!$omp enddo nowait


!-------------------------------

!  Myo2-crosslinker forces:


!$omp do schedule(guided,64)

   do n=1,npair_myo2lk

      jm=pair_myo2lk(1,n)
      jl=pair_myo2lk(2,n)

      dx=xmyo2(jm)-xlk(jl)
      dy=ymyo2(jm)-ylk(jl)
      dz=zmyo2(jm)-zlk(jl)

      d2=dx*dx+dy*dy+dz*dz

      if(d2<r_on2)then

         dist=sqrt(d2)

         f=fvdwmax/dist

         dfx=f*dx!*invdist
         dfy=f*dy!*invdist
         dfz=f*dz!*invdist

!$omp atomic
         fxmyo2(jm)=fxmyo2(jm)+dfx
!$omp atomic
         fymyo2(jm)=fymyo2(jm)+dfy
!$omp atomic
         fzmyo2(jm)=fzmyo2(jm)+dfz

!$omp atomic
         fxlk(jl)=fxlk(jl)-dfx
!$omp atomic
         fylk(jl)=fylk(jl)-dfy
!$omp atomic
         fzlk(jl)=fzlk(jl)-dfz

      else if(d2<r_off2)then

         dist=sqrt(d2)

!         invdist=1/dist

         ratio=(r_off-dist)/(dist-rvdw)

         f=kvdw*ratio*ratio/dist!(r_off-dist)*(r_off-dist)/(dist-rvdw)/(dist-rvdw)

         dfx=f*dx!*invdist
         dfy=f*dy!*invdist
         dfz=f*dz!*invdist

!$omp atomic
         fxmyo2(jm)=fxmyo2(jm)+dfx
!$omp atomic
         fymyo2(jm)=fymyo2(jm)+dfy
!$omp atomic
         fzmyo2(jm)=fzmyo2(jm)+dfz

!$omp atomic
         fxlk(jl)=fxlk(jl)-dfx
!$omp atomic
         fylk(jl)=fylk(jl)-dfy
!$omp atomic
         fzlk(jl)=fzlk(jl)-dfz

      end if

   end do

!$omp enddo nowait


!-------------------------------


!$omp do schedule(guided,64)

!  actin-actin forces:

   do n=1,npair_ac2

      ja1=pair_ac2(1,n)
      ja2=pair_ac2(2,n)

!      if(atyp(ja1)==0.or.atyp(ja2)==0)then
!         cycle
!      end if

      dx=xfa(ja1)-xfa(ja2)
      dy=yfa(ja1)-yfa(ja2)
      dz=zfa(ja1)-zfa(ja2)

      d2=dx*dx+dy*dy+dz*dz

      if(d2<r_on2)then

         dist=sqrt(d2)

         f=fvdwmax/dist

         dfx=f*dx!*invdist
         dfy=f*dy!*invdist
         dfz=f*dz!*invdist

!$omp atomic
         fxfa(ja1)=fxfa(ja1)+dfx
!$omp atomic
         fyfa(ja1)=fyfa(ja1)+dfy
!$omp atomic
         fzfa(ja1)=fzfa(ja1)+dfz

!$omp atomic
         fxfa(ja2)=fxfa(ja2)-dfx
!$omp atomic
         fyfa(ja2)=fyfa(ja2)-dfy
!$omp atomic
         fzfa(ja2)=fzfa(ja2)-dfz

      else if(d2<r_off2)then

         dist=sqrt(d2)

!         invdist=1/dist

         ratio=(r_off-dist)/(dist-rvdw)

         f=kvdw*ratio*ratio/dist!(r_off-dist)*(r_off-dist)/(dist-rvdw)/(dist-rvdw)

         dfx=f*dx!*invdist
         dfy=f*dy!*invdist
         dfz=f*dz!*invdist

!$omp atomic
         fxfa(ja1)=fxfa(ja1)+dfx
!$omp atomic
         fyfa(ja1)=fyfa(ja1)+dfy
!$omp atomic
         fzfa(ja1)=fzfa(ja1)+dfz

!$omp atomic
         fxfa(ja2)=fxfa(ja2)-dfx
!$omp atomic
         fyfa(ja2)=fyfa(ja2)-dfy
!$omp atomic
         fzfa(ja2)=fzfa(ja2)-dfz

      end if

   end do

!$omp enddo nowait


!-------------------------------


!  Myp2-Myp2 forces:

!$omp do schedule(guided,64)

   do n=1,npair_myp2

      jm1=pair_myp2(1,n)
      jm2=pair_myp2(2,n)

      dx=xmyp2(jm1)-xmyp2(jm2)
      dy=ymyp2(jm1)-ymyp2(jm2)
      dz=zmyp2(jm1)-zmyp2(jm2)

      d2=dx*dx+dy*dy+dz*dz

      if(d2<r_on2)then

         dist=sqrt(d2)

         f=fvdwmax/dist

         dfx=f*dx!*invdist
         dfy=f*dy!*invdist
         dfz=f*dz!*invdist

!$omp atomic
         fxmyp2(jm1)=fxmyp2(jm1)+dfx
!$omp atomic
         fymyp2(jm1)=fymyp2(jm1)+dfy
!$omp atomic
         fzmyp2(jm1)=fzmyp2(jm1)+dfz

!$omp atomic
         fxmyp2(jm2)=fxmyp2(jm2)-dfx
!$omp atomic
         fymyp2(jm2)=fymyp2(jm2)-dfy
!$omp atomic
         fzmyp2(jm2)=fzmyp2(jm2)-dfz

      else if(d2<r_off2)then

         dist=sqrt(d2)

!         invdist=1/dist

         ratio=(r_off-dist)/(dist-rvdw)

         f=kvdw*ratio*ratio/dist!(r_off-dist)*(r_off-dist)/(dist-rvdw)/(dist-rvdw)

         dfx=f*dx!*invdist
         dfy=f*dy!*invdist
         dfz=f*dz!*invdist

!$omp atomic
         fxmyp2(jm1)=fxmyp2(jm1)+dfx
!$omp atomic
         fymyp2(jm1)=fymyp2(jm1)+dfy
!$omp atomic
         fzmyp2(jm1)=fzmyp2(jm1)+dfz

!$omp atomic
         fxmyp2(jm2)=fxmyp2(jm2)-dfx
!$omp atomic
         fymyp2(jm2)=fymyp2(jm2)-dfy
!$omp atomic
         fzmyp2(jm2)=fzmyp2(jm2)-dfz

      end if

   end do

!$omp enddo nowait



!-------------------------------

!  Myo2-Myo2 forces:


!$omp do schedule(guided,64)

   do n=1,npair_myo2

      jm1=pair_myo2(1,n)
      jm2=pair_myo2(2,n)

      dx=xmyo2(jm1)-xmyo2(jm2)
      dy=ymyo2(jm1)-ymyo2(jm2)
      dz=zmyo2(jm1)-zmyo2(jm2)

      d2=dx*dx+dy*dy+dz*dz

      if(d2<r_on2)then

         dist=sqrt(d2)

         f=fvdwmax/dist

         dfx=f*dx!*invdist
         dfy=f*dy!*invdist
         dfz=f*dz!*invdist

!$omp atomic
         fxmyo2(jm1)=fxmyo2(jm1)+dfx
!$omp atomic
         fymyo2(jm1)=fymyo2(jm1)+dfy
!$omp atomic
         fzmyo2(jm1)=fzmyo2(jm1)+dfz

!$omp atomic
         fxmyo2(jm2)=fxmyo2(jm2)-dfx
!$omp atomic
         fymyo2(jm2)=fymyo2(jm2)-dfy
!$omp atomic
         fzmyo2(jm2)=fzmyo2(jm2)-dfz

      else if(d2<r_off2)then

         dist=sqrt(d2)

         ratio=(r_off-dist)/(dist-rvdw)

         f=kvdw*ratio*ratio/dist!(r_off-dist)*(r_off-dist)/(dist-rvdw)/(dist-rvdw)

         dfx=f*dx!*invdist
         dfy=f*dy!*invdist
         dfz=f*dz!*invdist

!$omp atomic
         fxmyo2(jm1)=fxmyo2(jm1)+dfx
!$omp atomic
         fymyo2(jm1)=fymyo2(jm1)+dfy
!$omp atomic
         fzmyo2(jm1)=fzmyo2(jm1)+dfz

!$omp atomic
         fxmyo2(jm2)=fxmyo2(jm2)-dfx
!$omp atomic
         fymyo2(jm2)=fymyo2(jm2)-dfy
!$omp atomic
         fzmyo2(jm2)=fzmyo2(jm2)-dfz

      end if

   end do

!$omp enddo nowait


!-------------------------------

!  Myp2-Myo2 forces:

!$omp do schedule(guided,64)

   do n=1,npair_mymy

      jm=pair_mymy(1,n)

      jp=pair_mymy(2,n)

      dx=xmyo2(jm)-xmyp2(jp)
      dy=ymyo2(jm)-ymyp2(jp)
      dz=zmyo2(jm)-zmyp2(jp)

      d2=dx*dx+dy*dy+dz*dz

      if(d2<r_on2)then

         dist=sqrt(d2)

         f=fvdwmax/dist

         dfx=f*dx!*invdist
         dfy=f*dy!*invdist
         dfz=f*dz!*invdist

!$omp atomic
         fxmyo2(jm)=fxmyo2(jm)+dfx
!$omp atomic
         fymyo2(jm)=fymyo2(jm)+dfy
!$omp atomic
         fzmyo2(jm)=fzmyo2(jm)+dfz

!$omp atomic
         fxmyp2(jp)=fxmyp2(jp)-dfx
!$omp atomic
         fymyp2(jp)=fymyp2(jp)-dfy
!$omp atomic
         fzmyp2(jp)=fzmyp2(jp)-dfz


      else if(d2<r_off2)then

         dist=sqrt(d2)

         ratio=(r_off-dist)/(dist-rvdw)

         f=kvdw*ratio*ratio/dist!(r_off-dist)*(r_off-dist)/(dist-rvdw)/(dist-rvdw)

         dfx=f*dx!*invdist
         dfy=f*dy!*invdist
         dfz=f*dz!*invdist

!$omp atomic
         fxmyo2(jm)=fxmyo2(jm)+dfx
!$omp atomic
         fymyo2(jm)=fymyo2(jm)+dfy
!$omp atomic
         fzmyo2(jm)=fzmyo2(jm)+dfz
      
!$omp atomic
         fxmyp2(jp)=fxmyp2(jp)-dfx
!$omp atomic
         fymyp2(jp)=fymyp2(jp)-dfy
!$omp atomic
         fzmyp2(jp)=fzmyp2(jp)-dfz

      end if

   end do

!$omp enddo nowait


!-------------------------------

!  crosslinker-crosslinker forces:

!$omp do schedule(guided,64)

   do n=1,npair_lk2

      jl1=pair_lk2(1,n)
      jl2=pair_lk2(2,n)

      dx=xlk(jl1)-xlk(jl2)
      dy=ylk(jl1)-ylk(jl2)
      dz=zlk(jl1)-zlk(jl2)

      d2=dx*dx+dy*dy+dz*dz

      if(d2<r_on2)then

         dist=sqrt(d2)

         f=fvdwmax/dist

         dfx=f*dx!*invdist
         dfy=f*dy!*invdist
         dfz=f*dz!*invdist

!$omp atomic
         fxlk(jl1)=fxlk(jl1)+dfx
!$omp atomic
         fylk(jl1)=fylk(jl1)+dfy
!$omp atomic
         fzlk(jl1)=fzlk(jl1)+dfz

!$omp atomic
         fxlk(jl2)=fxlk(jl2)-dfx
!$omp atomic
         fylk(jl2)=fylk(jl2)-dfy
!$omp atomic
         fzlk(jl2)=fzlk(jl2)-dfz

     else if(d2<r_off2)then

         dist=sqrt(d2)

!         invdist=1/dist

         ratio=(r_off-dist)/(dist-rvdw)

         f=kvdw*ratio*ratio/dist!(r_off-dist)*(r_off-dist)/(dist-rvdw)/(dist-rvdw)

         dfx=f*dx!*invdist
         dfy=f*dy!*invdist
         dfz=f*dz!*invdist

!$omp atomic
         fxlk(jl1)=fxlk(jl1)+dfx
!$omp atomic
         fylk(jl1)=fylk(jl1)+dfy
!$omp atomic
         fzlk(jl1)=fzlk(jl1)+dfz

!$omp atomic
         fxlk(jl2)=fxlk(jl2)-dfx
!$omp atomic
         fylk(jl2)=fylk(jl2)-dfy
!$omp atomic
         fzlk(jl2)=fzlk(jl2)-dfz

      end if

   end do

!$omp enddo nowait




!  membrane-membrane forces:

!$omp do schedule(guided,64)

   do n=1,npair_mb2

      n1=pair_mb2(1,n)
      n2=pair_mb2(2,n)

      dx=xmemb(n1)-xmemb(n2)
      dy=ymemb(n1)-ymemb(n2)
      dz=zmemb(n1)-zmemb(n2)

      if(boundtyp(n1)==1.and.boundtyp(n2)==2) dx=dx+shift

      d2=dx*dx+dy*dy+dz*dz

      dist=sqrt(d2)

!      invdist=1/dist


!     tethering:

      if(dist>l_pair)then

         f=kpair*(dist-l_pair)*(dist-l_pair)/dist

         dfx=f*dx!*invdist
         dfy=f*dy!*invdist
         dfz=f*dz!*invdist

!$omp atomic
         fxmemb(n1)=fxmemb(n1)-dfx
!$omp atomic
         fymemb(n1)=fymemb(n1)-dfy
!$omp atomic
         fzmemb(n1)=fzmemb(n1)-dfz

!$omp atomic
         fxmemb(n2)=fxmemb(n2)+dfx
!$omp atomic
         fymemb(n2)=fymemb(n2)+dfy
!$omp atomic
         fzmemb(n2)=fzmemb(n2)+dfz


      elseif(dist<l_mem)then

         f=-kpair*(dist-l_mem)*(dist-l_mem)/dist

         dfx=f*dx!*invdist
         dfy=f*dy!*invdist
         dfz=f*dz!*invdist

!$omp atomic
         fxmemb(n1)=fxmemb(n1)-dfx
!$omp atomic
         fymemb(n1)=fymemb(n1)-dfy
!$omp atomic
         fzmemb(n1)=fzmemb(n1)-dfz

!$omp atomic
         fxmemb(n2)=fxmemb(n2)+dfx
!$omp atomic
         fymemb(n2)=fymemb(n2)+dfy
!$omp atomic
         fzmemb(n2)=fzmemb(n2)+dfz

      end if




!  to keep layer structure of the membrane:


      n3=pairpart(1,n)
      n4=pairpart(2,n)

      if(n3==0)then
         cycle
      end if


      dx34=xmemb(n3)-xmemb(n4)
      dy34=ymemb(n3)-ymemb(n4)
      dz34=zmemb(n3)-zmemb(n4)

      if(boundtyp(n3)==1.and.boundtyp(n4)==2) dx34=dx34+shift

      if(boundtyp(n3)==2.and.boundtyp(n4)==1) dx34=dx34-shift


      xn=dy*dz34-dz*dy34
      yn=dz*dx34-dx*dz34
      zn=dx*dy34-dy*dx34

      invdist=1.0d0/sqrt(xn*xn+yn*yn+zn*zn)

      xn=xn*invdist
      yn=yn*invdist
      zn=zn*invdist

      dx=xmemb(n1)-xmemb(n3)
      dy=ymemb(n1)-ymemb(n3)
      dz=zmemb(n1)-zmemb(n3)

      if(boundtyp(n1)==1.and.boundtyp(n3)==2) dx=dx+shift

      if(boundtyp(n1)==2.and.boundtyp(n3)==1) dx=dx-shift


      proj=dx*xn+dy*yn+dz*zn

      dx=proj*xn
      dy=proj*yn
      dz=proj*zn


      dfx=kmemb*dx
      dfy=kmemb*dy
      dfz=kmemb*dz









!$omp atomic
      fxmemb(n1)=fxmemb(n1)-dfx
!$omp atomic
      fymemb(n1)=fymemb(n1)-dfy
!$omp atomic
      fzmemb(n1)=fzmemb(n1)-dfz

!$omp atomic
      fxmemb(n2)=fxmemb(n2)-dfx
!$omp atomic
      fymemb(n2)=fymemb(n2)-dfy
!$omp atomic
      fzmemb(n2)=fzmemb(n2)-dfz

!$omp atomic
      fxmemb(n3)=fxmemb(n3)+dfx
!$omp atomic
      fymemb(n3)=fymemb(n3)+dfy
!$omp atomic
      fzmemb(n3)=fzmemb(n3)+dfz

!$omp atomic
      fxmemb(n4)=fxmemb(n4)+dfx
!$omp atomic
      fymemb(n4)=fymemb(n4)+dfy
!$omp atomic
      fzmemb(n4)=fzmemb(n4)+dfz

   end do


!$omp enddo nowait





!$omp end parallel


   end subroutine

!=========================================================

   subroutine my2memtether(myo2num,myo2body,my2mem,fxmemb,fymemb,fzmemb,fxmyo2,fymyo2,fzmyo2, &
                xmemb,ymemb,zmemb,xmyo2,ymyo2,zmyo2,kmy2mem,lmy2mem)

   implicit none

   integer,value::myo2num
   integer jmb,jmy,nm

   integer,allocatable,dimension(:),intent(in)::my2mem
   integer,allocatable,dimension(:,:),intent(in)::myo2body

   double precision,value::kmy2mem,lmy2mem
   double precision dx,dy,dz,f,fx,fy,fz,dist,dist2,l2on

   double precision,allocatable,dimension(:),intent(in)::xmemb,ymemb,zmemb
   double precision,allocatable,dimension(:)::fxmemb,fymemb,fzmemb
   double precision,allocatable,dimension(:),intent(in)::xmyo2,ymyo2,zmyo2
   double precision,allocatable,dimension(:)::fxmyo2,fymyo2,fzmyo2

   l2on=lmy2mem*lmy2mem

!$omp parallel &
!$omp default(none) &


!$omp private(jmb,jmy,nm) &

!$omp shared(my2mem,myo2body,kmy2mem,lmy2mem,l2on) &
!$omp private(dx,dy,dz,f,fx,fy,fz,dist,dist2) &

!$omp shared(myo2num,xmemb,ymemb,zmemb) &
!$omp shared(fxmemb,fymemb,fzmemb) &
!$omp shared(xmyo2,ymyo2,zmyo2,fxmyo2,fymyo2,fzmyo2)


!$omp do schedule (guided,32)

   do nm=1,myo2num

      jmb=my2mem(nm)

      jmy=myo2body(1,nm)

      dx=xmemb(jmb)-xmyo2(jmy)
      dy=ymemb(jmb)-ymyo2(jmy)
      dz=zmemb(jmb)-zmyo2(jmy)

      dist2=dx*dx+dy*dy+dz*dz

      if(dist2>l2on)then

         dist=sqrt(dist2)

         f=kmy2mem*(dist-lmy2mem)/dist


         fx=f*dx!*invdist
         fy=f*dy!*invdist
         fz=f*dz!*invdist

         fxmemb(jmb)=fxmemb(jmb)-fx
         fymemb(jmb)=fymemb(jmb)-fy
         fzmemb(jmb)=fzmemb(jmb)-fz

         fxmyo2(jmy)=fxmyo2(jmy)+fx
         fymyo2(jmy)=fymyo2(jmy)+fy
         fzmyo2(jmy)=fzmyo2(jmy)+fz

      end if

   end do

!$omp enddo nowait

!$omp end parallel

   end subroutine

!=========================================================

   subroutine sethoops(nwahoop,xwhoop,nmemb,xmemb,l_mem,mbhoop,mbhooplen)

   implicit none

   integer,value:: nmemb,nwahoop
   integer n,jh,jget

   integer, allocatable, dimension(:)::mbhooplen
   integer, allocatable, dimension(:,:)::mbhoop

   double precision,value::l_mem
   double precision dist,d

   double precision, allocatable,intent(in), dimension(:)::xwhoop,xmemb

   mbhooplen=0

   do n=1,nmemb

      dist=2*l_mem

      jget=0

      do jh=1,nwahoop

         d=abs(xmemb(n)-xwhoop(jh))

         if(d<dist)then
            dist=d
            jget=jh
         end if

      end do

      if(jget>0)then

         mbhooplen(jget)=mbhooplen(jget)+1

         mbhoop(mbhooplen(jget),jget)=n

      end if

   end do

   end subroutine

!=========================================================

   subroutine wallremod(nmemb,nxsol,nphisol,nwall,jmbsol,wthick,gap,xgap,dphisol,wallrate, &
                xmemb,xboundmin,xboundmax,rwall,rmembmax,ywall,zwall)



   implicit none

   integer,value::nmemb,nxsol,nphisol,nwall
   integer jp,jx,n
   integer,allocatable,intent(in), dimension(:,:)::jmbsol

   double precision,value::wthick,gap,xgap,dphisol
   double precision wallrate
   double precision dwallrate
   double precision, allocatable,intent(in), dimension(:)::xmemb
   double precision, allocatable, dimension(:)::xboundmin,xboundmax,xmin,xmax
   double precision, allocatable, dimension(:,:)::rwall,rmembmax,ywall,zwall

   do jx=1,nxsol

      do jp=1,nphisol

         if(rwall(jp,jx)>rmembmax(jp,jx)+wthick+gap)then

            dwallrate=0.1d0*gap

            wallrate=wallrate+dwallrate/nwall

            rwall(jp,jx)=rwall(jp,jx)-dwallrate

            ywall(jp,jx)=rwall(jp,jx)*cos(jp*dphisol)
            zwall(jp,jx)=rwall(jp,jx)*sin(jp*dphisol)

         end if

      end do

   end do


!  boundary remodeling:

   allocate(xmin(nphisol),xmax(nphisol))
   xmin=0.0d0; xmax=0.0d0

   do n=1,nmemb

      jp=jmbsol(2,n)

      if(xmin(jp)>xmemb(n))then
         xmin(jp)=xmemb(n)
      end if

      if(xmax(jp)<xmemb(n))then
         xmax(jp)=xmemb(n)
      end if

   end do

!   xmin=minval(xmemb(1:nmemb))
!   xmax=maxval(xmemb(1:nmemb))

   do n=1,nphisol

      if(xmin(n)>xboundmin(n)+xgap)then
         xboundmin(n)=xboundmin(n)+0.1d0
      end if

      if(xboundmax(n)>xmax(n)+xgap)then
         xboundmax(n)=xboundmax(n)-0.1d0
      end if

   end do

   deallocate(xmin,xmax)

   end subroutine

!=========================================================


end module
