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

   subroutine setjob(nthreads,jobid,mbrad,mwid,rrad,rwid,rthick,falen,fanum,crlknum1,crlknum2, &
              myonum,node_size,kmemb,tension1,tension2,crlkorient,j_dep,cofilin,jxldel,jmyoturn,jxlturnover,halftime)

   implicit none

   integer nthreads,jobid,fanum,falen,crlknum1,crlknum2,myonum,nnode
   integer tension1,tension2,crlkorient,j_dep,cofilin,jxldel,jmyoturn,jxlturnover

   character (len=100) chara
   double precision mbrad,mwid,rrad,rwid,rthick,kmemb,node_size,halftime

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

   rrad=mbrad-10.0d0

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

10 read(1,2)chara
   if(chara(1:4)/='MYON')then
      goto 10
   end if
   read(1,*)myonum

11 read(1,2)chara
   if(chara(1:4)/='CRLK')then
      goto 11
   end if
   read(1,*)crlknum1
   read(1,*)crlknum2

17 read(1,2)chara
   if(chara(1:4)/='NODE')then
      goto 17
   end if
   read(1,*)nnode
   read(1,*)node_size
!   read(1,*)pollard

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

37 read(1,*)chara

   if(chara(1:4)/='MYOT')then
      goto 37
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

48 read(1,*)chara

   if(chara(1:4)/='HALF')then
      goto 48
   end if

   read(1,*)halftime


   close(1)

   end subroutine

!=========================================================

   subroutine makering(nbondwal,nwall,nmemb,nnode,nfa1,nfa2,fanum1,fanum2,nmyo,nbondmyo,myonum,nlk, &
              nbondlk,mybodylen,crlknum1,crlknum2,nxsol,nphisol,nmb2node,astart,alen,a2mem, &
              my2node,filid,apos,lkstart,bondwal,bondmyo,mybody,myhead,bondlk,memb2node,dxsol, &
              dphisol,xboundmin,xboundmax,xwall,ywall,zwall,xnorwall,ynorwall,znorwall,xmemb,ymemb,zmemb, &
              xsurf,ysurf,zsurf,xfa,yfa,zfa,xmyo,ymyo,zmyo,xlk,ylk,zlk,xnode,ynode,znode)



   implicit none

   integer,value::nbondwal,nwall,nmemb,nnode,nfa1,nfa2,fanum1,fanum2,nmyo,nbondmyo,myonum
   integer,value::nlk,nbondlk,mybodylen,crlknum1,crlknum2,nxsol,nphisol,nmb2node
   integer nfa,nact1,nact2,fanum,nfamax1,nfamax2,natom,natom_a1,natom_a2,natom_a2m
   integer nbond,nbond_a1,nbond_a2,nbond_a2m,fanummax,fanummax1,fanummax2,crlknum,crlknummax,crlknumactive,nbond0

   integer jstart,na,nmyoturn
   integer jm,nm,j
   integer izero,ires,n,jatom
   integer j1,j2,j3,j4,j5,j6,j7,j8,nline,i,jx,jp
   integer nf,ja,jnode,jmemb

   integer,allocatable,dimension(:),intent(in)::astart,alen,a2mem,my2node,filid,apos,lkstart
   integer,allocatable,dimension(:,:),intent(in)::bondwal,bondmyo,mybody,myhead,bondlk,memb2node
   integer,allocatable,dimension(:,:)::bond,mytyp,apar,fa2myo,fa2lk
   integer,allocatable,dimension(:)::lktyp

   real zero

   double precision,value::dxsol,dphisol,xboundmin,xboundmax

   double precision charg,mass,w1,w2

   double precision,allocatable,dimension(:,:),intent(in)::xwall,ywall,zwall
   double precision,allocatable,dimension(:,:),intent(in)::xnorwall,ynorwall,znorwall
   double precision,allocatable,dimension(:),intent(in)::xmemb,ymemb,zmemb
   double precision,allocatable,dimension(:,:),intent(in)::xsurf,ysurf,zsurf
   double precision,allocatable,dimension(:),intent(in)::xfa,yfa,zfa,xmyo,ymyo,zmyo,xlk,ylk,zlk

   double precision,allocatable,dimension(:),intent(in)::xnode,ynode,znode

   character(8) tex,typ,res,segid,resno


   open(1,file='ring.psf')

20  FORMAT(I8,1X,A4,1X,A4,1X,A3,2X,A3,2X,A4,2X,F9.6,6X,F8.4,11X,I1)
21  FORMAT(I8,1X,A)

   open(2,file='startring.pdb')

42 FORMAT(A6,I5,2X,A3,1X,A3,1X,I5,4X,3F8.3,2F6.2,6X,A4)
52 FORMAT(A50)


!  number of bonds to start:
   nbond=nbondwal!+nbondfa+nbondmyo+nbondlk

!  number of beads to start:
   natom=nwall+nmemb+nnode!+2*nbondfa+nmyo+nlk+nnode

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


!  to visualize tethering of F-actin to membrane:

   fanum=fanum1+fanum2

   fanummax1=4*fanum1

   fanummax2=4*fanum2

   fanummax=fanummax1+fanummax2


   natom=natom+2*fanummax

   natom_a2m=natom

   nbond=nbond+fanummax

   nbond_a2m=nbond

!  myosins:

   natom=natom+nmyo

   nbond=nbond+nbondmyo


!  to visualize tethering of myosin to nodes:

   natom=natom+2*myonum
   nbond=nbond+myonum

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

   do jx=1,nxsol

      do jp=1,nphisol

         jatom=jatom+1

         write(1,20)jatom,segid,resno,res,typ,tex,charg,mass,izero

         write(2,42)tex,jatom,typ,res,ires,xwall(jp,jx)/10,ywall(jp,jx)/10,zwall(jp,jx)/10,w1,w2,segid

      end do

   end do


   bond(1:2,1:nbondwal)=bondwal(1:2,1:nbondwal)

   nbond=nbondwal

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

!      nbond=nbond+1

!      bond(1,nbond)=jatom
!      bond(2,nbond)=jatom+nmemb

   end do

!  for nodes:

   res='NOD'
   typ='NOD'
   segid='NOD'

   ires=3
   write(resno,'(i1)')ires

   do n=1,nnode

      jatom=jatom+1

      write(1,20)jatom,segid,resno,res,typ,tex,charg,mass,izero

      write(2,42)tex,jatom,typ,res,ires,xnode(n)/10,ynode(n)/10,znode(n)/10,w1,w2,segid

   end do



!  for F-actin

   res='FAC'
   typ='FAC'
!   segid='FAC'

   ires=4
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

!  counter-clockwise:

   segid='FA2'

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



!  binding of actin plus ends to membrane:

   segid='A2M'
   typ='A2M'
   res='A2M'
   ires=5

   write(resno,'(i1)')ires

   do nf=1,fanum

      jatom=jatom+1

      write(1,20)jatom,segid,resno,res,typ,tex,charg,mass,izero

      ja=astart(nf)

      write(2,42)tex,jatom,typ,res,ires,xfa(ja)/10,yfa(ja)/10,zfa(ja)/10,w1,w2,segid

      jatom=jatom+1

      write(1,20)jatom,segid,resno,res,typ,tex,charg,mass,izero

      jmemb=a2mem(nf)

      write(2,42)tex,jatom,typ,res,ires,xmemb(jmemb)/10,ymemb(jmemb)/10,zmemb(jmemb)/10,w1,w2,segid

      nbond=nbond+1

      bond(1,nbond)=jatom-1

      bond(2,nbond)=jatom

   end do


   if(nbond<nbond_a2m)then

      nbond0=nbond

      do n=nbond0+1,nbond_a2m

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






!  for myosin

   res='MYO'
   typ='MYO'
   segid='MYO'

   ires=6
   write(resno,'(i1)')ires


   n=0

   do nm=1,myonum

      segid='MYB'

      do j=1,mybodylen

         n=n+1

         jatom=jatom+1

         write(1,20)jatom,segid,resno,res,typ,tex,charg,mass,izero

         write(2,42)tex,jatom,typ,res,ires,xmyo(n)/10,ymyo(n)/10,zmyo(n)/10,w1,w2,segid

      end do

      segid='MYH'

      do j=1,2

         jatom=jatom+1

         n=n+1

         write(1,20)jatom,segid,resno,res,typ,tex,charg,mass,izero

         write(2,42)tex,jatom,typ,res,ires,xmyo(n)/10,ymyo(n)/10,zmyo(n)/10,w1,w2,segid

      end do

   end do

   do n=1,nbondmyo

      nbond=nbond+1

      bond(1,nbond)=bondmyo(1,n)+jatom-nmyo
      bond(2,nbond)=bondmyo(2,n)+jatom-nmyo

   end do


!  binding of myosin to nodes:

   segid='M2N'
   typ='M2N'
   res='M2N'
   ires=7

   write(resno,'(i1)')ires

   do nm=1,myonum

      jatom=jatom+1

      write(1,20)jatom,segid,resno,res,typ,tex,charg,mass,izero

      jm=mybody(1,nm)

      write(2,42)tex,jatom,typ,res,ires,xmyo(jm)/10,ymyo(jm)/10,zmyo(jm)/10,w1,w2,segid

      jatom=jatom+1

      write(1,20)jatom,segid,resno,res,typ,tex,charg,mass,izero

      jnode=my2node(nm)

      write(2,42)tex,jatom,typ,res,ires,xnode(jnode)/10,ynode(jnode)/10,znode(jnode)/10,w1,w2,segid

      nbond=nbond+1

      bond(1,nbond)=jatom-1

      bond(2,nbond)=jatom

   end do

!  for crosslinkers:

   res='CLK'
   typ='CLK'
   segid='CLK'

   ires=8
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

   write(2,52)tex

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


!   write(1)rwall(1:nphisol,1:nxsol)


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

   write(1)nnode
   write(1)xnode(1:nnode)
   write(1)ynode(1:nnode)
   write(1)znode(1:nnode)

   nfa=nfa1+nfa2

   write(1)nfa
   write(1)xfa(1:nfa)
   write(1)yfa(1:nfa)
   write(1)zfa(1:nfa)

   write(1)nmyo
   write(1)xmyo(1:nmyo)
   write(1)ymyo(1:nmyo)
   write(1)zmyo(1:nmyo)

   write(1)nlk
   write(1)xlk(1:nlk)
   write(1)ylk(1:nlk)
   write(1)zlk(1:nlk)

   close(1)

!---------------------------------------------

!  write configuration

   open(1,file='rconf000.inp')

   write(1,*)'VISUAL'

   write(1,*)natom,natom_a1,natom_a2,natom_a2m

   write(1,*)'CELLWALL'
   write(1,*)nwall,nxsol,nphisol
   write(1,*)

   write(1,*)'MEMBRANE'
   write(1,*)nmemb
   write(1,*)

   write(1,*)'NODES'
   write(1,*)nnode
   do n=1,nnode
      write(1,*)memb2node(1:nmb2node,n)
   end do
   write(1,*)


   write(1,*)'FACTIN'
   write(1,*)nact1,nact2,nfamax1,nfamax2,nfa1,nfa2,fanum1,fanum2,fanummax1,fanummax2
   write(1,*)filid(1:nfa)
!   write(1,*)a2mem(1:nfa)

   allocate(apar(2,nfa))
   apar=0

   apar(1,astart(1:fanum))=-1

   write(1,*)apar(1:2,1:nfa)
   deallocate(apar)

   write(1,*)apos(1:nfa)

!   allocate(fa1stbound(fanum))
!   fa1stbound(1:fanum)=alen(1:fanum)

!   write(1,*)fa1stbound(1:fanum)
!   deallocate(fa1stbound)


   write(1,*)alen(1:fanum)

   write(1,*)astart(1:fanum)

   write(1,*)a2mem(1:fanum)

   write(1,*)

   nmyoturn=0
   write(1,*)'MYOSIN'
   write(1,*)nmyo,myonum,nmyoturn

   allocate(mytyp(2,myonum),fa2myo(2,myonum))
   mytyp=1; fa2myo=0

   do n=1,myonum
      write(1,*)my2node(n),myhead(1:2,n),mybody(1:mybodylen,n)
      write(1,*)mytyp(1:2,n),fa2myo(1:2,n)
   end do
   write(1,*)

   deallocate(mytyp,fa2myo)


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

!   do n=1,crlknum
!      write(1,*)crlk(1:lklen,n)
!      write(1,*)fa2lk(1:2,n)
!   end do
!   write(1,*)

   deallocate(fa2lk,lktyp)


   close(1)

   open(1,file='restart.inp')

   write(1,*)0,0
   write(1,*)0,0.0
   close(1)

   end subroutine

!=========================================================
   subroutine writedcd(nframe,junit,natom,natom_a1,natom_a2,natom_a2m,nxsol,nphisol,nmemb,nnode, &
              fanum1,fanum2,nmyo,myonum,nlk,alen,astart,a2mem,my2node,mybody,xfa,yfa,zfa,xmyo, &
              ymyo,zmyo,xnode,ynode,znode,xlk,ylk,zlk,xmemb,ymemb,zmemb,xwall,ywall,zwall)



   implicit none

   integer nframe
   integer,value::junit,natom,natom_a1,natom_a2,natom_a2m,nxsol,nphisol,nmemb
   integer,value::nnode,fanum1,fanum2,nmyo,myonum,nlk
   integer n,j,ja,jw,jm,jp,jx,jstart,jnode,jmemb

   integer,allocatable,intent(in),dimension(:)::alen,astart,a2mem,my2node
   integer,allocatable,intent(in),dimension(:,:)::mybody

   real,allocatable,dimension(:)::xw,yw,zw
   double precision,allocatable,dimension(:),intent(in)::xfa,yfa,zfa

   double precision,allocatable,intent(in),dimension(:)::xmyo,ymyo,zmyo,xnode,ynode,znode
   double precision,allocatable,intent(in),dimension(:)::xlk,ylk,zlk,xmemb,ymemb,zmemb
   double precision,allocatable,intent(in),dimension(:,:)::xwall,ywall,zwall


   nframe=nframe+1

   allocate(xw(natom),yw(natom),zw(natom))

   jw=0

   do jx=1,nxsol

      do jp=1,nphisol

         jw=jw+1

         xw(jw)=0.1*xwall(jp,jx)
         yw(jw)=0.1*ywall(jp,jx)
         zw(jw)=0.1*zwall(jp,jx)

      end do

   end do

   xw(jw+1:jw+nmemb)=0.1*xmemb(1:nmemb)
   yw(jw+1:jw+nmemb)=0.1*ymemb(1:nmemb)
   zw(jw+1:jw+nmemb)=0.1*zmemb(1:nmemb)

   jw=jw+nmemb


   xw(jw+1:jw+nnode)=0.1*xnode(1:nnode)
   yw(jw+1:jw+nnode)=0.1*ynode(1:nnode)
   zw(jw+1:jw+nnode)=0.1*znode(1:nnode)

   jw=jw+nnode


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

      xw(jw+1:natom_a2)=0.0
      yw(jw+1:natom_a2)=0.0
      zw(jw+1:natom_a2)=0.0

      jw=natom_a2

   end if


!  tethering actin to membrane:

   do n=1,fanum1+fanum2

      jmemb=a2mem(n)

      ja=astart(n)

      jw=jw+1

      xw(jw)=0.1*xmemb(jmemb)
      yw(jw)=0.1*ymemb(jmemb)
      zw(jw)=0.1*zmemb(jmemb)

      jw=jw+1

      xw(jw)=0.1*xfa(ja)
      yw(jw)=0.1*yfa(ja)
      zw(jw)=0.1*zfa(ja)

   end do

   if(jw<natom_a2m)then

      xw(jw+1:natom_a2m)=0.0
      yw(jw+1:natom_a2m)=0.0
      zw(jw+1:natom_a2m)=0.0

      jw=natom_a2m


   end if


!  visualize myosin:


   xw(jw+1:jw+nmyo)=0.1*xmyo(1:nmyo)
   yw(jw+1:jw+nmyo)=0.1*ymyo(1:nmyo)
   zw(jw+1:jw+nmyo)=0.1*zmyo(1:nmyo)

   jw=jw+nmyo

!  tethering myosin to nodes:

   do n=1,myonum

      jnode=my2node(n)

      jm=mybody(1,n)

      jw=jw+1

      xw(jw)=0.1*xnode(jnode)
      yw(jw)=0.1*ynode(jnode)
      zw(jw)=0.1*znode(jnode)

      jw=jw+1

      xw(jw)=0.1*xmyo(jm)
      yw(jw)=0.1*ymyo(jm)
      zw(jw)=0.1*zmyo(jm)

   end do

!  visualize crosslinkers:

   xw(jw+1:jw+nlk)=0.1*xlk(1:nlk)
   yw(jw+1:jw+nlk)=0.1*ylk(1:nlk)
   zw(jw+1:jw+nlk)=0.1*zlk(1:nlk)

   write(junit)xw(1:natom)
   write(junit)yw(1:natom)
   write(junit)zw(1:natom)

   deallocate(xw,yw,zw)

   end subroutine

!=========================================================

   subroutine getinfo(nstart,jfile,time0,runtime,natom,natom_a1,natom_a2,natom_a2m, &
              nwall,nxsol,nphisol,nmemb,nnode,nact1,nact2,nfamax1,nfamax2,nfa1,nfa2, &
              fanum1,fanum2,fanummax1,fanummax2,nmyo,myonum,nlk,crlknum1,crlknum2,crlknummax,crlknumactive,nmyoturn)



   implicit none

   integer(kind=8)::time0
   integer nstart,jfile,natom,natom_a1,natom_a2,natom_a2m,nwall,nxsol,nphisol
   integer nmemb,nnode,nact1,nact2,nfamax1,nfamax2,nfa1,nfa2,fanum1,fanum2,fanummax1,fanummax2
   integer nmyo,myonum,nlk,crlknum1,crlknum2,crlknummax,crlknumactive,nmyoturn
   double precision runtime

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

21 read(1,*)chara

   if(chara(1:4)/='VISU')then
      goto 21
   end if
   read(1,*)natom,natom_a1,natom_a2,natom_a2m



2  read(1,*)chara

   if(chara(1:4)/='CELL')then
      goto 2
   end if

   read(1,*)nwall,nxsol,nphisol

12 read(1,*)chara

   if(chara(1:4)/='MEMB')then
      goto 12
   end if

   read(1,*)nmemb

13 read(1,*)chara

   if(chara(1:4)/='NODE')then
      goto 13
   end if
   read(1,*)nnode

3  read(1,*)chara

   if(chara(1:4)/='FACT')then
      goto 3
   end if

   read(1,*)nact1,nact2,nfamax1,nfamax2,nfa1,nfa2,fanum1,fanum2,fanummax1,fanummax2

!   nfa=nfa1+nfa2

!   fanum=fanum1+fanum2


5  read(1,*)chara

   if(chara(1:4)/='MYOS')then
      goto 5
   end if

   read(1,*)nmyo,myonum,nmyoturn


8  read(1,*)chara

   if(chara(1:4)/='CROS')then
      goto 8
   end if

   read(1,*)nlk,crlknum1,crlknum2,crlknummax,crlknumactive

!   crlknum=crlknum1+crlknum2

   close(1)

   end subroutine

!=========================================================

   subroutine ringin(nstart,nxsol,nphisol,nmemb,fanum,nfa,myonum,nmyo,nnode,mybodylen,crlknum, &
               nlk,nmb2node,astart,alen,a2mem,my2node,filid,apos,lkstart,lktyp,mybody,myhead,mytyp,apar, &
                memb2node,fa2myo,fa2lk,xboundmin,xboundmax,dxsol,dphisol,xwall,ywall,zwall, &
                 xnorwall,ynorwall,znorwall,xfa,yfa,zfa,xmyo,ymyo,zmyo,xnode,ynode,znode, &
                  xlk,ylk,zlk,xmemb,ymemb,zmemb,xsurf,ysurf,zsurf,xnorsurf,ynorsurf,znorsurf)



   implicit none

   integer,value::nstart,nxsol,nphisol,nmemb,fanum,nfa,myonum,nmyo,nnode
   integer,value::mybodylen,crlknum,nlk,nmb2node
   integer n,jread(10)

   integer,allocatable,dimension(:)::astart,alen,a2mem,my2node,filid,apos,lkstart,lktyp
   integer,allocatable,dimension(:,:)::mybody,myhead,mytyp,apar,memb2node
   integer,allocatable,dimension(:,:)::fa2myo,fa2lk

   double precision xboundmin,xboundmax,dxsol,dphisol
   double precision,allocatable,dimension(:,:)::xwall,ywall,zwall,xnorwall,ynorwall,znorwall
   double precision,allocatable,dimension(:)::xfa,yfa,zfa,xmyo,ymyo,zmyo,xnode,ynode,znode
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



   read(1)jread(1)
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

   read(1)jread(1)!nnode
   read(1)xnode(1:nnode)
   read(1)ynode(1:nnode)
   read(1)znode(1:nnode)


   read(1)jread(1)!nfa
   read(1)xfa(1:nfa)
   read(1)yfa(1:nfa)
   read(1)zfa(1:nfa)

   read(1)jread(1)!nmyo
   read(1)xmyo(1:nmyo)
   read(1)ymyo(1:nmyo)
   read(1)zmyo(1:nmyo)

   read(1)jread(1)!nlk
   read(1)xlk(1:nlk)
   read(1)ylk(1:nlk)
   read(1)zlk(1:nlk)

   close(1)

!  read configuration

   open(1,file=fileconf)

!12 read(1,*)chara

!   if(chara(1:4)/='MEMB')then
!      goto 12
!   end if
!   read(1,*)
!   read(1,*)xboundmin,xboundmax

13 read(1,*)chara

   if(chara(1:4)/='NODE')then
      goto 13
   end if

   read(1,*)jread(1)

   do n=1,nnode
      read(1,*)memb2node(1:nmb2node,n)
   end do


3  read(1,*)chara

   if(chara(1:4)/='FACT')then
      goto 3
   end if

   read(1,*)jread(1:10)!nact1,nact2,nfamax1,nfamax2,nfa1,nfa2,fanum1,fanum2,fanummax1,fanummax2

   read(1,*)filid(1:nfa)
!   read(1,*)a2mem(1:nfa)
   read(1,*)apar(1:2,1:nfa)
   read(1,*)apos(1:nfa)
!   read(1,*)fa1stbound(1:fanum)
   read(1,*)alen(1:fanum)
   read(1,*)astart(1:fanum)
   read(1,*)a2mem(1:fanum)

!   do n=1,fanum
!      read(1,*)alen(n),fa1stbound(n),afil(1:alen(n),n)
!   enddo


5  read(1,*)chara

   if(chara(1:4)/='MYOS')then
      goto 5
   end if

   read(1,*)jread(1:3)!nmyo,myonum
   do n=1,myonum
      read(1,*)my2node(n),myhead(1:2,n),mybody(1:mybodylen,n)
      read(1,*)mytyp(1:2,n),fa2myo(1:2,n)
   end do



8  read(1,*)chara

   if(chara(1:4)/='CROS')then
      goto 8
   end if

   read(1,*)jread(1:5)

   read(1,*)lkstart(1:crlknum)
   read(1,*)lktyp(1:crlknum)
   read(1,*)fa2lk(1,1:crlknum)
   read(1,*)fa2lk(2,1:crlknum)

!   do n=1,crlknum1+crlknum2
!      read(1,*)crlk(1:lklen,n)
!      read(1,*)fa2lk(1:2,n)
!   end do


   close(1)

   end subroutine

!=========================================================

   subroutine solidset(nxsol,nphisol,nmemb,nfa,nmyo,nlk,jmbsol,jfasol,jmysol,jlksol, &
               pi,delta,dxsol,dphisol,xmemb,ymemb,zmemb,xfa,yfa,zfa,xmyo,ymyo,zmyo, &
                xlk,ylk,zlk,xwall,ywall,zwall,xnorwall,ynorwall,znorwall, &
                 xsurf,ysurf,zsurf,xnorsurf,ynorsurf,znorsurf,jsursol)


   implicit none

   integer,value:: nxsol,nphisol,nmemb,nfa,nmyo,nlk
   integer n,jx,jp,jx0,jp0,j1,j2,jxget,jpget

   integer,allocatable,dimension(:,:)::jmbsol,jfasol,jmysol,jlksol,jsursol

   double precision,value:: pi,delta,dxsol,dphisol
   double precision xmin,dxsolby2,piby2,dphisolby2,twopi,phi,arg
   double precision dx,dy,dz,d2,dist2

   double precision,allocatable,intent(in),dimension(:)::xmemb,ymemb,zmemb,xfa,yfa,zfa
   double precision,allocatable,intent(in),dimension(:)::xmyo,ymyo,zmyo,xlk,ylk,zlk
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
!$omp shared(nmyo,xmyo,jmysol,ymyo,zmyo) &

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

   do n=1,nmyo

      jx=1+(xmyo(n)-xmin)/dxsol

      if(xmyo(n)-xmin-(jx-1)*dxsol>dxsolby2)then
         jx=jx+1
      end if

      if(jx<1)then
         jx=1
      end if

      if(jx>nxsol)then
         jx=nxsol
      end if

      jx0=jx

      if(abs(ymyo(n))<delta)then
         if(zmyo(n)>0.0d0)then
            phi=piby2
         else
            phi=-piby2
         end if

      else

         arg=zmyo(n)/ymyo(n)

         phi=atan(arg)

         if(ymyo(n)<0.0d0)then
            phi=phi+pi
         end if

         if(arg<0.0d0.and.ymyo(n)>0.0d0)then
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

!      jmysol(2,n)=jp

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

            if(ymyo(n)*ysurf(jp,jx)+zmyo(n)*zsurf(jp,jx)<0.0) cycle

            dx=xmyo(n)-xsurf(jp,jx)
            dy=ymyo(n)-ysurf(jp,jx)
            dz=zmyo(n)-zsurf(jp,jx)

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

      jmysol(1,n)=jxget

      jmysol(2,n)=jpget

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

   subroutine ringout(nstart,natom,natom_a1,natom_a2,natom_a2m,nwall,nxsol,nphisol,nmemb,nnode, &
               nact1,nact2,nfamax1,nfamax2,nfa1,nfa2,fanummax1,fanummax2,fanum1,fanum2,nmyo,nmyoturn, &
                myonum,mybodylen,nlk,crlknum1,crlknum2,crlknummax,crlknumactive,nmb2node,astart,alen,a2mem,my2node, &
                 filid,apos,lkstart,lktyp,mybody,myhead,mytyp,fa2myo,fa2lk,apar,memb2node, &
                  dxsol,dphisol,xboundmin,xboundmax,xfa,yfa,zfa,xmyo,ymyo,zmyo,xnode,ynode, &
                   znode,xlk,ylk,zlk,xmemb,ymemb,zmemb,xsurf,ysurf,zsurf,xnorsurf,ynorsurf, &
                    znorsurf,xwall,ywall,zwall,xnorwall,ynorwall,znorwall)




   implicit none

   integer nstart
   integer,value::natom,natom_a1,natom_a2,natom_a2m,nwall,nxsol,nphisol,nmemb,nnode,nact1,nact2
   integer,value::nfamax1,nfamax2,nfa1,nfa2,fanummax1,fanummax2,fanum1,fanum2
   integer,value::nmyo,myonum,mybodylen,nmyoturn
   integer,value::nlk,crlknum1,crlknum2,crlknummax,crlknumactive,nmb2node

   integer n,crlknum,nfa,fanum

   integer,allocatable,intent(in),dimension(:)::astart,alen,a2mem,my2node,filid,apos,lkstart,lktyp
   integer,allocatable,intent(in),dimension(:,:)::mybody,myhead,mytyp,fa2myo,fa2lk,apar,memb2node

   double precision,value::dxsol,dphisol,xboundmin,xboundmax
   double precision,allocatable,intent(in),dimension(:)::xfa,yfa,zfa,xmyo,ymyo,zmyo,xnode,ynode,znode
   double precision,allocatable,intent(in),dimension(:)::xlk,ylk,zlk,xmemb,ymemb,zmemb
   double precision,allocatable,intent(in),dimension(:,:)::xsurf,ysurf,zsurf,xnorsurf,ynorsurf,znorsurf

   double precision,allocatable,intent(in),dimension(:,:)::xwall,ywall,zwall,xnorwall,ynorwall,znorwall

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

!  read coordinates

   open(1,file=filecoor,form='unformatted')

   write(1)dxsol,dphisol

   write(1)xwall(1:nphisol,1:nxsol)
   write(1)ywall(1:nphisol,1:nxsol)
   write(1)zwall(1:nphisol,1:nxsol)
   write(1)xnorwall(1:nphisol,1:nxsol)
   write(1)ynorwall(1:nphisol,1:nxsol)
   write(1)znorwall(1:nphisol,1:nxsol)

!   write(1)rwall(1:nphisol,1:nxsol)


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


   write(1)nnode
   write(1)xnode(1:nnode)
   write(1)ynode(1:nnode)
   write(1)znode(1:nnode)

   nfa=nfa1+nfa2

   write(1)nfa
   write(1)xfa(1:nfa)
   write(1)yfa(1:nfa)
   write(1)zfa(1:nfa)

   write(1)nmyo
   write(1)xmyo(1:nmyo)
   write(1)ymyo(1:nmyo)
   write(1)zmyo(1:nmyo)

   write(1)nlk
   write(1)xlk(1:nlk)
   write(1)ylk(1:nlk)
   write(1)zlk(1:nlk)

   close(1)

!  write configuration

   open(1,file=fileconf)

   write(1,*)'VISUAL'

   write(1,*)natom,natom_a1,natom_a2,natom_a2m

   write(1,*)'CELLWALL'
   write(1,*)nwall,nxsol,nphisol
   write(1,*)

   write(1,*)'MEMBRANE'
   write(1,*)nmemb
   write(1,*)

   write(1,*)'NODES'

   write(1,*)nnode
   do n=1,nnode
      write(1,*)memb2node(1:nmb2node,n)
   end do
   write(1,*)

   fanum=fanum1+fanum2

   write(1,*)'FACTIN'
   write(1,*)nact1,nact2,nfamax1,nfamax2,nfa1,nfa2,fanum1,fanum2,fanummax1,fanummax2
   write(1,*)filid(1:nfa)
!   write(1,*)a2mem(1:nfa)

   write(1,*)apar(1:2,1:nfa)

   write(1,*)apos(1:nfa)

!   write(1,*)fa1stbound(1:fanum)

   write(1,*)alen(1:fanum)

   write(1,*)astart(1:fanum)

   write(1,*)a2mem(1:fanum)


   write(1,*)'MYOSIN configuration'
   write(1,*)nmyo,myonum,nmyoturn
   do n=1,myonum
      write(1,*)my2node(n),myhead(1:2,n),mybody(1:mybodylen,n)
      write(1,*)mytyp(1:2,n),fa2myo(1:2,n)
   end do
   write(1,*)


   crlknum=crlknum1+crlknum2

   write(1,*)'CROSSLINKER'
   write(1,*)nlk,crlknum1,crlknum2,crlknummax,crlknumactive

   write(1,*)lkstart(1:crlknum)
   write(1,*)lktyp(1:crlknum)
   write(1,*)fa2lk(1,1:crlknum)
   write(1,*)fa2lk(2,1:crlknum)

!   do n=1,crlknum1+crlknum2
!      write(1,*)crlk(1:lklen,n)
!      write(1,*)fa2lk(1:2,n)
!   end do
!   write(1,*)


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

   subroutine neighbors(neinum5,mynei5,lknei5,xfa,yfa,zfa,astart,filid,alen,mybodylen, &
              mybody,myonum,xmyo,ymyo,zmyo,myhead,mynei,neinum,crlknum,lkstart,lklen,xlk,ylk,zlk,lknei)

   implicit none

   integer,value:: myonum,neinum,crlknum,lklen,neinum5,mybodylen
   integer nm,jmh,jmyo,jnei,ja,nl,jlk,jn,j1,j2,nf,jnode

   integer,allocatable,intent(in),dimension(:)::astart,filid,alen,lkstart
   integer,allocatable,intent(in),dimension(:,:)::myhead,mybody
   integer,allocatable,dimension(:,:,:)::mynei,lknei
   integer,allocatable,dimension(:,:,:),intent(in)::mynei5,lknei5
   double precision dx,dy,dz,d2,d2max,dx1,dy1,dz1,dx2,dy2,dz2

   double precision,allocatable,intent(in),dimension(:)::xfa,yfa,zfa,xmyo,ymyo,zmyo,xlk,ylk,zlk

   d2max=15.0d0*15.0d0
   mynei=0
   lknei=0

!$omp parallel &
!$omp default(none) &
!$omp private(nm,jmh,jmyo,jnei,jn,ja,dx,dy,dz,d2,dx1,dy1,dz1,dx2,dy2,dz2,jnode,j1,j2,nf) &
!$omp shared(myonum,myhead,xfa,yfa,zfa,xmyo,ymyo,zmyo,filid,astart) &
!$omp shared(mynei,neinum,mynei5,neinum5,d2max,mybodylen,mybody,alen) &

!$omp private(nl,jlk) &
!$omp shared(lkstart,crlknum,lklen,xlk,ylk,zlk,lknei,lknei5)

!$omp do schedule(guided,32)

   do nm=1,myonum


!      jnode=my2node(nm)

      j1=mybody(1,nm)

      j2=mybody(mybodylen,nm)

      dx1=xmyo(j1)-xmyo(j2)
      dy1=ymyo(j1)-ymyo(j2)
      dz1=zmyo(j1)-zmyo(j2)

      do jmh=1,2

         jmyo=myhead(jmh,nm)


         jnei=0

         do jn=1,neinum5

            if(mynei5(jn,jmh,nm)==0)then
               exit
            end if

            ja=mynei5(jn,jmh,nm)

            nf=filid(ja)

!            if(a2node(nf)==jnode.or.ja==astart(nf))then
!               cycle
!            end if


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


            dx=xfa(ja)-xmyo(jmyo)
            dy=yfa(ja)-ymyo(jmyo)
            dz=zfa(ja)-zmyo(jmyo)

            d2=dx*dx+dy*dy+dz*dz

            if(d2<d2max)then
               jnei=jnei+1
               mynei(jnei,jmh,nm)=ja
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

   subroutine dcdheader(junit,jfile,ntotal)

   implicit none

   character coor*4,zero*1,charid1*1,charid2*2,charid3*3
   character (len=80) string1,string2,filedcd
   integer ifirst,nframe,nfreq,zeros5(5),peroff,zeros7(7),two,twentyfour,ntot
   integer,value:: junit,jfile,ntotal
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
   WRITE(JUNIT)NTOTAL

   end subroutine

!=========================================================

   subroutine myocycle(natp,myonum,neinum,nrand,mybodylen,jrand,apar,apos,filid, &
              mytyp,fa2myo,myhead,mybody,mynei,p1_hydr,p2_bind,p3_pi_rele,p4_adp_rele, &
              p5_ubind,dt,invmylen2,lbind,rands,xmyo,ymyo,zmyo,xfa,yfa,zfa)

   implicit none

   integer,value::myonum,neinum,nrand,mybodylen
   integer natp,jrand

   integer nm,jh,jtem,omp_get_thread_num,jmyo,ncan,jn,ja,jf,j,jexit,j1,j2,seed,iseed
   integer,allocatable,dimension(:),intent(in)::apos,filid
   integer,allocatable,dimension(:,:)::mytyp,fa2myo,apar
   integer,allocatable,dimension(:,:),intent(in)::myhead,mybody
   integer,allocatable,dimension(:,:,:),intent(in)::mynei
   integer acan(neinum)!,jtem0(nthreads)

   real(kind=8)::d2max,dx,dy,dz,d2,psum,rtem,rat0,normal,ratio,pact,lbind2
   real(kind=8)::pcan(neinum)
   real(kind=8),value::p1_hydr,p2_bind,p3_pi_rele,p4_adp_rele,p5_ubind,dt,invmylen2,lbind
   real(kind=8),allocatable,dimension(:),intent(in)::rands,xmyo,ymyo,zmyo,xfa,yfa,zfa


!     changing myosin head status:

!   do tid=1,nthreads
!      jtem0(tid)=jrand+2*myonum+tid-1
!   end do

   rat0=0.8

   normal=1.0d0/(1.0d0-rat0)

   normal=normal*normal

   lbind2=lbind*lbind

   CALL SYSTEM_CLOCK(COUNT=seed)

!$omp parallel &
!$omp default(none) &
!$omp private(nm,jh,jtem,jmyo,ncan,jn,ja,jf,j,jexit,acan,pcan,j1,j2,iseed) &
!$omp private(d2max,dx,dy,dz,d2,psum,rtem,ratio,pact) &
!$omp shared(myonum,neinum,rat0,normal,nrand,seed) &
!$omp shared(jrand,apos,filid,mytyp,fa2myo,apar,myhead,mynei) &
!$omp shared(p1_hydr,p2_bind,p3_pi_rele,p4_adp_rele,p5_ubind,dt,lbind2) &
!$omp shared(rands,xmyo,ymyo,zmyo,xfa,yfa,zfa,invmylen2,mybody,mybodylen) &
!$omp reduction(+:natp)


   iseed=seed+omp_get_thread_num( )

!$omp do schedule (guided,32)

   do nm=1,myonum

      do jh=1,2

         jtem=(nm-1)*2+jh+jrand

         if(jtem<=nrand)then

            rtem=rands(jtem)
         else

            rtem=r8_uniform_01(iseed)

!            rtem=0.0d0
         end if


         if(mytyp(jh,nm)==1)then

            j1=mybody(1,nm)

            j2=mybody(mybodylen,nm)

            dx=xmyo(j1)-xmyo(j2)
            dy=ymyo(j1)-ymyo(j2)
            dz=zmyo(j1)-zmyo(j2)

            d2=dx*dx+dy*dy+dz*dz

            ratio=d2*invmylen2

            if(ratio>rat0)then
               pact=normal*(ratio-rat0)**2
            else
               pact=0.0d0
            endif


            if(pact*p1_hydr*dt>rtem)then
               mytyp(jh,nm)=2
!               natp=natp+1

            end if

         else if(mytyp(jh,nm)==2)then

            if(p2_bind*dt>rtem)then

!              max binding distance square:

!               d2max=100.0d0
               d2max=15.0d0*15.0d0

               jmyo=myhead(jh,nm)

               ncan=0

               do jn=1,neinum

                  if(mynei(jn,jh,nm)==0)then
                     exit
                  end if

                  ja=mynei(jn,jh,nm)

                  if(apar(1,ja)/=0)then
                     cycle
                  end if

                  jf=filid(ja)

!                  if(apos(ja)<=fa1stbound(jf).and.tension2==1)then
!                     cycle
!                  end if

                  dx=xfa(ja)-xmyo(jmyo)
                  dy=yfa(ja)-ymyo(jmyo)
                  dz=zfa(ja)-zmyo(jmyo)

                  d2=dx*dx+dy*dy+dz*dz

                  if(d2<d2max)then

                     ncan=ncan+1

                     pcan(ncan)=lbind2*(d2max-d2)/d2/(d2max-lbind2)!1.0d0-d2/d2max

                     pcan(ncan)=min(1.0d0,pcan(ncan))

                     acan(ncan)=ja

                  end if

               end do

               if(ncan>0)then

                  !omp critical

!                  tid=omp_get_thread_num()+1
!                  jtem0(tid)=jtem0(tid)+nthreads
!                  jtem=jtem0(tid)!jrand+4*myonum+1

                  !omp end critical

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

                           mytyp(jh,nm)=3
                           fa2myo(jh,nm)=ja
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

         else if(mytyp(jh,nm)==3)then

            if(p3_pi_rele*dt>rtem)then
               mytyp(jh,nm)=4
               natp=natp+1
            end if

         else if(mytyp(jh,nm)==4)then

            if(p4_adp_rele*dt>rtem)then
               mytyp(jh,nm)=5
            end if

         else

            if(p5_ubind*dt>rtem)then
               mytyp(jh,nm)=1
               ja=fa2myo(jh,nm)
               apar(1:2,ja)=0
               fa2myo(jh,nm)=0
            end if

         end if

      end do

   end do

!$omp end do nowait
!$omp end parallel

   jrand=jrand+2*myonum

   end subroutine
!=========================================================

   subroutine crlkcycle(jrand,nrand,crlknum,crlknum1,fanum,lklen,neinum,fa2lk,apar, &
              lkstart,lktyp,filid,lknei,plk_ubind1,plk_ubind2, &
              dt,plk_bind,lbind,rands,xlk,ylk,zlk,xfa,yfa,zfa)

   implicit none

   integer nl,jl,jtem,omp_get_thread_num,ja,jlk,filid0,ncan,jn,j,jexit,seed,iseed
   integer jrand
   integer,value::nrand,crlknum,crlknum1,fanum,lklen,neinum
   integer,allocatable,dimension(:,:)::fa2lk,apar
   integer,allocatable,dimension(:),intent(in)::filid,lkstart!,lktyp
   integer,allocatable,dimension(:)::lktyp!fa1stbound
   integer,allocatable,dimension(:,:,:),intent(in)::lknei
   integer acan(neinum)!,jtem0(nthreads)

   real(kind=8)::rtem,plk_off,d2max,dx,dy,dz,d2,psum,lbind2
   real(kind=8),value::plk_ubind1,plk_ubind2,dt,plk_bind,lbind
   real(kind=8),allocatable,dimension(:),intent(in)::rands,xlk,ylk,zlk,xfa,yfa,zfa
   real(kind=8) pcan(neinum)

!     update binding status of crosslinkers to actin

!   jchange=0

!   do tid=1,nthreads
!      jtem0(tid)=jrand+2*crlknum+tid-1
!   end do


   CALL SYSTEM_CLOCK(COUNT=seed)

   lbind2=lbind*lbind


!$omp parallel &
!$omp default(none) &
!$omp private(nl,jl,jtem,ja,jlk,filid0,ncan,jn,j,jexit,iseed) &
!$omp private(acan,pcan,rtem,plk_off,d2max,dx,dy,dz,d2,psum) &
!$omp shared(jrand,nrand,crlknum,crlknum1,lklen,neinum,seed) &
!$omp shared(fa2lk,apar,lkstart,lktyp,filid,lknei) &
!$omp shared(plk_ubind1,plk_ubind2,dt,plk_bind,lbind2,rands,xlk,ylk,zlk,xfa,yfa,zfa) 
!omp reduction(+:jchange)

   iseed=seed+omp_get_thread_num( )

!$omp do schedule (guided, 32)


   do nl=1,crlknum

      if(lktyp(nl)==0) cycle

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

               fa2lk(jl,nl)=0

!-------------------------------------
! no need to update fa1stbound
!              update 1st crosslinked bead on a filament:

!               jf=filid(ja)

!               jstart=astart(jf)


!               if(apos(ja)==fa1stbound(jf))then

!                  jchange=jchange+1

                  !omp critical

!                  fa1stbound(jf)=-1

                  !omp end critical

!               end if


!               fa2lk(jl,nl)=0

!              update first crosslinked bead on the former crosslinking partner filament:

!               if(jl==1)then
!                  ja=fa2lk(2,nl)
!               else
!                  ja=fa2lk(1,nl)
!               end if

!               if(ja>0)then

!                  jf=filid(ja)


!                  if(apos(ja)==fa1stbound(jf))then

!                     jchange=jchange+1

                     !omp critical

!                     fa1stbound(jf)=-1

                     !omp end critical


!                  end if

!               end if
!--------------------------------------------------------

            end if

!        crosslink binds to filament:

         else

            if(plk_bind*dt>rtem)then

!              max binding distance square:

!               d2max=100.0d0
               d2max=15.0d0*15.0d0

               if(jl==1)then

                  jlk=lkstart(nl) !crlk(1,nl)

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

                     acan(ncan)=ja

                  end if

               end do

               if(ncan==0)then
                  cycle
               end if

               !omp critical
!               tid=omp_get_thread_num()+1
!               jtem0(tid)=jtem0(tid)+nthreads!jrand+2*crlknum+1
!               jtem=jtem0(tid)

               !omp end critical

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


!----------------------------------------
! no need to update fa1stbound

!                        if(jl==1)then
!                           ja1=fa2lk(2,nl)
!                        else
!                           ja1=fa2lk(1,nl)
!                        end if


!                        if(ja1>0)then

!                           jf=filid(ja)

!                           if(apos(ja)<fa1stbound(jf))then
!                              fa1stbound(jf)=apos(ja)
!                           end if

!                           jf=filid(ja1)

!                           if(apos(ja1)<fa1stbound(jf))then
!                              fa1stbound(jf)=apos(ja1)
!                           end if


!                        end if
!------------------------------------------

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

   jrand=jrand+2*crlknum

!   if(jchange==0)then
!      return
!   end if


!omp parallel &
!omp default(none) &
!omp private(jf,jstart,j,nl) &
!omp shared(fanum,fa1stbound,astart,alen,apar,fa2lk)
!omp do schedule (guided, 32)

!   do jf=1,fanum

!      if(fa1stbound(jf)>0)then
!         cycle
!      end if

!      fa1stbound(jf)=alen(jf)

!      jstart=astart(jf)


!      do j=1,alen(jf)

!         if(apar(2,jstart+j-1)<0)then

!            nl=-apar(2,jstart+j-1)

!            if(fa2lk(1,nl)>0.and.fa2lk(2,nl)>0)then
!               fa1stbound(jf)=j
!               exit
!            end if

!         end if

!      end do

!   end do

!omp end do nowait
!omp end parallel 


   end subroutine

!=========================================================

   subroutine tethering(nnode,nmb2node,memb2node,fanum,astart,a2mem,myonum,mybody,my2node, &
               fxnode,fynode,fznode,fxmemb,fymemb,fzmemb,fxfa,fyfa,fzfa,fxmyo,fymyo,fzmyo, &
                xnode,ynode,znode,xmemb,ymemb,zmemb,xfa,yfa,zfa,xmyo,ymyo,zmyo,knode,lnode,lnode2, &
                 ktether,ltether,ltether2,npairnode_mm,pairnode_mm)

   implicit none

   integer,value::nnode,fanum,myonum,nmb2node,npairnode_mm
   integer jnode,jmb,ja,jmy,j,nf,nm,n1,n2,nm1,nm2,jm1,jm2,j1,j2,n

   integer,allocatable,dimension(:),intent(in)::a2mem,my2node,astart
   integer,allocatable,dimension(:,:),intent(in)::memb2node,mybody,pairnode_mm

   double precision,value::ktether,ltether,ltether2,knode,lnode,lnode2
   double precision dx,dy,dz,f,fx,fy,fz,dist,dist2,radnode,rad
   double precision xno,yno,zno,xmb,ymb,zmb

   double precision l2max,lmax
   double precision d2min,dxmin,dymin,dzmin,dx0,dy0,dz0,dx1,dy1,dz1,dx2,dy2,dz2


   double precision,allocatable,dimension(:),intent(in)::xnode,ynode,znode,xmemb,ymemb,zmemb
   double precision,allocatable,dimension(:)::fxnode,fynode,fznode,fxmemb,fymemb,fzmemb
   double precision,allocatable,dimension(:),intent(in)::xfa,yfa,zfa,xmyo,ymyo,zmyo
   double precision,allocatable,dimension(:)::fxfa,fyfa,fzfa,fxmyo,fymyo,fzmyo

   fxnode=0.0d0
   fynode=0.0d0
   fznode=0.0d0

   lmax=10.0d0

   l2max=lmax*lmax


!$omp parallel &
!$omp default(none) &


!$omp shared(nnode,fanum,myonum) &
!$omp private(jnode,jmb,ja,jmy,j,nf,nm) &

!$omp shared(a2mem,my2node,astart,nmb2node,memb2node,mybody,ktether,ltether,ltether2,knode,lnode,lnode2) &
!$omp private(dx,dy,dz,f,fx,fy,fz,dist,dist2,radnode,rad,xno,yno,zno,xmb,ymb,zmb) &

!$omp shared(xnode,ynode,znode,xmemb,ymemb,zmemb) &
!$omp shared(fxnode,fynode,fznode,fxmemb,fymemb,fzmemb) &
!$omp shared(xfa,yfa,zfa,xmyo,ymyo,zmyo,fxfa,fyfa,fzfa,fxmyo,fymyo,fzmyo) &

!$omp shared(npairnode_mm,pairnode_mm,l2max,lmax) &
!$omp private(n1,n2,nm1,nm2,jm1,jm2,j1,j2,n) &
!$omp private(d2min,dxmin,dymin,dzmin,dx0,dy0,dz0,dx1,dy1,dz1,dx2,dy2,dz2)


!  nodes tethered to membrane:

!$omp do

   do jnode=1,nnode

      xno=xnode(jnode)
      yno=ynode(jnode)
      zno=znode(jnode)

      radnode=sqrt(yno*yno+zno*zno)

      do j=1,2!nmb2node

         jmb=memb2node(j,jnode)

         if(jmb==0)then
            exit
         end if

         xmb=xmemb(jmb)
         ymb=ymemb(jmb)
         zmb=zmemb(jmb)

         dx=xno-xmb
         dy=yno-ymb
         dz=zno-zmb

         dist2=dx*dx+dy*dy+dz*dz

         if(dist2>lnode2)then

            dist=sqrt(dist2)

            f=knode*(dist-lnode)/dist


            fx=f*dx!*invdist
            fy=f*dy!*invdist
            fz=f*dz!*invdist

            fxnode(jnode)=fxnode(jnode)-fx
            fynode(jnode)=fynode(jnode)-fy
            fznode(jnode)=fznode(jnode)-fz

            fxmemb(jmb)=fxmemb(jmb)+fx
            fymemb(jmb)=fymemb(jmb)+fy
            fzmemb(jmb)=fzmemb(jmb)+fz

         end if

         rad=sqrt(ymb*ymb+zmb*zmb)

         f=knode*(radnode-rad)/radnode

         fy=f*yno
         fz=f*zno

         fynode(jnode)=fynode(jnode)-fy
         fznode(jnode)=fznode(jnode)-fz

         fymemb(jmb)=fymemb(jmb)+fy
         fzmemb(jmb)=fzmemb(jmb)+fz

      end do
   end do

!$omp enddo nowait


!  actin tethered to membrane:

!$omp do schedule (guided,32)

   do nf=1,fanum

      jmb=a2mem(nf)

      ja=astart(nf)

      dx=xmemb(jmb)-xfa(ja)
      dy=ymemb(jmb)-yfa(ja)
      dz=zmemb(jmb)-zfa(ja)

      dist2=dx*dx+dy*dy+dz*dz

      if(dist2>ltether2)then

         dist=sqrt(dist2)

         f=ktether*(dist-ltether)/dist


         fx=f*dx!*invdist
         fy=f*dy!*invdist
         fz=f*dz!*invdist

         fxmemb(jmb)=fxmemb(jmb)-fx
         fymemb(jmb)=fymemb(jmb)-fy
         fzmemb(jmb)=fzmemb(jmb)-fz

         fxfa(ja)=fxfa(ja)+fx
         fyfa(ja)=fyfa(ja)+fy
         fzfa(ja)=fzfa(ja)+fz

      end if

   end do

!$omp enddo nowait


!  myosin tethered to nodes:

!$omp do schedule (guided,32)

   do nm=1,myonum

      jnode=my2node(nm)

      jmy=mybody(1,nm)

      dx=xnode(jnode)-xmyo(jmy)
      dy=ynode(jnode)-ymyo(jmy)
      dz=znode(jnode)-zmyo(jmy)

      dist2=dx*dx+dy*dy+dz*dz

      if(dist2>lnode2)then

         dist=sqrt(dist2)

         f=knode*(dist-lnode)/dist


         fx=f*dx!*invdist
         fy=f*dy!*invdist
         fz=f*dz!*invdist

!$omp atomic
         fxnode(jnode)=fxnode(jnode)-fx
!$omp atomic
         fynode(jnode)=fynode(jnode)-fy
!$omp atomic
         fznode(jnode)=fznode(jnode)-fz

         fxmyo(jmy)=fxmyo(jmy)+fx
         fymyo(jmy)=fymyo(jmy)+fy
         fzmyo(jmy)=fzmyo(jmy)+fz

      end if

   end do

!$omp enddo nowait


!--------------------------------------

!  to prevent tethers from slipping past each other:


!$omp do

   do n=1,npairnode_mm

      nm1=pairnode_mm(1,n)

      nm2=pairnode_mm(2,n)

      jm1=mybody(1,nm1)

      jm2=mybody(1,nm2)

      n1=my2node(nm1)

      n2=my2node(nm2)

      dx1=0.1d0*(xmyo(jm1)-xnode(n1))
      dy1=0.1d0*(ymyo(jm1)-ynode(n1))
      dz1=0.1d0*(zmyo(jm1)-znode(n1))

      dx2=0.1d0*(xmyo(jm2)-xnode(n2))
      dy2=0.1d0*(ymyo(jm2)-ynode(n2))
      dz2=0.1d0*(zmyo(jm2)-znode(n2))

      d2min=l2max+1.0d0

      dx0=xnode(n1)-xnode(n2)
      dy0=ynode(n1)-ynode(n2)
      dz0=znode(n1)-znode(n2)

      do j1=1,11

         do j2=1,11

            dx=dx0+(j1-1)*dx1-(j2-1)*dx2
            dy=dy0+(j1-1)*dy1-(j2-1)*dy2
            dz=dz0+(j1-1)*dz1-(j2-1)*dz2

            dist2=dx*dx+dy*dy+dz*dz

            if(dist2<d2min)then

               d2min=dist2

               dxmin=dx
               dymin=dy
               dzmin=dz

            end if

         end do

      end do

      if(d2min<l2max) then

         dist=sqrt(d2min)

         f=10*(lmax-dist)/(d2min+1.0d0)/dist

!         f=1.0d0/(d2min+1.0d0)

         fx=f*dxmin
         fy=f*dymin
         fz=f*dzmin

!$omp atomic
         fxnode(n1)=fxnode(n1)+fx
!$omp atomic
         fynode(n1)=fynode(n1)+fy
!$omp atomic
         fznode(n1)=fznode(n1)+fz

!$omp atomic
         fxnode(n2)=fxnode(n2)-fx
!$omp atomic
         fynode(n2)=fynode(n2)-fy
!$omp atomic
         fznode(n2)=fznode(n2)-fz

!$omp atomic
         fxmyo(jm1)=fxmyo(jm1)+fx
!$omp atomic
         fymyo(jm1)=fymyo(jm1)+fy
!$omp atomic
         fzmyo(jm1)=fzmyo(jm1)+fz

!$omp atomic
         fxmyo(jm2)=fxmyo(jm2)-fx
!$omp atomic
         fymyo(jm2)=fymyo(jm2)-fy
!$omp atomic
         fzmyo(jm2)=fzmyo(jm2)-fz


      end if

   end do

!$omp enddo nowait


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


         n2=n2+1!jstart+jf-1  ! used to be afil(jf,nf)

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


         n3=n3+1!jstart+jf  ! used to be afil(jf+1,nf)

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

   subroutine myoforce(jforce,fxbond,fybond,fzbond,fxangl,fyangl,fzangl, &
                       myonum,myhead,mytyp,mybody,mybodylen,xmyo,ymyo,zmyo,k_a,l_mh,l_mb, &
                       kmh_thet,kmb_thet,thet_mh1,thet_mh2,thet_mb,delta,invdelta,beta)

   implicit none

   integer,value:: jforce,myonum,mybodylen
   integer nm,jm,n1,n2,n3

   integer,allocatable,intent(in),dimension(:,:)::myhead,mybody,mytyp

   double precision,value::k_a,l_mh,l_mb,kmh_thet,kmb_thet,thet_mh1,thet_mh2,thet_mb,delta,invdelta,beta
   double precision dx,dy,dz,dx1,dy1,dz1,dx3,dy3,dz3,invdist1,invdist3,thet0
   double precision f0,cos_t,cos_t0,thet,drep,dfx1,dfy1,dfz1,dfx3,dfy3,dfz3,f,fx,fy,fz
   double precision x1,y1,z1,x2,y2,z2,x3,y3,z3,dist

   double precision,allocatable,intent(in),dimension(:)::xmyo,ymyo,zmyo
   double precision,allocatable,dimension(:)::fxangl,fyangl,fzangl,fxbond,fybond,fzbond



!$omp parallel &
!$omp default(none) &
!$omp private(nm,jm,n1,n2,dx,dy,dz,dist,f,x1,y1,z1,x2,y2,z2) &
!$omp shared(myonum,mybodylen,myhead,mybody) &
!$omp shared(k_a,l_mh,l_mb) &
!$omp shared(xmyo,ymyo,zmyo,fxbond,fybond,fzbond)

!$omp do schedule(guided,32)

   do nm=1,myonum

!     force on the body:

      n2=mybody(1,nm)

      fxbond(n2)=0.0d0
      fybond(n2)=0.0d0
      fzbond(n2)=0.0d0

      x2=xmyo(n2)
      y2=ymyo(n2)
      z2=zmyo(n2)

      do jm=2,mybodylen

         x1=x2
         y1=y2
         z1=z2

         n1=n2


         n2=n2+1!mybody(jm,nm)

         x2=xmyo(n2)
         y2=ymyo(n2)
         z2=zmyo(n2)


         dx=x1-x2
         dy=y1-y2
         dz=z1-z2

         dist=sqrt(dx*dx+dy*dy+dz*dz)

!         invdist=1.0d0/dist

         f=k_a*(dist-l_mb)/dist

         fxbond(n2)=f*dx!*invdist
         fybond(n2)=f*dy!*invdist
         fzbond(n2)=f*dz!*invdist

         fxbond(n1)=fxbond(n1)-fxbond(n2)
         fybond(n1)=fybond(n1)-fybond(n2)
         fzbond(n1)=fzbond(n1)-fzbond(n2)

!         fxbond(n2)=fxbond(n2)+fx
!         fybond(n2)=fybond(n2)+fy
!         fzbond(n2)=fzbond(n2)+fz

      end do

!     force on the head:

      n2=mybody(mybodylen,nm)

      do jm=1,2

         n1=n2+jm!myhead(jm,nm)

!         n2=mybody(mybodylen,nm)


         dx=xmyo(n1)-xmyo(n2)
         dy=ymyo(n1)-ymyo(n2)
         dz=zmyo(n1)-zmyo(n2)

         dist=sqrt(dx*dx+dy*dy+dz*dz)

!         invdist=1.0d0/dist

         f=-k_a*(dist-l_mh)/dist

         fxbond(n1)=f*dx!*invdist
         fybond(n1)=f*dy!*invdist
         fzbond(n1)=f*dz!*invdist

         fxbond(n2)=fxbond(n2)-fxbond(n1)
         fybond(n2)=fybond(n2)-fybond(n1)
         fzbond(n2)=fzbond(n2)-fzbond(n1)

      end do


!      do jm=3,4

!         n1=myhead(jm,nm)

!         n2=mybody(mybodylen,nm)


!         dx=xmyo(n1)-xmyo(n2)
!         dy=ymyo(n1)-ymyo(n2)
!         dz=zmyo(n1)-zmyo(n2)

!         dist=sqrt(dx*dx+dy*dy+dz*dz)


!         f=-k_a*(dist-l_mh)/dist

!         fxbond(n1)=f*dx!*invdist
!         fybond(n1)=f*dy!*invdist
!         fzbond(n1)=f*dz!*invdist

!         fxbond(n2)=fxbond(n2)-fxbond(n1)
!         fybond(n2)=fybond(n2)-fybond(n1)
!         fzbond(n2)=fzbond(n2)-fzbond(n1)

!      end do

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
!$omp shared(myonum,mybodylen,myhead,mybody,mytyp) &
!$omp shared(kmh_thet,kmb_thet,thet_mh1,thet_mh2,thet_mb,delta,invdelta,beta) &
!$omp shared(xmyo,ymyo,zmyo,fxangl,fyangl,fzangl)

!$omp do schedule(guided,32)


   do nm=1,myonum

      n2=mybody(1,nm)

      fxangl(n2)=0.0d0
      fyangl(n2)=0.0d0
      fzangl(n2)=0.0d0

      n3=mybody(2,nm)

      fxangl(n3)=0.0d0
      fyangl(n3)=0.0d0
      fzangl(n3)=0.0d0

      x3=xmyo(n3)
      y3=ymyo(n3)
      z3=zmyo(n3)

      dx3=x3-xmyo(n2)
      dy3=y3-ymyo(n2)
      dz3=z3-zmyo(n2)

      invdist3=1.0d0/sqrt(dx3*dx3+dy3*dy3+dz3*dz3)

      do jm=2,mybodylen-1


         dx1=-dx3
         dy1=-dy3
         dz1=-dz3

         invdist1=invdist3

         x2=x3
         y2=y3
         z2=z3

         n1=n2

         n2=n3


         n3=n3+1!mybody(jm+1,nm)

         x3=xmyo(n3)
         y3=ymyo(n3)
         z3=zmyo(n3)

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

         fxangl(n1)=fxangl(n1)+dfx1

!        force on n1 along y:

         dy=dy1+delta

         drep=sqrt(dx1*dx1+dy*dy+dz1*dz1)

         cos_t=(dx1*dx3+dy*dy3+dz1*dz3)/drep*invdist3


         dfy1=f0*(cos_t-cos_t0)

         fyangl(n1)=fyangl(n1)+dfy1

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

         fxangl(n3)=dfx3

!        force on n3 along y:

         dy=dy3+delta

         drep=sqrt(dx3*dx3+dy*dy+dz3*dz3)

         cos_t=(dx1*dx3+dy1*dy+dz1*dz3)/drep*invdist1

         dfy3=f0*(cos_t-cos_t0)

         fyangl(n3)=dfy3

!        force on n3 along z:

         dz=dz3+delta

         drep=sqrt(dx3*dx3+dy3*dy3+dz*dz)

         cos_t=(dx1*dx3+dy1*dy3+dz1*dz)/drep*invdist1

         dfz3=f0*(cos_t-cos_t0)

         fzangl(n3)=dfz3

!        forces on n2:

         fxangl(n2)=fxangl(n2)-dfx1-dfx3
         fyangl(n2)=fyangl(n2)-dfy1-dfy3
         fzangl(n2)=fzangl(n2)-dfz1-dfz3

      end do

!     force on the head:

      n2=mybody(mybodylen,nm)
      n3=n2-1

      do jm=1,2

         n1=n2+jm!myhead(jm,nm)


!         n2=mybody(1,nm)

!         n3=mybody(2,nm)


         dx1=xmyo(n1)-xmyo(n2)
         dy1=ymyo(n1)-ymyo(n2)
         dz1=zmyo(n1)-zmyo(n2)

         invdist1=1.0d0/sqrt(dx1*dx1+dy1*dy1+dz1*dz1)


         if(mytyp(jm,nm)==2.or.mytyp(jm,nm)==3)then
            thet0=thet_mh2
         else
            thet0=thet_mh1
         end if

         dx3=xmyo(n3)-xmyo(n2)
         dy3=ymyo(n3)-ymyo(n2)
         dz3=zmyo(n3)-zmyo(n2)

         invdist3=1.0d0/sqrt(dx3*dx3+dy3*dy3+dz3*dz3)

         cos_t0=(dx1*dx3+dy1*dy3+dz1*dz3)*invdist1*invdist3

         thet=acos((1.0d0-beta)*cos_t0)

         f0=kmh_thet*(thet-thet0)/sin(thet)*invdelta

!        force on n1 along x:

         dx=dx1+delta

         drep=sqrt(dx*dx+dy1*dy1+dz1*dz1)

         cos_t=(dx*dx3+dy1*dy3+dz1*dz3)/drep*invdist3

         dfx1=f0*(cos_t-cos_t0)

         fxangl(n1)=dfx1

!        force on n1 along y:

         dy=dy1+delta

         drep=sqrt(dx1*dx1+dy*dy+dz1*dz1)

         cos_t=(dx1*dx3+dy*dy3+dz1*dz3)/drep*invdist3

         dfy1=f0*(cos_t-cos_t0)

         fyangl(n1)=dfy1

!        force on n1 along z:

         dz=dz1+delta

         drep=sqrt(dx1*dx1+dy1*dy1+dz*dz)

         cos_t=(dx1*dx3+dy1*dy3+dz*dz3)/drep*invdist3

         dfz1=f0*(cos_t-cos_t0)

         fzangl(n1)=dfz1

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



!      do jm=3,4

!         n1=myhead(jm,nm)

!         n2=mybody(mybodylen,nm)

!         n3=mybody(mybodylen-1,nm)


!         dx1=xmyo(n1)-xmyo(n2)
!         dy1=ymyo(n1)-ymyo(n2)
!         dz1=zmyo(n1)-zmyo(n2)

!         invdist1=1.0d0/sqrt(dx1*dx1+dy1*dy1+dz1*dz1)


!         if(mytyp(jm,nm)==2.or.mytyp(jm,nm)==3)then
!            thet0=thet_mh2
!         else
!            thet0=thet_mh1
!         end if

!         dx3=xmyo(n3)-xmyo(n2)
!         dy3=ymyo(n3)-ymyo(n2)
!         dz3=zmyo(n3)-zmyo(n2)

!         invdist3=1.0d0/sqrt(dx3*dx3+dy3*dy3+dz3*dz3)

!         cos_t0=(dx1*dx3+dy1*dy3+dz1*dz3)*invdist1*invdist3

!         thet=acos((1.0d0-beta)*cos_t0)

!         f0=kmh_thet*(thet-thet0)/sin(thet)*invdelta

!        force on n1 along x:

!         dx=dx1+delta

!         drep=sqrt(dx*dx+dy1*dy1+dz1*dz1)

!         cos_t=(dx*dx3+dy1*dy3+dz1*dz3)/drep*invdist3

!         dfx1=f0*(cos_t-cos_t0)

!         fxangl(n1)=dfx1

!        force on n1 along y:

!         dy=dy1+delta

!         drep=sqrt(dx1*dx1+dy*dy+dz1*dz1)

!         cos_t=(dx1*dx3+dy*dy3+dz1*dz3)/drep*invdist3

!         dfy1=f0*(cos_t-cos_t0)

!         fyangl(n1)=dfy1

!        force on n1 along z:

!         dz=dz1+delta

!         drep=sqrt(dx1*dx1+dy1*dy1+dz*dz)

!         cos_t=(dx1*dx3+dy1*dy3+dz*dz3)/drep*invdist3

!         dfz1=f0*(cos_t-cos_t0)

!         fzangl(n1)=dfz1

!        force on n3 along x:

!         dx=dx3+delta

!         drep=sqrt(dx*dx+dy3*dy3+dz3*dz3)

!         cos_t=(dx1*dx+dy1*dy3+dz1*dz3)/drep*invdist1

!         dfx3=f0*(cos_t-cos_t0)

!         fxangl(n3)=fxangl(n3)+dfx3

!        force on n3 along y:

!         dy=dy3+delta

!         drep=sqrt(dx3*dx3+dy*dy+dz3*dz3)

!         cos_t=(dx1*dx3+dy1*dy+dz1*dz3)/drep*invdist1

!         dfy3=f0*(cos_t-cos_t0)

!         fyangl(n3)=fyangl(n3)+dfy3

!        force on n3 along z:

!         dz=dz3+delta

!         drep=sqrt(dx3*dx3+dy3*dy3+dz*dz)

!         cos_t=(dx1*dx3+dy1*dy3+dz1*dz)/drep*invdist1

!         dfz3=f0*(cos_t-cos_t0)

!         fzangl(n3)=fzangl(n3)+dfz3

!        forces on n2:

!         fxangl(n2)=fxangl(n2)-dfx1-dfx3
!         fyangl(n2)=fyangl(n2)-dfy1-dfy3
!         fzangl(n2)=fzangl(n2)-dfz1-dfz3

!      end do



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


   subroutine fabindmyo(fxfa,fyfa,fzfa,xfa,yfa,zfa,fxmyo,fymyo,fzmyo,xmyo,ymyo,zmyo, &
                     myhead,fa2myo,myonum,kbind,lbind,invl_a)


   implicit none

   integer,value::myonum
   integer nm,jm,ja,jmyo

   integer,allocatable,intent(in),dimension(:,:)::myhead,fa2myo

   double precision,value::kbind,lbind,invl_a
   double precision dx,dy,dz,dist,f,fx,fy,fz,dx1,dy1,dz1

   double precision,allocatable,intent(in),dimension(:)::xfa,yfa,zfa,xmyo,ymyo,zmyo
   double precision,allocatable,dimension(:)::fxfa,fyfa,fzfa,fxmyo,fymyo,fzmyo

!integer na

!na=643

!$omp parallel &
!$omp default(none) &
!$omp private(nm,jm,ja,jmyo,dx,dy,dz,dx1,dy1,dz1,dist,f,fx,fy,fz) &
!$omp shared(myonum,myhead,fa2myo,kbind,lbind,invl_a) &
!$omp shared(xfa,yfa,zfa,xmyo,ymyo,zmyo,fxfa,fyfa,fzfa,fxmyo,fymyo,fzmyo)

!$omp do schedule(guided,32)

   do nm=1,myonum

      do jm=1,2

         ja=fa2myo(jm,nm)

         if(ja==0)then
            cycle
         end if

         jmyo=myhead(jm,nm)

         dx=xfa(ja)-xmyo(jmyo)
         dy=yfa(ja)-ymyo(jmyo)
         dz=zfa(ja)-zmyo(jmyo)

         dist=sqrt(dx*dx+dy*dy+dz*dz)

         f=kbind*(dist-lbind)/dist

!         invdist=1.0d0/dist

         fx=f*dx!*invdist
         fy=f*dy!*invdist
         fz=f*dz!*invdist

         fxfa(ja)=fxfa(ja)-fx
         fyfa(ja)=fyfa(ja)-fy
         fzfa(ja)=fzfa(ja)-fz

         fxmyo(jmyo)=fxmyo(jmyo)+fx
         fymyo(jmyo)=fymyo(jmyo)+fy
         fzmyo(jmyo)=fzmyo(jmyo)+fz

!if(ja==na)then
!print*,'myosin',jmyo
!print*,fxfa(ja),fyfa(ja),fzfa(ja)
!end if

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

         fxmyo(jmyo)=fxmyo(jmyo)+fx
         fymyo(jmyo)=fymyo(jmyo)+fy
         fzmyo(jmyo)=fzmyo(jmyo)+fz

!if(ja==na)then
!print*,'angle'
!print*,fxfa(ja),fyfa(ja),fzfa(ja)
!print*,dx,dy,dz
!print*,dx1,dy1,dz1
!end if


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

   integer,allocatable,intent(in),dimension(:,:)::fa2lk
  integer,allocatable,intent(in),dimension(:)::lkstart
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

   subroutine rforceset(jrforce,nrforce,nmemb,nfa,nmyo,nlk,pi,rxmemb,rymemb,rzmemb, &
              rxfa,ryfa,rzfa,rxmyo,rymyo,rzmyo,rxlk,rylk,rzlk,k_scale,k_scale_lk)

   implicit none

   integer,value::nrforce,nmemb,nfa,nmyo,nlk
   integer jrforce,j,nmax,n,nh1,nh2,omp_get_thread_num

   double precision,value::pi,k_scale,k_scale_lk
   double precision,allocatable,dimension(:,:)::rxfa,ryfa,rzfa,rxmyo,rymyo,rzmyo,rxlk,rylk,rzlk
   double precision,allocatable,dimension(:,:)::rxmemb,rymemb,rzmemb
   double precision,allocatable,dimension(:)::rh1,rh2,rx,ry,rz
!   double precision k,kmb
   integer seed,iseed


   nmax=nmemb+nfa+nmyo+nlk
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
!$omp shared(nrforce,nmax,nh1,nh2,nmemb,nfa,nmyo,nlk,pi,k_scale,k_scale_lk) &
!$omp shared(rxfa,ryfa,rzfa,rxmyo,rymyo,rzmyo,rxlk,rylk,rzlk,rxmemb,rymemb,rzmemb,seed) 

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

      rxmyo(1:nmyo,j)=k_scale_lk*rx(n+1:n+nmyo)
      rymyo(1:nmyo,j)=k_scale_lk*ry(n+1:n+nmyo)
      rzmyo(1:nmyo,j)=k_scale_lk*rz(n+1:n+nmyo)

      n=n+nmyo

      rxlk(1:nlk,j)=k_scale_lk*rx(n+1:n+nlk)
      rylk(1:nlk,j)=k_scale_lk*ry(n+1:n+nlk)
      rzlk(1:nlk,j)=k_scale_lk*rz(n+1:n+nlk)


   end do


!$omp enddo nowait
!$omp end parallel



   deallocate(rh1,rh2)


   end subroutine

!=========================================================

   subroutine newactin(nmemb,myonum,crlknum,falen,nmono,jdir,nxsol,nphisol,nnode,nmb2node, &
               fanum,fanum1,fanum2,nfa,nfa1,nfa2,jsursol,memb2node,apar,astart,alen, &
                fadist,apos,filid,a2mem,jfasol,fa2myo,fa2lk,l_a,l_mem, &
                 pi,ltether,xrmin,xrmax,xfa,yfa,zfa,xmemb,ymemb,zmemb, &
                  xnorsurf,ynorsurf,znorsurf,xsurf,ysurf,zsurf)

   implicit none

   integer,value::nmemb,myonum,crlknum,falen,nmono,jdir,nxsol,nphisol,nnode,nmb2node
   integer fanum,fanum1,fanum2,nfa,nfa1,nfa2

   integer length,na,nnew,nm,jm,nl,jl,jx0,jp0,j,jx,jp,jxget,jpget,j1,j2,n

   integer,allocatable,dimension(:,:),intent(in)::jsursol,memb2node

   integer,allocatable,dimension(:)::apar(:,:)
   integer,allocatable,dimension(:)::astart,alen,fadist,apos,filid,a2mem
   integer,allocatable,dimension(:,:)::jfasol,fa2myo,fa2lk

!   integer,allocatable,dimension(:,:)::apartem
   integer,allocatable,dimension(:)::mark
!   integer,allocatable,dimension(:,:)::jfasoltem

   real(kind=8),value::l_a,l_mem,pi,ltether,xrmin,xrmax

!   real(kind=8)::distmin,distmax

   real(kind=8)::d2max,dshift,r,xn0,yn0,zn0,distance,x0,y0,z0,phi,rad,dphi
   real(kind=8)::dist2,d2,dx,dy,dz,xn,yn,zn,proj,proj0

   real(kind=8),allocatable,dimension(:)::xfa,yfa,zfa
   real(kind=8),allocatable,dimension(:),intent(in)::xmemb,ymemb,zmemb
   real(kind=8),allocatable,dimension(:,:),intent(in)::xnorsurf,ynorsurf,znorsurf,xsurf,ysurf,zsurf
   real(kind=8),allocatable,dimension(:)::xfatem,yfatem,zfatem

   d2max=1e10!25*l_mem*l_mem!ltether*ltether

!   l_memby2=0.5d0*l_mem

!   dshift=0.1d0*l_a

!   xmin=-rwid/2
!   xmax=-xmin

!  defining ring boundary:

!   xrmin=1e10
!   xrmax=-xrmin

!   distmin=1e10
!   distmax=0.0d0

!   do na=1,nfa,10

!      xrmin=min(xrmin,xfa(na))

!      xrmax=max(xrmax,xfa(na))

!      jx=jfasol(1,na)

!      jp=jfasol(2,na)

!      xn=xnorsurf(jp,jx)
!      yn=ynorsurf(jp,jx)
!      zn=znorsurf(jp,jx)

!      dx=xfa(na)-xsurf(jp,jx)
!      dy=yfa(na)-ysurf(jp,jx)
!      dz=zfa(na)-zsurf(jp,jx)

!      proj=dx*xn+dy*yn+dz*zn

!      distmin=min(distmin,proj)

!      distmax=max(distmax,proj)

!   end do


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

         xfa(n+length)=xfa(n)
         yfa(n+length)=yfa(n)
         zfa(n+length)=zfa(n)

      end do


!      allocate(xfatem(nfa2),yfatem(nfa2),zfatem(nfa2))

!      xfatem(1:nfa2)=xfa(nfa1+1:nfa)
!      yfatem(1:nfa2)=yfa(nfa1+1:nfa)
!      zfatem(1:nfa2)=zfa(nfa1+1:nfa)

!      allocate(apartem(2,nfa2))
!      apartem(1:2,1:nfa2)=apar(1:2,nfa1+1:nfa)

!      allocate(fa1stboundtem(fanum2))
!      fa1stboundtem(1:fanum2)=fa1stbound(fanum1+1:fanum)

      do n=fanum,fanum1+1,-1

         astart(n+1)=astart(n)+length

         alen(n+1)=alen(n)

         a2mem(n+1)=a2mem(n)

      end do


!      allocate(alentem(fanum2))
!      alentem(1:fanum2)=alen(fanum1+1:fanum)

!      allocate(jfasoltem(2,nfa2),fadisttem(nfa2))
!      jfasoltem(1:2,1:nfa2)=jfasol(1:2,nfa1+1:nfa)
!      fadisttem(1:nfa2)=fadist(nfa1+1:nfa)

!      allocate(apostem(nfa2))
!      apostem(1:nfa2)=apos(nfa1+1:nfa)

!      do n=nfa,nfa1+1,-1

!         filid(n+length)=filid(n)+1

!      end do

!      filid(nfa1+1+length:nfa+length)=filid(nfa1+1:nfa)+1

!      allocate(a2memtem(nfa2))
!      a2memtem(1:fanum2)=a2mem(fanum1+1:fanum)

      na=nfa1

      nnew=fanum1+1

      do nm=1,myonum

         do jm=1,2

            if(fa2myo(jm,nm)>nfa1)then

               fa2myo(jm,nm)=fa2myo(jm,nm)+length

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

!   fa1stbound(nnew)=length

   fadist(na+1:na+length)=1

!   a2mem(na+1:na+length)=0

!  pick a membrane bead as a reference:

   allocate(mark(nmemb))

   mark=0

   mark(a2mem(1:fanum))=1

   do n=1,nnode
      do j=1,2!nmb2node
         if(memb2node(j,n)>0) mark(memb2node(j,n))=1
      end do
   end do


11 call random_number(r)

   jm=nmemb*r+1

   if(mark(jm)==1) goto 11

   if(xmemb(jm)<xrmin.or.xmemb(jm)>xrmax) goto 11

   jx0=jsursol(1,jm)

   jp0=jsursol(2,jm)

   xn0=xnorsurf(jp0,jx0)
   yn0=ynorsurf(jp0,jx0)
   zn0=znorsurf(jp0,jx0)

   deallocate(mark)

!  tethering F-actin to membrane:

   a2mem(nnew)=jm

!  pick a distance away from the membrane:

   call random_number(r)

   distance=(ltether-l_mem)*r+l_mem

!   distance=ringthick*r+l_mem
!   distance=(distmax-distmin)*r+distmin



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


   x0=xsurf(jp0,jx0)+distance*xn0
   y0=ysurf(jp0,jx0)+distance*yn0
   z0=zsurf(jp0,jx0)+distance*zn0



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

!      xfa(nfa1+1:nfa)=xfatem(1:nfa2)
!      yfa(nfa1+1:nfa)=yfatem(1:nfa2)
!      zfa(nfa1+1:nfa)=zfatem(1:nfa2)

!      apar(1:2,nfa1+1:nfa)=apartem(1:2,1:nfa2)

!      fa1stbound(fanum1+1:fanum)=fa1stboundtem(1:fanum2)

!      alen(fanum1+1:fanum)=alentem(1:fanum2)

!      jfasol(1:2,nfa1+1:nfa)=jfasoltem(1:2,1:nfa2)

!      fadist(nfa1+1:nfa)=fadisttem(1:nfa2)

!      apos(nfa1+1:nfa)=apostem(1:nfa2)

!      a2mem(fanum1+1:fanum)=a2memtem(1:fanum2)


!      deallocate(xfatem,yfatem,zfatem,apartem,alentem,jfasoltem,fadisttem,apostem,a2memtem)


   else

      fanum=nnew

      fanum2=fanum2+1

      nfa2=nfa2+length

      nfa=na

   end if




   end subroutine

!=========================================================


   subroutine newactin1(nmono,nxsol,nphisol,falen,nnode,nmb2node,fanum,nfa,jsursol,memb2node,apar,jfasol,astart,alen, &
              fadist,apos,filid,a2node,l_a,l_mem,rwid,ringthick,pi,xfa,yfa,zfa, &
              xmemb,ymemb,zmemb,xnorsurf,ynorsurf,znorsurf,xsurf,ysurf,zsurf)

   implicit none

   integer,value::nmono,nxsol,nphisol,falen,nnode,nmb2node
   integer fanum,nfa

   integer length,na,nnew,jm,jmax,jx0,jp0,j,jx,jp,jxget,jpget,j1,j2,n,jnode,jdir

   integer,allocatable,dimension(:,:),intent(in)::jsursol,memb2node

   integer,allocatable,dimension(:,:)::apar,jfasol
   integer,allocatable,dimension(:)::astart,alen,fadist,apos,filid,a2node

   integer,allocatable,dimension(:)::ncount

   real(kind=8),value::l_a,l_mem,rwid,ringthick,pi

   real(kind=8)::d2max,l_memby2,dshift,xmin,xmax,r,xn0,yn0,zn0,distance,x0,y0,z0,phi,rad,dphi
   real(kind=8)::dist2,d2,dx,dy,dz,xn,yn,zn,proj,proj0,psum

   real(kind=8),allocatable,dimension(:)::xfa,yfa,zfa
   real(kind=8),allocatable,dimension(:),intent(in)::xmemb,ymemb,zmemb
   real(kind=8),allocatable,dimension(:,:),intent(in)::xnorsurf,ynorsurf,znorsurf,xsurf,ysurf,zsurf

   real(kind=8),allocatable,dimension(:)::prob


   d2max=1e4!36*l_mem*l_mem

   l_memby2=0.5d0*l_mem

   dshift=0.1d0*l_a

   xmin=-rwid/2
   xmax=-xmin

!  pick a length for filament:

   call random_number(r)

   length=falen+(r-0.5)*(falen/2)

   if(length>nmono)then
      length=nmono
   end if

   na=nfa

   nnew=fanum+1

   alen(nnew)=length

   astart(nnew)=na+1

   filid(na+1:na+length)=nnew

   apar(1:2,na+1:na+length)=0
   apar(1,na+1)=-1

   fadist(na+1:na+length)=1

!  pick a mother node:

   allocate(prob(nnode),ncount(nnode))

   ncount=0

   do n=1,fanum-1

      jnode=a2node(n)

      ncount(jnode)=ncount(jnode)+1

   end do

   prob(1:nnode)=1.0d0/ncount(1:nnode)

   psum=sum(prob(1:nnode))

   prob(1:nnode)=prob(1:nnode)/psum


   call random_number(r)

   psum=0.0d0

   do n=1,nnode

      psum=psum+prob(n)

      if(psum>r)then

         jnode=n

         exit

      end if

   end do

   a2node(nnew)=jnode

!  pick a membrane bead as a reference:

   do j=1,nmb2node

      if(memb2node(j,jnode)==0)then
         jmax=j-1
         exit
      end if
   end do

   call random_number(r)

   j=jmax*r+1

   jm=memb2node(j,jnode)

   jx0=jsursol(1,jm)

   jp0=jsursol(2,jm)

   xn0=xnorsurf(jp0,jx0)
   yn0=ynorsurf(jp0,jx0)
   zn0=znorsurf(jp0,jx0)

!  pick a distance away from the membrane:

   call random_number(r)

   distance=ringthick*r+l_mem

   x0=xsurf(jp0,jx0)+distance*xn0
   y0=ysurf(jp0,jx0)+distance*yn0
   z0=zsurf(jp0,jx0)+distance*zn0

!  pick direction of elongation

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

!     update coordinates for the next bead:


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

      if(distance-proj0>l_memby2)then

         x0=x0+xn0*dshift
         y0=y0+yn0*dshift
         z0=z0+zn0*dshift

      elseif(proj0-distance>l_memby2)then

         x0=x0-xn0*dshift
         y0=y0-yn0*dshift
         z0=z0-zn0*dshift

      end if

   end do

   fanum=nnew

   nfa=na

   end subroutine




!=========================================================

   subroutine randforce(nmemb,fxmemb,fymemb,fzmemb, &
               nfa,fxfa,fyfa,fzfa,nmyo,fxmyo,fymyo,fzmyo,nlk,fxlk,fylk,fzlk, &
                jrforce,rxmemb,rymemb,rzmemb, &
                 rxfa,ryfa,rzfa,rxmyo,rymyo,rzmyo,rxlk,rylk,rzlk)

   implicit none

   integer,value::nmemb,nfa,nmyo,nlk
   integer jrforce

   double precision,allocatable,intent(in),dimension(:,:)::rxfa,ryfa,rzfa,rxmyo,rymyo,rzmyo,rxlk,rylk,rzlk
   double precision,allocatable,intent(in),dimension(:,:)::rxmemb,rymemb,rzmemb

   double precision,allocatable,dimension(:)::fxfa,fyfa,fzfa,fxmyo,fymyo,fzmyo,fxlk,fylk,fzlk
   double precision,allocatable,dimension(:)::fxmemb,fymemb,fzmemb

!   double precision,allocatable,dimension(:)::rh1,rh2,rx,ry,rz



   jrforce=jrforce+1

   fxmemb(1:nmemb)=rxmemb(1:nmemb,jrforce)
   fymemb(1:nmemb)=rymemb(1:nmemb,jrforce)
   fzmemb(1:nmemb)=rzmemb(1:nmemb,jrforce)


   fxfa(1:nfa)=rxfa(1:nfa,jrforce)
   fyfa(1:nfa)=ryfa(1:nfa,jrforce)
   fzfa(1:nfa)=rzfa(1:nfa,jrforce)

   fxmyo(1:nmyo)=rxmyo(1:nmyo,jrforce)
   fymyo(1:nmyo)=rymyo(1:nmyo,jrforce)
   fzmyo(1:nmyo)=rzmyo(1:nmyo,jrforce)

   fxlk(1:nlk)=rxlk(1:nlk,jrforce)
   fylk(1:nlk)=rylk(1:nlk,jrforce)
   fzlk(1:nlk)=rzlk(1:nlk,jrforce)

   end subroutine

!=========================================================

   subroutine depoly(nfa,nfa1,nfa2,fanum,fanum1,fanum2,apar,apos,filid,fadist,jfasol,astart,alen,a2mem, &
              xfa,yfa,zfa,myonum,mytyp,fa2myo,crlknum,fa2lk,p_dep,dt,tension1,l_a)


   implicit none

   integer,value::myonum,crlknum,tension1
   integer nfa,nfa1,nfa2,fanum,fanum1,fanum2

   integer nnew,check,jcheck,nfil,newfanum1

   integer n,j,ja,j0,jstop,length,nm,jh,nl,jl,jstart,jcontinue

   integer,allocatable,dimension(:)::astart,alen,a2mem,apos,filid,fadist,map
   integer,allocatable,dimension(:,:)::fa2myo,fa2lk,mytyp,jfasol,apar

   double precision,value::p_dep,dt,l_a
   double precision half,dist0,dist,dx,dy,dz,dcheck,dchange,prob,r

   double precision,allocatable,dimension(:)::xfa,yfa,zfa

   check=0

   half=0.5d0*l_a

!$omp parallel &
!$omp default(none) &
!$omp private(n,jstop,jstart,length,j,ja,r,jh,nm,jl,nl,jcontinue,j0) &
!$omp shared(fanum,alen,astart,tension1,apar,p_dep,dt,fa2myo,mytyp,fa2lk) &
!$omp private(jcheck,dx,dy,dz,dist0,dist,dcheck,dchange,prob) &
!$omp shared(l_a,half,xfa,yfa,zfa) &

!$omp reduction(+:check)

!$omp do schedule(guided,32)

   do n=1,fanum



      length=alen(n)


      jstart=astart(n)

      jstop=jstart+length-1

      if(tension1==1)then

         if(apar(1,jstop)>0.or.apar(1,jstop-1)>0)then
            cycle
         end if


      end if


      jcontinue=0

!     action of cofilin:


      if(length>20)then

         j0=jstart+length/2

         dx=xfa(j0)-xfa(jstart)
         dy=yfa(j0)-yfa(jstart)
         dz=zfa(j0)-zfa(jstart)

         dist0=sqrt(dx*dx+dy*dy+dz*dz)

         dcheck=l_a

         do ja=j0+1,jstop-1

            dx=xfa(ja)-xfa(jstart)
            dy=yfa(ja)-yfa(jstart)
            dz=zfa(ja)-zfa(jstart)

            dist=sqrt(dx*dx+dy*dy+dz*dz)

            dchange=dist-dist0

            dist0=dist

            if(dchange<dcheck.and.apar(1,ja)==0.and.apar(2,ja)==0)then
               dcheck=dchange
               jcheck=ja
            end if

         end do


         if(dcheck<half)then

            prob=1.0d0-dcheck/l_a

            call random_number(r)

            if(prob>r)then

               jcontinue=1

               jstart=jcheck

            end if

         end if

      end if

!     spontaneous depolymerization:

      if(jcontinue==0)then

         call random_number(r)

         if(p_dep*dt>r)then

            jcontinue=1

            if(length>5)then

               jstart=jstop

            end if

         end if

      end if


      if(jcontinue==1)then

         do ja=jstart,jstop

            alen(n)=alen(n)-1

            if(apar(1,ja)>0)then


!              check myosin binding:

               if(apar(2,ja)>0)then

                  jh=apar(1,ja)

                  nm=apar(2,ja)

                  fa2myo(jh,nm)=0

                  mytyp(jh,nm)=1


!              check crosslinker binding:

               else

                  jl=apar(1,ja)

                  nl=-apar(2,ja)

                  fa2lk(jl,nl)=0


               end if

            end if

!            apar(1,ja)=-1
!            apar(2,ja)=0

!            a2mem(ja)=0

            check=check+1


         end do

      end if




   end do

!$omp end do nowait
!$omp end parallel



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

      a2mem(nfil)=a2mem(n)

      do j=1,alen(n)

         ja=jstart+j-1

         nnew=nnew+1

         map(ja)=nnew

         apar(1:2,nnew)=apar(1:2,ja)

         apos(nnew)=apos(ja)


         filid(nnew)=nfil

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


   do n=1,myonum

      do j=1,2

         if(fa2myo(j,n)>0)then
            fa2myo(j,n)=map(fa2myo(j,n))
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

   subroutine myoturnover(nmyoturn,fanum,nnode,nmb2node,nmemb,mybodylen,nxsol,nphisol,a2mem,mytyp,jmysol,memb2node, &
               mybody,myhead,jsursol,jfasol,mydist,l_mb,l_mem,pi,pmyturn,dt,dphisol,xrmin,xrmax, &
                lnode,lnode2,xmyo,ymyo,zmyo,xnode,ynode,znode,xmemb,ymemb,zmemb, &
                 xnorsurf,ynorsurf,znorsurf,xsurf,ysurf,zsurf)

   implicit none

   integer,value::fanum,nnode,nmb2node,nmemb,mybodylen,nxsol,nphisol

   integer nmyoturn

   integer n,jb,jh,nm,nm1,nm2,jx0,jp0,jm,jx,jp,j,jxget,jpget,j1,j2,nmax

   integer,allocatable,dimension(:),intent(in)::a2mem

   integer,allocatable,dimension(:,:)::mytyp,jmysol,memb2node
   integer,allocatable,dimension(:,:),intent(in)::mybody,myhead,jsursol,jfasol
   integer,allocatable,dimension(:)::mydist

   integer,allocatable,dimension(:)::mark

   real(kind=8),value::l_mb,l_mem,pi,pmyturn,dt,dphisol,xrmin,xrmax,lnode,lnode2

   real(kind=8)::dshift,r,xn0,yn0,zn0,x0,y0,z0,phi,rad
   real(kind=8)::dist2,d2,dx,dy,dz,xn,yn,zn,proj,proj0,xmb,ymb,zmb
   real(kind=8)::dx1,dy1,dz1,dx2,dy2,dz2,dx3,dy3,dz3,dx4,dy4,dz4,cosp,sinp
   real(kind=8)::xdir,ydir,zdir,invdist,distance

   real(kind=8),allocatable,dimension(:)::xmyo,ymyo,zmyo,xnode,ynode,znode
   real(kind=8),allocatable,dimension(:),intent(in)::xmemb,ymemb,zmemb
   real(kind=8),allocatable,dimension(:,:),intent(in)::xnorsurf,ynorsurf,znorsurf,xsurf,ysurf,zsurf



   allocate(mark(nmemb))
   mark=0

   do n=1,nnode
      do j=1,2!nmb2node
         if(memb2node(j,n)>0) mark(memb2node(j,n))=1
      end do
   end do

   mark(a2mem(1:fanum))=1

!$omp parallel &
!$omp default(none) &
!$omp private(n,jb,jh,nm,nm1,nm2,jx0,jp0,jm,jx,jp,j,jxget,jpget,j1,j2,nmax) &
!$omp private(dshift,r,xn0,yn0,zn0,x0,y0,z0,phi,rad) &
!$omp private(dist2,d2,dx,dy,dz,xn,yn,zn,proj,proj0,xmb,ymb,zmb) &
!$omp private(dx1,dy1,dz1,dx2,dy2,dz2,dx3,dy3,dz3,dx4,dy4,dz4,cosp,sinp) &
!$omp private(xdir,ydir,zdir,invdist,distance) &
!$omp shared(nnode,nmb2node,nmemb,mybodylen,nxsol,nphisol) &
!$omp shared(mytyp,jmysol,memb2node,mybody,myhead,jsursol,jfasol,mydist,mark) &
!$omp shared(l_mb,l_mem,pi,pmyturn,dt,dphisol,xrmin,xrmax,lnode,lnode2) &
!$omp shared(xmyo,ymyo,zmyo,xnode,ynode,znode,xmemb,ymemb,zmemb) &
!$omp shared(xnorsurf,ynorsurf,znorsurf,xsurf,ysurf,zsurf) &

!$omp reduction(+:nmyoturn)

!$omp do

   do n=1,nnode

      nm1=2*n-1
      nm2=nm1+1

      if(mytyp(1,nm1)/=2.or.mytyp(2,nm1)/=2.or.mytyp(1,nm2)/=2.or.mytyp(2,nm2)/=2) cycle

      call random_number(r)

      if(pmyturn*dt<r) cycle

      nmyoturn=nmyoturn+2

      mytyp(1:2,nm1)=1

      mytyp(1:2,nm2)=1

!     preserve heads configuration:

      jb=mybody(mybodylen,nm1)

      jh=myhead(1,nm1)

      dx1=xmyo(jh)-xmyo(jb)
      dy1=ymyo(jh)-ymyo(jb)
      dz1=zmyo(jh)-zmyo(jb)

      jh=myhead(2,nm1)

      dx2=xmyo(jh)-xmyo(jb)
      dy2=ymyo(jh)-ymyo(jb)
      dz2=zmyo(jh)-zmyo(jb)

      jb=mybody(mybodylen,nm2)

      jh=myhead(1,nm2)

      dx3=xmyo(jh)-xmyo(jb)
      dy3=ymyo(jh)-ymyo(jb)
      dz3=zmyo(jh)-zmyo(jb)

      jh=myhead(2,nm2)

      dx4=xmyo(jh)-xmyo(jb)
      dy4=ymyo(jh)-ymyo(jb)
      dz4=zmyo(jh)-zmyo(jb)

      memb2node(1:nmb2node,n)=0

!     pick a tethering membrane bead:

12    call random_number(r)

      jm=nmemb*r+1

      if(xmemb(jm)<xrmin.or.xmemb(jm)>xrmax) goto 12


!$omp atomic
      mark(jm)=mark(jm)+1

      if(mark(jm)>1) goto 12

      memb2node(1,n)=jm

      jx0=jsursol(1,jm)

      jp0=jsursol(2,jm)

      xn0=xnorsurf(jp0,jx0)
      yn0=ynorsurf(jp0,jx0)
      zn0=znorsurf(jp0,jx0)

      xmb=xmemb(jm)
      ymb=ymemb(jm)
      zmb=zmemb(jm)

!     angular shift of myosin:

      jm=mybody(1,nm1)

      jp=jmysol(2,jm)

      phi=(jp0-jp)*dphisol

      cosp=cos(phi)

      sinp=sin(phi)

!     update distance signal:

      mydist(jm:jm+mybodylen+1)=1

      jm=mybody(1,nm2)

      mydist(jm:jm+mybodylen+1)=1

!     new position of node:

      xnode(n)=xmb
      ynode(n)=ymb
      znode(n)=zmb

      jm=1

      do nm=1,nmemb

         if(jm==2)then
            exit
         end if

         if(mark(nm)==1) cycle

         dx=xmb-xmemb(nm)
         dy=ymb-ymemb(nm)
         dz=zmb-zmemb(nm)

         if(dx*dx+dy*dy+dz*dz<lnode2)then
            jm=jm+1
            memb2node(jm,n)=nm
            mark(nm)=1

            xnode(n)=xnode(n)+xmemb(nm)
            ynode(n)=ynode(n)+ymemb(nm)
            znode(n)=znode(n)+zmemb(nm)

         end if

      end do

      nmax=jm

      xnode(n)=xnode(n)/jm
      ynode(n)=ynode(n)/jm
      znode(n)=znode(n)/jm

!     positions of the first myosin:

!     of the first bead:

      call random_number(r)

      distance=lnode*(1.0d0+r)*0.5d0

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

      do j=1,mybodylen

         jm=mybody(j,nm1)

         xmyo(jm)=x0
         ymyo(jm)=y0
         zmyo(jm)=z0

         jmysol(1,jm)=jx0

         jmysol(2,jm)=jp0

         if(j==mybodylen)then

            jmysol(1,myhead(1:2,nm1))=jx0

            jmysol(2,myhead(1:2,nm1))=jp0

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
            print*,'error in jxget for Myo1',jx0,jp0
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

      jb=mybody(mybodylen,nm1)

      x0=xmyo(jb)
      y0=ymyo(jb)
      z0=zmyo(jb)

      jh=myhead(1,nm1)

      xmyo(jh)=x0+dx1

      ymyo(jh)=y0+dy1*cosp+dz1*sinp

      zmyo(jh)=z0+dz1*cosp-dy1*sinp

      jh=myhead(2,nm1)

      xmyo(jh)=x0+dx2

      ymyo(jh)=y0+dy2*cosp+dz2*sinp

      zmyo(jh)=z0+dz2*cosp-dy2*sinp


!     positions of the second myosin:

      call random_number(r)

      jm=nmax*r+1

      jm=memb2node(jm,n)

      jx0=jsursol(1,jm)

      jp0=jsursol(2,jm)

      xn0=xnorsurf(jp0,jx0)
      yn0=ynorsurf(jp0,jx0)
      zn0=znorsurf(jp0,jx0)

      xmb=xmemb(jm)
      ymb=ymemb(jm)
      zmb=zmemb(jm)

!     position of the first bead:

      call random_number(r)

      distance=lnode*(1.0d0+r)*0.5d0

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

      do j=1,mybodylen

         jm=mybody(j,nm2)

         xmyo(jm)=x0
         ymyo(jm)=y0
         zmyo(jm)=z0

         jmysol(1,jm)=jx0

         jmysol(2,jm)=jp0

         if(j==mybodylen)then

            jmysol(1,myhead(1:2,nm2))=jx0

            jmysol(2,myhead(1:2,nm2))=jp0

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

      jb=mybody(mybodylen,nm2)

      x0=xmyo(jb)
      y0=ymyo(jb)
      z0=zmyo(jb)

      jh=myhead(1,nm2)

      xmyo(jh)=x0+dx3

      ymyo(jh)=y0+dy3*cosp+dz3*sinp

      zmyo(jh)=z0+dz3*cosp-dy3*sinp

      jh=myhead(2,nm2)

      xmyo(jh)=x0+dx4

      ymyo(jh)=y0+dy4*cosp+dz4*sinp

      zmyo(jh)=z0+dz4*cosp-dy4*sinp

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



   subroutine rmxlinker(crlknum,crlknummax,crlknumactive,fa2lk,lktyp,plk_remove,dt)

   implicit none
   integer,value::crlknum,crlknummax
   integer crlknumactive
   integer nl
   integer,allocatable,dimension(:,:),intent(in)::fa2lk
   integer,allocatable,dimension(:)::lktyp

   real(kind=8),value::plk_remove,dt
   real(kind=8)::r,prob

   if(crlknumactive<=crlknummax) return

   prob=plk_remove*dt

!$omp parallel &
!$omp default(none) &
!$omp private(nl,r) &
!$omp shared(crlknum,fa2lk,prob,lktyp) &
!$omp reduction(+:crlknumactive)
!$omp do

   do nl=1,crlknum

      if(fa2lk(1,nl)==0.and.fa2lk(2,nl)==0)then

         call random_number(r)

         if(prob>r)then

            lktyp(nl)=0

            crlknumactive=crlknumactive-1

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
                  end if
 
                  exit
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

   subroutine allpairs(nmemb,xmemb,ymemb,zmemb,neinum,mynei,lknei, &
               nfa,xfa,yfa,zfa,nmyo,xmyo,ymyo,zmyo,nlk,xlk,ylk,zlk,fanum,astart,alen, &
                myonum,myhead,mybody,mybodylen,crlknum,lkstart,lklen, &
                 npair_myac,pair_myac,npair_lkac,pair_lkac,npair_mylk,pair_mylk, &
                  npair_ac2,pair_ac2,npair_my2,pair_my2,npair_lk2,pair_lk2, &
                   npair_mb2,pair_mb2,r_off,cos_t_2)

   implicit none

   integer,value::nmemb,nfa,nmyo,nlk,fanum,myonum,crlknum,lklen,mybodylen,neinum
   integer npair_myac,npair_lkac,npair_mylk,npair_ac2,npair_my2,npair_lk2
   integer npair_mb2
   integer ja,ja1,ja2,jm,jm1,jm2,jl,jl1,jl2,nf1,nf2,jf1,jf2,nm1,nm2,jb1,jb2,jnei
   integer jh1,jh2,nl1,nl2,jc1,jc2,jmb,jo,ji,j,n,jstart1,jstart2

   integer,allocatable,intent(in),dimension(:)::alen,astart,lkstart
   integer,allocatable,intent(in),dimension(:,:)::mybody,myhead
   integer,allocatable,dimension(:,:)::pair_myac,pair_lkac,pair_mylk,pair_ac2,pair_my2,pair_lk2
   integer,allocatable,dimension(:,:)::pair_mb2
   integer,allocatable,dimension(:)::pairtyp
   integer,allocatable,dimension(:,:,:)::mynei,lknei

   double precision,value::r_off,cos_t_2
   double precision dx,dx2,rad1_2,rad2_2,dy,dz,d2,d2max,dmax,dist,cos_2

   double precision,allocatable,intent(in),dimension(:)::xfa,yfa,zfa,xmyo,ymyo,zmyo,xlk,ylk,zlk
   double precision,allocatable,intent(in),dimension(:)::xmemb,ymemb,zmemb
   double precision,allocatable,dimension(:)::pairlen
   integer,allocatable,dimension(:)::pairnum
   integer,allocatable,dimension(:,:)::pairtem
   integer nmax,nnei




!   d2max=2*r_off*r_off

   d2max=2*r_off

   nmax=max(nfa,nmyo,nlk,nmemb)

   nnei=500

   allocate(pairtem(nnei,nmax),pairnum(nmax))


!  actin-myosin pairs:

   mynei=0

!$omp parallel &
!$omp default(none) &
!$omp private(ja,n,jm,nm1,jb1,jh1,jnei,dx,dy,dz,d2,jl) &
!$omp shared(nfa,myonum,mybodylen,mybody,myhead,xfa,yfa,zfa,d2max,pairtem,pairnum) &
!$omp shared(xmyo,ymyo,zmyo,neinum,mynei,nnei)
!$omp do schedule(guided,64)



   do nm1=1,myonum

      do jb1=1,mybodylen

         jm=mybody(jb1,nm1)

         n=0

         do ja=1,nfa

            dx=xfa(ja)-xmyo(jm)
            dy=yfa(ja)-ymyo(jm)
            dz=zfa(ja)-zmyo(jm)

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

         jm=myhead(jh1,nm1)

         n=0

         jnei=0

         do ja=1,nfa

            dx=xfa(ja)-xmyo(jm)
            dy=yfa(ja)-ymyo(jm)
            dz=zfa(ja)-zmyo(jm)

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
                  mynei(jnei,jh1,nm1)=ja
               end if

            end if

         end do

         pairnum(jm)=n

      end do

   end do

!$omp end do nowait
!$omp end parallel

   npair_myac=0

   do jm=1,nmyo

      if(pairnum(jm)==0)then
         cycle
      end if

      do j=1,pairnum(jm)

         npair_myac=npair_myac+1

         pair_myac(1,npair_myac)=jm
         pair_myac(2,npair_myac)=pairtem(j,jm)

      end do

   end do

!print*,'actin-myosin',npair_myac

!-------------------------------------

!     myosin-crosslinker pairs:



!$omp parallel &
!$omp default(none) &
!$omp private(n,jm,dx,dy,dz,d2,jl) &
!$omp shared(nmyo,d2max,pairtem,nlk,pairnum) &
!$omp shared(xmyo,ymyo,zmyo,xlk,ylk,zlk,nnei)
!$omp do schedule(guided,64)

   do jm=1,nmyo



      n=0

      do jl=1,nlk

         dx=xmyo(jm)-xlk(jl)
         dy=ymyo(jm)-ylk(jl)
         dz=zmyo(jm)-zlk(jl)

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


   npair_mylk=0


   do jm=1,nmyo

      if(pairnum(jm)==0)then
         cycle
      end if

      do j=1,pairnum(jm)!nnei


         npair_mylk=npair_mylk+1

         pair_mylk(1,npair_mylk)=jm
         pair_mylk(2,npair_mylk)=pairtem(j,jm)

      end do

   end do


!print*,'myosin-crlk',npair_mylk

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
!  myosin-myosin pairs:



!$omp parallel &
!$omp default(none) &
!$omp private(nm1,jb1,jm1,n,nm2,jb2,jm2,dx,dy,dz,d2,jh1,jh2) &
!$omp shared(myonum,d2max,pairtem,mybody,myhead,mybodylen,pairnum) &
!$omp shared(xmyo,ymyo,zmyo,nnei)
!$omp do schedule(guided,32)

   do nm1=1,myonum-1

      do jb1=1,mybodylen

         jm1=mybody(jb1,nm1)

         n=0


         do nm2=nm1+1,myonum

!           body-body:

            do jb2=1,mybodylen

               jm2=mybody(jb2,nm2)

               dx=xmyo(jm1)-xmyo(jm2)
               dy=ymyo(jm1)-ymyo(jm2)
               dz=zmyo(jm1)-zmyo(jm2)

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

            do jh2=1,2

               jm2=myhead(jh2,nm2)

               dx=xmyo(jm1)-xmyo(jm2)
               dy=ymyo(jm1)-ymyo(jm2)
               dz=zmyo(jm1)-zmyo(jm2)

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

      do jh1=1,2

         jm1=myhead(jh1,nm1)


         n=0

         do nm2=nm1+1,myonum

!           head-body:

            do jb2=1,mybodylen

               jm2=mybody(jb2,nm2)

               dx=xmyo(jm1)-xmyo(jm2)
               dy=ymyo(jm1)-ymyo(jm2)
               dz=zmyo(jm1)-zmyo(jm2)

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

            do jh2=1,2

               jm2=myhead(jh2,nm2)

               dx=xmyo(jm1)-xmyo(jm2)
               dy=ymyo(jm1)-ymyo(jm2)
               dz=zmyo(jm1)-zmyo(jm2)

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


   npair_my2=0


   do jm=1,nmyo-2-mybodylen ! the last myosin beads don't have pairs

      if(pairnum(jm)==0)then
         cycle
      end if


      do j=1,pairnum(jm)!nnei

         npair_my2=npair_my2+1

         pair_my2(1,npair_my2)=jm

         pair_my2(2,npair_my2)=pairtem(j,jm)


      end do

   end do


!print*,'myo-myo',npair_my2

!--------------------------------------------------

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

   subroutine setpair(nmemb,xmemb,ymemb,zmemb,xfa,yfa,zfa,xmyo,ymyo,zmyo,xlk,ylk,zlk, &
               npair5_myac,pair5_myac,npair5_lkac,pair5_lkac,npair5_mylk,pair5_mylk, &
                npair_myac,pair_myac,npair_lkac,pair_lkac,npair_mylk,pair_mylk, &
                 npair5_ac2,pair5_ac2,npair5_my2,pair5_my2,npair5_lk2,pair5_lk2, &
                  npair_ac2,pair_ac2,npair_my2,pair_my2,npair_lk2,pair_lk2, &
                   npair5_mb2,pair5_mb2,npair_mb2,pair_mb2,r_off,l_pair,thet2by2, &
                    pairpart,boundtyp,l_mem,xboundmin,xboundmax,shift, &
                     nnode,xnode,ynode,znode,npairnode_mm,pairnode_mm)



   implicit none

   integer,value::nmemb,npair5_myac,npair5_lkac,npair5_mylk,npair5_ac2,npair5_my2,npair5_lk2,npair5_mb2

   integer,value::nnode

   integer npair_myac,npair_lkac,npair_mylk,npair_ac2,npair_my2,npair_lk2,npair_mb2

   integer npairnode_mm


   integer ja,ja1,ja2,jm,jm1,jm2,jl,jl1,jl2
   integer jmb,jo,n
   integer np1,np2,jb,jc,jd,nnew,jcycle
   integer nmax

   integer n1,n2,n3,n4,jexit,j1,j2,j3,j4,m1,m2,jp,np,nm1,nm2

   integer,allocatable,intent(in),dimension(:,:)::pair5_myac,pair5_lkac,pair5_mylk,pair5_ac2
   integer,allocatable,intent(in),dimension(:,:)::pair5_mb2,pair5_my2,pair5_lk2
   integer,allocatable,dimension(:,:)::pair_myac,pair_lkac,pair_mylk,pair_ac2,pair_my2,pair_lk2
   integer,allocatable,dimension(:,:)::pair_mb2,pairpart,pairnode_mm
   integer,allocatable,dimension(:)::boundtyp
   integer,allocatable,dimension(:)::pairtyp,mark,partlen

   integer,allocatable,dimension(:,:)::partner


   double precision,value::r_off,l_pair,thet2by2,l_mem,xboundmin,xboundmax,shift
   double precision dx,dx2,arg,rad1,rad2,dy,dz,d2,d2max,dist,xc,yc,zc,xsum,ysum,zsum,xmin,xmax
   double precision xl1,yl1,zl1,xl2,yl2,zl2,xp1,yp1,zp1,xp2,yp2,zp2,xp3,yp3,zp3

   double precision,allocatable,intent(in),dimension(:)::xfa,yfa,zfa,xmyo,ymyo,zmyo,xlk,ylk,zlk
   double precision,allocatable,intent(in),dimension(:)::xmemb,ymemb,zmemb,xnode,ynode,znode
   double precision,allocatable,dimension(:)::pairlen,lentem


   DOUBLE PRECISION EX1,EY1,EZ1,EX2,EY2,EZ2,TX,TY,TZ,PX,PY,PZ,QX,QY,QZ
   DOUBLE PRECISION DET,INVDET,T,U,V,delta,low,up,up1
!--------------------------

!   integer zero,nbond,ires,jrepeat,ncheck
!   integer j1,j2,j3,j4,j5,j6,j7,j8,nline,i

!   integer,allocatable,dimension(:,:)::bond,pairtem1,pairtem2


!   double precision charg,mass,w1,w2


!   character(8) tex,typ,res,segid,resno


   d2max=2*r_off*r_off


   nmax=max(npair5_myac,npair5_lkac,npair5_mylk,npair5_ac2,npair5_my2,npair5_lk2,npair5_mb2)



   allocate(mark(nmax))

   mark(1:npair5_myac)=0

!$omp parallel &
!$omp default(none) &
!$omp private(n,ja,jm,dx,dy,dz,d2) &
!$omp shared(npair5_myac,pair5_myac,xfa,yfa,zfa,xmyo,ymyo,zmyo,d2max,mark)
!$omp do

   do n=1,npair5_myac

      jm=pair5_myac(1,n)

      ja=pair5_myac(2,n)


      dx=xfa(ja)-xmyo(jm)
      dy=yfa(ja)-ymyo(jm)
      dz=zfa(ja)-zmyo(jm)

      d2=dx*dx+dy*dy+dz*dz

      if(d2<d2max)then
         mark(n)=1
      end if

   end do

!$omp end do nowait
!$omp end parallel


   npair_myac=0


   do n=1,npair5_myac

      if(mark(n)==1)then

         npair_myac=npair_myac+1

         pair_myac(1:2,npair_myac)=pair5_myac(1:2,n)

      end if

   end do

!-------------------------------------

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

   mark(1:npair5_mylk)=0


!$omp parallel &
!$omp default(none) &
!$omp private(n,jm,jl,dx,dy,dz,d2) &
!$omp shared(npair5_mylk,pair5_mylk,xmyo,ymyo,zmyo,xlk,ylk,zlk,d2max,mark)
!$omp do

   do n=1,npair5_mylk

      jm=pair5_mylk(1,n)

      jl=pair5_mylk(2,n)


      dx=xlk(jl)-xmyo(jm)
      dy=ylk(jl)-ymyo(jm)
      dz=zlk(jl)-zmyo(jm)

      d2=dx*dx+dy*dy+dz*dz

      if(d2<d2max)then
         mark(n)=1
      end if

   end do


!$omp end do nowait
!$omp end parallel

   npair_mylk=0


   do n=1,npair5_mylk

      if(mark(n)==1)then

         npair_mylk=npair_mylk+1

         pair_mylk(1:2,npair_mylk)=pair5_mylk(1:2,n)

      end if

   end do


!--------------------------------------


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


   mark(1:npair5_my2)=0

!$omp parallel &
!$omp default(none) &
!$omp private(n,jm1,jm2,dx,dy,dz,d2) &
!$omp shared(npair5_my2,pair5_my2,xmyo,ymyo,zmyo,d2max,mark)
!$omp do

   do n=1,npair5_my2

      jm1=pair5_my2(1,n)

      jm2=pair5_my2(2,n)

      dx=xmyo(jm1)-xmyo(jm2)
      dy=ymyo(jm1)-ymyo(jm2)
      dz=zmyo(jm1)-zmyo(jm2)

      d2=dx*dx+dy*dy+dz*dz

      if(d2<d2max)then
         mark(n)=1
      end if

   end do

!$omp end do nowait
!$omp end parallel

   npair_my2=0


   do n=1,npair5_my2

      if(mark(n)==1)then

         npair_my2=npair_my2+1

         pair_my2(1:2,npair_my2)=pair5_my2(1:2,n)

      end if

   end do


!--------------------------------------


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



   deallocate(pairtyp,pairlen)


!-------------------------------

!  list of tetrahedrons

!  partners of the pairs in tetrahedrons:

   allocate(partner(20,nmemb),partlen(nmemb))

   partlen=0

   do n=1,npair_mb2

      n1=pair_mb2(1,n)

      n2=pair_mb2(2,n)

      partlen(n1)=partlen(n1)+1

      partner(partlen(n1),n1)=n2

      partlen(n2)=partlen(n2)+1

      partner(partlen(n2),n2)=n1

   end do

   pairpart=0

   do n=1,npair_mb2

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


!-----------------------------------------------------------------


!  to prevent tethers from slipping past each other:

   d2=400.0d0

   npairnode_mm=0

   do n1=1,nnode-1

      do n2=n1+1,nnode

         dx=xnode(n1)-xnode(n2)
         dy=ynode(n1)-ynode(n2)
         dz=znode(n1)-znode(n2)

         if(dx*dx+dy*dy+dz*dz>d2) cycle

         do nm1=2*n1-1,2*n1

            do nm2=2*n2-1,2*n2

               npairnode_mm=npairnode_mm+1

               pairnode_mm(1,npairnode_mm)=nm1

               pairnode_mm(2,npairnode_mm)=nm2

            end do

         end do

      end do

   end do



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

   subroutine solidupdate(nxsol,nphisol,nmemb,nfa,nmyo,nlk,jmbsol,jfasol,jmysol,jlksol, &
               pi,delta,dxsol,dphisol,xmemb,ymemb,zmemb,xfa,yfa,zfa,xmyo,ymyo,zmyo, &
                xlk,ylk,zlk,xwall,ywall,zwall,xnorwall,ynorwall,znorwall, &
                 xsurf,ysurf,zsurf,xnorsurf,ynorsurf,znorsurf,jsursol)



   implicit none

   integer,value:: nxsol,nphisol,nmemb,nfa,nmyo,nlk
   integer n,jx,jp,jx0,jp0,j1,j2,jxget,jpget

   integer,allocatable,dimension(:,:)::jmbsol,jfasol,jmysol,jlksol,jsursol

   double precision,value:: pi,delta,dxsol,dphisol
   double precision xmin,dxsolby2,piby2,dphisolby2,twopi,phi,arg
   double precision dx,dy,dz,d2,dist2

   double precision,allocatable,intent(in),dimension(:)::xmemb,ymemb,zmemb,xfa,yfa,zfa
   double precision,allocatable,intent(in),dimension(:)::xmyo,ymyo,zmyo,xlk,ylk,zlk
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
!$omp shared(nmyo,xmyo,jmysol,ymyo,zmyo) &

!$omp shared(nlk,xlk,jlksol,ylk,zlk)


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

!$omp do schedule(guided,64)

   do n=1,nmyo

      jx0=jmysol(1,n)

      jp0=jmysol(2,n)

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

            dx=xmyo(n)-xsurf(jp,jx)
            dy=ymyo(n)-ysurf(jp,jx)
            dz=zmyo(n)-zsurf(jp,jx)

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

      jmysol(1,n)=jxget

      jmysol(2,n)=jpget

   end do

!$omp enddo nowait

!----------------------------

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


   subroutine solids(nxsol,nphisol,dxsol,dphisol,nmemb,xmemb,ymemb,zmemb,jmbsol, &
                     nfa,xfa,yfa,zfa,jfasol,nmyo,xmyo,ymyo,zmyo,jmysol, &
                     nlk,xlk,ylk,zlk,jlksol,pi,delta)

   implicit none

   integer,value:: nxsol,nphisol,nmemb,nfa,nmyo,nlk
   integer n,jx,jp

   integer,allocatable,dimension(:,:)::jmbsol,jfasol,jmysol,jlksol

   double precision,value:: pi,delta,dxsol,dphisol
   double precision xmin,dxsolby2,piby2,dphisolby2,twopi,phi,arg

   double precision,allocatable,intent(in),dimension(:)::xmemb,ymemb,zmemb,xfa,yfa,zfa
   double precision,allocatable,intent(in),dimension(:)::xmyo,ymyo,zmyo,xlk,ylk,zlk


   xmin=-(nxsol-1)/2*dxsol

   dxsolby2=0.5d0*dxsol

   piby2=0.5d0*pi

   dphisolby2=0.5d0*dphisol

   twopi=2*pi

!$omp parallel &
!$omp default(none) &
!$omp private(n,jx,phi,arg,jp) &
!$omp shared(nmemb,xmemb,xmin,dxsol,dxsolby2,nxsol,jmbsol,ymemb,delta,zmemb,piby2) &
!$omp shared(dphisolby2,nphisol,pi,twopi,dphisol)
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

      jmbsol(1,n)=jx

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

      jmbsol(2,n)=jp

   end do

!$omp end do nowait
!$omp end parallel


!open(1,file='check2.inp')

!do n=1,nmemb
!write(1,*)n,jmbsol(1:2,n)
!end do

!close(1)
!stop


!$omp parallel &
!$omp default(none) &
!$omp private(n,jx,phi,arg,jp) &
!$omp shared(nfa,xfa,xmin,dxsol,dxsolby2,nxsol,jfasol,yfa,delta,zfa,piby2) &
!$omp shared(dphisolby2,nphisol,pi,twopi,dphisol)

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

      jfasol(1,n)=jx

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

      jfasol(2,n)=jp

   end do

!$omp end do nowait
!$omp end parallel


!$omp parallel &
!$omp default(none) &
!$omp private(n,jx,phi,arg,jp) &
!$omp shared(nmyo,xmyo,xmin,dxsol,dxsolby2,nxsol,jmysol,ymyo,delta,zmyo,piby2) &
!$omp shared(dphisolby2,nphisol,pi,twopi,dphisol)
!$omp do schedule(guided,64)

   do n=1,nmyo

      jx=1+(xmyo(n)-xmin)/dxsol

      if(xmyo(n)-xmin-(jx-1)*dxsol>dxsolby2)then
         jx=jx+1
      end if

      if(jx<1)then
         jx=1
      end if

      if(jx>nxsol)then
         jx=nxsol
      end if

      jmysol(1,n)=jx

      if(abs(ymyo(n))<delta)then
         if(zmyo(n)>0.0d0)then
            phi=piby2
         else
            phi=-piby2
         end if

      else

         arg=zmyo(n)/ymyo(n)

         phi=atan(arg)

         if(ymyo(n)<0.0d0)then
            phi=phi+pi
         end if

         if(arg<0.0d0.and.ymyo(n)>0.0d0)then
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

      jmysol(2,n)=jp

   end do

!$omp end do nowait
!$omp end parallel


!$omp parallel &
!$omp default(none) &
!$omp private(n,jx,phi,arg,jp) &
!$omp shared(nlk,xlk,xmin,dxsol,dxsolby2,nxsol,jlksol,ylk,delta,zlk,piby2) &
!$omp shared(dphisolby2,nphisol,pi,twopi,dphisol)
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

      jlksol(1,n)=jx

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

      jlksol(2,n)=jp

   end do

!$omp end do nowait
!$omp end parallel

   end subroutine

!=========================================================

   subroutine solidrads(nthreads,rwall,rmemb,rmembmax,wthick,nxsol,nphisol,nmemb,ymemb,zmemb,jmbsol)

   implicit none

   integer,value::nxsol,nphisol,nmemb,nthreads

   integer jp,jx,n,omp_get_thread_num,tid,nsum

   double precision,value::wthick
   double precision rad,rsum

   integer,allocatable,intent(in),dimension(:,:)::jmbsol
   integer,allocatable,dimension(:,:,:)::ncount

   double precision,allocatable,intent(in),dimension(:)::ymemb,zmemb

   double precision,allocatable,dimension(:,:)::rwall,rmemb,rmembmax
   double precision,allocatable,dimension(:,:,:)::radtemp,radtempmax

   allocate(ncount(nthreads,nphisol,nxsol))
   allocate(radtemp(nthreads,nphisol,nxsol))
   allocate(radtempmax(nthreads,nphisol,nxsol))

   ncount=0


   radtemp=0.0d0
   radtempmax=0.0d0

!$omp parallel &
!$omp default(none) &
!$omp private(n,jx,jp,rad,tid) &
!$omp shared(nmemb,jmbsol,ncount,ymemb,zmemb,radtemp,radtempmax)
!$omp do schedule(guided,64)

   do n=1,nmemb

      jx=jmbsol(1,n)

      jp=jmbsol(2,n)

      tid=omp_get_thread_num()+1
!tid=1
      ncount(tid,jp,jx)=ncount(tid,jp,jx)+1

      rad=sqrt(ymemb(n)*ymemb(n)+zmemb(n)*zmemb(n))

      radtemp(tid,jp,jx)=radtemp(tid,jp,jx)+rad

      if(rad>radtempmax(tid,jp,jx))then
         radtempmax(tid,jp,jx)=rad
      end if


   end do

!$omp enddo nowait
!$omp end parallel

!open(1,file='check1.inp')

!do jx=1,nxsol

!      do jp=1,nphisol

!         write(1,*),jx,jp,sum(ncount(1:nthreads,jp,jx)),sum(radtemp(1:nthreads,jp,jx))

!      end do

!   end do

!close(1)

!stop

   do jx=1,nxsol

      do jp=1,nphisol

         nsum=sum(ncount(1:nthreads,jp,jx))

         if(nsum==0)then

            rmemb(jp,jx)=rwall(jp,jx)-wthick
            rmembmax(jp,jx)=rmemb(jp,jx)

         else

            rsum=sum(radtemp(1:nthreads,jp,jx))

            rmemb(jp,jx)=rsum/nsum

            rmembmax(jp,jx)=maxval(radtempmax(1:nthreads,jp,jx))

         end if

      end do

   end do

   deallocate(ncount,radtemp,radtempmax)

   end subroutine

!=========================================================

   subroutine constraints(nfa,nmemb,nmyo,nlk,jupdate,nxsol,nphisol,nthreads, &
               fadist,mydist,lkdist,jmbsol,jsursol,jfasol,jmysol,jlksol,nsurf, &
                kwall,lsqueez,wthick,k_mem,l_mem,xwall,ywall,zwall,xnorwall, &
                 ynorwall,znorwall,xmemb,ymemb,zmemb,xboundmin,xboundmax,xrmin,xrmax,xfa,yfa,zfa, &
                  xmyo,ymyo,zmyo,xlk,ylk,zlk,xsurf,ysurf,zsurf,xnorsurf,ynorsurf, &
                   znorsurf,fxmemb,fymemb,fzmemb,fxmyo,fymyo,fzmyo,fxlk,fylk,fzlk, &
                    fxfa,fyfa,fzfa)



   implicit none

   integer,value:: nfa,nmemb,nmyo,nlk,jupdate,nxsol,nphisol,nthreads
   integer n,jx,jp,jm,tid,omp_get_thread_num

!   integer,allocatable,intent(in),dimension(:)::a2mem
   integer,allocatable,dimension(:)::fadist,mydist,lkdist
   integer,allocatable,intent(in),dimension(:,:)::jmbsol,jsursol,jfasol,jmysol,jlksol,nsurf

   double precision,value::kwall,lsqueez,wthick,k_mem,l_mem,xboundmin,xboundmax,xrmin,xrmax

   double precision dx,dy,dz,dist,f,d2,rad,l_memby2,dfx,dfy,dfz
   double precision xn,yn,zn

   double precision,allocatable,intent(in),dimension(:,:)::xwall,ywall,zwall,xnorwall,ynorwall,znorwall
   double precision,allocatable,intent(in),dimension(:)::xmemb,ymemb,zmemb
   double precision,allocatable,intent(in),dimension(:)::xfa,yfa,zfa
   double precision,allocatable,intent(in),dimension(:)::xmyo,ymyo,zmyo,xlk,ylk,zlk
   double precision,allocatable,intent(in),dimension(:,:)::xsurf,ysurf,zsurf,xnorsurf,ynorsurf,znorsurf
   double precision,allocatable,dimension(:)::fxmemb,fymemb,fzmemb,fxmyo,fymyo,fzmyo,fxlk,fylk,fzlk
   double precision,allocatable,dimension(:)::fxfa,fyfa,fzfa

   double precision,allocatable,dimension(:,:,:)::ftem
   double precision,allocatable,dimension(:,:)::fxsurf,fysurf,fzsurf


   l_memby2=0.5d0*l_mem

!$omp parallel &
!$omp default(none) &
!$omp private(n,jx,jp,xn,yn,zn,dx,dy,dz,dist,f) &
!$omp shared(nmemb,jmbsol,xwall,ywall,zwall,xnorwall,ynorwall,znorwall) &
!$omp shared(kwall,lsqueez,wthick) &
!$omp shared(xmemb,ymemb,zmemb,fxmemb,fymemb,fzmemb) &
!$omp shared(xboundmin,xboundmax) 

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

!  interacting between actin and membrane:

   allocate(ftem(nthreads,nphisol,nxsol))
   ftem=0.0d0

!   l2max=ltether*ltether


!$omp parallel &
!$omp default(none) &
!$omp private(n,jx,jp,jm,dx,dy,dz,xn,yn,zn,dist,f,tid,d2,dfx,dfy,dfz) &
!$omp shared(nfa,jfasol,fadist,jupdate,ftem) &
!$omp shared(k_mem,l_mem,l_memby2) &
!$omp shared(xmemb,ymemb,zmemb,fxmemb,fymemb,fzmemb,xboundmin,xboundmax,xrmin,xrmax) &
!$omp shared(xfa,yfa,zfa,fxfa,fyfa,fzfa) &
!$omp shared(xsurf,ysurf,zsurf,xnorsurf,ynorsurf,znorsurf) &

!$omp shared(nmyo,jmysol,mydist) &
!$omp shared(xmyo,ymyo,zmyo,fxmyo,fymyo,fzmyo) &

!$omp shared(nlk,jlksol,lkdist) &
!$omp shared(xlk,ylk,zlk,fxlk,fylk,fzlk)



!$omp do schedule(guided,64)

   do n=1,nfa


      if(fadist(n)==1.or.jupdate==1)then

         jx=jfasol(1,n)

         jp=jfasol(2,n)

         dx=xfa(n)-xsurf(jp,jx)
         dy=yfa(n)-ysurf(jp,jx)
         dz=zfa(n)-zsurf(jp,jx)

         xn=xnorsurf(jp,jx)
         yn=ynorsurf(jp,jx)
         zn=znorsurf(jp,jx)

         dist=dx*xn+dy*yn+dz*zn

         if(dist<l_mem)then

            f=k_mem*(l_mem-dist)

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

!$omp do schedule(guided,64)


   do n=1,nmyo


      jp=jmysol(2,n)

      if(mydist(n)==1.or.jupdate==1)then

         jx=jmysol(1,n)

         dx=xmyo(n)-xsurf(jp,jx)
         dy=ymyo(n)-ysurf(jp,jx)
         dz=zmyo(n)-zsurf(jp,jx)

         xn=xnorsurf(jp,jx)
         yn=ynorsurf(jp,jx)
         zn=znorsurf(jp,jx)


         dist=dx*xn+dy*yn+dz*zn

         if(dist<l_mem)then

            f=k_mem*(l_mem-dist)

            fxmyo(n)=fxmyo(n)+f*xn
            fymyo(n)=fymyo(n)+f*yn
            fzmyo(n)=fzmyo(n)+f*zn

            mydist(n)=1

            tid=omp_get_thread_num()+1

            ftem(tid,jp,jx)=ftem(tid,jp,jx)-f

         else

            if(dist<l_mem+l_memby2)then
               mydist(n)=1
            else
               mydist(n)=0
            end if

         end if

      end if

!     blocking from the x boundaries:


      if(xmyo(n)<xrmin)then
         fxmyo(n)=fxmyo(n)+xrmin-xmyo(n)
      else if(xmyo(n)>xrmax)then
         fxmyo(n)=fxmyo(n)+xrmax-xmyo(n)
      end if

   end do

!$omp enddo nowait



!---------------------------------


!$omp do schedule(guided,32)

   do n=1,nlk

      jp=jlksol(2,n)

      if(lkdist(n)==1.or.jupdate==1)then

         jx=jlksol(1,n)

         dx=xlk(n)-xsurf(jp,jx)
         dy=ylk(n)-ysurf(jp,jx)
         dz=zlk(n)-zsurf(jp,jx)

         xn=xnorsurf(jp,jx)
         yn=ynorsurf(jp,jx)
         zn=znorsurf(jp,jx)


         dist=dx*xn+dy*yn+dz*zn

         if(dist<l_mem)then

            f=k_mem*(l_mem-dist)

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

   subroutine nonbond(npair_myac,pair_myac,npair_lkac,pair_lkac,npair_mylk,pair_mylk, &
               npair_ac2,pair_ac2,npair_my2,pair_my2,npair_lk2,pair_lk2, &
                npair_mb2,pair_mb2,pairpart,boundtyp, &
                 xfa,yfa,zfa,fxfa,fyfa,fzfa, &
                  xmyo,ymyo,zmyo,fxmyo,fymyo,fzmyo,&
                   xlk,ylk,zlk,fxlk,fylk,fzlk, &
                    xmemb,ymemb,zmemb,fxmemb,fymemb,fzmemb, &
                     rvdw,kvdw,r_off,r_on2,r_off2,fvdwmax, &
                      kpair,l_pair,l_mem,kmemb,shift)


   implicit none

   integer,value::npair_myac,npair_lkac,npair_mylk,npair_ac2,npair_my2,npair_lk2
   integer,value::npair_mb2

   integer n,ja,ja1,ja2,jm,jm1,jm2,jl,jl1,jl2,n1,n2,n3,n4

   integer,allocatable,intent(in),dimension(:,:)::pair_myac,pair_lkac,pair_mylk,pair_ac2,pair_my2,pair_lk2
   integer,allocatable,intent(in),dimension(:,:)::pair_mb2,pairpart

   integer,allocatable,intent(in),dimension(:)::boundtyp

   double precision,value::rvdw,kvdw,r_off,r_on2,r_off2,fvdwmax

   double precision,value::kpair,l_pair,l_mem,kmemb,shift

   double precision dx,dy,dz,d2,dist,invdist,f,dfx,dfy,dfz,ratio
   double precision invd2,dx34,dy34,dz34,proj,xn,yn,zn

   double precision,allocatable,intent(in),dimension(:)::xfa,yfa,zfa,xmyo,ymyo,zmyo,xlk,ylk,zlk
   double precision,allocatable,intent(in),dimension(:)::xmemb,ymemb,zmemb
   double precision,allocatable,dimension(:)::fxmemb,fymemb,fzmemb
   double precision,allocatable,dimension(:)::fxfa,fyfa,fzfa,fxmyo,fymyo,fzmyo,fxlk,fylk,fzlk

!   halfbound=0.5d0*shift


!$omp parallel &
!$omp default(none) &

!$omp shared(npair_myac,npair_lkac,npair_mylk,npair_ac2,npair_my2,npair_lk2,npair_mb2) &
!$omp private(n,ja,ja1,ja2,jm,jm1,jm2,jl,jl1,jl2) &
!$omp shared(pair_myac,pair_lkac,pair_mylk,pair_ac2,pair_my2,pair_lk2,pair_mb2) &
!$omp shared(rvdw,kvdw,r_off,r_on2,r_off2,fvdwmax) &
!$omp shared(kpair,l_pair,l_mem) &
!$omp private(dx,dy,dz,d2,dist,invdist,f,dfx,dfy,dfz,ratio) &
!$omp shared(xfa,yfa,zfa,xmyo,ymyo,zmyo,xlk,ylk,zlk) &
!$omp shared(xmemb,ymemb,zmemb,fxmemb,fymemb,fzmemb) &
!$omp shared(fxfa,fyfa,fzfa,fxmyo,fymyo,fzmyo,fxlk,fylk,fzlk) &
!$omp private(n1,n2,n3,n4,invd2,dx34,dy34,dz34,proj,xn,yn,zn) &
!$omp shared(pairpart,kmemb,boundtyp,shift)


!  actin-myosin forces:

!$omp do schedule(guided,64)

   do n=1,npair_myac

      jm=pair_myac(1,n)
      ja=pair_myac(2,n)


      dx=xfa(ja)-xmyo(jm)
      dy=yfa(ja)-ymyo(jm)
      dz=zfa(ja)-zmyo(jm)

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
         fxmyo(jm)=fxmyo(jm)-dfx
!$omp atomic
         fymyo(jm)=fymyo(jm)-dfy
!$omp atomic
         fzmyo(jm)=fzmyo(jm)-dfz

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
         fxmyo(jm)=fxmyo(jm)-dfx
!$omp atomic
         fymyo(jm)=fymyo(jm)-dfy
!$omp atomic
         fzmyo(jm)=fzmyo(jm)-dfz

      end if

   end do

!$omp enddo nowait




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




!  myosin-crosslinker forces:

!$omp do schedule(guided,64)

   do n=1,npair_mylk

      jm=pair_mylk(1,n)
      jl=pair_mylk(2,n)

      dx=xmyo(jm)-xlk(jl)
      dy=ymyo(jm)-ylk(jl)
      dz=zmyo(jm)-zlk(jl)

      d2=dx*dx+dy*dy+dz*dz

      if(d2<r_on2)then

         dist=sqrt(d2)

         f=fvdwmax/dist

         dfx=f*dx!*invdist
         dfy=f*dy!*invdist
         dfz=f*dz!*invdist

!$omp atomic
         fxmyo(jm)=fxmyo(jm)+dfx
!$omp atomic
         fymyo(jm)=fymyo(jm)+dfy
!$omp atomic
         fzmyo(jm)=fzmyo(jm)+dfz

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
         fxmyo(jm)=fxmyo(jm)+dfx
!$omp atomic
         fymyo(jm)=fymyo(jm)+dfy
!$omp atomic
         fzmyo(jm)=fzmyo(jm)+dfz

!$omp atomic
         fxlk(jl)=fxlk(jl)-dfx
!$omp atomic
         fylk(jl)=fylk(jl)-dfy
!$omp atomic
         fzlk(jl)=fzlk(jl)-dfz

      end if

   end do

!$omp enddo nowait




!$omp do schedule(guided,64)

!  actin-actin forces:

   do n=1,npair_ac2

      ja1=pair_ac2(1,n)
      ja2=pair_ac2(2,n)


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




!  myosin-myosin forces:

!$omp do schedule(guided,64)

   do n=1,npair_my2

      jm1=pair_my2(1,n)
      jm2=pair_my2(2,n)

      dx=xmyo(jm1)-xmyo(jm2)
      dy=ymyo(jm1)-ymyo(jm2)
      dz=zmyo(jm1)-zmyo(jm2)

      d2=dx*dx+dy*dy+dz*dz

      if(d2<r_on2)then

         dist=sqrt(d2)

         f=fvdwmax/dist

         dfx=f*dx!*invdist
         dfy=f*dy!*invdist
         dfz=f*dz!*invdist

!$omp atomic
         fxmyo(jm1)=fxmyo(jm1)+dfx
!$omp atomic
         fymyo(jm1)=fymyo(jm1)+dfy
!$omp atomic
         fzmyo(jm1)=fzmyo(jm1)+dfz

!$omp atomic
         fxmyo(jm2)=fxmyo(jm2)-dfx
!$omp atomic
         fymyo(jm2)=fymyo(jm2)-dfy
!$omp atomic
         fzmyo(jm2)=fzmyo(jm2)-dfz

      else if(d2<r_off2)then

         dist=sqrt(d2)

!         invdist=1/dist

         ratio=(r_off-dist)/(dist-rvdw)

         f=kvdw*ratio*ratio/dist!(r_off-dist)*(r_off-dist)/(dist-rvdw)/(dist-rvdw)

         dfx=f*dx!*invdist
         dfy=f*dy!*invdist
         dfz=f*dz!*invdist

!$omp atomic
         fxmyo(jm1)=fxmyo(jm1)+dfx
!$omp atomic
         fymyo(jm1)=fymyo(jm1)+dfy
!$omp atomic
         fzmyo(jm1)=fzmyo(jm1)+dfz

!$omp atomic
         fxmyo(jm2)=fxmyo(jm2)-dfx
!$omp atomic
         fymyo(jm2)=fymyo(jm2)-dfy
!$omp atomic
         fzmyo(jm2)=fzmyo(jm2)-dfz

      end if

   end do

!$omp enddo nowait





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

   subroutine measures(junit2,nfa,crlknumactive,natp,jsignal,nmemb,nfamax1,nfamax2,crlknummax, &
              wallrate,aconc1,aconc2,lkconc,thet2by2,cos_t_2,printtime,runtime, &
              pi,l_mem,xfa,yfa,zfa,ymemb,zmemb)

   integer,value::junit2,nfa,crlknumactive,natp,jsignal,nmemb
   integer nfamax1,nfamax2,crlknummax

   integer ncount(4),n,j
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

   subroutine wallremod(nmemb,xmemb,rmembmax,xboundmin,xboundmax,xgap, &
                        nxsol,nphisol,rwall,ywall,zwall,wthick,gap,dphisol)



   implicit none

   integer,value::nmemb,nxsol,nphisol
   integer jp,jx

   double precision,value::wthick,gap,xgap,dphisol
   double precision xboundmin,xboundmax,xmin,xmax
   double precision, allocatable,intent(in), dimension(:)::xmemb
   double precision, allocatable, dimension(:,:)::rwall,rmembmax,ywall,zwall

   do jx=1,nxsol

      do jp=1,nphisol

         if(rwall(jp,jx)>rmembmax(jp,jx)+wthick+gap)then

            rwall(jp,jx)=rwall(jp,jx)-0.1d0*gap

            ywall(jp,jx)=rwall(jp,jx)*cos(jp*dphisol)
            zwall(jp,jx)=rwall(jp,jx)*sin(jp*dphisol)

         end if

      end do

   end do


!  boundary remodeling:

   xmin=minval(xmemb(1:nmemb))
   xmax=maxval(xmemb(1:nmemb))

   if(xmin>xboundmin+xgap)then
      xboundmin=xboundmin+0.1d0
   end if

   if(xboundmax>xmax+xgap)then
      xboundmax=xboundmax-0.1d0
   end if

   end subroutine

!=========================================================


end module
