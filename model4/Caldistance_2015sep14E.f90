   program distance

   implicit none

   character zero*1,charid1*1,charid2*2,charid3*3
   character (len=64) filecoor,fileconf,chara
   character (len=64) filedist,fileangl

   integer nrun

   integer nstart,nxsol,nphisol,nmemb,fanum,nfa,nwalrad,nwalsurf,nwall,nnode
   integer jread(1:20),nact1,nact2,nact,nfamax1,nfamax2,nfa1,nfa2,fanum1,fanum2

   integer,allocatable,dimension(:)::filid,apos,alen,astart,a2mem,fa1stbound
   integer,allocatable,dimension(:,:)::apar


   double precision dxsol,dphisol,xboundmin,xboundmax

   double precision,allocatable,dimension(:,:)::xwall
   real,allocatable,dimension(:)::xwalsurf,xwalrep
   double precision,allocatable,dimension(:)::xmemb,ymemb,zmemb,xfa,yfa,zfa,xnode
   double precision,allocatable,dimension(:,:)::xsurf,ysurf,zsurf,xnorsurf,ynorsurf,znorsurf

   integer n,jx,jp,jx0,jp0,j1,j2,jxget,jpget

   double precision pi,delta,gam,phisum
   double precision xmin,dxsolby2,piby2,dphisolby2,twopi,phi,arg
   double precision dx,dy,dz,d2,dist2,dist,xnor,ynor,znor

   integer j,jmax,jangmax,nsumangl

   integer,allocatable,dimension(:)::ncount,nangle,endfil


   print*,'which run?'

   read(*,*)nrun

   if(nrun<10)then
      write(charid1,'(i1)')nrun
      filedist='distance_run'//charid1//'.txt'
      fileangl='angle_run'//charid1//'.txt'

   else
      write(charid2,'(i2)')nrun
      filedist='distance_run'//charid2//'.txt'
      fileangl='angle_run'//charid2//'.txt'

   end if


   print*,'pick configuration #:'

   read(*,*)nstart

!--------------------------------------------------------

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

   open(1,file=filecoor,form='unformatted')

   open(51,file=fileconf)



2  read(51,*)chara

   if(chara(1:4)/='CELL')then
      goto 2
   end if

   read(51,*)nwall,nxsol,nphisol!,nwalsurf,nwalrad


   allocate(xwall(nphisol,nxsol))!,xwalsurf(nwalsurf),xwalrep(nwalrad))


   read(1)dxsol,dphisol

   read(1)xwall(1:nphisol,1:nxsol)
   read(1)xwall(1:nphisol,1:nxsol)
   read(1)xwall(1:nphisol,1:nxsol)
   read(1)xwall(1:nphisol,1:nxsol)
   read(1)xwall(1:nphisol,1:nxsol)
   read(1)xwall(1:nphisol,1:nxsol)

!   read(1)xwalsurf(1:nwalsurf)
!   read(1)xwalsurf(1:nwalsurf)
!   read(1)xwalsurf(1:nwalsurf)

!   read(1)xwalrep(1:nwalrad)
!   read(1)xwalrep(1:nwalrad)
!   read(1)xwalrep(1:nwalrad)
!   read(1)xwalrep(1:nwalrad)
!   read(1)xwalrep(1:nwalrad)
!   read(1)xwalrep(1:nwalrad)

12 read(51,*)chara

   if(chara(1:4)/='MEMB')then
      goto 12
   end if

   read(51,*)nmemb


   allocate(xmemb(nmemb),ymemb(nmemb),zmemb(nmemb))

   allocate(xsurf(nphisol,nxsol),ysurf(nphisol,nxsol),zsurf(nphisol,nxsol))
   allocate(xnorsurf(nphisol,nxsol),ynorsurf(nphisol,nxsol),znorsurf(nphisol,nxsol))


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

13 read(51,*)chara

   if(chara(1:4)/='NODE')then
      goto 13
   end if
   read(51,*)nnode

   allocate(xnode(nnode))

   read(1)jread(1)
   read(1)xnode(1:nnode)
   read(1)xnode(1:nnode)
   read(1)xnode(1:nnode)


3  read(51,*)chara

   if(chara(1:4)/='FACT')then
      goto 3
   end if

   read(51,*)nact1,nact2,nfamax1,nfamax2,nfa1,nfa2,fanum1,fanum2,jread(1:2)!,ntether

   nact=nact1+nact2

   nfa=nfa1+nfa2

   nfa=nfa1+nfa2

   fanum=fanum1+fanum2

   allocate(filid(nfa),a2mem(nfa),apar(2,nfa),apos(nfa),fa1stbound(fanum),alen(fanum),astart(fanum))

   read(51,*)filid(1:nfa)
   read(51,*)apar(1:2,1:nfa)
   read(51,*)apos(1:nfa)
   read(51,*)alen(1:fanum)
   read(51,*)astart(1:fanum)

   allocate(xfa(nact),yfa(nact),zfa(nact))

   read(1)jread(1)!nfa
   read(1)xfa(1:nfa)
   read(1)yfa(1:nfa)
   read(1)zfa(1:nfa)

   close(1)

   close(51)

   print*,'nfa',nfa   
!--------------------------------------------------------

   allocate(ncount(1000),nangle(1000),endfil(nfa))

   endfil=0

   endfil(astart(1:fanum))=1


   ncount=0

   nangle=0

   jmax=0

   jangmax=0

   nsumangl=0

   phisum=0.0

   pi=3.141592653589793239D0

   delta=0.000001d0

   xmin=-(nxsol-1)/2*dxsol

   dxsolby2=0.5d0*dxsol

   piby2=0.5d0*pi

   dphisolby2=0.5d0*dphisol

   twopi=2*pi

   gam=1.0d0

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

            if(dist2>d2.and.arg>0.0d0)then

               dist2=d2

               dist=arg

               jxget=jx

               jpget=jp

               xnor=xnorsurf(jp,jx)
               ynor=ynorsurf(jp,jx)
               znor=znorsurf(jp,jx)


            end if

         end do

      end do

      if(jxget==0)then
         print*,'error: could not get indices',n,jx0,jp0
!         stop
         cycle
      end if

      j=dist/gam+1

      ncount(j)=ncount(j)+1

      if(jmax<j) jmax=j


      if(mod(n,2)==1.and.endfil(n)==0)then

         dx=xfa(n)-xfa(n-1)
         dy=yfa(n)-yfa(n-1)
         dz=zfa(n)-zfa(n-1)

         arg=dx*xnor+dy*ynor+dz*znor

         arg=abs(arg)/(sqrt(dx*dx+dy*dy+dz*dz))/(1.0d0+delta)

         phi=90.0-acos(arg)*180.0d0/pi

         phisum=phisum+phi

         j=phi+1

         nangle(j)=nangle(j)+1

         if(jangmax<j) jangmax=j

         nsumangl=nsumangl+1

      end if

   end do

   open(1,file=filedist)

   write(1,*)'distance for configuration #',nstart
   write(*,*)'distance for configuration #',nstart

   print*,'jmax',jmax

   do j=1,jmax

      dx=j*gam

      write(1,*)dx,1.0*ncount(j)/nfa

   end do

   close(1)

   open(1,file=fileangl)

   write(1,*)'angle for configuration #',nstart
   write(*,*)'angle for configuration #',nstart

   print*,'jangmax',jangmax

   print*,'average angle =',phisum/nsumangl,'deg'

   write(1,*)'average angle =',phisum/nsumangl,'deg'

   do j=1,jangmax

      write(1,*)1.0*(j-1),1.0*nangle(j)/nsumangl

   end do

   close(1)



   end 
