module declare

   implicit none

   integer(kind=8)::count_rate,cr,timestart,time0,timerun
   real(kind=8)::rate


   integer(kind=8)::nstep,jstep,jskip

   integer jobid,fanum,falen,myonum,mybodylen,nmemb,nfa,nbondfa,nmyo,nbondmyo,maxfalen
   integer nf,nb,jx,ja,ja1,ja2,jf,n,length,jdir,nm,jm,n1,n2,nstart,jfile,junit,junit1,junit2
   integer j,nl,jl,nltem,fanummax,nwalrad,nwalrep,nwalsurf,nwalsurf_p,nwalsurf_x,nlayer
   integer nrand,jrand,jreset,neinum,neinum5,ntotal,nbinds,nframe,jforce1,jforce2,jforce3
   integer jprint,nrforce,jrforce,crlknum1,crlknum2,crlknum,nlk,lklen,nbondlk,nnode,nabound,ntether,ntethertem
   integer days,hours,mins,secs,nthreads,jf1,jf2,j1,j2
   integer npair5_myac,npair5_lkac,npair5_mylk,npair5_ac2,npair5_my2,npair5_lk2,npair5_mb2
   integer npair_myac,npair_lkac,npair_mylk,npair_ac2,npair_my2,npair_lk2
   integer npair_acmb,npair_mymb,npair_lkmb,npair_mb2,nnei
   integer nhoop,nwall,nw,jp,nmeho,nbondwal,jm0,npick,nphi
   integer nxsol,nphisol,jh,jupdate,tension1,tension2,crlkorient,j_dep,cofilin,jmyoturn,jxldel,jxlturnover
   integer natom,natom_a1,natom_a2,natom_a2m,nact1,nact2,nfamax1,nfamax2,fanummax1,fanummax2
   integer jstart,nfahalf,nfafull,astart0,nfa1,nfa2,fanum1,fanum2,jbreak,nact,nmono
   integer natp,jsignal,crlknummax,crlknumactive,nturnover,nmyoturn

   integer,allocatable,dimension(:)::astart,alen,adir,mydir,filid,a2mem,mbhooplen,waltyp,boundtyp
   integer,allocatable,dimension(:,:)::bondfa,bondmyo,mybody,myhead,mytyp,bondwal,mbhoop
   integer,allocatable,dimension(:,:)::fa2myo,fa2lk,bondlk,apar,wposid
   integer,allocatable,dimension(:,:,:)::mynei,lknei,mynei5,lknei5
   integer,allocatable,dimension(:,:)::pair5_myac,pair5_lkac,pair5_mylk,pair5_ac2
   integer,allocatable,dimension(:,:)::pair5_my2,pair5_lk2,pair5_mb2
   integer,allocatable,dimension(:,:)::pair_myac,pair_lkac,pair_mylk,pair_ac2,pair_my2,pair_lk2
   integer,allocatable,dimension(:,:)::pair_acmb,pair_mymb,pair_lkmb,pair_mb2,pairpart
   integer,allocatable,dimension(:,:)::jfasol,jmysol,jlksol,jmbsol,jsursol,nsurf
   integer,allocatable,dimension(:)::fadist,mydist,lkdist,fa1stbound,apos,lktyp,lkstart,mark
!   integer,allocatable,dimension(:)::fa1stboundtem,alentem,astarttem

   double precision walrad,mbrad,rrad,rwid,rthick,mbrad2,invmbrad,mwid,kwall,lsqueez,wallbead
   double precision l_lk1,l_lk2,l_lk,k_lk,klk_thet,thet_lk,phi_node_size,l_mem,k_mem,k_orient,kblock
   double precision pi,k_a,l_a,invl_a,ka_thet,thet_fa,kmh_thet,kmb_thet,thet_mh1,thet_mh2,thet_mb,l_mh,l_mb
   double precision xh,yh,zh,xgap,prob,xboundmin,xboundmax,xrmin,xrmax,shift
   double precision ltether,ktether,pres,fmemb,fturgor,df,mbradnew,p_tether
   double precision rvdw,kvdw,r_on,r_off,r_on2,r_off2,fvdwmax,krad
   double precision kpair,l_pair,rpair,thetpair,thet2by2,cos_t_2,kmemb!_thet,thet_memb
   double precision rvdw_mb,r_on_mb,r_off_mb,r_on2_mb,r_off2_mb,fvdwmax_mb
   double precision dphi,x0,y0,z0,xdel,r,phi0,phi,thick,rad,dx,dy,dz,dytem,dztem,dist
   double precision dx1,dy1,dz1,dx2,dy2,dz2,dx3,dy3,dz3,dx4,dy4,dz4
   double precision sinp,cosp,rxy,rxz,ryy,ryz,rzy,rzz
   double precision xmin,xmax,wthick,d2,dist2,fmax2,runtime,oldtime,printtime,dt,gap
   double precision dxsol,dphisol,wallrate,ringthick,aconc1,aconc2,lkconc,invmylen2

   double precision kbind,lbind,p1_hydr,p1_hydr0,p2_bind,p3_pi_rele,p4_adp_rele,p5_ubind,p_dep
   double precision delta,invdelta,beta,gam,fmaxfa2,fmaxmyo2,fmaxlk2,plk_bind,fmaxmb2
   double precision plk_ubind1,plk_ubind2,plk_off,k_scale,k_scale_lk,plk_remove,pmyturn,plk_turn

   double precision halftime

   real,allocatable,dimension(:)::xwalsurf,ywalsurf,zwalsurf
   real,allocatable,dimension(:)::xwalrad,ywalrad,zwalrad,xwalrep,ywalrep,zwalrep

   double precision,allocatable,dimension(:)::xfa,yfa,zfa,xmyo,ymyo,zmyo
   double precision,allocatable,dimension(:)::rwhoop,xwhoop
   double precision,allocatable,dimension(:)::xmemb,ymemb,zmemb
   double precision,allocatable,dimension(:,:)::xwall,ywall,zwall,rwall,rmemb,rmembmax
   double precision,allocatable,dimension(:)::xlk,ylk,zlk
   double precision,allocatable,dimension(:)::fxfa,fyfa,fzfa,fxfarep,fyfarep,fzfarep
   double precision,allocatable,dimension(:)::fxanglfa,fyanglfa,fzanglfa
   double precision,allocatable,dimension(:)::fxmyo,fymyo,fzmyo
   double precision,allocatable,dimension(:)::fxmyorep,fymyorep,fzmyorep
   double precision,allocatable,dimension(:)::fxanglmyo,fyanglmyo,fzanglmyo
   double precision,allocatable,dimension(:)::fxlk,fylk,fzlk,fxlkrep,fylkrep,fzlkrep
   double precision,allocatable,dimension(:)::fxangllk,fyangllk,fzangllk
   double precision,allocatable,dimension(:)::rands,phi_node
   double precision,allocatable,dimension(:,:)::rxfa,ryfa,rzfa,rxmyo,rymyo,rzmyo,rxlk,rylk,rzlk
   double precision,allocatable,dimension(:,:)::rxmemb,rymemb,rzmemb
!   double precision,allocatable,dimension(:)::fymemfa,fzmemfa,fymemmyo,fzmemmyo,fymemlk,fzmemlk
   double precision,allocatable,dimension(:)::fxmemb,fymemb,fzmemb
!   double precision,allocatable,dimension(:)::fxoll_vdw,fyoll_vdw,fzoll_vdw
!   double precision,allocatable,dimension(:)::fxill_vdw,fyill_vdw,fzill_vdw
   double precision,allocatable,dimension(:)::fxmemb_vdw,fymemb_vdw,fzmemb_vdw

   double precision,allocatable,dimension(:)::fxfa_vdw,fyfa_vdw,fzfa_vdw
   double precision,allocatable,dimension(:)::fxmyo_vdw,fymyo_vdw,fzmyo_vdw
   double precision,allocatable,dimension(:)::fxlk_vdw,fylk_vdw,fzlk_vdw
   double precision,allocatable,dimension(:)::fxmemb_cst,fymemb_cst,fzmemb_cst
   double precision,allocatable,dimension(:)::fxfa_cst,fyfa_cst,fzfa_cst
   double precision,allocatable,dimension(:)::fxmyo_cst,fymyo_cst,fzmyo_cst
   double precision,allocatable,dimension(:)::fxlk_cst,fylk_cst,fzlk_cst
   double precision,allocatable,dimension(:,:)::rmbsol,rwlsol

   double precision,allocatable,dimension(:,:)::xnorwall,ynorwall,znorwall,xnorsurf,ynorsurf,znorsurf
   double precision,allocatable,dimension(:,:)::xsurf,ysurf,zsurf
   real fa1,fa2,fa3,fa4,fa5

   character(len=91)::a91

   real timecheck,myorate

end module
