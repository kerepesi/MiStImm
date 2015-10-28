/**********************************************************************/
/*MiStImm - Microscopic Stochastic Immun System Simulation Software*/

/*Tamas Szabados with students Tamas Kiss, Kristof Horompoly, Endre Szecsei, Csaba Kerepesi*/

/*Comment: Introducing helper T cells, MHC class II, peptide
universe, interleukins, B-T interactions, B cell maturation, antibodies.*/
/**********************************************************************/

/* Compile in Windows: Dev-C++ 4.9.9.2.*/
/* Compile in Linux (with gcc 4.4.6): gcc immune2_jul18_2013.c -lm -o immune2_jul18_2013 */


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define NSUBINTM 16
#define TARGMAX 1000
#define NTYPESMX 100	/* maximum number of cell types */

     typedef struct black *pntrbl;
     /* B cells */
    struct black {
     short bmatur;      /*level of maturation; 0:native B cell,
						1-2: intermediate states, 3: memory cell, 4: plasma cell */
     short bli,blj;     /* subrectangle indices */
     short blx,bly,blr; /* coordinates of B cell receptor and affinity radius */
     float blt1,blt2;   /* birth and death times */
     short bagd;        /* > 0: B cell-antigen distance dd; -1: free;
                           -2: reproducing - not in use; 
                           -dd-3: attached, but nothing can happen */
     short bpepx,bpepy; /*coordinates of B cell and antibody peptide*/
     long nab1;         /*number of antibodies of this type*/          
     short mhc2x[3];    /* MHCII proteins x coord. */   
     short mhc2y[3];    /* MHCII proteins y coord. */
     short mhc_stress[3];
	 short mhcth[3];	/*type of last MHC-Th interaction*/
     float touch_t[3];	 
     short stress; 		/* 0: stress stae of every MHCII are 0; 
						1: there is at least one MHCII with stress state one */
	 short firstkill;	
     float last_il_act;	 
     pntrbl pbl;
    };

    typedef struct thelper *pntrth;
    /* T helper cells */
    struct thelper {
     short tcrx,tcry;     /*coordinates of T cell receptor*/
     float tht1,tht2;     /* birth and death times */
     short thpd;          /* T helper-peptide distance; -1: free,
                                -2: reproducing - not in use */
     short stress;        /* 0: default, 1: IL1 was received, 
                          2: contacted an excited MHCII */
	 short tmatur;		/* 0: default; 1: survived positive and negative selection in the Thymus*/
     float last_il_act;
     pntrth pth;
    };

    typedef struct white *pntrw;
    /* self antigens */
    struct white {
     short wturnon;    /* 0: not existing, 1: existing */
     short xw,yw;      /* coordinates in shape space */
     long nw,nwkilled; /* no. of antigens, killed antigens */
     short wpepx,wpepy;   /* coordinates of peptide */
     float t0w,tauw0,tauw,thnw,etanw; /* initial time, ave. divison time, par's */
     pntrw pnw;
    };

    typedef struct red *pntrr;
    /* non-self antigens */
    struct red {
     short rturnon;    /* 0: not existing, 1: existing */
     short xr,yr;      /* coordinates in shape space */
     long nr,nrkilled; /* no. of antigens, resp., killed antigens */
     short rpepx,rpepy;   /*coordinates of peptide*/
     float t0r,taur0,taur,thnr,etanr; /* initial time, ave. divison time, par's */
     pntrr pnr;
    };

    typedef struct actual *pntra;
    /* events */
    struct actual {
     short etype;  /*the type of the event*/
     void *content; /*the pointer to the actual object*/
     float lambda;     /*the lambda parameter of the object*/
     pntra pact;      /*pointer to the next actual record*/
    };

/* constants */
    float small = 1.0e-10;
    float large = 1.0e5;
    short perturb = 3;

/* pointers */
    pntrw pw,pw1,pwold,pw0;
    pntrr pr,pr1,prold;
    pntrbl emptlbl;
    pntrth emptlth,thlast;
    pntra actfirst,actnew,awhich,emptlast;

short
    comptype,srchtype,rndtype,writtype,outfreq,screenout,mediumreprod,
    lrestore,limmunst,lbradius,
    xmax,xmax2,delx,delx2,rcell,rblgr,nsubint,blkill,
    nm,nb,r0,r0s,rminnew,rminsprd,drwidth,nwtypes,nrtypes,errcode,
    nma,nmres,idum,
    thrad,pxmax,pxmax2,nth,
    rminth,rmaxth,epepx,epepy,rminb;

long
    nrmax,nwa,nra,nwres,nrres,
    nbdied,nbkilled,nbres1,nbres2,previt,prevlit,nab,nnil,nnil1,
    nthdied,nthkild,th_stress,bl_stress;

float
    t,tnew,tstart,tmax,timmst,tbirth,weakr,t0r1,tres,pmem,pmut,
    tlifeb0,tlifmem0,taum0,taubm0,taub0,tauba0,taubr0,bcell_tau_stress_control0,tkill0,tbdiv0,
    tlifeb,tlifmem,taum,taubm,taub,tauba,taubr,bcell_tau_stress_control,tkill,tbdiv,
    thd,thr,thccb,thwd,thnm,thnmb,thparalize1,thparalize2,thba,thil,
    thwdth,thpth,thccth,thdth,thilth,thnmth,thnmpos,thnmneg,
    etad,etar,etaccb,etawd,etanm,etanmb,etaparalize,etaba,etail,etapselt,etanselt,
    etawdth,etapth,etaccth,etadth,etailth,etanmth,tcell_tcrit_il1,tcell_tau_stress_control,
    crnew,crspread,dummy,
    taubil0,taudil0,tauab0,taubab0,taudab0,tauthm0,tauthymus0,thccthn1,tauth0,tcell_tau_stress_control0,
    etaccthn1,tauthr0,tlifeth0,
    taubil,taudil,tauab,taubab,taudab,tauthm,tauthymus,tauth,tau_sel_bcell0,tau_sel_bcell,
    tauthr,tlifeth,
	kth0,kth1,kth2,kb0,kb1,kb2;

/* arrays */
    pntrbl target[TARGMAX];
    pntrw targetw[NTYPESMX];
    pntrr targetr[NTYPESMX];
    float cdf[TARGMAX];
    pntrbl blk[NSUBINTM][NSUBINTM];
    short target_ind[TARGMAX];
	
    FILE *inda,*outd,*out2;
    char *indas,*outds,*out2s;

/* time handling */
    char *timestr;
    time_t tvari1,tvari2;
    double tdif;
   
   float tau_stress; 
   float dring; 
   float tcrit_il;
   float tcrit_stress;
   float tau_prodil1;/*lambda pmeter of B cell producing IL1*/
   float taudil1;/*death rate of IL1*/
   float tauil1;/*attack rate of IL1*/
   float tauil;/*attack rate of IL2*/
   float negselp; /*negative selection probability in the thymus*/
   float posselp; /*positive selection probability in the thymus*/
   float bcellselp;
   
   float sreprod_crit,tau_prodil10,taudil10,tauil10,tauil0;
   
   
/**********************************************************************/
/*rounding real to integer*/

 long roundspec(float x)
/**********************************************************************/

{
 return(floor(x+0.5));
}/*roundspec*/;


/**********************************************************************/
/*writes error messages*/

 void writeerr(void)
/**********************************************************************/

{
 
 if (errcode == 1) {
	 printf(" Error %4d%s%s\n",errcode,":   ","Inner pointer error");
	 fprintf(outd,"\n");
	 fprintf(outd," Error %4d%s%s\n",errcode,":   ","Inner pointer error");
	 fprintf(outd,"\n");
	 } 
 if (errcode == 2) {
	 printf(" Error %4d%s%s\n",errcode,":   "," Input data error (nsubint > NSUBINTM) ");
	 fprintf(outd,"\n");
	 fprintf(outd," Error %4d%s%s\n",errcode,":   "," Input data error (nsubint > NSUBINTM) ");
	 fprintf(outd,"\n");
	 }	 
 if (errcode == 3) {
	 printf(" Error %4d%s%s\n",errcode,":   "," Input data error (nwtypes or nrtypes > ntypesmax) ");
	 fprintf(outd,"\n");
	 fprintf(outd," Error %4d%s%s\n",errcode,":   "," Input data error (nwtypes or nrtypes > ntypesmax) ");
	 fprintf(outd,"\n");
	 } 
 if (errcode == 5) {
	 printf(" Error %4d%s%s\n",errcode,":   "," Not enough memory (for array target) ");
	 fprintf(outd,"\n");
	 fprintf(outd," Error %4d%s%s\n",errcode,":   "," Not enough memory (for array target) ");
	 fprintf(outd,"\n");
	 } 
 if (errcode == 6) {
	 printf(" Error %4d%s%s\n",errcode,":   "," Too many red cells (nr > nrmax) ");
	 fprintf(outd,"\n");
	 fprintf(outd," Error %4d%s%s\n",errcode,":   "," Too many red cells (nr > nrmax) ");
	 fprintf(outd,"\n");
	 } 
 if (errcode == 7) {
	 printf(" Error %4d%s%s\n",errcode,":   "," Not enough memory (for array actual) ");
	 fprintf(outd,"\n");
	 fprintf(outd," Error %4d%s%s\n",errcode,":   "," Not enough memory (for array actual) ");
	 fprintf(outd,"\n");
	 }
 if (errcode == 8) {
	 printf(" Error %4d%s%s\n",errcode,":   "," Not enough memory (for array black) ");
	 fprintf(outd,"\n");
	 fprintf(outd," Error %4d%s%s\n",errcode,":   "," Not enough memory (for array black) ");
	 fprintf(outd,"\n");
	 }	 
 if (errcode == 10) {
	 printf(" Error %4d%s%s\n",errcode,":   "," Trying to handle a non-existing actual ");
	 fprintf(outd,"\n");
	 fprintf(outd," Error %4d%s%s\n",errcode,":   "," Trying to handle a non-existing actual ");
	 fprintf(outd,"\n");
	 }	 
 if (errcode == 11) {
	 printf(" Error %4d%s%s\n",errcode,":   "," Not enough memory (for array white) ");
	 fprintf(outd,"\n");
	 fprintf(outd," Error %4d%s%s\n",errcode,":   "," Not enough memory (for array white) ");
	 fprintf(outd,"\n");
	 }
	 
 if (errcode == 12) {
	 printf(" Error %4d%s%s\n",errcode,":   "," Not enough memory (for array red) ");
	 fprintf(outd,"\n");
	 fprintf(outd," Error %4d%s%s\n",errcode,":   "," Not enough memory (for array red) ");
	 fprintf(outd,"\n");
	 } 
 if (errcode == 14) {
	 printf(" Error %4d%s%s\n",errcode,":   "," Trying to handle a non-existing black (14)");
	 fprintf(outd,"\n");
	 fprintf(outd," Error %4d%s%s\n",errcode,":   "," Trying to handle a non-existing black (14) ");
	 fprintf(outd,"\n");
	 }	 
 if (errcode == 18) {
	 printf(" Error %4d%s%s\n",errcode,":   "," Trying to handle a non-existing T helper ");
	 fprintf(outd,"\n");
	 fprintf(outd," Error %4d%s%s\n",errcode,":   "," Trying to handle a non-existing T helper ");
	 fprintf(outd,"\n");
	 } 
 if (errcode == 19) {
	 printf(" Error %4d%s%s\n",errcode,":   "," Unsuccessful search ");
	 fprintf(outd,"\n");
	 fprintf(outd," Error %4d%s%s\n",errcode,":   "," Unsuccessful search ");
	 fprintf(outd,"\n");
	 }
if (errcode == 20) {
	 printf(" Error %4d%s%s\n",errcode,":   "," Trying to handle a non-existing black (20) ");			/*Temporary */
	 fprintf(outd,"\n");
	 fprintf(outd," Error %4d%s%s\n",errcode,":   "," Trying to handle a non-existing black (20) ");
	 fprintf(outd,"\n");
	 }
if (errcode == 21) {
	 printf(" Error %4d%s%s\n",errcode,":   "," Trying to handle a non-existing black (21) ");			/*Temporary */
	 fprintf(outd,"\n");
	 fprintf(outd," Error %4d%s%s\n",errcode,":   "," Trying to handle a non-existing black (21) ");
	 fprintf(outd,"\n");
	 }
if (errcode == 22) {
	 printf(" Error %4d%s%s\n",errcode,":   "," Trying to handle a non-existing black (22) ");			/*Temporary */
	 fprintf(outd,"\n");
	 fprintf(outd," Error %4d%s%s\n",errcode,":   "," Trying to handle a non-existing black (22) ");
	 fprintf(outd,"\n");
	 }
if (errcode == 23) {
	 printf(" Error %4d%s%s\n",errcode,":   "," Trying to handle a non-existing black (23) ");			/*Temporary */
	 fprintf(outd,"\n");
	 fprintf(outd," Error %4d%s%s\n",errcode,":   "," Trying to handle a non-existing black (23) ");
	 fprintf(outd,"\n");
	 }
if (errcode == 24) {
	 printf(" Error %4d%s%s\n",errcode,":   "," Trying to handle a non-existing black (24) ");			/*Temporary */
	 fprintf(outd,"\n");
	 fprintf(outd," Error %4d%s%s\n",errcode,":   "," Trying to handle a non-existing black (24) ");
	 fprintf(outd,"\n");
	 }
if (errcode == 101) {
	 printf(" Error %4d%s%s\n",errcode,":   "," Trying to handle a non-existing actual (101) ");			/*Temporary */
	 fprintf(outd,"\n");
	 fprintf(outd," Error %4d%s%s\n",errcode,":   "," Trying to handle a non-existing actual (101) ");
	 fprintf(outd,"\n");
	 }
	 

}/*writeerr*/;


/**********************************************************************/
/*sets the parameters*/

 void setpar(void)
/**********************************************************************/

{
 float coef;
 short i,j;


 out2=fopen(out2s,"w");

 
 actfirst=NULL;
 actnew=NULL;
 awhich=NULL;
 emptlast=NULL;
 emptlbl=NULL;
 emptlth=NULL;

 coef=-1.0;
 tlifeb=coef/tlifeb0;
 tlifmem=coef/tlifmem0;
 taum=coef/taum0;
 taubm=coef/taubm0;
 taub=coef/(taub0);
 tauba=coef/tauba0;
 taubr=coef/taubr0;
 
 taubil=coef/taubil0;
 taudil=coef/taudil0;
 tauab=coef/tauab0;
 taubab=coef/taubab0;
 taudab=coef/taudab0;
 tauthm=coef/tauthm0;
 tauthymus=coef/tauthymus0;
 tau_sel_bcell=coef/tau_sel_bcell0;
 tauth=coef/tauth0;
 tauthr=coef/tauthr0;
 tlifeth=coef/tlifeth0;
 tkill=coef/tkill0;
 tbdiv=coef/tbdiv0;
 
 tau_prodil1=coef/tau_prodil10;
 taudil1=coef/taudil10;
 tauil1=coef/tauil10;
 tauil=coef/tauil0;
 bcell_tau_stress_control=coef/bcell_tau_stress_control0;
 tcell_tau_stress_control=coef/tcell_tau_stress_control0;

 pw=pw1;

 while (pw != NULL )
 {
   pw->tauw=coef/(pw->tauw0);
   pw=pw->pnw;
 }/*while*/;

 pr=pr1;

 while (pr != NULL )
 {
   pr->taur=coef/(pr->taur0);
   pr=pr->pnr;
 }/*while*/;

 /*initiate black and T*/
 /**********************/

 thlast=NULL;

 for (i=1; i<= nsubint; i++) 
 {
  for (j=1; j<= nsubint; j++) 
   blk[i][j]=NULL;
 }/*i*/;

 previt=-1;
 prevlit=-outfreq;


  tstart=0.0;

  nb=0;
  nth=0;
  nab=0;
  nnil=0;
  nnil1=0;
  nbdied=0;
  nbkilled=0;
  nthdied=0;
  nthkild=0;
  th_stress=0;
  bl_stress=0;


}/*setpar*/;


/**********************************************************************/
/*stores the output*/

 void storeout(short first)
/**********************************************************************/

{
    long nt,it;
    short lout;

 if (first)
 {

  fprintf(outd,"\n");
  fprintf(outd,"The numbers of white, red, black, and bone marrow cells %s\n",
               "during the simulation, ");
  fprintf(outd,"died black and T helper cells: coordinates, %s\n",
               "(radius), times of birth and death ");
  fprintf(outd,"\n");

 }/*if*/;

 it=roundspec(t);

 if (it == previt)
  goto soend;

 previt=it;

  if (screenout == 1) {	
	printf("t=%6.0f",t);
  }/*if*/

 lout=(((it % outfreq) == 0) || (it-prevlit >= outfreq));

 if (lout)
 {
  fprintf(outd,"t=%6.0f",t);
  prevlit=it;
 }

 nwa=0;
 pw=pw1;

 while (pw != NULL )
 {

   if (pw->wturnon)
    nt=pw->nw;
   else
    nt=0;

   nwa=nwa+nt;
   pw=pw->pnw;
 
 }/*while*/;

 nra=0;
 pr=pr1;

 while (pr != NULL )
 {

   if (pr->rturnon)
    nt=pr->nr;
   else
    nt=0;

   nra=nra+nt;
   pr=pr->pnr;

 }/*while*/;

 if (t < timmst)
  nma=0;
 else
  nma=nm;
	
 if (screenout == 1) {	
	printf("%s%7ld%s%7ld%s%4d%s%7ld%s%4d%s%7ld%s%4d%s%7ld%s%7ld%s%7ld\n"," nW=",
        nwa," nR=",nra," nB=",nb," nAb=",nab," nTh=",nth," nIL=",nnil," nM=",
        nma," nIL1",nnil1," th_stress",th_stress," bl_stress",bl_stress);
 }
		
 if (lout)
 {

  fprintf(outd,"%s%7ld%s%7ld%s%4d%s%7ld%s%4d%s%7ld%s%4d\n"," nW=",nwa," nR=",
               nra," nB=",nb," nAb=",nab," nTh=",nth," nIL=",nnil," nM=",nma);
 }/*if*/;
 
 soend:;

}/*storeout*/;


/**********************************************************************/
/*includes an entry in the "actual" list*/

 void include(short lnew, void *lbwr, short et, float lam)
/**********************************************************************/

{

 pntra actold,ai;

 if (lnew)
 {
  actold=actnew;

  if (emptlast == NULL)
  {

   actnew=(struct actual *) malloc(sizeof(struct actual));

   if (actnew == NULL)
   {
    errcode=7;
    goto iend;
   }
  }
  else
  {
   actnew=emptlast;
   emptlast=emptlast->pact;
  }/*if*/;

  if (actfirst == NULL)
   actfirst=actnew;

  if (actold != NULL)
   actold->pact=actnew;

   actnew->etype=et;
   actnew->lambda=lam;
   actnew->content=lbwr;
   actnew->pact=NULL;

 }
 else
 {  /*modification only*/

  ai=actfirst;

  while (ai != NULL )
  {

    if ((ai->content == lbwr) && (ai->etype == et))
    {
     ai->lambda=lam;
     goto iend;
    }
    else
     ai=ai->pact;

  }/*while*/;

  errcode=101;

 }/*if*/;

 iend:;

}/*include*/;


/**********************************************************************/
/*removes an entry from the "actual" list*/

 void remov(void *lbwr, short et)
/**********************************************************************/

{

 pntra actold,ai;

 /*searching the object among the actuals*/

 actold=NULL;
 ai=actfirst;

 while (ai != NULL )
 {

   if ((ai->content == lbwr) && (ai->etype == et))
   {
    ai->etype=0;
    ai->lambda=0.0;
    ai->content=NULL;

    if (ai == actfirst)
     actfirst=ai->pact;
    else
     actold->pact=ai->pact;

    ai->pact=emptlast;
    emptlast=ai;

    if (ai == actnew)
     actnew=actold;

    goto rend;

   }/*if*/;

   actold=ai;
   ai=ai->pact;

 }/*while*/;

 errcode=10;

 rend:;

}/*remov*/;


/**********************************************************************/
/*generates a random number uniformly distributed on [0, 1)*/

 float randr(void)
/**********************************************************************/

{

 float rn;

 rn=rand()/(RAND_MAX+1.0);

 return rn;

} /*randr*/;


/**********************************************************************/
/*generates a random number uniformly distributed on the integers 
  {0,...,n-1}*/

 long randi(long n)
/**********************************************************************/

{

 long rn;

 rn=floor(n*randr());

 return rn;

} /*randi*/;


/**********************************************************************/
/*generates an exponential random number*/

 float randexp(float coef)
/**********************************************************************/

{

 float rn;

 rn=log(1.0-randr())/coef;

 if (rn > large)
  rn=large;

 return rn;

} /*randexp*/;


/**********************************************************************/
/*computes the distance*/

 short d(short x1,short y1,short x2,short y2)
/**********************************************************************/

{

 short dd,xd,yd;

 xd=abs(x2-x1);
 yd=abs(y2-y1);
 if (yd >= xd)
  dd=yd;
 else
  dd=xd;

 return dd;

} /*d*/;


/**********************************************************************/
/*generates a uniform random number on an "annulus"*/

 void unif(short xc,short yc,short rr,short way,short xm,short xm2, 
           short *xx,short *yy)
/**********************************************************************/

{

 short x,y,r1,r2;

 if (way == 1)
 {
  r1=0;
  r2=rr;
 }
 else
 {
  r1=rr-drwidth;
  r2=rr+drwidth;
 }

 lback:;

 x=roundspec((2*r2)*randr());
 y=roundspec((2*r2)*randr());
 *xx=xc-r2+x;
 *yy=yc-r2+y;

 if (d(xc,yc,*xx,*yy) < r1)
  goto lback;

 if ((*xx >= xm) || (*xx < 0) || (*yy >= xm2) || (*yy < -xm2))
  goto lback;

} /*unif*/;


/**********************************************************************/
/*computes a gate shaped function*/

 float gate(long cc,float aa,float bb,float et)
/**********************************************************************/

{

 float gat,ce;

 if ((cc < 0) || ((cc == 0) && (fabs(aa) >= small)))
  gat=1.0;
 else
 {
  if (cc == 0)
   gat=1.0;
  else
  {
   ce=exp(et*log(cc));

   if (fabs(aa) < small)
    gat=bb/(bb+ce);
   else
    gat=(ce*bb)/((aa+ce)*(bb+ce));

  }
 }

 return gat;

}/*gate*/;


/**********************************************************************/
/*computes a gate shaped function of a real argument*/

 float gater(float cc,float aa,float bb,float et)
/**********************************************************************/

{

 float gat,ce;

 if ((cc < 0.0) || ((fabs(cc) < small) && (fabs(aa) >= small)))
  gat=1.0;
 else
 {
  if (fabs(cc) < small)
   gat=1.0;
  else
  {
   ce=exp(et*log(cc));

   if (fabs(aa) < small)
    gat=bb/(bb+ce);
   else
    gat=(ce*bb)/((aa+ce)*(bb+ce));
  }
 }

 return gat;

}/*gater*/;


/**********************************************************************/
/*computes the bounds of the affected cells*/

 void bounds(short x,short y,short r,
             short *i1,short *i2,short *j1,short *j2)
/**********************************************************************/

{

 *i1=1+((x-r) / delx);
 *i2=1+((x+r) / delx);
 *j1=1+((xmax2+y-r) / delx);
 *j2=1+((xmax2+y+r) / delx);

 if (*i1 < 1)
  *i1=1;

 if (*i2 > nsubint)
  *i2=nsubint;

 if (*j1 < 1)
  *j1=1;

 if (*j2 > nsubint)
  *j2=nsubint;

}/*bounds*/;


/**********************************************************************/
/*computes the weighted number of objects of a B cell or antibody*/

 float scountob(pntrbl l)
/**********************************************************************/

{

    short xx,yy,rr,i,j,i1,i2,j1,j2,dblack,dwhite,dred;
    long c1;
    float we,th,et;

    pntrbl k;

 xx=l->blx;
 yy=-(l->bly);
 rr=l->blr;


 bounds(xx,yy,rr,&i1,&i2,&j1,&j2);

 we=0.0;

 for ( i=i1; i<= i2; i++)  
 {
  for ( j=j1; j<= j2; j++) 
  {
   k=blk[i][j];

   while (k != NULL )
   {

     dblack=d(k->blx,k->bly,xx,yy);

     if ((dblack <= rr) && (k->blt1 <= t))
     {

       c1=l->nab1;

       if ((l->blt2 > t) && (l->bmatur <= 3) && (k != l))
        c1=c1+1;                                    /*B cell is available*/

       if (c1 > 0)
       {
          we=we+c1*gate(dblack,0.0,thwd,etawd);
       }/*if*/;
     }/*if*/;

     k=k->pbl;

   }/*while*/;
  }/*j*/;
 }/*i*/;

 /*position of the whites*/

 pw=pw1;

 while (pw != NULL )
 {

   dwhite=d(xx,yy,pw->xw,pw->yw);

   if ((dwhite <= rr) && (pw->wturnon))
   {

     we=we+(pw->nw)*gate(dwhite,0.0,thwd,etawd);

   }/*if*/;

   pw=pw->pnw;

 }/*while*/;

 /*position of the reds*/

 pr=pr1;

 while (pr != NULL )
 {

   dred=d(xx,yy,pr->xr,pr->yr);

   if ((dred <= rr) && (pr->rturnon))
   {

     we=we+(pr->nr)*gate(dred,0.0,thwd,etawd);

   }/*if*/;

   pr=pr->pnr;

 }/*while*/;

 return we;

}/*scountob*/;


/**********************************************************************/
/*computes the number of B cells and antibodies in a neighborhood*/

 long countbl(pntrbl l,short i1,short i2,short j1,short j2,
              short xx,short yy,short rr,short how,short *nbty)
/**********************************************************************/

{

    short i,j,dblack;
    long c,c1,nb1;
    float et,th;

    pntrbl k;

 c=0;
 nb1=0;



 if (how <= 1)
  cdf[0]=0.0;

 for ( i=i1; i<= i2; i++) 
 {
  for ( j=j1; j<= j2; j++)  
  {
   k=blk[i][j];

   while (k != NULL )
   {

     dblack=d(k->blx,k->bly,xx,yy);

     if ((dblack <= rr) && (k->blt1 <= t))
     {

      if (how == 3)         /*counting the number of B cells only*/
      {
       if (k->blt2 > t)
        nb1=nb1+1;
      }
      else
      {

       c1=k->nab1;

       if ((k->blt2 > t) && (k->bmatur <= 3) && ((k != l) || (how != 1)))
        c1=c1+1;                                    /*B cell is available*/

       if (c1 > 0)
       {
        c=c+1;
        nb1=nb1+c1;

        if (how <= 1)
        {

         if (c > TARGMAX)
         {
          errcode=5;
          goto cend;
         }

         target[c]=k;

          cdf[c]=cdf[c-1]+c1*gate(dblack,0.0,thwd,etawd);

        }/*if*/;
       }/*if*/;
      }/*if*/;
     }/*if*/;

     k=k->pbl;

   }/*while*/;
  }/*j*/;
 }/*i*/;

 cend:;

 *nbty=c;
 return nb1;

}/*countbl*/;


/**********************************************************************/
/*computes the number of T cells in a neighborhood*/

 long countth(short xx,short yy,short dr)
/**********************************************************************/

{

 short dthelper,c;
 pntrth k;
 
 c=0;
 
  if(thlast==NULL)
  goto rend;              
  k=thlast;
   

  while (k != NULL )
  {
     dthelper=d(k->tcrx,k->tcry,xx,yy);

     if ((dthelper <= dr) && (k->tht1 <= t) && (k->tht2 > t))
     {
		c=c+1;
	 }
     k=k->pth;
  }/*while*/
	
	rend:;
	
	return c;
		 
}/*countth*/

/**********************************************************************/
/*computes the number of objects in a neighborhood*/

 long countobj(pntrbl l,short xx,short yy,short rr,short how, 
               short *nbty,short *nwty,short *nrty, 
               float *weigb,float *weigw,float *weigr)
/**********************************************************************/

{
    short i1,i2,j1,j2,ii,c,dwhite,dred;
    long counto,nb1,nw1,nr1;
    float th,et;

 nw1=0;
 nr1=0;

 bounds(xx,yy,rr,&i1,&i2,&j1,&j2);

 nb1=countbl(l,i1,i2,j1,j2,xx,yy,rr,how,&c);

 if (errcode > 0)
  goto coend;

 *nbty=c;
 *weigb=cdf[c];
 ii=0;

 /*position of the whites*/

 pw=pw1;

 while (pw != NULL )
 {

   dwhite=d(xx,yy,pw->xw,pw->yw);

   if ((dwhite <= rr) && (pw->wturnon))
   {
    nw1=nw1+pw->nw;
    c=c+1;
    ii=ii+1;
    targetw[ii]=pw;

     cdf[c]=cdf[c-1]+(pw->nw)*gate(dwhite,0.0,thwd,etawd);

   }/*if*/;

   pw=pw->pnw;

 }/*while*/;

 *nwty=ii;
 *weigw=cdf[c];
 ii=0;

 /*position of the reds*/

 pr=pr1;

 while (pr != NULL )
 {

   dred=d(xx,yy,pr->xr,pr->yr);

   if ((dred <= rr) && (pr->rturnon))
   {
    nr1=nr1+pr->nr;
    c=c+1;
    ii=ii+1;
    targetr[ii]=pr;

     cdf[c]=cdf[c-1]+(pr->nr)*gate(dred,0.0,thwd,etawd);

   }/*if*/;

   pr=pr->pnr;

 }/*while*/;

 *nrty=ii;
 *weigr=cdf[c];

 coend:;

 counto=nb1+nw1+nr1;

 return counto;

}/*countobj*/;

 
/**********************************************************************/
 
 float scounpep(pntrth l)
/**********************************************************************/

{

    short i,j,xx,yy,dpep;
    float we,th,et;

    pntrbl k;

 we=0.0;

  xx=l->tcrx;
  yy=l->tcry;

 for ( i=1; i<= nsubint; i++) 
 {
  for ( j=1; j<= nsubint; j++) 
  {

   k=blk[i][j];

   while (k != NULL )
   {

     int l;
     
     for(l=0;l<3;l++)
     {
		if ((k->mhc2x[l]==epepx) && (k->mhc2y[l]==epepy))
		{
		goto epepend;
		}
		
        dpep=d(k->mhc2x[l],k->mhc2y[l],xx,-yy);

        if ((dpep <= thrad) && (k->blt1 <= t) && 
        (k->blt2 >= t) && (k->bmatur <= 3))
        {                               /*the peptide is available*/

			we=we+gate(dpep,0.0,thwdth,etawdth);

        }/*if*/;
		epepend:;
		
	}
        
     k=k->pbl;

   }/*while*/;

  }/*j*/;
 }/*i*/;

 return we;

}/*scounpep*/;


/**********************************************************************/
/*the number of MHC class II bound peptides in a neighborhood*/

short countpep(short xx,short yy,float *weig)
/**********************************************************************/

{
    short i,j,dpep,c;
    float th,et;

    pntrbl k;

 c=0;
 cdf[0]=0.0;


 for ( i=1; i<= nsubint; i++) 
 {
  for ( j=1; j<= nsubint; j++) 
  {

   k=blk[i][j];

   while (k != NULL )
   {

     int l;
     
     for(l=0;l<3;l++)
     {
	 if ((k->mhc2x[l]==epepx) && (k->mhc2y[l]==epepy))
		{
		goto epepend2;
		}
		
        dpep=d(k->mhc2x[l],k->mhc2y[l],xx,yy);

        if ((dpep <= thrad) && (k->blt1 <= t) && 
        (k->blt2 >= t) && (k->bmatur <= 3))
        {                               /*the peptide is available*/
            c=c+1;

            if (c > TARGMAX)
            {
               errcode=5;
               goto cend;
            }

            target[c]=k;
            target_ind[c]=l;

               cdf[c]=cdf[c-1]+gate(dpep,0.0,thwdth,etawdth);

        }/*if*/;
		epepend2:;
     }  
     k=k->pbl;

   }/*while*/;

  }/*j*/;
 }/*i*/;

 cend:;

 *weig=cdf[c];
 return c;

}/*countpep*/;


/**********************************************************************/
/*computes the effective radius of the black offsprings*/

 short rho(short d,short rr)
/**********************************************************************/

{
 short rh;
 float b;
 
 switch (lbradius) {

  case 0: case 1: rh=r0; break;

  case 2: rh=roundspec((r0-rminnew)*randr())+rminnew; break;

  case 3: rh=roundspec((rr-rminnew)*randr())+rminnew; break;

  default:
   rh=roundspec(crnew*d)+rminnew; break;

 }/*switch*/;

 return rh;

}/*rho*/;


/**********************************************************************/
/*computes the radius of the cell where the black offsprings are spread*/

 short rspread(short d)
/**********************************************************************/

{

 short rs;

 switch (lbradius) {

  case 0: case 2: case 3: rs=r0s; break;

  default:
   rs=roundspec(crspread*d)+rminsprd; break;

 }/*switch*/;

 return rs;

}/*rspread*/;


/**********************************************************************/
/*modifies the time parameter of free B cells*/

 void tbfreem(void)
/**********************************************************************/

{

/*	float weigb,weigw,weig,thpara2;
    short i,j,nw1,xc,yc,rr,nbty,nwty,nrty;
	long n2;


    pntrbl l;

 for ( i=1; i<= nsubint; i++) 
 {
  for ( j=1; j<= nsubint; j++) 
  {
   l=blk[i][j];
   while (l != NULL )
   {
     if ((l->blt1 <= t) && (l->blt2 >= t) && (l->bmatur <= 3))
     {                                           
      
*/
	  /*the objectives are in the conjugate cells*/
/*
	    xc=l->blx;
		yc=l->bly;
		rr=l->blr;
*/      
	   /*n2=countobj(l,xc,-yc,rr,1,&nbty,&nwty,&nrty,&weigb,&weigw,&weigr);*/
	   
/*	   weig=scountob(l);
	    
		thpara2=exp(etaparalize*log(thparalize2));
		
			include(0,l,24,taub*gate(weig,0.0,thpara2,etaparalize)); */ /* modifying B cell action process*/
	    
			/*fprintf(out2,"%8.2f%20.2f%20.2f\n",t,gate(weig,0.0,thpara,etaparalize));*/
			
 /*    };

     l=l->pbl;

   };

  };
 };

*/
 
}/*tbfreem*/;


/**********************************************************************/
/*modifies the time parameter of free T helper cells*/

 void tthfreem(void)
/**********************************************************************/

{
 /*   short i,j;
    float weig;

	
	
    pntrth l;

 l=thlast;

 while (l != NULL )
 {

     if ((l->tht1 <= t) && (l->tht2 >= t))      
                                                 
     {
*/
      /*the objectives are in the conjugate cells*/

 /*     weig=scounpep(l);
*/
	  
	  /*fprintf(out2,"%f\n","weig=",weig);*/
	  
/*      include(0,l,25,weig*tauth);

     };

     l=l->pth;

 };
*/
 
}/*tthfreem*/;


/**********************************************************************/
/*a black is born*/

 void blborn(short ldelay, pntrbl l)
/**********************************************************************/

{
    short xx,yy,rr,bmat,lam;
    float weigr;
 
   nb=nb+1;
   xx=l->blx;
   yy=l->bly;
   rr=l->blr;
   bmat=l->bmatur;

 if (ldelay != 0)
 {
	l->blt1=t;
	remov(l,0);	 /*removing the birth from the actual list*/
 }/*if*/;

 if (bmat==0) 
 {
	 include(1,l,28,tau_sel_bcell);		 /*bone_marrow_filter*/
 }
 
 if (bmat <= 3)
 {
	if(nnil>0)							
	{
		lam=nnil*nb*tauil;
		include(0,NULL,22,lam);				/*modifying il2 actions (il2_act)*/
	}
	/*weigr=scountob(l);*/				  /*the objectives are in the conjugate cells*/
	if (bmat >= 1) {
		include(1,l,24,taub);           /*normal actions (attack_bl) */
	}
	/*include(1,l,23,bcell_tau_stress_control);    moved to the first succesful action */
 }
 else
	if (bmat == 4)                   /*plasma cell*/
	{
	   include(1,l,9,taubab);		/*plasma cell action (plasmab)*/
	}
	else
	   errcode= 14;

 if (errcode > 0)
  goto rend;

 if (bmat == 3)          		 /*memory cell is born*/
  include(1,l,2,tlifmem);		 /*bldies*/
 else
  include(1,l,2,tlifeb);		 /*bldies*/

 if (errcode > 0)
  goto rend;

 /*tbfreem();*/			 /*modifying the time parameters of free B and Th cells*/
 /*tthfreem(); */

 rend:;

}/*blborn*/;


/**********************************************************************/
/*creates or prepares the birth of a new black*/

 void newblack(short inew,short jnew,short xnew,short ynew,short rnew,
               short bmnew,short pxnew[3],short pynew[3],short ldelay,float rt)
/**********************************************************************/

{
 pntrbl lold,lnew;
 int k;

 lold=blk[inew][jnew];

 if (emptlbl == NULL)
 {
  lnew=(struct black *) malloc(sizeof(struct black));
  if (lnew == NULL)
  {
   errcode=8;
   goto nbend;
  }
 }
 else
 {
  lnew=emptlbl;
  emptlbl=emptlbl->pbl;
 }/*if*/;

 blk[inew][jnew]=lnew;
  lnew->bli=inew;
  lnew->blj=jnew;
  lnew->blx=xnew;
  lnew->bly=ynew;
  lnew->blr=rnew;
  lnew->bmatur=bmnew;
  lnew->bagd=-1;
  lnew->bpepx=xnew;
  lnew->bpepy=ynew;
  lnew->nab1=0;
  lnew->last_il_act=0;
  lnew->stress=0;
  lnew->firstkill=0;
  
  for ( k=0;k < 3; k++)
  {
     lnew->mhc_stress[k]=0;
     lnew->touch_t[k]=-1;
	 lnew->mhcth[k]=-1;
     lnew->mhc2x[k]=pxnew[k];
     lnew->mhc2y[k]=pynew[k];          
  }
  
  if (ldelay == 0)
   lnew->blt1=t;
  else
   lnew->blt1=tmax+1;

  lnew->blt2=tmax+1;
  lnew->pbl=lold;

 if (ldelay == 0)	 /*if (ldelay = 0) a native black is to be born, otherwise an offspring*/
	blborn(0,lnew);
 else
  include(1,lnew,0,rt);		/*blborn*/

 nbend:;
}/*newblack*/;

/**********************************************************************/
/*a black is removed*/

 void blremove(pntrbl l,short i,short j)
/**********************************************************************/

{
 pntrbl pbli,pblold;

 if (l == blk[i][j])
  blk[i][j]=l->pbl;
 else
 {
  pblold=blk[i][j];
  while (pblold != NULL )
  {
   pbli=pblold->pbl;
   if (l == pbli)
   {
    pblold->pbl=l->pbl;
    goto lfound;
   }
   else
    pblold=pbli;
  }/*while*/;

  errcode=20;
  goto lrout;

 }/*if*/;

 lfound:;

  l->pbl=emptlbl;
  emptlbl=l;
  l->bli=0;
  l->blj=0;
  l->blx=0;
  l->bly=0;
  l->blr=0;
  l->blt1=0.0;
  l->blt2=0.0;
  l->bmatur=0;
  l->bagd=0;
  l->bpepx=0;
  l->bpepy=0;
  l->nab1=0;
  l->stress=0;
  l->last_il_act=0;
  l->firstkill=0;
  
  int ll;
  for( ll=0;ll < 3;ll++)
  {
     l->mhc2x[ll]=0;
     l->mhc2y[ll]=0;
     l->mhc_stress[ll]=0;
     l->touch_t[ll]=0;
	 l->mhcth[ll]=-1;
  }

 lrout:;
}/*blremove*/;


/**********************************************************************/
/*a black dies*/

 void bldies(pntrbl l,short plus)
/**********************************************************************/

{
    short i,j,xx,yy,rr,m,bmat;
    float t1,t2,lam;
    short lab;

	lab=0;

  if ((l->blt1 > t) || (l->blt2 < t) || (l->bmatur > 4))		 /*checking the existence of the black*/
  {
   errcode=21;
   goto lout;
  }

  i=l->bli;
  j=l->blj;
  xx=l->blx;
  yy=l->bly;
  rr=l->blr;
  t1=l->blt1;
  bmat=l->bmatur;

  if (l->nab1 > 0)     /*B cell dies, but antibodies remain*/
  {
   l->bmatur=5;
   lab=1;
  }

  if (plus == 1)
  {
   l->blt2=t;
   nbdied=nbdied+1;
  }
  else
  {
   l->blt2=-t;
   nbkilled=nbkilled+1;
  }

  t2=l->blt2;

 nb=nb-1;

 if (bmat != 4)			 /*deleting from the actual list: life and action*/
 { 
  lam=nb*nnil*tauil;
   if(nnil>0)
   {
    include(0,NULL,22,lam);		/*modifying il2 actions*/        
   }
   if (bmat>=1) 
   {
		remov(l,24);                         /*removing normal action (attack_bl)*/  
   }
   if (bmat==0) 
   {
		remov(l,28);                         /*removing bone_marrow_filter*/  
   }	
   if (l->firstkill==1) 
   {
		remov(l,23);                         /*removing stress_control_bcell*/           
   }  
   if (l->stress > 0)
   { 
     remov(l,20);                       	/*removing il1 production (blprodil)*/           
     bl_stress--;
	 l->stress=0;	 
   }
 }        

 remov(l,2);		/*bldies*/

 if (bmat == 4)     /*if it is a plasma cell to stop Ab production*/
  remov(l,9);

 if (writtype >= 1)		 /*writing the black record to the output file*/
  fprintf(outd,"black: %4d%4d%6d%6d%6d%3d%6d%9ld%12.1f%12.1f\n",i,j,xx,yy,rr,bmat,l->stress,l->nab1,t1,t2);

 ldel:;

 if (1-lab)		 /*moving the black record from the used list to the empty list if (no antibodies remain*/
  blremove(l,i,j);

 /*tbfreem();*/	 	/*modifying the time parameters of free B and Th cells*/
 /*tthfreem();*/

 lout:;

}/*bldies*/;

/**********************************************************************/
/*a plasma cell creates an amount of antibodies*/

 void plasmab(pntrbl l)
/**********************************************************************/

{
 float lam,lamd,weigr;

  nab=nab+1;
  l->nab1=l->nab1+1;

  /*the objectives are in the conjugate cells*/

  weigr=scountob(l);

  lam=weigr*tauab*(l->nab1);
  lamd=taudab*(l->nab1);

  if (l->nab1 == 1)
  {
   include(1,l,10,lamd);
   include(1,l,11,lam);
  }
  else
  {
   include(0,l,10,lamd);
   include(0,l,11,lam);
  }/*if*/;

 /*modifying the time parameters of free B cells*/

 tbfreem();

}/*plasmab*/;


/**********************************************************************/
/*an amount of antibodies dies*/

 void abdies(pntrbl l)
/**********************************************************************/

{
 float lam,lamd,weigr;

  nab=nab-1;
  l->nab1=l->nab1-1;

  if (l->nab1 == 0)
  {
   remov(l,10);
   remov(l,11);

   /*removing the black completely*/

   if (l->bmatur == 5)
    blremove(l,l->bli,l->blj);

  }
  else
  {
   /*the objectives are in the conjugate cells*/

   weigr=scountob(l);

   lam=weigr*tauab*(l->nab1);
   lamd=taudab*(l->nab1);
   include(0,l,10,lamd);
   include(0,l,11,lam);
  }/*if*/;

 /*modifying the time parameters of free B cells*/

 tbfreem();

}/*abdies*/;


/**********************************************************************/
/*B cell or antibodies kill cells*/

 void killer(short borab,pntrbl l,short xc,short yc,short rr,short nbty,
             short nwty,short nrty,float weigb,float weigw,float weigr,
             short *res,short *dd,short *pepx,short *pepy)
/**********************************************************************/

{
    pntrbl lt;
    pntrw lw;
    pntrr lr;
    short ii;
    float rn,rtime,rt;
	*res=0;
	rn=weigr*randr();		/*random selection of a target*/
	if (rn <= weigb)
    {
		/*a black-black reaction*/
		for ( ii=1; ii<= nbty; ii++) 
		{
			if (rn <= cdf[ii])
			goto loutb;
		}/*ii*/;
		loutb:;
		lt=target[ii];
		*dd=d(xc,-yc,lt->blx,lt->bly);
		*pepx=lt->bpepx;
		*pepy=lt->bpepy;
			*res=2;
			/*B cell or antibody?*/
			if (lt->nab1 == 0)                         /*B cell only*/
			{
				if ((l == lt) && (borab == 1))
				{
				errcode=22;
				goto loend;
				}
				if (blkill == 1)
					bldies(lt,-1);
				else
				{
					*res=1;
					*pepx=epepx;
					*pepy=epepy;
				}
			}
		else
      {

       if ((lt->blt2 > t) && ((lt != l) || (borab != 1)))
       {                                   /*B cell plus antibodies*/
        rn=(1.0+(lt->nab1))*randr();

        if (rn < 1)                         /*the B cell is chosen*/
        {
         if ((blkill == 1) && (l != lt))
          bldies(lt,-1);
         else
         {
			*res=1;
			*pepx=epepx;
			*pepy=epepy;
         }
        }
        else                                   /*antibody is chosen*/
         abdies(lt);

       }
       else                         /*Otherwise only antibodies are here*/
        abdies(lt);

      }/*if*/;
   }
   else
   {
    if (rn <= weigw)
    {


     for ( ii=nbty+1; ii<= nbty+nwty; ii++) 
     {
      if (rn <= cdf[ii])
       goto loutw;
     }/*ii*/;

     loutw:;

     lw=targetw[ii-nbty];

      *dd=d(xc,-yc,lw->xw,lw->yw);
      *pepx=lw->wpepx;
      *pepy=lw->wpepy;

       *res=2;
       lw->nw=(lw->nw)-1;

       if (lw->nw == 0)  /*remove from the actual list*/
        remov(lw,3);
       else
       {
        rtime=(lw->tauw)*((lw->nw)*gate(lw->nw,0.0,lw->thnw,lw->etanw));
        include(0,lw,3,rtime);
       }

       lw->nwkilled=lw->nwkilled+1;

       /*modifying the time parameters of free B cells*/

       tbfreem();

    }
    else
    {

    /*a red dies*/

     for ( ii=nbty+nwty+1; ii<= nbty+nwty+nrty; ii++) 
     {
      if (rn <= cdf[ii])
       goto loutr;
     }/*ii*/;

     loutr:;

     lr=targetr[ii-nbty-nwty];

      *dd=d(xc,-yc,lr->xr,lr->yr);
      
      *pepx=lr->rpepx;
      *pepy=lr->rpepy;
      *res=2;  
	   
       lr->nr=(lr->nr)-1;  
	   		
		if ((0 <= lr->nr) && (lr->nr < 10))  {

		fprintf(out2,"%s%4d%8.2f\n"," Number of red cells: ",lr->nr,t);
					
		}/*if*/

       if (lr->nr == 0)  /*remove from the actual list*/
        remov(lr,4);
       else
       {
        rtime=(lr->taur)*((lr->nr)*gate(lr->nr,0.0,lr->thnr,lr->etanr));
        include(0,lr,4,rtime);
       }
       lr->nrkilled=lr->nrkilled+1;
       /*modifying the time parameters of free B cells*/
       tbfreem();
    }
   }/*if*/;

  if ((*res == 2) && (l->firstkill==0) && (mediumreprod)) {
		rt=bcell_tau_stress_control;  					/* *(nm*gate(nm,0.0,thnmb,etanmb));*/
		include(1,l,23,rt); 
		l->firstkill=1;
  }

   loend:;

}/*killer*/;


/**********************************************************************/
/*an amount of antibodies acts*/

 void abacts(pntrbl l)
/**********************************************************************/

{
    short xc,yc,rr,res,dd,bmat,nab1t,nbty,nwty,nrty,pepx,pepy;
    long n2;
    float lam,lamd,weigb,weigw,weigr,rn,pba;
 /*checking the existence of the antibodies*/
  if ((l->blt1 > t) || (l->bmatur <= 3))
  {
   errcode=23;
   goto abend;
  }
  xc=l->blx;
  yc=l->bly;
  rr=l->blr;
  dd=l->bagd;
  bmat=l->bmatur;
  nab1t=l->nab1;

 /*the objectives are in the conjugate cells*/
 n2=countobj(l,xc,-yc,rr,0,&nbty,&nwty,&nrty,&weigb,&weigw,&weigr);

 if (errcode > 0)
  goto abend;

 if (n2 == 0)      /*there exists no object*/
  pba=0.0;
 else
 {
  pba=(1.0-gater(weigr,0.0,thba,etaba))/**gate(rr,0.0,thr,etar)*/;   /*!!!!!*/
 }/*if*/;

 rn=randr();

 if (rn >= pba)      /*nothing happens*/
 {
  lam=weigr*tauab*nab1t;
  include(0,l,11,lam);
  goto abend;
 }
                          /*the antibodies are being attached to an object
                          ************************************************/

 killer(0,l,xc,yc,rr,nbty,nwty,nrty,weigb,weigw,weigr,&res,&dd,&pepx,&pepy);

 if (res <= 1)      /*!!! there was no actual killing, no antibody was used*/
  goto abend;

 nab1t=l->nab1;

 if (nab1t == 0)      /*it could kill its own last antibody too!*/
 {
  /*modifying the time parameters of free B cells*/

  tbfreem();

  goto abend;
 }

 nab=nab-1;
 nab1t=nab1t-1;
 l->nab1=nab1t;

 if (nab1t == 0)
 {
  remov(l,10);
  remov(l,11);

  /*removing the black completely*/

  if (bmat == 5)
   blremove(l,l->bli,l->blj);

 }
 else
 {
  lam=weigr*tauab*nab1t;
  lamd=taudab*nab1t;
  include(0,l,10,lamd);
  include(0,l,11,lam);
 }/*if*/;

 /*modifying the time parameters of free B cells*/

 tbfreem();

 abend:;

}/*abacts*/;


/**********************************************************************/
/*reproduction process of an activated black*/

 void reprodbl(short reprod_type, pntrbl l)
/**********************************************************************/

{
   short xc,yc,rr,bmat,i1,i2,j1,j2,res,dd,nbty,nwty,nrty,
    rnew,rsprd,xnew,ynew,inew,jnew,bmnew,ndumm;
    long n0,n2,conc;
    float thpara1,thpara2,weig,weigb,weigw,weigr,rn,pba,pbr,thb,weakconc;
    int k;
    short ind;
    short pepx[3];
    short pepy[3];
    
    xc=l->blx;
    yc=l->bly;
    rr=l->blr;
    dd=l->bagd;
    bmat=l->bmatur;
      
    if (rcell >= xmax)
     n0=nb;
    else
    {
      bounds(xc,yc,rcell,&i1,&i2,&j1,&j2);
      n0=countbl(l,i1,i2,j1,j2,xc,yc,rcell,3,&ndumm);
	  
      if (errcode > 0)
         goto rend;
     }/*if*/;

	 	bounds(xc,yc,rr,&i1,&i2,&j1,&j2);
		conc=scountob(l)-countbl(l,i1,i2,j1,j2,xc,yc,rr,3,&ndumm);
		thpara1=exp(etaparalize*log(thparalize1));
		thpara2=exp(etaparalize*log(thparalize2));	
		thb=exp(etaccb*log(thccb*nm)); 
		weakconc=nm;
   
   	  if (reprod_type == 0) {
			pbr=kb0*gate(dd,0.0,rminb+thd,etad)*(1-gate(dd,0.0,rminb-thd,etad))*gate(rr,0.0,thr,etar)*gate(n0,0.0,thb,etaccb)*(1-gate(conc,0.0,weakconc,etaparalize));
			
			fprintf(out2,"B cell weak reprod try: \t t=%7.2f \t bmatur=%3d \t gdd1=%3.2f \t gdd2=%3.2f \t grr=%3.2f \t gn0=%3.2f \t gconc2=%3.2f \t pbr=%4.3f \n", t, bmat, gate(dd,0.0,rminb+thd,etad), (1-gate(dd,0.0,rminb-thd,etad)), gate(rr,0.0,thr,etar), gate(n0,0.0,thb,etaccb), 1-gate(conc,0.0,weakconc,etaparalize), pbr);
			
			}
	  else if (reprod_type == 1) {
			pbr=kb1*gate(dd,0.0,thd,etad)*gate(rr,0.0,thr,etar)*gate(n0,0.0,thb,etaccb)*gate(conc,0.0,thpara2,etaparalize)*(1-gate(conc,0.0,thpara1,etaparalize)); 
			
			fprintf(out2,"B cell medium reprod try: \t t=%7.2f \t bmatur=%3d \t gdd=%3.2f \t grr=%3.2f \t gn0=%3.2f \t gconc1=%4.3f \t gconc2=%3.2f \t pbr=%4.3f \n", t, bmat, gate(dd,0.0,thd,etad), gate(rr,0.0,thr,etar), gate(n0,0.0,thb,etaccb), gate(conc,0.0,thpara2,etaparalize), 1-gate(conc,0.0,thpara1,etaparalize), pbr);	
	  }
	  else if (reprod_type == 2) {	
			pbr=kb2*gate(dd,0.0,thd,etad)*gate(rr,0.0,thr,etar)*gate(n0,0.0,thb,etaccb)*gate(conc,0.0,thpara2,etaparalize)*(1-gate(conc,0.0,thpara1,etaparalize));
	
			fprintf(out2,"B cell strong reprod try: \t t=%7.2f \t bmatur=%3d \t gdd=%3.2f \t grr=%3.2f \t gn0=%3.2f \t gconc1=%4.3f \t gconc2=%3.2f \t pbr=%4.3f \n", t, bmat, gate(dd,0.0,thd,etad), gate(rr,0.0,thr,etar), gate(n0,0.0,thb,etaccb), gate(conc,0.0,thpara2,etaparalize), 1-gate(conc,0.0,thpara1,etaparalize), pbr);	
	  }

     rn=randr();
	 
     if(rn<pbr) 
     {  
	    if ((reprod_type == 2) || ((reprod_type == 1) && (l->bmatur == 1)))		/*|| (l->bmatur == 2)) to remove?*/
		{
			bmnew=l->bmatur+1;     		/*if bmatur = 3 the memory cell becomes a plasma cell immediately*/
			if (bmnew == 3)    			/*the new B cell becomes a plasma or memory cell*/
			   {
				if (randr() > pmem)  	/*becomes plasma cell*/
				 bmnew=4;
			   }
		}
		else
		{
			bmnew=l->bmatur;
        }

        if (errcode > 0)
           goto rend;

        if (randr() < pmut)				/*mutation*/
        {                                  
           rnew=rho(dd,rr);
           rsprd=rspread(dd);
           unif(xc,yc,rsprd,srchtype,xmax,xmax2,&xnew,&ynew);
        }
        else							/*no mutation, only perturbation for the better visualization*/
        {                           
           xnew=xc+(randi(2*perturb+1)-perturb);
           ynew=yc+(randi(2*perturb+1)-perturb);
           rnew=rr;
        }

        inew=1+(xnew / delx);
        jnew=1+((xmax2+ynew) / delx);
        
        for(k=0;k<3;k++)
        {
           ind=randi(3);   
           pepx[k]=epepx;
           pepy[k]=epepy;                      
        }
		
		if (writtype == 2) 
		{ 

				if (reprod_type == 0 ) 
				{
					fprintf(out2,"%25s","B cell weak reprod:");

				}
				else if (reprod_type == 1) 
				{
					fprintf(out2,"%25s","B cell medium reprod:");
					
				}
				else if (reprod_type == 2) 
				{
					fprintf(out2,"%25s","B cell strong reprod:");
				
				}
	
				fprintf(out2,"%10.2f%8s%2d%8d%8d%8d%8d%6s%8d%10s%8d%8s%8d%5d%5d%5d%5d%5d%5d \n",t,"bmatur:",bmnew,xc,yc,xnew,ynew,"rnew:",rnew,"dd old:",dd,"rr:",rr,l->mhc2x[0],l->mhc2y[0],l->mhc2x[1],l->mhc2y[1],l->mhc2x[2],l->mhc2y[2]);
	
		}/*if*/
			
		newblack(inew,jnew,xnew,ynew,rnew,bmnew,pepx,pepy,1,taubr); 

     }
	 
	
    
  rend:;
  
};


/**********************************************************************/
/*a white is born*/

 void wborn(pntrw lw)
/**********************************************************************/

{
 float rtime;

  if (lw->nw > 0)
  {

    lw->nw=(lw->nw)+1;
    rtime=(lw->tauw)*((lw->nw)*gate(lw->nw,0.0,lw->thnw,lw->etanw));

    include(0,lw,3,rtime);

  }/*if*/;

 /*modifying the time parameters of free B cells*/

 tbfreem();

}/*wborn*/;


/**********************************************************************/
/*a red is born*/

 void redborn(pntrr lr)
/**********************************************************************/

{
 float rtime;
 
  if (lr->nr > 0)
  {
    lr->nr=(lr->nr)+1;
	
    if (lr->nr > nrmax)
    {
     errcode=6;
     goto rbend;
    }
    rtime=(lr->taur)*((lr->nr)*gate(lr->nr,0.0,lr->thnr,lr->etanr));
    include(0,lr,4,rtime);
  }

 /*tbfreem();*/		 /*modifying the time parameters of free B cells*/

 rbend:;
}/*redborn*/;


/**********************************************************************/
/*a bone marrow cell is born*/

 void bmborn(void)
/**********************************************************************/

{
 float rtime;
  short j;
 pntrth k;
 

  nm=nm+1;

  /*modifying the next birth of a bone marrow cell*/

  rtime=taum*(nm*gate(nm,0.0,thnm,etanm));

  include(0,NULL,7,rtime);

  if (errcode > 0)
   goto bmbend;

  /*modifying the native black production by bone marrow cells*/

  rtime=taubm*(nm*gate(nm,0.0,thnmb,etanmb));

  include(0,NULL,5,rtime);

  if (errcode > 0)
   goto bmbend;

  /*modifying the native T helper production by bone marrow cells*/

  rtime=tauthm*(nm*gate(nm,0.0,thnmth,etanmth));

  include(0,NULL,8,rtime);

  if (errcode > 0)
   goto bmbend;

  bmbend:;

}/*bmborn*/;


/**********************************************************************/
/*the bone marrow produces a black*/

 void bmprodb(void)
/**********************************************************************/

{
 short xnew,ynew,inew,jnew;
 int l,count,k;
 short mhc2x[3],mhc2y[3];
 
  unif(xmax2,0,xmax2,1,xmax,xmax2,&xnew,&ynew);   /*uniform generation on the whole configuration space*/
  inew=1+(xnew / delx);
  jnew=1+((xmax2+ynew) / delx);
  
  for(l=0;l<3;l++)
  {
    if(nwtypes > 1)
    {
       k=randi(nwtypes)+1;   
       pw=pw1;
       count=1;
       while (pw != NULL && count < k)
       {
          count++;
          pw=pw->pnw;
       }
    } 
    else
    {
       pw=pw1;
    }
	mhc2x[l]=epepx;			/*empty peptide*/
    mhc2y[l]=epepy;		
  } 
  newblack(inew,jnew,xnew,ynew,r0,0,mhc2x,mhc2y,0,0.0);
}/*bmprodb*/;

/**********************************************************************/
/*a T helper is born*/

 void thborn(short ldelay,pntrth l)
/**********************************************************************/

{
    short xx,yy,dthw;
    float weig,lam2,rt;
    short lnorm;
	

  nth=nth+1;
  
  if(nnil1 > 0)/*modifying il1 actions*/
  {
    lam2=nnil1*nth*tauil1;
    include(0,NULL,21,lam2);
  }

  xx=l->tcrx;
  yy=l->tcry;

  if (ldelay != 0)
  {
   l->tht1=t;

   /*removing the birth from the actual list*/

   remov(l,14);

  }/*if*/;

 /*finally, the T helper shall die*/

 include(1,l,15,tlifeth);

 if (errcode > 0)
  goto rend;

 /*preparation for ( the first action*/

 /*the objectives are in the conjugate cells*/

 /*weig=scounpep(l);*/

 include(1,l,25,tauth);/*include normal action*/  /*removed weig* factor*/
 include(1,l,26,tauthymus);

 if (errcode > 0)
  goto rend;

 if (ldelay != 0)
  goto rend;         /*only newly born T helper cells should be tested*/

 
 rend:;

}/*thborn*/;



/**********************************************************************/
/*creates or prepares the birth of a new T helper cell*/

 void newth(short xnew,short ynew,short ldelay,float rt,short *thmake)
/**********************************************************************/

{
    pntrth lold,lnew;
    float pth,gvp,gvn,rn0;
    short dthw,dmin;

 *thmake=0;

 lold=thlast;

 if (emptlth == NULL)
 {
  lnew=(struct thelper *) malloc(sizeof(struct thelper));  

  if (lnew == NULL)
  {
   errcode=15;
   goto nthend;
  }
 }
 else
 {
  lnew=emptlth;
  emptlth=emptlth->pth;
 }/*if*/;

 thlast=lnew;

  lnew->tcrx=xnew;
  lnew->tcry=ynew;
  
  if (ldelay == 0)
   lnew->tht1=t;
  else
   lnew->tht1=tmax+1;

  lnew->tht2=tmax+1;
  lnew->thpd=-1;
  lnew->pth=lold;
  lnew->stress=0;
  lnew->last_il_act=-tcrit_il-1;
  lnew->tmatur=0;

 /*if (ldelay = 0) a native T helper is to be born, otherwise an offspring*/

 if (ldelay == 0)

  thborn(0,lnew);

 else

  include(1,lnew,14,rt);

 nthend:;

}/*newth*/;


/**********************************************************************/
/*the bone marrow produces a T helper cell*/

 void bmprodth(void)
/**********************************************************************/

{
    short xnew,ynew;
    short thmake;

  thmake=1;

  /*uniform generation on the whole peptide space*/

  while (thmake )
  {

   unif(pxmax2,0,pxmax2,1,pxmax,pxmax2,&xnew,&ynew);

   newth(xnew,ynew,0,0.0,&thmake);

  }/*while*/;

}/*bmprodth*/;


/**********************************************************************/
/*a T helper dies*/

 void thdies(pntrth l,short plus)
/**********************************************************************/

{
    short xx,yy,dthw,m,tp;
    float t1,t2,lam,lamd,lam2;
    short lnorm;
    pntrth pthpr,pthi,pthold;


  nth=nth-1;
  
  if(l->stress>0)
  {
    th_stress--;
    remov(l,17);           /*remove thprodil2*/
	remov(l,27);
    l->stress=0;
  }

  if (nnil1 > 0)/*modifying il1 actions*/
  {
    lam2=nnil1*nth*tauil1;   
    include(0,NULL,21,lam2);
  }

 /*checking the existence of the T helper*/

  if ((l->tht1 > t) || (l->tht2 < t))
  {
   errcode=18;
   goto lout; 
  }
  
  xx=l->tcrx;
  yy=l->tcry;
  tp=l->thpd;
  t1=l->tht1;
  pthpr=l->pth;

  if (plus == 1)
  {
   l->tht2=t;
   nthdied=nthdied+1;
  }
  else
  {
   l->tht2=-t;
   nthkild=nthkild+1;
  }

  t2=l->tht2;

 /*deleting from the actual list: life and action*/

 remov(l,25);/*removing normal action*/
 /*remov(l,16);*/
 remov(l,15);

 if  (l->tmatur==0) { 
	remov(l,26);
 }
 

 /*writing the T helper record to the output file*/

if (writtype >= 1)
 fprintf(outd,"Thelper: %6d%6d%6d%12.1f%12.1f\n",xx,yy,l->thpd,t1,t2);

 /*moving the T helper record from the used list to the empty list*/

 if (l == thlast)
  thlast=pthpr;
 else
 {
  pthold=thlast;

  while (pthold != NULL )
  {
   pthi=pthold->pth;

   if (l == pthi)
   {
    pthold->pth=l->pth;
    goto lfound;
   }
   else
    pthold=pthi;

  }/*while*/;

  errcode=18;
  goto lout;


 }/*if*/;

 lfound:;

  l->pth=emptlth;
  emptlth=l;
  l->tcrx=0;
  l->tcry=0;
  l->tht1=0.0;
  l->tht2=0.0;
  l->thpd=0;
  l->last_il_act=0;
  l->tmatur=0;

 lout:;

}/*thdies*/;


/**********************************************************************/
/*a T helper cell creates an amount of interleukins*/

 void thprodil(pntrth l)
/**********************************************************************/

{
 float lam;
 float lam2;

  /*modifying the death rate*/

  nnil=nnil+1;
  lam=taudil*nnil;
  lam2=nnil*nb*tauil;

  if (nnil == 1)
  {
    include(1,NULL,18,lam);
    include(1,NULL,22,lam2);
  }
  else
  {
    include(0,NULL,18,lam);
    include(0,NULL,22,lam2);
  }
}/*thprodil*/;


/**********************************************************************/
/*an amount of interleukins dies*/

 void il2dies(void)
/**********************************************************************/

{
  float lam,lam2;

  nnil=nnil-1;
  lam=taudil*nnil;
  lam2=tauil*nnil*nb;

  if (nnil == 0)
  {
   remov(NULL,18);
   remov(NULL,22);
  }
  else
  { 
    include(0,NULL,18,lam);
    include(0,NULL,22,lam2);/*TODO*/
  }

}/*il2dies*/;



/**********************************************************************/
/*reproduction of T helper cells*/

 void reprodth(short reprod_type, pntrth l)
/**********************************************************************/

{
   pntrbl lt;
pntrw lw;
pntrr lr;

short ii,xc,yc,i1,i2,j1,j2,dd,nbty,n0,n2,xnew,ynew;
float weig,rn,rn1,ptha,pthr,thth,ththn1,rtime,n1;
short thdum;

xc=l->tcrx;
yc=l->tcry;
dd=l->thpd;
n2=countpep(xc,-yc,&weig);

 if (errcode > 0)
  goto rend;
  	
	n0=countth(xc,yc,xmax);	 /*the total number of T helper cells*/ 
	n1=countth(xc,yc,dring);		 /*the number of T helper cells in a neigbour*/

	  if (reprod_type == 0) {		
			/*if ((dd<rmaxth) && (rminth<dd))*/								/*weak reprod*/
			
			thth=exp(etaccth*log(thccth*nm));
			ththn1=exp(etaccthn1*log(thccthn1*nm));
			pthr=kth0*gate(n0,0.0,thth,etaccth)*gate(n1,0.0,ththn1,etaccthn1); 

	  } 
	  else if (reprod_type == 1) {										/*medium reprod*/
			thth=exp(etaccth*log(thccth*nm));
			ththn1=exp(etaccthn1*log(thccthn1*nm));
			pthr=kth1*gate(n0,0.0,thth,etaccth)*gate(n1,0.0,ththn1,etaccthn1);
	  
			fprintf(out2,"T cell medium reprod try: \t t=%7.2f \t tmatur=%3d \t gn0=%3.2f \t gn1=%3.2f \t pthr=%4.3f \n", 	t, l->tmatur, gate(n0,0.0,thth,etaccth), gate(n1,0.0,ththn1,etaccthn1), pthr);
	
	  }
	  else if (reprod_type == 2) {										/*strong reprod*/
			thth=exp(etaccth*log(thccth*nm));
			ththn1=exp(etaccthn1*log(thccthn1*nm));
			pthr=kth2*gate(dd,0.0,thdth,etadth)*gate(n0,0.0,thth,etaccth)*gate(n1,0.0,ththn1,etaccthn1);
			
			fprintf(out2,"T cell strong reprod try: \t t=%7.2f \t tmatur=%3d \t gdd=%3.2f \t gn0=%3.2f \t gn1=%3.2f \t pthr=%4.3f \n", t, l->tmatur, gate(dd,0.0,thdth,etadth), gate(n0,0.0,thth,etaccth), gate(n1,0.0,ththn1,etaccthn1), pthr);
			
			
	  }  

  rn=randr();
 
      if ( rn < pthr)
      {                   
        /* include(0,l,16,tauthr);/*????????????????????*/

         if (errcode > 0)
            goto rend;

         xnew=xc+(randi(2*perturb+1)-perturb);
         ynew=yc+(randi(2*perturb+1)-perturb);
		 
		 
		if (writtype == 2) { 
			if (reprod_type == 0) {
				fprintf(out2,"%25s","T cell weak reprod:");
				fprintf(out2,"%10.2f%8s%2d%8d%8d%8d%8d%5s%8.3f%5s%8.3f\n",t,"tmatur:",l->tmatur,xc,yc,xnew,ynew,"gn0:",gate(n0,0.0,thth,etaccth),"gn1:",gate(n1,0.0,thth,etaccth));
			}
			else if (reprod_type == 1) {
				fprintf(out2,"%25s","T cell medium reprod:");
				fprintf(out2,"%10.2f%8s%2d%8d%8d%8d%8d%5s%8.3f%5s%8.3f\n",t,"tmatur:",l->tmatur,xc,yc,xnew,ynew,"gn0:",gate(n0,0.0,thth,etaccth),"gn1:",gate(n1,0.0,thth,etaccth));
			}
			else if (reprod_type == 2 ) {
				fprintf(out2,"%25s","T cell strong reprod:");
				fprintf(out2,"%10.2f%8s%2d%8d%8d%8d%8d%5s%8.3f%5s%8.3f\n",t,"tmatur:",l->tmatur,xc,yc,xnew,ynew,"gn0:",gate(n0,0.0,thth,etaccth),"gn1:",gate(n1,0.0,thth,etaccth));
			}
	
		
		}/*if*/
		 
				
        newth(xnew,ynew,1,tauthr,&thdum);
		 
      }
      else
      {                 
         l->thpd=-1;
      };
   
 rend:;
}/*reprodth*/;


/**********************************************************************/
/*Positive and negative selection in the Thymus*/

 void thymus(pntrth l)
/**********************************************************************/

{
  short xx,yy,dthw,neg,pos;
  short lnorm;
 
  lnorm=0;
  neg=0;
  pos=0;
  
  xx=l->tcrx;
  yy=l->tcry;
  
  pw=pw1;/*pw0?*/

  while (pw != NULL )
  {
   dthw=d(xx,-yy,pw->wpepx,pw->wpepy);

   if (dthw <= rminth)
   {
		neg=1;
	
		if (randr()<negselp)
		{
			neg=2;
		}

    goto rend;

   }/*if*/;

   if (dthw <= rmaxth)
   {
	lnorm=1;
   }
	
   pw=pw->pnw;

  }/*while*/;

  if (1-lnorm)
  {
	pos=1;
	
	if (randr()<posselp)
	{
	pos=2;
	}
	
  }/*if*/;

  rend:;
 
	if (neg==2) 
	{
		thdies(l,3);
	}
	else 
	{
		if (pos==2) 
		{
			thdies(l,2);
		}
		else 
		{ 
			if ((neg==0) && (pos==0) && (1-comptype))		/*Regulatory T-cells only in the ERS model*/ 
			{
				l->tmatur=2;
			}
			else
			{
				l->tmatur=1;
			}
			remov(l,26);			/*removing thymus for this T-cell*/
		}
	}
	
 }/*thymus*/
 
/**********************************************************************/
/*Negative selection of B cells in the bone marrow*/

 void bone_marrow_filter (pntrbl l)
/**********************************************************************/

{
  short xx,yy,dbw,neg;
 
  neg=0;

  xx=l->blx;
  yy=l->bly;

  pw=pw1; /*pw0?*/

  while (pw != NULL )
  {
	   dbw=d(xx,-yy,pw->xw,pw->yw);

	   if (dbw <= rminb)
	   {
			neg=1;					/*the B cell is near to at least one self cell*/
			if (randr()<bcellselp)
			{
				neg=2;				/*the B cell is filtered*/
			}
		goto rend;
	   }	
	   pw=pw->pnw;
  }

  rend:;
 
	if (neg==2) 
	{
		bldies(l,-1);	
	}
	else 
	{ 
		l->bmatur=1;
	    include(1,l,24,taub);           /*normal actions*/
		remov(l,28);		/*remove bone_marrow_filter*/
	}			

 }/*bone_marrow_filter*/


/**********************************************************************/
void stress_control_bcell(pntrbl l) {
/**********************************************************************/
 int i;
 if (1-comptype)
 {
	if ((l->bmatur<=3) && (l->bmatur>=1)) {
		for (i=0;i<3;i++) {
			if ((l->mhc2x[i]!=epepx) && (l->mhc2y[i]!=epepy)) {
				if(t-l->touch_t[i] > tcrit_stress) {							
					l->mhc_stress[i]=1;
					if (l->stress==0) {
						l->stress=1;
						include(1,l,20,tau_prodil1);	   /*Il1 production starts*/                    
						bl_stress++;	    
					}
					if (t-l->last_il_act <= tcrit_il) {
						l->mhcth[i]=1;				
					}
				} else {
					l->mhc_stress[i]=0;
					if ((l->stress == 1) &&  (l->mhc_stress[0]==0) &&  (l->mhc_stress[1]==0) &&  (l->mhc_stress[2]==0)) {
						l->stress=0;
						remov(l,20);
						bl_stress--;
					}
				}	
			}
		}
	 }
 }
}

/**********************************************************************/
void stress_control_tcell(pntrth l) {
/**********************************************************************/
 
 if (mediumreprod) {
	if ((t-l->last_il_act > tcell_tcrit_il1) && (l->stress==1)) {
			th_stress--;
			l->stress=0;
			remov(l,27); 	/*remove stress_control_tcell process*/
			remov(l,17);	/*remove il2 production*/
		}
	}	
 }

/**********************************************************************/
void blprodil(pntrbl l) {
/**********************************************************************/
 
 float lam;
 float lam2;
  /*modifying the death rate and il1 attack rate*/
  nnil1=nnil1+1;
  lam=taudil1*nnil1;
  lam2=nnil1*nth*tauil1;
  if (nnil1 == 1) {
    include(1,NULL,19,lam);		/*il1dies*/
    include(1,NULL,21,lam2);	/*il1_act*/
  } else {
    include(0,NULL,19,lam);
    include(0,NULL,21,lam2);  
  }
}


/**********************************************************************/
void il1dies(void)
/**********************************************************************/
{
  float lam,lam2;

  nnil1=nnil1-1;
  lam=taudil1*nnil1;
  lam2=nnil1*nth*tauil1;/*modifying the death rate and il1 attack rate*/

  if (nnil1 == 0)
  {
    remov(NULL,19);
    remov(NULL,21);
  }
  else
  {
    include(0,NULL,19,lam);
    include(0,NULL,21,lam2);
  }
}


/**********************************************************************/
void il1_act(void)
/**********************************************************************/
{
 int i,j;
 pntrth th_it;
 float rt;
  
  if(thlast==NULL)
    goto lout;              
  th_it=thlast;
  
  i=randi(nth)+1;/*choosing TH cell*/
  j=1;
  
  while ( j<=i)
  {
  
	if ((th_it->tht1 <= t) && ( th_it->tht2 >= t))
	{
		   j++; 	
	}
	if (j<=i) 
		{
		th_it=th_it->pth;
		} 
  } 
  
  if ((th_it->tht1 > t) || (th_it->tht2 < t))
  {
   errcode=18;
   goto lout;
  }
  
  il1dies();		/*The used IL1 dies*/
  
  th_it->last_il_act=t;
   
	if ((th_it->stress==0) && (th_it->tmatur==1))			/*if Th cell is unstressed, it begins to produce il2, and state changes to 1*/
    {
       th_it->stress=1;
       include(1,th_it,17,taubil); /*il2 production starts*/
	   rt=tcell_tau_stress_control;  /**(nm*gate(nm,0.0,thnmth,etanmth));*/
	   include(1,th_it,27,rt); 		/*stress control process starts*/
       th_stress++;
     } 
	 
  lout:;
}
/*some il1 dies*/

/**********************************************************************/
void il2_act(void)
/**********************************************************************/
{
 int i,k,ii,jj;
  pntrbl bl_it;
  
  ii=1;
  jj=1;
  
  i=randi(nb)+1;/*choosing a B cell*/
  
  k=1;
  
for (ii=1; ii<= nsubint; ii++) 
 {
  for (jj=1; jj<= nsubint; jj++) 
  {
   bl_it=blk[ii][jj];

   while (bl_it != NULL )
   {
    if(k==i)
     goto fend;
   
    bl_it=bl_it->pbl;
    k=k+1;
   
   }/*while*/;
  }/*j*/;
 }/*i*/;
   
  fend:;
  
  if ((bl_it != NULL) &&  (bl_it->blt1 <= t) && (bl_it->blt2 >= t) && (bl_it->bmatur <= 3) && (bl_it->bmatur >= 1)) 
  {  
		bl_it->last_il_act=t;
		il2dies(); 		/*The used IL2 dies*/
  }

}

/**********************************************************************/
/*B cell attempts to kill*/
void attack_bl(pntrbl l)
/**********************************************************************/
{
  short xc,yc,rr,bmat,i1,i2,j1,j2,res,dd,nbty,nwty,nrty,
  rnew,rsprd,xnew,ynew,inew,jnew,bmnew,ndumm,pepx,pepy,ind,mhctha;
  long n0,n2;
  float weigb,weigw,weigr,rn,pba,pbr,thb;

  if ((l->blt1 > t) || (l->blt2 < t) || (l->bmatur > 3))
  {
   errcode=24;
   goto rend;
  }

  xc=l->blx;
  yc=l->bly;
  rr=l->blr;
  dd=l->bagd;
  bmat=l->bmatur;

  if (bmat==0) 	
  {
  goto rend;
  }
  
  n2=countobj(l,xc,-yc,rr,1,&nbty,&nwty,&nrty,&weigb,&weigw,&weigr); 
  /*counting objects*/

 if (errcode > 0)
	goto rend;
 if (n2 == 0)     
  pba=0.0;
 else
 {
	pba=(1.0-gater(weigr,0.0,thba,etaba)) ; /*probability of killing*/
 };
 rn=randr();
 if (rn < pba)
 {                
  killer(1,l,xc,yc,rr,nbty,nwty,nrty,weigb,weigw,weigr,&res,&dd,&pepx,&pepy);
  switch (res) {
    case 0: /*l->bagd=-dd-3; break; unsuccesfull kill*/
    case 1: l->bagd=dd; break;/*attached but nothing happens*/
    case 2:/*killing*/
    {
        l->bagd=dd;
		ind=0;
		mhctha=-1;
		while (ind<=2)
		{
			if (mhctha < l->mhcth[ind])
			{
				mhctha=l->mhcth[ind];
			}
			ind=ind+1;
		}	
		if ((mhctha==0) && (1-comptype))
		{
			reprodbl(mhctha,l);
		}
		else if ((mhctha==1) && (t-l->last_il_act<=tcrit_il) && (1-comptype))
		{
			reprodbl(mhctha,l);
			l->mhc_stress[ind]=1;
			if (l->stress==0) 
			{
				l->stress=1;
				include(1,l,20,tau_prodil1);	                       
				bl_stress++;	    
			}
		}
		else if  (mhctha==2) 				/* removed: && (t-l->last_il_act<=tcrit_il) ) */
		{
			reprodbl(mhctha,l);
			
			if (1-comptype) 
			{
				l->mhc_stress[ind]=2;
				
				if (l->stress==0) 
				{
					l->stress=2;
					include(1,l,20,tau_prodil1);	                       
					bl_stress++;	    
				}
				else if (l->stress==1) 
				{
					l->stress=2;
				}
			}
		}

		ind=0;

		while ((l->touch_t[ind]>=0) &&  (ind<=2))	/*searching the first unused mhc*/
		{
			ind=ind+1;
		}/*while*/
		
		if (ind==3)					/*there was no unused mhc*/
		{
		 ind=randi(3);
		 if ((t-l->touch_t[ind] < tcrit_stress) && (mhctha != 2)) 
			{
				l->mhc_stress[ind]=0;		/*set mhc stress zero*/
				if ((l->stress == 1) &&  (l->mhc_stress[0]==0) &&  (l->mhc_stress[1]==0) &&  (l->mhc_stress[2]==0)) 
				{
					l->stress=0;
					remov(l,20);
					bl_stress--;
				}	
			}
		}
        
		l->mhc2x[ind]=pepx;
        l->mhc2y[ind]=pepy;

		l->touch_t[ind]=t;			/*starting point of the critical time*/
		      

       tthfreem();
	   
       break;                    
    
    }
   };
 } 
 else 
 {
/*   include(0,l,24,taub);*/
  /* l->bagd=-1; B cell becomes free*/
 }


 rend:;
}


/**********************************************************************/
void attack_th(pntrth l)
/**********************************************************************/
{
     pntrbl lt;
    pntrw lw;
    pntrr lr;

    short ii,xc,yc,i1,i2,j1,j2,dd,nbty,n0,n2,xnew,ynew,ind;
    float weig,rn,rn1,ptha,pthr,thth,rtime;
    short thdum;

 /*checking the existence of the T helper*/

  if ((l->tht1 > t) || (l->tht2 < t))
  {
   errcode=18;
   goto rend;
  }

  xc=l->tcrx;
  yc=l->tcry;
  dd=l->thpd;
  
 /*the objectives are in the conjugate cells*/

 n2=countpep(xc,-yc,&weig);

 if (errcode > 0)
  goto rend;

 if (n2 == 0)                /*there exists no object*/
  ptha=0.0;
 else
 {
  ptha=1.0-gater(weig,0.0,thpth,etapth);/*TODO ezt pontosítani*/
 }/*if*/;

 rn=randr();
 
 if (rn < ptha)
 {                    /*the T helper cell is being attached to an object
                       **************************************************/
   /*random selection of a target*/

   rn1=weig*randr();

   for ( ii=1; ii<= n2; ii++) 
   {
    if (rn1 <= cdf[ii])
     goto loutb;
   }/*ii*/;

   loutb:;
  
   lt=target[ii];
   
   ind=target_ind[ii];/*mhcII index*/
   
   dd=d(xc,-yc,lt->mhc2x[ind],lt->mhc2y[ind]);
    
   l->thpd=dd;
		
	if ((dd<sreprod_crit) && (lt->bmatur<=3) && (lt->bmatur>=1) && (l->tmatur==1))	/*strong T and B cell reproduction*/
	{
		lt->mhcth[ind]=2;
		reprodth(2,l);
		if (1-comptype)			/*only for ERS model*/
		{
			if (l->stress==0) 
			{
				include(1,l,17,taubil);
				include(1,l,27,tcell_tau_stress_control);
				th_stress++;	
			}
			l->stress=2;
		}
	}
	else 
	{
		if ((weakr) && (lt->mhc_stress[ind] == 0) && (l->tmatur==2) && (lt->bmatur>=1) && (dd<rmaxth) && (rminth<dd))	/*weak T and B cell reproduction*/	
		{	
			lt->mhcth[ind]=0;
			reprodth(0,l);
		}
		else 
		{
			if ((mediumreprod) && (l->tmatur==1) && (l->stress==1) && (lt->bmatur>=1))
			{
				reprodth(1,l);				/*Th cell medium reproduction*/
			}
		}
	}

	lt->touch_t[ind]=t;
 }
 else
 
 {
     l->thpd=-1;/*th cell becomes free*/
 }
rend:;
}/*attack_th*/;


/**********************************************************************/
/*something must happen at a specific instant during simulation*/

 void action(void)
/**********************************************************************/

{
    void *whbwr;
    short et;

 if (awhich != NULL)
 {

   whbwr=awhich->content;
   et=awhich->etype;

  switch (et) {
 
    case 0: blborn(1,whbwr);  /*black is born at the periphery*/
            break;  
/*    case 1: reprodbl(2,whbwr); */ /*reproduction  of black*/
            break;  
    case 2: bldies(whbwr,1);  /*black dies because its lifespan terminates*/
            break;  
    case 3: wborn(whbwr);     /*white is born*/
            break;  
    case 4: redborn(whbwr);   /*red is born*/
            break;  
    case 5: bmprodb();        /*bone marrow creates a native black*/
            break;  
    case 7: bmborn();         /*bone marrow cell is born*/
            break;  
    case 8: bmprodth();       /*bone marrow creates a T helper cell*/
            break;  
    case 9: plasmab(whbwr);  /*a plasma cell creates an amount of antibodies*/
            break;  
   case 10: abdies(whbwr);    /*an amount of antibodies dies*/
            break;  
   case 11: abacts(whbwr);    /*an amount of antibodies acts*/
            break;   
   case 14: thborn(1,whbwr);  /*T helper is born at the periphery*/
            break;  
   case 15: thdies(whbwr,1); /*T helper dies because its lifespan terminates*/
            break;  
 /*  case 16: reprodth(2,whbwr); */ /*reproduction of T helper cell*/
            break;  
   case 17: thprodil(whbwr); /*a T helper creates some interleukin2*/
            break;  
   case 18: il2dies();         /*some interleukin2  die*/
            break;  
   case 19: il1dies();        /*some interleukin1  die*/
            break;
   case 20: blprodil(whbwr);  /*a B cell creates some interleukin1*/
            break;
   case 21: il1_act();          /*il1 action*/
            break;
   case 22: il2_act();          /*il2 action*/
            break;
   case 23: stress_control_bcell(whbwr);   /*a B cell checks its stress state*/
            break;
   case 24: attack_bl(whbwr);   /* B cell action*/
            break;
   case 25: attack_th(whbwr);    /*T helper action*/
            break;
   case 26: thymus(whbwr);    	/*Positive and negative selection in the thymus*/
            break;
   case 27: stress_control_tcell(whbwr);    	/*a T cell checks its stress state*/
            break;
   case 28: bone_marrow_filter(whbwr);    	/*Negative selection of B cells in the bone marrow*/
            break;
  default:
            errcode=1;
            break;  
  }/*switch*/;

 }/*if*/;

 
}/*action*/;


/**********************************************************************/
/*setting initial times*/

 void inittime(void)
/**********************************************************************/

{
    short i,j,k,x,y;
    float rt;


if ((t >= timmst) && (limmunst == 1))
{

 limmunst=0;

/*time generation for the birth of bone marrow*/
/**********************************************/

 rt=taum*(nm*gate(nm,0.0,thnm,etanm));

 include(1,NULL,7,rt);

 if (errcode > 0)
  goto itend;


/*time generation for the action of bone marrow*/
/***********************************************/

 rt=taubm*(nm*gate(nm,0.0,thnmb,etanmb));

 include(1,NULL,5,rt);        /*B cells*/

 if (errcode > 0)
  goto itend;

 rt=tauthm*(nm*gate(nm,0.0,thnmth,etanmth));

 include(1,NULL,8,rt);         /*T helper cells*/

 if (errcode > 0)
  goto itend;

}/*if*/;


/*time generation for the whites*/
/********************************/

 pw=pw1;

 while (pw != NULL )
 {

   if ((t >= pw->t0w) && (pw->wturnon == 0))
   {
    pw->wturnon=1;
    rt=(pw->tauw)*((pw->nw)*gate(pw->nw,0.0,pw->thnw,pw->etanw));

    include(1,pw,3,rt);

    if (errcode > 0)
     goto itend;
   }/*if*/;

   pw=pw->pnw;

 }/*while*/;


/*time generation for the reds*/
/********************************/

 pr=pr1;

 while (pr != NULL )
 {

   if ((t >= pr->t0r) && (pr->rturnon == 0))
   {
    pr->rturnon=1;
    rt=(pr->taur)*((pr->nr)*gate(pr->nr,0.0,pr->thnr,pr->etanr));

    include(1,pr,4,rt);

    if (errcode > 0)
     goto itend;
   }/*if*/;

   pr=pr->pnr;

 }/*while*/;

 itend:;

}/*inittime*/;


/**********************************************************************/
/*determining the next event*/

 void nextevnt(void)
/**********************************************************************/

{
    float sl,randy;

    pntra ai;

/*the time instant and the object of the next event*/

 awhich=NULL;
 sl=0.0;
 ai=actfirst;

 while (ai != NULL )
 {

   sl=sl+ai->lambda;
   ai=ai->pact;

 }/*while*/;

 if (sl != 0.0)
  t=t+randexp(sl);

 randy=randr()*sl;
 sl=0.0;
 ai=actfirst;

 while (ai != NULL )
 {

   sl=sl+ai->lambda;

   if (sl < randy)
   {
    awhich=ai;
    goto levnt;
   }

   ai=ai->pact;

 }/*while*/;

 levnt:;

}/*nextevnt*/;


/**********************************************************************/
/*simulation process*/

 void simulate(void)
/**********************************************************************/

{

 limmunst=1;

 setpar();

 if (errcode > 0)
  goto simend;

/*loop for time and actions*/

 t=tstart;

 storeout(1);

 while (t <= tmax ) {

  if (errcode > 0)
   goto simend;

  inittime();

  if (errcode > 0)
   goto simend;
   
  if (actfirst == NULL) {
  
    t=t+1.0;  
	
  } 
  else {

   nextevnt();

   if (errcode > 0)
    goto simend;

   action();

   if (errcode > 0)
    goto simend;

  }
  
  storeout(0);
  

 }/*t*/;
 

 simend:;


}/*simulate*/;


/**********************************************************************/
/*writes headlines to screen and inputs data*/

 void datainp(void)
/**********************************************************************/

{
 short i;

 printf("A STOCHASTIC MODEL OF THE IMMUNE SYSTEM \n");

 if ((inda=fopen(indas,"r")) == NULL) {perror(indas); exit(1);};

 fscanf(inda,"%d%d%d%d%d%d\n",&comptype,&srchtype,&rndtype,&writtype,&outfreq,&lbradius);
 fscanf(inda,"%d%d%d%f%ld%d\n",&screenout,&lrestore,&mediumreprod,&weakr,&nrmax,&nm);
 fscanf(inda,"%f%f%f%d%d%d%d\n",&timmst,&tbirth,&tmax,&xmax,&nsubint,&r0,&r0s);
 fscanf(inda,"%d%d%f%f%f%d\n",&rcell,&rblgr,&tlifeb0,&tlifmem0,&pmem,&blkill);
 fscanf(inda,"%f%d%f%d%d%f\n",&crnew,&rminnew,&crspread,&rminsprd,&drwidth,&tau_sel_bcell0);
 fscanf(inda,"%d%d%d%d%d%d%d\n",&thrad,&pxmax,&rminth,&rmaxth,&epepx,&epepy,&rminb);
 fscanf(inda,"%f%f%f%f%f%f%f\n",&pmut,&taum0,&taubm0,&taub0,&tauba0,&taubr0,&bcell_tau_stress_control0);
 fscanf(inda,"%f%f%f%f%f%f\n",&taubil0,&taudil0,&tauab0,&taubab0,&taudab0,&tauthm0);
 fscanf(inda,"%f%f%f%f%f%f\n",&tauthymus0,&thccthn1,&tauth0,&etaccthn1,&tauthr0,&tlifeth0);
 fscanf(inda,"%f%f%f%f%f\n",&sreprod_crit,&dring,&tcrit_il,&tcrit_stress,&tau_prodil10);         
 fscanf(inda,"%f%f%f%f%f%f\n",&taudil10,&tauil10,&tauil0,&negselp,&posselp,&bcellselp);                
 fscanf(inda,"%f%f%f%f%f%f%f\n",&thnm,&thnmb,&thparalize1,&thparalize2,&etanm,&etanmb,&etaparalize);
 fscanf(inda,"%f%f%f%f%f%f\n",&thd,&thr,&thccb,&thwd,&tkill0,&tbdiv0);
 fscanf(inda,"%f%f%f%f%f%f\n",&etad,&etar,&etaccb,&etawd,&etapselt,&etanselt);
 fscanf(inda,"%f%f%f%f%f\n",&thba,&thil,&thwdth,&thpth,&thccth);
 fscanf(inda,"%f%f%f%f%f\n",&etaba,&etail,&etawdth,&etapth,&etaccth);
 fscanf(inda,"%f%f%f%f%f\n",&thdth,&thilth,&thnmth,&thnmpos,&thnmneg);
 fscanf(inda,"%f%f%f%f%f\n",&etadth,&etailth,&etanmth,&tcell_tcrit_il1,&tcell_tau_stress_control0);
 fscanf(inda,"%f%f%f%f%f%f\n",&kth0,&kth1,&kth2,&kb0,&kb1,&kb2);
 
 fscanf(inda,"%d%d\n",&nwtypes,&nrtypes);

 /*storing the empty MHC II among the white*/

 pw=(struct white *) malloc(sizeof(struct white));

 if (pw == NULL)
 {
  errcode=11;
  goto inpend;
 }

 pw0=pw;
 pwold=pw;

  pw->nw=0;
  pw->xw=0;
  pw->yw=0;
  pw->wpepx=epepx;
  pw->wpepy=epepy;
  pw->t0w=timmst;
  pw->tauw0=0.0;
  pw->tauw=0.0;
  pw->thnw=0.0;
  pw->etanw=0.0;
  pw->nwkilled=0;
  pw->wturnon=1;

 /*reading the whites*/

 for ( i=1; i<= nwtypes; i++) 
 {
  pw=(struct white *) malloc(sizeof(struct white));

  if (pw == NULL)
  {
   errcode=11;
   goto inpend;
  }

  if (i == 1)
   pw1=pw;

  pwold->pnw=pw;
  pwold=pw;

   fscanf(inda,"%ld%d%d%d%d%f%f%f%f\n",&(pw->nw),&(pw->xw),&(pw->yw),
               &(pw->wpepx),&(pw->wpepy),&(pw->t0w),&(pw->tauw0),&(pw->thnw),
               &(pw->etanw));
   pw->nwkilled=0;

 }/*i*/;

 if (nwtypes == 0)
 {
  pw=NULL;
  pw1=NULL;
  pwold=NULL;
 }
 else
  pw->pnw=NULL;

 /*reading the reds*/

 for ( i=1; i<= nrtypes; i++) 
 {
  pr=(struct red *) malloc(sizeof(struct red));

  if (pr == NULL)
  {
   errcode=12;
   goto inpend;
  }

  if (i == 1)
  {
   pr1=pr;
   prold=pr;
  }
  else
  {
   prold->pnr=pr;
   prold=pr;
  }

   fscanf(inda,"%ld%d%d%d%d%f%f%f%f\n",&(pr->nr),&(pr->xr),&(pr->yr),
               &(pr->rpepx),&(pr->rpepy),&(pr->t0r),&(pr->taur0),&(pr->thnr),
               &(pr->etanr));
   pr->nrkilled=0;

 }/*i*/;

 if (nrtypes == 0)
 {
  pr=NULL;
  pr1=NULL;
  prold=NULL;
 }
 else
  pr->pnr=NULL;

 inpend:;

 fclose(inda);

}/*datainp*/;


/**********************************************************************/
/*writes headlines and input data to screen/to file*/

 void datawrit(void)
/**********************************************************************/

{
 short i;

  if ((outd=fopen(outds,"w")) == NULL) {perror(outds); exit(1);}

  fprintf(outd,"\n");
  tvari1=time(NULL);
  timestr=ctime(&tvari1);
  if (timestr == NULL)
   perror(outds);
  fprintf(outd,"Date and time: %s\n",timestr); 
  fprintf(outd,"\n");

  fprintf(outd,"A STOCHASTIC MODEL OF THE IMMUNE SYSTEM \n");
  fprintf(outd,"\n");


 fprintf(outd,"%s%d%s%d%s%d%s%d%s%d%s%d\n"," comptype= ",comptype,
              " srchtype= ",srchtype," rndtype= ",rndtype," writtype= ",
              writtype," outfreq= ",outfreq," lbradius= ",lbradius);
 fprintf(outd,"%s%d%s%d%s%d%s%5.0f%s%ld%s%d\n"," screenout= ",screenout,
              " lrestore= ",lrestore," mediumreprod= ",mediumreprod," weakr= ",
              weakr," nrmax= ",nrmax," nm= ",nm);
 fprintf(outd,"%s%4.0f%s%5.0f%s%6.0f%s%d%s%d%s%d%s%d\n"," timmst= ",timmst,
              " tbirth= ",tbirth," tmax= ",tmax," xmax= ",xmax," nsubint= ",
              nsubint," r0= ",r0," r0s= ",r0s);
 fprintf(outd,"%s%d%s%d%s%5.1f%s%5.1f%s%5.3f%s%d\n"," rcell= ",rcell,
              " rblgr= ",rblgr," tlifeb0= ",tlifeb0," tlifmem0= ",
              tlifmem0," pmem= ",pmem," blkill= ",blkill);
 fprintf(outd,"%s%5.1f%s%d%s%5.1f%s%d%s%d%s%5.2f\n"," crnew= ",crnew,
              " rminnew= ",rminnew," crspread= ",crspread," rminsprd= ",
              rminsprd," drwidth= ",drwidth," tau_sel_bcell0= ",tau_sel_bcell0);
 fprintf(outd,"%s%d%s%d%s%d%s%d%s%d%s%d%s%d\n"," thrad= ",thrad," pxmax= ",pxmax,
              " rminth= ",rminth," rmaxth= ",rmaxth," epepx= ",epepx,
              " epepy= ",epepy, " rminb= ",rminb);
 fprintf(outd,"%s%5.3f%s%5.1f%s%5.1f%s%5.2f%s%5.2f%s%5.2f%s%5.1f\n"," pmut= ",
              pmut," taum0= ",taum0," taubm0= ",taubm0," taub0= ",
              taub0," tauba0= ",tauba0," taubr0= ",taubr0," bcell_tau_stress_control0= ",bcell_tau_stress_control0);
 fprintf(outd,"%s%5.3f%s%5.3f%s%5.2f%s%5.3f%s%5.3f%s%5.1f\n"," taubil0= ",
              taubil0," taudil0= ",taudil0," tauab0= ",tauab0," taubab0= ",taubab0,
              " taudab0= ",taudab0," tauthm0= ",tauthm0);
 fprintf(outd,"%s%5.1f%s%5.1f%s%5.2f%s%5.2f%s%5.2f%s%5.1f\n"," tauthymus0= ",
              tauthymus0," thccthn1= ",thccthn1," tauth0= ",tauth0," etaccthn1= ",
              etaccthn1," tauthr0= ",tauthr0," tlifeth0= ",tlifeth0);
 fprintf(outd,"%s%5.2f%s%5.2f%s%5.2f%s%5.2f%s%5.2f\n"," sreprod_crit= ",sreprod_crit,
              " dring= ",dring," tcrit_il= ",tcrit_il," tcrit_stress= ",
              tcrit_stress," tau_prodil10= ",tau_prodil10);
 fprintf(outd,"%s%5.2f%s%5.2f%s%5.2f%s%5.2f%s%5.2f%s%5.2f\n"," taudil10= ",taudil10," tauil10= ",
              tauil10," tauil0= ",tauil0," negselp= ",negselp," posselp= ",posselp, " bcellselp= ",bcellselp);            
 fprintf(outd,"%s%5.1f%s%5.1f%s%5.1f%s%5.1f%s%3.1f%s%3.1f%s%3.1f\n"," thnm= ",
              thnm," thnmb= ",thnmb," thparalize1= ",thparalize1," thparalize2= ",thparalize2," etanm= ",etanm,
              " etanmb= ",etanmb," etaparalize= ",etaparalize);
 fprintf(outd,"%s%5.1f%s%5.1f%s%5.1f%s%5.1f%s%5.2f%s%5.2f\n"," thd= ",thd,
              " thr= ",thr," thccb= ",thccb," thwd= ",thwd," tkill0= ",
              tkill0," tbdiv0= ",tbdiv0);
 fprintf(outd,"%s%3.1f%s%3.1f%s%3.1f%s%3.1f%s%3.1f%s%3.1f\n"," etad= ",
              etad," etar= ",etar," etaccb= ",etaccb," etawd= ",etawd,
              " etapselt= ",etapselt," etanselt= ",etanselt);
 fprintf(outd,"%s%5.1f%s%5.1f%s%5.1f%s%5.2f%s%5.3f\n"," thba= ",thba,
              " thil= ",thil," thwdth= ",thwdth," thpth= ",thpth,
              " thccth= ",thccth);
 fprintf(outd,"%s%3.1f%s%3.1f%s%3.1f%s%3.1f%s%3.1f\n"," etaba= ",
              etaba," etail= ",etail," etawdth= ",etawdth," etapth= ",
              etapth," etaccth= ",etaccth);
 fprintf(outd,"%s%5.1f%s%5.1f%s%5.1f%s%5.1f%s%5.1f\n"," thdth= ",
              thdth," thilth= ",thilth," thnmth= ",
              thnmth," thnmpos= ",thnmpos," thnmneg= ",thnmneg);
 fprintf(outd,"%s%3.1f%s%3.1f%s%3.1f%s%3.1f%s%3.1f\n"," etadth= ",
              etadth," etailth= ",etailth," etanmth= ",etanmth," tcell_tcrit_il1= ",
              tcell_tcrit_il1," tcell_tau_stress_control0= ",tcell_tau_stress_control0);
 fprintf(outd,"%s%3.1f%s%3.1f%s%3.1f%s%3.1f%s%3.1f%s%3.1f\n"," kth0= ",
              kth0," kth1= ",kth1," kth2= ",kth2," kb0= ",
              kb0," kb1= ",kb1," kb2= ",kb2);
			  

 fprintf(outd,"%s%d%s%d\n"," nwtypes= ",nwtypes," nrtypes= ",nrtypes);

 pw=pw1;

 for ( i=1; i<= nwtypes; i++) 
 {

   fprintf(outd,"%s%4ld%s%4d%s%4d%s%4d%s%4d%s%4.0f%s%5.1f%s%7.0f%s%3.1f\n",
                "nw=",pw->nw," xw=",pw->xw," yw=",pw->yw," wpx=",
                pw->wpepx," wpy=",pw->wpepy," t0w=",pw->t0w," tw=",
                pw->tauw0," thnw=",pw->thnw," enw=",pw->etanw);
   pw=pw->pnw;

 }/*i*/;

 pr=pr1;

 for ( i=1; i<= nrtypes; i++) 
 {

   fprintf(outd,"%s%4ld%s%4d%s%4d%s%4d%s%4d%s%4.0f%s%5.1f%s%7.0f%s%3.1f\n",
                "nr=",pr->nr," xr=",pr->xr," yr=",pr->yr," rpx=",
                pr->rpepx," rpy=",pr->rpepy," t0r=",pr->t0r," tr=",
                pr->taur0," thnr=",pr->thnr," enr=",pr->etanr);
   pr=pr->pnr;

 }/*i*/;

}/*datawrit*/;


/**********************************************************************/
/*initializing*/

 void initpar(void)
/**********************************************************************/

{
 short i;

 xmax2=xmax / 2;
 delx=xmax / nsubint;
 delx2=delx / 2;
 pxmax2=pxmax / 2;


 thnm=exp(etanm*log(thnm));
 thnmb=exp(etanmb*log(thnmb));


 thd=exp(etad*log(thd));
 thr=exp(etar*log(thr));
 thwd=exp(etawd*log(thwd));

 thba=exp(etaba*log(thba));
 thil=exp(etail*log(thil));
 thwdth=exp(etawdth*log(thwdth));
 thpth=exp(etapth*log(thpth));
 thdth=exp(etadth*log(thdth));
 thilth=exp(etailth*log(thilth));
 thnmth=exp(etanmth*log(thnmth));

 /*thccb and thccth are multipliers of nm !*/

 pw=pw1;

 for ( i=1; i<= nwtypes; i++) 
 {

   pw->wturnon=0;
   pw->thnw=exp((pw->etanw)*log(pw->thnw));
   pw=pw->pnw;

 }/*i*/;

 pr=pr1;

 for ( i=1; i<= nrtypes; i++) 
 {

   pr->rturnon=0;
   pr->thnr=exp((pr->etanr)*log(pr->thnr));
   pr=pr->pnr;

 }/*i*/;

}/*init*/;


/**********************************************************************/
/*writes output data to file*/

 void outfile(void)
/**********************************************************************/

{
    pntrbl l;
    pntrth lt;
    pntra ei;

    short i,j,k;
    short et;
    float lam;

	
 fprintf(outd,"\n");
 fprintf(outd," The numbers of dead cells after the simulation: \n");
 fprintf(outd,"\n");
 fprintf(outd,"%s%6ld%s%6ld%s%6ld%s%6ld"," nbdied= ",nbdied,"  nbkilled= ",
              nbkilled,"  nthdied= ",nthdied,"  nthkilled= ",nthkild);
 fprintf(outd,"\n");

 
 pw=pw1;

 for ( i=1; i<= nwtypes; i++) 
 {

   fprintf(outd,"%4d%s%6ld\n",i,")   nwkilled= ",pw->nwkilled);
   pw=pw->pnw;

 }/*i*/;

 fprintf(outd,"\n");

 pr=pr1;

 for ( i=1; i<= nrtypes; i++) 
 {

   fprintf(outd,"%4d%s%6ld\n",i,")   nrkilled= ",pr->nrkilled);
   pr=pr->pnr;

 }/*i*/;

 fprintf(outd,"\n");

 if (writtype >= 1)
 {
  fprintf(outd,"\n");
  fprintf(outd," The list of black cells at the end: \n");
  fprintf(outd,
  " (indices, coordinates, radius, matur., Ag-d, #Ab, times of birth and death) \n");
  fprintf(outd,"\n");

  for ( i=1; i<= nsubint; i++) 
  {
   for ( j=1; j<= nsubint; j++) 
   {

    k=1;
    l=blk[i][j];

    while (l != NULL )
    {
     
      fprintf(outd,"%4d%4d%6d%s%6d%6d%6d%3d%6d%9ld%12.1f%12.1f\n",l->bli,
                   l->blj,k,"  ",l->blx,l->bly,l->blr,
                   l->bmatur,l->bagd,l->nab1,l->blt1,l->blt2);
		
	  if (writtype == 2) {
		fprintf(out2,"%4d%4d%6d%s%6d%6d%6d%3d%6d%9ld%12.1f%12.1f\n",l->bli,			
                   l->blj,k,"  ",l->blx,l->bly,l->blr,
                   l->bmatur,l->bagd,l->nab1,l->blt1,l->blt2);
	
		fprintf(out2,"%4d%4d%4d%s%4d%4d%4d%s%4d%4d%4d\n",							
					l->mhc2x[0],l->mhc2y[0],l->mhc_stress[0],"  ",
					l->mhc2x[1],l->mhc2y[1],l->mhc_stress[1],"  ",
					l->mhc2x[2],l->mhc2y[2],l->mhc_stress[2]); 
	  }/*if*/
	  
      l=l->pbl;
      k=k+1;
 
    }/*while*/;
   }/*j*/;
  }/*i*/;
  
  fprintf(outd,"\n");
  fprintf(outd,"\n");
  fprintf(outd," The list of black cells in the green domain at the end: \n");
  fprintf(outd,"\n");


  fprintf(outd,"\n");
  fprintf(outd," The list of T helper cells at the end: \n");
  fprintf(outd," (cell coordinates, MHC II+peptide dist., times of birth and death) \n");
  fprintf(outd,"\n");

    lt=thlast;

    while (lt != NULL )
    {
   
     fprintf(outd,"%6d%6d%6d%12.1f%12.1f\n",lt->tcrx,lt->tcry,lt->thpd,
                  lt->tht1,lt->tht2);
     lt=lt->pth;
     
    }/*while*/;
  
  fprintf(outd,"\n");
  fprintf(outd,"\n");
  fprintf(outd," The list of T helper cells under positive sel. at the end: \n");
  fprintf(outd,"\n");

  fprintf(outd,"\n");
  fprintf(outd,"\n");
  fprintf(outd," The list of T helper cells under negative sel. at the end: \n");
  fprintf(outd,"\n");
  
 }/*if*/;


 fprintf(outd,"\n");

}/*outfile*/;


/**********************************************************************/
/*main program*/
/**********************************************************************/

int main( long argc, char *argv[] )

{


 char buffer1[20];	
 char buffer2[20];	

 errcode=0;
 indas = "indat1";



 datainp();
 
 if (argc == 2) 			/*1 plus the number of arguments -> for multiple run*/
	rndtype=atoi(argv[1]);	/* command line argument converted to integer */

 if (errcode > 0)
		goto termnt; 

	
		
 sprintf(buffer1, "outdat-%d", time(NULL));		/* unique file name */
		
 outds = buffer1;

 sprintf(buffer2, "out2-%d", time(NULL));		/* unique file name */
		
 out2s = buffer2;
 
 
 datawrit();

 printf("\n");

 initpar();

 if (comptype)  /*comptype=1 if CRS model*/
 {
 rminth=rmaxth;
 rmaxth=xmax;		/*switching off the positive selection*/
 } 

 if (nsubint > NSUBINTM)
 {
  errcode=2;
  goto termnt;
 }

 if ((nwtypes > NTYPESMX) || (nrtypes > NTYPESMX))
 {
  errcode=3;
  goto termnt;
 }

 if (rndtype == 0)
 {
  srand(time(NULL) % 37);
  rndtype=100;
 }

 for ( idum=1; idum<= rndtype; idum++) 
  dummy=randr();

  simulate();
  
  
  outfile();

 termnt:;
 

 if (errcode != 0)		
  writeerr();

 fprintf(outd,"\n");
 tvari2=time(NULL);
 tdif=difftime(tvari2,tvari1);
 timestr=ctime(&tvari2);
 if (timestr == NULL)
  perror(outds);
 fprintf(outd,"Termination: %s%s%le\n",timestr," Running time in seconds: ",
              tdif);
 fprintf(outd,"\n");
 fclose(outd);

 return 0;

}
