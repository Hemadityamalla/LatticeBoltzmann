#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#ifndef __APPLE__
#include <malloc.h>
#endif

/* Profilo teorico v_x(y) = 1/2 * (f / nu) * y * (L_y - y) */
/* v_x(max) = 1/8 * (f / nu) * L_y^2 */

/*************************************************************
 *	                                                         *
 *   y                                                       *
 *                                                           *
 *   ^                                                       *
 *   |                                                       *
 *   |                                                       *
 *   |                                                       *
 *   |                                                       *
 *   +----------------------------> x                        *
 *                                                           *
 *************************************************************/

/*********************************************
 *                                           *
 *        x  y                               *
 * 0    (+0,+0)         6    2    5          *
 * 1    (+1,+0)          \   |   /           *
 * 2    (+0,+1)           \  |  /            *
 * 3    (-1,+0)            \ | /             *
 * 4    ( 0,-1)      3 <---- 0 ----> 1       *
 * 5    (+1,+1)            / | \             *
 * 6    (-1,+1)           /  |  \            *
 * 7    (-1,-1)          /   |   \           *
 * 8    (+1,-1)         7    4    8          * 
 *                                           *
 *********************************************/
/*
 * 
 * \nu = c_s^2 (tau -1/2) 
 * D2Q9 c_s^2 = 1/3.
 *
 */

typedef struct {
  double p[9];
} pop;

typedef struct {
  double vx, vy;
} velocity;

#define ff 1e-6
#define tau 9.34465

//#define KOLMO


#define MAX_STEP 1000000

#define NX 20
#define NY 32
#define Uo 1.0 //Only used for kolmogorov


#define cs2  (1.0 /  3.0)
#define cs22 (2.0 *  cs2)
#define cssq (2.0 /  9.0)

#define rt0  (4.0 /  9.0)
#define rt1  (1.0 /  9.0)
#define rt2  (1.0 / 36.0)

/* WARNING vx is rho.vx */
#define vx(a) (a.p[1]+a.p[5]+a.p[8]-a.p[3]-a.p[6]-a.p[7])
#define vy(a) (a.p[2]+a.p[5]+a.p[6]-a.p[4]-a.p[7]-a.p[8])

#define  m(a) (a.p[0]+a.p[1]+a.p[2]+a.p[3]+a.p[4]+a.p[5]+a.p[6]+a.p[7]+a.p[8])

#define delta 0.1

/*
double ff = 8.0*(tau-1.0) * v_max 
*/

/* *************** */
pop p[NY+2][NX+2],p_eq[NY+2][NX+2];

velocity v[NY+2][NX+2];
velocity vold[NY+2][NX+2];
int dirFlag = 1;

pop pzero;

void bc() 
{
  int pp;
  int x, y;
  int xm, xp;

  for (x=0; x<NX+2; x++) {
    for (pp=0; pp<9; pp++) {
      p[   0][x].p[pp] = 0.0;
      p[NY+1][x].p[pp] = 0.0;
    }
  }

  for (y=0; y<NY+2; y++) {
    for (pp=0; pp<9; pp++) {
      p[y][   0].p[pp] = 0.0;
      p[y][NX+1].p[pp] = 0.0;
    }
  }

#ifdef KOLMO

  /* Periodic bc along y */
  for (x=1; x<NX+1; x++) {
    p[0][   x] = p[NY][x];
    p[NY+1][x] = p[1][x];
  }

  /* Periodic bc along x */
  for (y=0; y<NY+2; y++) {
    p[y][   0] = p[y][NX];
    p[y][NX+1] = p[y][ 1];
  }

#else

  for (x=1; x<NX+1; x++) {
    xm = x-1;
    xp = x+1;

    /* At bottom */
    p[0][x ].p[2] = p[1][x].p[4]; 
    p[0][xm].p[5] = p[1][x].p[7]; 
    p[0][xp].p[6] = p[1][x].p[8]; 

    /* At Top */
    p[NY+1][x ].p[4] = p[NY][x].p[2]; 
    p[NY+1][xm].p[8] = p[NY][x].p[6]; 
    p[NY+1][xp].p[7] = p[NY][x].p[5]; 
    
  }

  /* Periodic bc along x */
  for (y=1; y<NY+1; y++) {
    p[y][   0] = p[y][NX];
    p[y][NX+1] = p[y][ 1];
  }

#endif

}

void displace() 
{
  int x, y, pp;
  pop buffer[NY+2][NX+2];
  int xm, xp, ym, yp;

  for (y=1; y<NY+1; y++) {
    ym = y-1;
    yp = y+1;
    for (x=1; x<NX+1; x++) {
      xm = x-1;
      xp = x+1;

      buffer[y][x].p[0] = p[y ][x ].p[0]; 
      buffer[y][x].p[1] = p[y ][xm].p[1]; 
      buffer[y][x].p[5] = p[ym][xm].p[5]; 
      buffer[y][x].p[2] = p[ym][x ].p[2]; 
      buffer[y][x].p[6] = p[ym][xp].p[6]; 
      buffer[y][x].p[3] = p[y ][xp].p[3]; 
      buffer[y][x].p[7] = p[yp][xp].p[7]; 
      buffer[y][x].p[4] = p[yp][x ].p[4]; 
      buffer[y][x].p[8] = p[yp][xm].p[8]; 
    }  
  }

  for (x=1; x<NX+1; x++) {
    for (y=1; y<NY+1; y++) {
      for (pp=0; pp<9; pp++){ 
	p[y][x].p[pp] = buffer[y][x].p[pp];
      }
    }
  }

}

void collide() 
{
  int x, y, pp;
  double u, v;
  double usq, vsq, u2, v2;
  double sumsq, sumsq2;
  double invtau, rho;
  double ui, vi, uv;
  //pop p_eq;

  invtau = 1.0 / tau;
  
  for (y=1; y<NY+1; y++) {
    for (x=1; x<NX+1; x++) {

      /*rho = m(p[y][x]);

      u = vx(p[y][x]) / rho;
      v = vy(p[y][x]) / rho;
      usq = u * u;
      vsq = v * v;

      sumsq  = (usq + vsq) / cs22;
      sumsq2 = sumsq * (1.0 - cs2) / cs2;
      u2 = usq / cssq;
      v2 = vsq / cssq;

      ui = u / cs2;
      vi = v / cs2;
      uv = ui * vi;

      p_eq.p[0] = rho * rt0 * (1.0 - sumsq);
      p_eq.p[1] = rho * rt1 * (1.0 - sumsq + u2 + ui);
      p_eq.p[2] = rho * rt1 * (1.0 - sumsq + v2 + vi);
      p_eq.p[3] = rho * rt1 * (1.0 - sumsq + u2 - ui);
      p_eq.p[4] = rho * rt1 * (1.0 - sumsq + v2 - vi);
      p_eq.p[5] = rho * rt2 * (1.0 + sumsq2 + ui + vi + uv);
      p_eq.p[6] = rho * rt2 * (1.0 + sumsq2 - ui + vi - uv);
      p_eq.p[7] = rho * rt2 * (1.0 + sumsq2 - ui - vi + uv);
      p_eq.p[8] = rho * rt2 * (1.0 + sumsq2 + ui - vi - uv);*/

      for (pp=0; pp<9; pp++)
	p[y][x].p[pp] = p[y][x].p[pp] - invtau * (p[y][x].p[pp] - p_eq[y][x].p[pp]);
    }  
  }
}

void compute_feq()
{

	double u, v;
  	double usq, vsq, u2, v2;
  	double sumsq, sumsq2;
  	double rho;
  	double ui, vi, uv;
	
	
  
  for (int y=1; y<NY+1; y++) {
    for (int x=1; x<NX+1; x++) {

      rho = m(p[y][x]);

      u = (vx(p[y][x]) / rho) ;
      v = (vy(p[y][x]) / rho) ;
      usq = u * u;
      vsq = v * v;

      sumsq  = (usq + vsq) / cs22;
      sumsq2 = sumsq * (1.0 - cs2) / cs2;
      u2 = usq / cssq;
      v2 = vsq / cssq;

      ui = u / cs2;
      vi = v / cs2;
      uv = ui * vi;

      p_eq[y][x].p[0] = rho * rt0 * (1.0 - sumsq);
      p_eq[y][x].p[1] = rho * rt1 * (1.0 - sumsq + u2 + ui);
      p_eq[y][x].p[2] = rho * rt1 * (1.0 - sumsq + v2 + vi);
      p_eq[y][x].p[3] = rho * rt1 * (1.0 - sumsq + u2 - ui);
      p_eq[y][x].p[4] = rho * rt1 * (1.0 - sumsq + v2 - vi);
      p_eq[y][x].p[5] = rho * rt2 * (1.0 + sumsq2 + ui + vi + uv);
      p_eq[y][x].p[6] = rho * rt2 * (1.0 + sumsq2 - ui + vi - uv);
      p_eq[y][x].p[7] = rho * rt2 * (1.0 + sumsq2 - ui - vi + uv);
      p_eq[y][x].p[8] = rho * rt2 * (1.0 + sumsq2 + ui - vi - uv);
	}
	}	

}






void init() 
{
  int x, y, pp;
  double rho, pi=3.141592;
  double u, v, usq, vsq, u2, v2;
  double ui, vi, sumsq, sumsq2, uv;
  
  u = 0.0; /* vx */
  v = 0.0; /* vy */
    
#ifdef KOLMO

  rho = 1.0;
  for (y=1; y<NY+1; y++){
    for (x=1; x<NX+1; x++) {

      u = Uo*sin(2.0*pi*y/(double) NY);
 
      usq = u * u;
      vsq = v * v;

      sumsq  = (usq + vsq) / cs22;
      sumsq2 = sumsq * (1.0 - cs2) / cs2;
      u2 = usq / cssq;
      v2 = vsq / cssq;

      ui = u / cs2;
      vi = v / cs2;
      uv = ui * vi;

      p[y][x].p[0] = rho * rt0 * (1.0 - sumsq);
      p[y][x].p[1] = rho * rt1 * (1.0 - sumsq + u2 + ui);
      p[y][x].p[2] = rho * rt1 * (1.0 - sumsq + v2 + vi);
      p[y][x].p[3] = rho * rt1 * (1.0 - sumsq + u2 - ui);
      p[y][x].p[4] = rho * rt1 * (1.0 - sumsq + v2 - vi);
      p[y][x].p[5] = rho * rt2 * (1.0 + sumsq2 + ui + vi + uv);
      p[y][x].p[6] = rho * rt2 * (1.0 + sumsq2 - ui + vi - uv);
      p[y][x].p[7] = rho * rt2 * (1.0 + sumsq2 - ui - vi + uv);
      p[y][x].p[8] = rho * rt2 * (1.0 + sumsq2 + ui - vi - uv);
    }
  }

#else

  for (pp = 0; pp < 9; pp++)
    pzero.p[pp] = 1.0/9.0;

  for (y=0; y<NY+2; y++)
    for (x=0; x<NX+2; x++) 
	p[y][x] = pzero;

#endif

}

void force() 
{
  double ff_true;
  int x, y;

  ff_true = ff / 6.0;
  for (y=1; y<NY+1; y++)
    for (x=1; x<NX+1; x++) {
      p[y][x].p[1] += ff_true;
      p[y][x].p[5] += ff_true;
      p[y][x].p[8] += ff_true;
      p[y][x].p[3] -= ff_true;
      p[y][x].p[6] -= ff_true;
      p[y][x].p[7] -= ff_true;
    }
}

void write_pop(int tstep)
{
  FILE *fout;
  char fname[128], OutDir[128];
  int x, y, pp, check;

  /* create output directory */

#ifdef KOLMO
  sprintf(OutDir, "Production_KOLMOGOROV");
  check=mkdir(OutDir, 0755);
  if(!check){
    fprintf(stderr,"Directory Production_KOLMOGOROV created. \n");
  }else{
  	if (dirFlag == 1){
    	fprintf(stderr,"Directory Production_KOLMOGOROV already exsist!\n");
    	dirFlag = 0;
    }
  }
#else
  sprintf(OutDir, "Production_POISEUILLE");
  check=mkdir(OutDir, 0755);
  if(!check){
    fprintf(stderr,"Directory Production_POISEUILLE created. \n");
  }else{
  	if (dirFlag == 1){
    	fprintf(stderr,"Directory Production_POISEUILLE already exsist!\n");
    	dirFlag = 0;
    }
  }
#endif


  /* Here dumps the populations */
#ifdef DUMP_POP
  for (pp=0; pp<9; pp++) {
    sprintf(fname,"%s/pop.%d.%d",OutDir,tstep,pp);
    fout = fopen(fname,"w");
    for (y=1; y<NY+1; y++) {
      for (x=1; x<NX+1; x++) 
	fprintf(fout,"%d %d %g %g\n", x, y,
		p[y][x].p[pp],
		p[y][x].p[pp]);
      fprintf(fout,"\n");
    }
    fclose(fout);
  }
#endif

  /* Here dumps the velocity field */
  sprintf(fname,"%s/vel.%d",OutDir,tstep);
  fout = fopen(fname,"w");
  for (y=1; y<NY+1; y++) {
    for (x=1; x<NX+1; x++) 
      fprintf(fout,"%d %d %g %g\n", x, y,
	      vx(p[y][x])/m(p[y][x]),
	      vy(p[y][x])/m(p[y][x]));
    fprintf(fout,"\n");
  }
  fclose(fout);
  
      /*Here dumps the velprofile field */
  sprintf(fname,"%s/velprofile.%d",OutDir,tstep);
  fout = fopen(fname,"w");
  
    for (y=1; y<NY+1; y++){
      double tvx,tvy;
      tvx = vx(p[y][15])/m(p[y][15]);
      tvy = vy(p[y][15])/m(p[y][15]);  
      fprintf(fout,"%d %d %g\n", 15, y, tvx);
	      }
    fprintf(fout,"\n");
  fclose(fout);
  #ifdef KOLMO
  /*Here dumps the velprofile in x direction field */
  sprintf(fname,"%s/velprofilekolmo.%d",OutDir,tstep);
  fout = fopen(fname,"w");
  int yval = 10;
    for (x=1; x<NX+1; x++){
      double tvx,tvy;
      tvx = vx(p[yval][x])/m(p[yval][x]);
      tvy = vy(p[yval][x])/m(p[yval][x]);  
      fprintf(fout,"%d %d %g\n \n", x, yval, tvx);
	      }
    fprintf(fout,"\n");
  fclose(fout);
  #endif
  
      /* Here dumps the density field */
  sprintf(fname,"%s/density.%d",OutDir,tstep);
  fout = fopen(fname,"w");
  for (y=1; y<NY+1; y++) {
    for (x=1; x<NX+1; x++){
      //double tvx,tvy;
      //tvx = (vx(p[y][x])/m(p[y][x]) + 0.5*force[0]);
      //tvy = (vy(p[y][x])/m(p[y][x]) + 0.5*force[1]);  
      fprintf(fout,"%d %d %g\n", x, y, m(p[y][x]));
	      }
    fprintf(fout,"\n");
  }
  fclose(fout);
  
  
}

int main(int argc, char** argv)
{
  FILE *ferr;
  FILE *fmass;
  FILE *velAmp;
  
  int i;

  int x,y;
  double error;
  double tmpx, tmpy;
  double tmpvelAmp = 0.0;
  double mass;

  init();

  ferr  = fopen("error.dat","w");
  fmass = fopen("mass.dat","w");
  velAmp = fopen("velAmp.dat","w");
  for (i=0; i<MAX_STEP; i++) {
    /* fprintf(stderr,"time step %d\n",i); */
    //printf("Time Step: %d \n",i);
    bc();    
    
    for (y=1; y<NY+1; y++) 
      for (x=1; x<NX+1; x++) {
	v[y][x].vx = vx(p[y][x]);
	v[y][x].vy = vy(p[y][x]);
      }
    
    if (i%1 == 0) {
		error = 0.0;
      	for (y=1; y<NY+1; y++) 
			for (x=1; x<NX+1; x++) {
	  			tmpx = (v[y][x].vx - vold[y][x].vx);///vold[y][x].vx;
	  			tmpy = (v[y][x].vy - vold[y][x].vy);///vold[y][x].vy;
	  			tmpvelAmp = 0.0;
	  			/*if (tmpvelAmp > v[y][x].vx)
	  			{
	  				tmpvelAmp = v[y][x].vx;
	  			}*/
	  
	  			error += (tmpx * tmpx + tmpy * tmpy);
			}
	  fprintf(velAmp,"%d %g\n",i,v[32][10]);
	  fflush(velAmp);
      fprintf(ferr,"%d %g\n",i,error);
      fflush(ferr);

      if ((error < 10e-20) && (i!=0)) {
        fprintf(stderr,"Run termalized\n");
	exit(1);
      }
      
      write_pop(i);    

      mass = 0.0;
      for (y=1; y<NY+1; y++) 
	for (x=1; x<NX+1; x++) 
	  mass += m(p[y][x]);
      fprintf(fmass,"%d %g\n",i,mass);  
      fflush(fmass);
    }

    displace(); /* */
    bc();
    compute_feq();
    collide();

#ifndef KOLMO
    force(); 
#endif

    for (y=1; y<NY+1; y++) 
      for (x=1; x<NX+1; x++) {
	vold[y][x].vx = v[y][x].vx;
	vold[y][x].vy = v[y][x].vy;
      }

  } /* for i */
  printf("Max Iterations Done! \n");
  fclose(ferr);
  fclose(fmass);
  fclose(velAmp);
  exit(1);
}
