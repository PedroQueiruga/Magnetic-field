//Physics.cc

// System Headers:
#include <stdlib.h>

// Local Headers:
#include "Physics.h"
#include "Mathematics.h"

// Constructor
Physics::Physics ( )       
{
  
}

// Destructor
Physics::~Physics ( )
{
        
}

//Lorentz force
float Physics::FLorentz(double F[3], double q, double v[3], double B[3])
{
	
	Mathematics * Math = new Mathematics();
	double prod[3];
	Math->VecProd(prod,v,B);
	
	for(int i=0;i<=2;i++) F[i]=q*prod[i];
	
	
  return 0;
}
//Colocar a saida primeiro / t:Intervalo de tempo

//Velocity
float Physics::Velocity(double v[3],double v0[3],double a[3],double t)
{
  Mathematics * Math = new Mathematics();
  int i,k;
  double v2[3];
  
  //printf(" v0=(%.6e,%.6e,%.6e) a=(%.3e,%.3e,%.3e) t=%.3e \n",v0[0],v0[1],v0[2],a[0],a[1],a[2],t);
  
  for(i=0;i<=2;i++) {
	  v[i]=v0[i]+(a[i]*t);
  } 

  double vv=Math->Mod(v); 
  //printf(" v=(%.6e,%.6e,%.6e), modv=%.12e\n",v[0],v[1],v[2],vv);
   
  Math->Norm(v2,v); 
  //printf(" v/vv=(%.6e,%.6e,%.6e), mod(v/vv)=%.12e \n",v[0]/vv,v[1]/vv,v[2]/vv,sqrt((v[0]*v[0]+v[1]*v[1]+v[2]*v[2])/(vv*vv)));
  //printf(" v2  =(%.6e,%.6e,%.6e), mod(v2)=%.12e \n",v2[0],v2[1],v2[2],sqrt(v2[0]*v2[0]+v2[1]*v2[1]+v2[2]*v2[2]));
  
  double v0n=Math->Mod(v0); 
		
  for(k=0;k<=2;k++) v[k]=v2[k]*v0n;
  //printf(" v=(%.9e,%.9e,%.9e) \n\n",v[0],v[1],v[2]);
 
  return 0;
}

//Particle trajectory
float Physics::Trajectory(double r[3],double r0[3],double v0[3],double a[3], double t)
{
	
	
  for(int i=0;i<=2;i++) r[i]=r0[i]+v0[i]*t+0.5*(a[i]*t*t);
  
  
  return 0;
}

//Angular velocity

double Physics::Angvelocity(double w,double q, double b, double m)
{
	Mathematics * Math = new Mathematics();
	//double b;
	//b=Math->Mod(B);
	w=(q*b)/m;
	
	
  return w;
}

//Moviment period

double Physics::period(double T,double w)
{	
	
	T=(2*pi)/w;
	
  return T;
}

//Larmor Radius

double Physics::Larmor(double rL, double m, double q, double vmod, double bmod){
	
	Mathematics * Math = new Mathematics();
	//double b,vmod2;
	
	//b=Math->Mod(B);
	//vmod2=Math->Mod(v);	
	rL=(m*vmod)/(q*bmod);
	
	return rL;
}

//Pitch angle

double Physics::Pitch(double vper[3],double vpar[3]){
	
	Mathematics * Math = new Mathematics();
	double vmodper,vmodpar, angle;
	
	vmodpar=Math->Mod(vpar);
	vmodper=Math->Mod(vper);
	
	angle=atan(vmodper/vmodpar);
	
	/*printf("vmodpar=%.3e ",vmodpar);
	printf("vmodper=%.3e ",vmodper);
	
	*/
	
	//printf("angle=%.3e \n",angle);
	
	return angle;

}

//Variable field

void Physics::VarB(double B[3],double B0[3],double r[3],double rL){
	
	Mathematics * Math = new Mathematics();
	int i,k;
	double rmod,r1,rb;//r1 fator para diminuir o rmod
	double modB,modB0; 
	
	r1=0.0001*rL;
	rb=100000;
	
	rmod=Math->Mod(r);
	modB0=Math->Mod(B0);	
	
	//printf("\n  r=%.3e, B0=%.3e, ",rmod,modB0);
		
	for(k=0;k<=2;k++) B[k]=B0[k];		
	if(rmod>1){
		for(k=0;k<=2;k++) B[k]=((rb*B0[k]*r1)/rmod);
	}
	modB=Math->Mod(B);	
	//printf("B=%.9e \n",modB);	
	
//dipole field
}


void Physics::Dip(double BD[3], double r[3]){
	
	//BE:Campo dipolar em coordenadas esféricas, E:Conversão do r para esféricas, BD:Campo dipolar em coordenadas cartesianas
	
	Mathematics * Math = new Mathematics();
	int k;
	double E[3],BE[3],mi0=12.5663706143592*pow(10.0,-7),mi=8.00*pow(10.0,22);
	double r1=7*pow(10.0,6);
	
	Math->Coordenada(E,r);
	
//	for(k=0;k<=2;k++) E[k]=vec1[k];
	//printf("r=%.6e phi=%.6e theta=%.6e\n",E[0],E[1],E[2]);
	
	//calculando o campo dipolar em coordenadas esféricas:
	
	BE[0]=(mi0*mi*(cos(E[2])))/(2*pi*E[0]*E[0]*E[0]);
	
	BE[1]=0;
	
	BE[2]=(mi0*mi*(sin(E[2])))/(4*pi*E[0]*E[0]*E[0]);
	
	//printf("BEX=%.6e BEY=%.6e BEZ=%.6e\n",BE[0],BE[1],BE[2]);
	
	//Voltando o campo dipolar para coordenadas cartesianas:
	
	BD[0]=(BE[2]*cos(E[2])+BE[0]*sin(E[2]))*cos(E[2]);
	
	BD[1]=(BE[2]*cos(E[2])+BE[0]*sin(E[2]))*sin(E[2]);
	
	BD[2]= BE[0]*cos(E[2])-BE[2]*sin(E[2]);
	
	//printf("BDX=%.6e BDY=%.6e BDZ=%.6e\n",BD[0],BD[1],BD[2]);
	
}
	


