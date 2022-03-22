// Prog.cc
// This Program is

#include "Main.h"
int main (int par_exec1, char* pars_exec[])
{    
//Objects creation
	Mathematics * Math = new Mathematics();
    Physics * Phys = new Physics();
    
	FILE * outfile;
	outfile = fopen("Positions.dat","w");
	
	FILE * outfile2;
	outfile2 = fopen("larmor.dat","w");

    printf("*******PROPAGATION OF COSMIC RAYS*******\n\n");
       
    v0[0]=c/3,v0[1]=0,v0[2]=c/30;
    B0[0]=0,B0[1]=0,B0[2]=6e-7;
    
    modv=Math->Mod(v0);	
	printf("v0=%.3e m/s = %.3f*c \n",modv,modv/c);
	modB=Math->Mod(B0);
	printf("B0=%.3e T = %.3e G\n",modB,modB/1E-4);
	
	rL0=Phys->Larmor(rL,m,q,modv,modB);
	printf("rL=%.3e m = %.3e UA\n\n",rL0,rL0/UA);

    x_max=O[0]+2*rL0,y_max=O[1]+2*rL0,z_max=O[2]+10*rL0;
    x_min=O[0]-2*rL0,y_min=O[1]-2*rL0,z_min=O[2]-0*rL0;

//	x_max=r0[0]+100*rL0,y_max=r0[1]+100*rL0,z_max=r0[2]+100*rL0;
//    x_min=r0[0]-100*rL0,y_min=r0[1]-100*rL0,z_min=r0[2]-0*rL0;
	
	r[3]=z_min;
		
	printf("UA=%.3e \n",UA);

    fprintf(outfile,"%.3e %.3e %.3e %.3e %.3e %.3e\n",x_min/UA,y_min/UA,z_min/UA,x_max/UA,y_max/UA,z_max/UA); //Give the max and minimum values
	
		
    w=Phys->Angvelocity(w,q,modB,m);
    T=Phys->period(T,w);
    t=T/1000; 
    
	printf("w=%e rad/s\n",w);	
    printf("T=%e s\n",T);
    printf("t=%e s\n",t);
    
	int NP=7*int(T/t); //Modificar n° de pontos
	printf("NP=%d \n\n",NP);	
	
//	LOOP PRINCIPAL

	for(i=1;i<=10000;i++){
		
//		if(i<=10)	 printf("i=%2d\n v=(%.9e,%.9e,%.9e) ",i,v0[0],v0[1],v0[2]);		
		//if(i<=3)	 printf("i=%2d\n v=(%.6e,%.6e,%.6e) r=(%.3e,%.3e,%.3e) \n ",i,v0[0],v0[1],v0[2],r0[0],r0[1],r0[2]);

		modr=Math->Mod(r0);	
		//if(i<=30) printf("r0x=%.6e r0y=%.6e r0z=%.6e\n",r0[0],r0[1],r0[2]);		
		//if(i<=30)	 printf("r0=%.3e \n",modr);
		
	//	Phys->Dip(BD,r0);
		
	//	if(i<=30) printf("BDX=%.6e BDY=%.6e BDZ=%.6e\n",BD[0],BD[1],BD[2]);
		
	//	Phys->VarB(B,B0,r0,rL0);
		
		modB=Math->Mod(B0);
		
		Math->Projections(vpar,vper,v0,B0);
		
		if(i<=30) printf("modB=%.6e\n",modB);
		
		//if(i<=30) printf("vperX=%.6e vperY=%.6e vpervZ=%.6e\n",vper[0],vper[1],vper[2]);
		//if(i<=30) printf("vparX=%.6e vparY=%.6e vparZ=%.6e\n",vpar[0],vpar[1],vpar[2]);
		
		modvpar=Math->Mod(vpar);
		modvper=Math->Mod(vper);
		
		//if(i<=30)	printf("modvper=%.6e, modvpar=%.6e \n",modvper,modvpar); 
		
		rL=Phys->Larmor(rL,m,q,modvper,modB);		
		//if(i<=30)	printf("B=%.16e, rL=%.16e \n",modB,rL);
		w=Phys->Angvelocity(w,q,modB,m);
		//vol=(w*t*i)/(2*pi);
		T=Phys->period(T,w);
		
		if(i<=30) printf("T=%.6e, w=%.6e, rL=%.6e \n",T,w,rL);
		
		fprintf(outfile2,"%.3e %.3e %.3e\n",t*i,rL,T); //tempo, raio de larmor, voltas
		
		angle=Phys->Pitch(vper,vpar);		
		if(i<=30) printf("angle=%.16e \n",angle);
		
		Phys->FLorentz(F,q,vper,B0); //Chama a função, passa os valores com os nomes que eu utilizo.

		for(k=0;k<=2;k++) {
			a[k]=F[k]/m;
		}
		
		if(i<=3)	printf("F=(%.3e,%.3e,%.3e) a=(%.3e,%.3e,%.3e) v0=(%.3e,%.3e,%.3e)\n",F[0],F[1],F[2],a[0],a[1],a[2],v0[0],v0[1],v0[2]);
		
		Phys->Velocity(v,vper,a,t); 
		
		//if(i<=30) printf("vX=%.6e vY=%.6e vZ=%.6e\n",v[0],v[1],v[2]);
		
		modv=Math->Mod(v); 
		
		if(i<=30) printf(" modv=%.6e \n",modv);
//		if(i<=3)	 printf(" modv=%.6e \n",modv);
//		if(i<=3)	 printf(" v=(%.3e,%.3e,%.3e) \n",v[0],v[1],v[2]);
		
		//printf("modvpar=%.7e \n",modvpar);
	
		for(k=0;k<=2;k++) v[k]=v[k]+vpar[k];
		
//		if(i<=30) printf("vX=%.6e vY=%.6e vZ=%.6e\n",v[0],v[1],v[2]);
		
		//if(i<=30) printf("vj=%.6e \n",v);
		
		/*double vnorm[3];
		
		Math->Norm(vnorm,v); 
	
  		double v0n=Math->Mod(v0); 
		
  		for(k=0;k<=2;k++) v[k]=vnorm[k]*v0n;*/
  		
		Phys->Trajectory(r,r0,v0,a,t);
		
//		if(i<=30) printf("rX=%.6e rY=%.6e rZ=%.6e\n",r[0],r[1],r[2]);
	
//  	fprintf(outfile,"%.3e %.3e %.3e \n", r[0],r[1],r[2]);
		fprintf(outfile,"%.3e %.3e %.3e \n",r[0]/UA,r[1]/UA,r[2]/UA);
		
		/*
		if(r[0]<x_min){
			x_min=r[0];
		}
		if(r[1]<y_min){
			y_min=r[1];
		}
		if(r[2]<z_min){
			z_min=r[2];
		}
		if(r[0]>x_max){
			x_max=r[0];
		}
		if(r[1]>y_max){
			y_max=r[1];
		}
		if(r[2]>z_max){
			z_max=r[2];
		}
		*/
		
		for(k=0;k<=2;k++) {
			v0[k]=v[k];
			r0[k]=r[k];	
		}
				
	}//fim do loop principal
	
//	printf("...\ni=%2d\n v=(%.3e,%.3e,%.3e) r=(%.3e,%.3e,%.3e) modv=%.16e \n",i,v[0],v[1],v[2],r[0],r[1],r[2],modv);
//	rewind(outfile);
//	fprintf(outfile,"%.3e %.3e %.3e %.3e %.3e %.3e\n",x_min/UA,y_min/UA,z_min/UA,x_max/UA,y_max/UA,z_max/UA); //Give the max and minimum values
	

	fclose(outfile2);
	fclose(outfile);
    printf("\n\n*******END OF RUN*******\n\n");
 	
	system("pause");
	return 0;
}
