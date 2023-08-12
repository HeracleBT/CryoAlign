#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include "mrc.h"

extern CMD cmd;

int chkcmdline( int argc, char **argv,CMD *cmd){
        char **p;
	int num=0;
	void errmsg();
        if(argc < 2){
		errmsg();
                return(FALSE);
        }
	//default values
        p=argv;
	cmd->Nthr=2;
	cmd->dreso=16.00;
	cmd->ssize=7.0;
	cmd->th1=0.00;
	cmd->Mode=0;

        while (--argc){
         p++;
         if(**p == '-'){
          	switch(*(*p+1)){
		 case 'i':
			strcpy(cmd->filename,*(++p));
                	--argc; break;
		 case 'a':
			strcpy(cmd->file,*(++p));
                	--argc; break;
		 case 'c':
			cmd->Nthr=atoi(*(++p)); 
			--argc; break;
		 case 't':
			cmd->th1=atof(*(++p)); 
			--argc; break;
		 case 'g':
			cmd->dreso=atof(*(++p)); 
			--argc; break;
		 case 's':
                        cmd->ssize=atof(*(++p));
                        --argc; break;
		 case 'V':
                        cmd->Mode=1;
                        break;
		 default: 
		  	fprintf(stderr,"No such a option: %s\n",*p+1); 
			errmsg(); 
			return(FALSE); 
		 break;
	  }
	 }
        }
        return(TRUE);
}

void errmsg(){
	puts("Usage: Vectorize -a [MAP.mrc]  [(option)]");

	puts("---Options---");
	printf("-t [float] : Threshold of density map1 def=%.3f\n",0.00);
	printf("-g [float] : Bandwidth of the Gaussian filter\n");
        printf("             def=16.0, sigma = 0.5*[float]\n");
	printf("-s [float] : Sampling voxel spacing def=7.0\n");
	printf("-c [int  ] : Number of cores for threads def=%d\n",2);
	printf("-V 	   : Vector coordinate information def=false\n");

}

int permuted_index(int ordermode,int count, unsigned nc, unsigned nr, unsigned ns) {
        unsigned ic, ir, is;
        unsigned long ncr, q;

        ncr = nc*nr;
        is = count / ncr;
        q = count - is*ncr;
        ir = q / nc;
        ic = q - ir*nc;

        switch(ordermode) {
                case 1:
                        return ic+ir*nc+is*nc*nr;
                case 2:
                        return ic+is*nc+ir*nc*ns;
                case 3:
                        return ir+ic*nr+is*nr*nc;
                case 4:
                        return is+ic*ns+ir*ns*nc;
                case 5:
                        return ir+is*nr+ic*nr*ns;
                case 6:
                        return is+ir*ns+ic*ns*nr;
                default:
                        exit(0);
        }
}


bool readmrc(MRC *mrc, char *filename){//读取文件,调用permuted_index
 FILE *fpin;
 int ibuf,res,mode;
 int mapr,mapc,maps;
 int xdim,ydim,zdim;
 int nsymbt;
 float fbuf,*tmp_dens;

 if((fpin=fopen(filename,"rb")) == NULL){ 
  fprintf(stderr,"Can't open %s\n",filename); 
  return(true); 
 }
 
 res=fread(&(mrc->xdim),sizeof(int),1,fpin);
 res=fread(&(mrc->ydim),sizeof(int),1,fpin);
 res=fread(&(mrc->zdim),sizeof(int),1,fpin);
 xdim=mrc->xdim;
 ydim=mrc->ydim;
 zdim=mrc->zdim;
 //printf("#reading %s\n",filename);
 //printf("#XYZ dim: %d %d %d\n",mrc->xdim,mrc->ydim,mrc->zdim);

 res=fread(&(mode),sizeof(int),1,fpin);
 //ignore
 res=fread(&(mrc->ncstart),sizeof(int),1,fpin);
 //printf("#NCSTART %d\n",mrc->ncstart);
 res=fread(&(mrc->nrstart),sizeof(int),1,fpin);
 //printf("#NRSTART %d\n",mrc->nrstart);
 res=fread(&(mrc->nsstart),sizeof(int),1,fpin);
 //printf("#NSSTART %d\n",mrc->nsstart);

 //mxyz
 res=fread(&(mrc->mx),sizeof(int),1,fpin);
 res=fread(&(mrc->my),sizeof(int),1,fpin);
 res=fread(&(mrc->mz),sizeof(int),1,fpin);
 //printf("#MXYZ: %d %d %d\n",mrc->mx,mrc->my,mrc->mz);

 //xyz len
 res=fread(&(mrc->xlen),sizeof(float),1,fpin);
 res=fread(&(mrc->ylen),sizeof(float),1,fpin);
 res=fread(&(mrc->zlen),sizeof(float),1,fpin);
 //printf("#LenXYZ: %f %f %f\n",mrc->xlen,mrc->ylen,mrc->zlen);


 //abg
 res=fread(&(mrc->alpha),sizeof(float),1,fpin);
 res=fread(&(mrc->beta),sizeof(float),1,fpin);
 res=fread(&(mrc->gamma),sizeof(float),1,fpin);
 //printf("#abg: %f %f %f\n",mrc->alpha,mrc->beta,mrc->gamma);

 //map crs
 res=fread(&(mrc->mapc),sizeof(int),1,fpin);
 res=fread(&(mrc->mapr),sizeof(int),1,fpin);
 res=fread(&(mrc->maps),sizeof(int),1,fpin);
 mapc=mrc->mapc;
 mapr=mrc->mapr;
 maps=mrc->maps;
 //printf("#crs: %d %d %d\n",mrc->mapc,mrc->mapr,mrc->maps);

 res=fread(&(mrc->dmin),sizeof(float),1,fpin);
 res=fread(&(mrc->dmax),sizeof(float),1,fpin);
 res=fread(&(mrc->dmean),sizeof(float),1,fpin);
 res=fread(&(mrc->ispg),sizeof(int),1,fpin);
 //printf("#dmax,dmin, dmean, ispg: %f %f %f %d\n",mrc->dmax,mrc->dmin,mrc->dmean,mrc->ispg);

 //93-96
 res=fread(&(mrc->nsymbt),sizeof(int),1,fpin);
 nsymbt=mrc->nsymbt;

 //97-196 Extra
 fseek(fpin,4*25,SEEK_CUR);

 //197-208 ORIGIN
 res=fread(&(mrc->orgxyz),sizeof(float),3,fpin);
 //printf("#orgXYZ: %f %f %f\n",mrc->orgxyz[0],mrc->orgxyz[1],mrc->orgxyz[2]);
 
 //ignore MAP 209-212
 char text[4];
 res=fread(text,sizeof(char)*4,1,fpin);

 if(strncmp(text,"MAP",3)){
  printf("Format Error!!!\n");
  return true;
 }
 char machst[4];
 bool swap=false;
 res=fread(machst,sizeof(char)*4,1,fpin);

 if(machst[0] == 0x44){
  swap = false;
  //printf("#little-endian mode\n");
 }else if(machst[0] == 0x11){
  swap = true;
  //printf("big-endian mode\n");
 }

 if(swap){
  printf("WARNING THIS IS BIG-ENDIAN FILES!!\n");
  return true;
 }

 int ordermode=0;
 	if(mapc==1 && mapr==2 && maps==3) {
                ordermode = 1;
        }
        else if(mapc==1 && mapr==3 && maps==2) {
                ordermode = 2;
        }
        else if(mapc==2 && mapr==1 && maps==3) {
                ordermode = 3;
        }
        else if(mapc==2 && mapr==3 && maps==1) {
                ordermode = 4;
        }
        else if(mapc==3 && mapr==1 && maps==2) {
                ordermode = 5;
        }
        else if(mapc==3 && mapr==2 && maps==1) {
                ordermode = 6;
        }
        else if(ordermode == 0) {
         //printf("Input file gives malformed dimension ordering.");
         return true;
        }
  //printf("#Order Mode= %d\n",ordermode);
  mrc->NumVoxels = mrc->xdim*mrc->ydim*mrc->zdim;
 //printf("#Nvoxels= %d\n",mrc->NumVoxels);

 if((mrc->dens=(float*)malloc(sizeof(float)*mrc->NumVoxels))==NULL)
  return true;

 //fin.ignore(4*(256-54)+nsymbt);
 fseek(fpin,4*(256-54)+nsymbt,SEEK_CUR);
 
 	switch(mode) {
 	 case 0: // char - converted to float, testing for signed-ness
	   printf("Cannot read mode 0 mrc file\n");
           return true;
	   //break;
         case 1: // 16-bit float
	  //printf("#Reading 16-bit mrc file\n");
          for(int i = 0; i<mrc->NumVoxels; ++i)
	   res=fread(&(mrc->dens[permuted_index(ordermode,i,xdim,ydim,zdim)]),2,1,fpin);
           break;
         case 2: // 32-bit float
	  //printf("#Reading 32-bit mrc file\n");
          for(int i = 0; i<mrc->NumVoxels; ++i){
	   res=fread(&(mrc->dens[permuted_index(ordermode,i,xdim,ydim,zdim)]),sizeof(float),1,fpin);
	  }
          break;
         default:
          printf("Unknown floating-point mode specified.");
	  return true;
        }

 	mrc->widthx = mrc->xlen / (double) mrc->mx;
        mrc->widthy = mrc->ylen / (double) mrc->my;
        mrc->widthz = mrc->zlen / (double) mrc->mz;

	if(fabs(mrc->widthx - mrc->widthy)>0.000001 || 
	 fabs(mrc->widthx -  mrc->widthz)>0.000001 ||
	 fabs(mrc->widthy -  mrc->widthz)>0.000001){

	 printf("#ERROR: grid sizes are different %f %f %f\n",
	  mrc->widthx,mrc->widthy,mrc->widthz);
	 printf("PLEASE USE CUBIC MRC MAP DATA\n");
	 return true;
	}

	int nx,ny,nz;
	switch(ordermode){
		case 1:
			nx = xdim; ny = ydim; nz = zdim;
			break;
		case 2:
			nx = xdim; ny = zdim; nz = ydim;
			break;
		case 3:
			nx = ydim; ny = xdim; nz = zdim;
			break;
		case 4:
			nx = zdim; ny = xdim; nz = ydim;
			break;
		case 5:
			nx = ydim; ny = zdim; nz = xdim;
			break;
		case 6:
			nx = zdim; ny = ydim; nz = xdim;
			break;
		default:
			printf("Input file gives malformed dimension ordering.");
			return true;
	}

 mrc->xdim = nx;
 mrc->ydim = ny;
 mrc->zdim = nz;

 //printf("#XYZ dim: %d %d %d\n",mrc->xdim,mrc->ydim,mrc->zdim);
 //printf("#XYZ dim: %d \n",mrc->xdim*mrc->ydim*mrc->zdim);

 fclose(fpin);



 if(mrc->ncstart!=0||mrc->nrstart!=0||mrc->nsstart!=0){
  mrc->orgxyz[0]=mrc->orgxyz[0]+mrc->ncstart*mrc->widthx;
  mrc->orgxyz[1]=mrc->orgxyz[1]+mrc->nrstart*mrc->widthy;
  mrc->orgxyz[2]=mrc->orgxyz[2]+mrc->nsstart*mrc->widthz;
 }



 return false;
}


bool fastVEC(MRC *m,MRC *M){
 int i,ind;
 int xydim=m->xdim*m->ydim;
 int Ndata=M->xdim*M->ydim*M->zdim;
 //malloc

 if((M->vec=(double **)malloc(sizeof(double *)*Ndata))==NULL)
  return true;
 if((M->dens=(float *)malloc(sizeof(double)*Ndata))==NULL)
  return true;
 for(i=0;i<Ndata;i++){
  if((M->vec[i]=(double *)malloc(sizeof(double)*3))==NULL)
  return true;
 }

 if((M->xyz=(int **)malloc(sizeof(int *)*Ndata))==NULL)
  return true;
 for(i=0;i<Ndata;i++){
  if((M->xyz[i]=(int *)malloc(sizeof(int)*3))==NULL)
  return true;
 }


 //puts("#Start VEC");
 //Setup Filter
 //Gaussian kernel dreso=window size
 double dreso=cmd.dreso;//def=16
 double gstep=m->widthx;//def=7.0
 double fs=(dreso/gstep)*0.5;
 fs=fs*fs;
 double fsiv=1.000/fs;
 double fmaxd=(dreso/gstep)*2.0;
//  printf("#maxd= %f\n",fmaxd);

 double dsum=0;
 int Nact=0;

 #pragma omp parallel for reduction(+:Nact) schedule(dynamic,5)
 for(int x=0;x<M->xdim;x++){
  double rx,ry,rz,d2;
 for(int y=0;y<M->ydim;y++){
 for(int z=0;z<M->zdim;z++){

  int stp[3],endp[3],ind2,ind;
  double pos[3],pos2[3],ori[3];
  double tmpcd[3];
  double v,dtotal,rd;


  //real cd -> m-cd
  pos[0]=(x*M->widthx+M->orgxyz[0]-m->orgxyz[0])/m->widthx;
  pos[1]=(y*M->widthx+M->orgxyz[1]-m->orgxyz[1])/m->widthx;
  pos[2]=(z*M->widthx+M->orgxyz[2]-m->orgxyz[2])/m->widthx;

  ind=M->xdim*M->ydim*z+M->xdim*y+x;
  //printf("%d %d %d %d\n",x,y,z,ind);
  //printf("%f %f %f %d\n",pos[0],pos[1],pos[2],ind);

  //check density
  if(pos[0]<0||pos[1]<0||pos[2]<0||
     pos[0]>=m->xdim||pos[1]>=m->ydim||pos[2]>=m->zdim){
   M->dens[ind]=0;
   M->vec[ind][0]=M->vec[ind][1]=M->vec[ind][2]=0.00;
   continue;
  }





  int ind0=m->xdim*m->ydim*(int)pos[2]+m->xdim*(int)pos[1]+(int)pos[0];
  if(m->dens[ind0]==0){
   M->dens[ind]=0;
   M->vec[ind][0]=M->vec[ind][1]=M->vec[ind][2]=0.00;
   continue;
  }

  ori[0]=pos[0];
  ori[1]=pos[1];
  ori[2]=pos[2];
   //Start Point
   stp[0]=(int)(pos[0]-fmaxd);
   stp[1]=(int)(pos[1]-fmaxd);
   stp[2]=(int)(pos[2]-fmaxd);

   if(stp[0]<0)stp[0]=0;
   if(stp[1]<0)stp[1]=0;
   if(stp[2]<0)stp[2]=0;


   endp[0]=(int)(pos[0]+fmaxd+1);
   endp[1]=(int)(pos[1]+fmaxd+1);
   endp[2]=(int)(pos[2]+fmaxd+1);

   if(endp[0]>=m->xdim) endp[0]=m->xdim;
   if(endp[1]>=m->ydim) endp[1]=m->ydim;
   if(endp[2]>=m->zdim) endp[2]=m->zdim;

   dtotal=0;
   pos2[0]=pos2[1]=pos2[2]=0;
   for(int xp=stp[0];xp<endp[0];xp++){
    rx=(double)xp-pos[0];
    rx=rx*rx;
   for(int yp=stp[1];yp<endp[1];yp++){
    ry=(double)yp-pos[1];
    ry=ry*ry;
   for(int zp=stp[2];zp<endp[2];zp++){
    rz=(double)zp-pos[2];
    rz=rz*rz;
    d2=rx+ry+rz;
    //d=exp(-1.50*fr*fsiv)*amap(ii,jj,kk)
    ind2=xydim*zp+m->xdim*yp+xp;
    v=exp(-1.50*d2*fsiv)*m->dens[ind2];
    dtotal+=v;
    //if(v>0)
    //printf("d %f %d %d %d\n",v,xp,yp,zp);
    pos2[0]+=v*(double)xp;
    pos2[1]+=v*(double)yp;
    pos2[2]+=v*(double)zp;
   }}}
   M->dens[ind]=dtotal;
   if(dtotal==0.00){
    M->vec[ind][0]=M->vec[ind][1]=M->vec[ind][2]=0.00;
    continue;
   }
   rd=1.00/dtotal;
   pos2[0]*=rd;
   pos2[1]*=rd;
   pos2[2]*=rd;
   tmpcd[0]=pos2[0]-pos[0];
   tmpcd[1]=pos2[1]-pos[1];
   tmpcd[2]=pos2[2]-pos[2];

   double dvec=sqrt(tmpcd[0]*tmpcd[0]+tmpcd[1]*tmpcd[1]+tmpcd[2]*tmpcd[2]);
   if(dvec==0.00) dvec=1.00;//For low resolution maps!!
   double rdvec=1.000/dvec;
   M->vec[ind][0]=tmpcd[0]*rdvec;
   M->vec[ind][1]=tmpcd[1]*rdvec;
   M->vec[ind][2]=tmpcd[2]*rdvec;
   M->xyz[ind][0]=x;
   M->xyz[ind][1]=y;
   M->xyz[ind][2]=z;
   //printf("%d %d %d \n",x,y,z);
   //printf("%d %d %d %d\n",x,y,z,ind);
   //printf("%d %d %d\n",M->xyz[ind][0],M->xyz[ind][1],M->xyz[ind][2]);
   //printf("%f %f %f %f\n",M->vec[ind][0],M->vec[ind][1],M->vec[ind][2],M->dens[ind]);
   //printf("%f %f %f\n",M->vec[ind][0],M->vec[ind][1],M->vec[ind][2]);

   dsum+=dtotal;
   Nact++;

 }}}
 //puts("#End LDP");

  
 //Add Average & STD
 M->dsum=dsum;
 M->Nact=Nact;
 M->ave=dsum/(double)Nact;
 dsum=0;
 double dsum2=0;
 for(int i=0;i<M->xdim*M->ydim*M->zdim;i++)
  if(M->dens[i]>0){//Cross correlation
   dsum+=(M->dens[i])*(M->dens[i]);
   dsum2+=(M->dens[i]-M->ave)*(M->dens[i]-M->ave);
  }
 M->std_norm_ave=sqrt(dsum2);
 M->std=sqrt(dsum);
 printf("%f %f %f\n",M->ave,M->std,M->std_norm_ave);
 
 for(int x=0;x<M->xdim;x++){
 for(int y=0;y<M->ydim;y++){
 for(int z=0;z<M->zdim;z++){    
    ind=M->xdim*M->ydim*z+M->xdim*y+x;
    if(M->xyz[ind][0]!=0){
        printf("%d %d %d %d\n",ind,M->xyz[ind][0],M->xyz[ind][1],M->xyz[ind][2]);
        printf("%f %f %f %f\n",M->vec[ind][0],M->vec[ind][1],M->vec[ind][2],M->dens[ind]);
    }   
 
 }}}

 

 return false;
}



void SetUpVoxSize(MRC *m,MRC *M,double t,double ssize){
 int Ndata=m->xdim*m->ydim*m->zdim;
 int xydim=m->xdim*m->ydim;
 int xdim=m->xdim;

 double cent[3];
 double d2,dmax=0;

 cent[0]=m->xdim*0.5;
 cent[1]=m->ydim*0.5;
 cent[2]=m->zdim*0.5;

  //If t<0
 if(t<0){
        for(int x=0;x<m->xdim;x++){
        for(int y=0;y<m->ydim;y++){
        for(int z=0;z<m->zdim;z++){
         int ind=xydim*z+xdim*y+x;
         m->dens[ind]-=t;
        }}}
        t=0.00;
 }
 //---------------

 for(int x=0;x<m->xdim;x++){
 for(int y=0;y<m->ydim;y++){
 for(int z=0;z<m->zdim;z++){
  int ind=xydim*z+xdim*y+x;
  if(m->dens[ind]<t){
   m->dens[ind]=0.00;//filter
   continue;
  }

  d2=((double)x-cent[0])*((double)x-cent[0])
    +((double)y-cent[1])*((double)y-cent[1])
    +((double)z-cent[2])*((double)z-cent[2]);
  if(d2>dmax)
   dmax=d2;
 }}}
 //printf("#dmax= %f size=%d\n",sqrt(dmax)/m->widthx,(int)(2*sqrt(dmax)/m->widthx));

 //real coordinates 实坐标
 M->cent[0]=cent[0]*m->widthx+m->orgxyz[0];
 M->cent[1]=cent[1]*m->widthx+m->orgxyz[1];
 M->cent[2]=cent[2]*m->widthx+m->orgxyz[2];

 M->widthx=ssize;
 //printf("Widthx= %f\n",M->widthx);
 printf("%f\n",M->widthx);

 m->dmax=sqrt(dmax)*m->widthx;
 int tmp_size=2*sqrt(dmax)*m->widthx/M->widthx;

 int a=2;
 while(1){
  if(a>tmp_size)
   break;
  a*=2;
 }
 int b=3;
 while(1){
  if(b>tmp_size)
   break;
  b*=2;
 }
 if(b<a)
  a=b;
 b=9;
 while(1){
  if(b>tmp_size)
   break;
  b*=2;
 }
 if(b<a)
  a=b;

 M->xdim=M->ydim=M->zdim=a;
 M->orgxyz[0]=M->cent[0]-0.5*a*M->widthx;
 M->orgxyz[1]=M->cent[1]-0.5*a*M->widthx;
 M->orgxyz[2]=M->cent[2]-0.5*a*M->widthx;
 //printf("Nvox= %d*%d*%d\n",M->xdim,M->ydim,M->zdim);
 //printf("cent= %f %f %f\n",M->cent[0],M->cent[1],M->cent[2]);
 //printf("ori= %f %f %f\n",M->orgxyz[0],M->orgxyz[1],M->orgxyz[2]);

 printf("%d %d %d\n",M->xdim,M->ydim,M->zdim);
 printf("%f %f %f\n",M->cent[0],M->cent[1],M->cent[2]);
 printf("%f %f %f\n",M->orgxyz[0],M->orgxyz[1],M->orgxyz[2]);
}


void ShowVec(MRC *M){

 int i,ind;
 double tmp[3];
 for(int x=0;x<M->xdim;x++){
 for(int y=0;y<M->ydim;y++){
 for(int z=0;z<M->zdim;z++){
  ind=M->xdim*M->ydim*z+M->xdim*y+x;

  i=ind;
  if(M->dens[i]==0.00&&ind!=0) continue;
  tmp[0]=x*M->widthx+M->orgxyz[0];
  tmp[1]=y*M->widthx+M->orgxyz[1];
  tmp[2]=z*M->widthx+M->orgxyz[2];

  //printf("H %8.3f%8.3f%8.3f\n",tmp[0],tmp[1],tmp[2]);
  printf("H       %f        %f        %f\n",tmp[0],tmp[1],tmp[2]);
  


  tmp[0]=(x+M->vec[i][0])*M->widthx+M->orgxyz[0];
  tmp[1]=(y+M->vec[i][1])*M->widthx+M->orgxyz[1];
  tmp[2]=(z+M->vec[i][2])*M->widthx+M->orgxyz[2];
  
  //printf("H %8.3f%8.3f%8.3f\n",tmp[0],tmp[1],tmp[2]);
  printf("H       %f        %f        %f\n",tmp[0],tmp[1],tmp[2]);
 
  
 }}}
 
}








