#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <omp.h>
#include "mrc.h"

double gettimeofday_sec()
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + (double)tv.tv_usec*1e-6;
}

CMD cmd;

int main(int argc, char **argv)
{
 double t1=gettimeofday_sec();
 double t4;
 MRC mrc,mrc1;
 //Get cmd
 if(chkcmdline(argc,argv,&cmd)==FALSE)
  return(0);

 //Set threads
 if(cmd.Nthr < omp_get_num_procs()){
  omp_set_num_threads(cmd.Nthr);
 }else{
  omp_set_num_threads(omp_get_num_procs());
 }
 
 if(readmrc(&mrc1,cmd.file))
  return(0);

 //Setup 2^x size 

 MRC mrcN1;

 SetUpVoxSize(&mrc1,&mrcN1,cmd.th1,cmd.ssize);


 //Convert to vector
 if(fastVEC(&mrc1,&mrcN1))
  return(0);
 
 if(cmd.Mode==1){
   ShowVec(&mrcN1);
   }
  /*ShowVec(&mrcN1);
 t4=gettimeofday_sec();*/
 //printf("#FINISHED TOTAL TIME= %f\n",t4-t1);
 return 0;

}

