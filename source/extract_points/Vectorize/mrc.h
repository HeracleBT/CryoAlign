#define LIN 256
#define TRUE 0
#define FALSE -1

typedef struct{
	char filename[LIN];
	char file[LIN];
	int Nthr;
	double dreso;
	int Mode;
	double th1;
	double ssize;
} CMD;

typedef struct{
	char filename[LIN];
	int xdim,ydim,zdim;//mrc information
	int ncstart,nrstart,nsstart;//
	int mx,my,mz;//
	float xlen,ylen,zlen;//
	float alpha,beta,gamma;//
	int mapc,mapr,maps;//
	float dmin,dmax,dmean;//
	int ispg;//
	int nsymbt;//
	float orgxyz[3];
	int NumVoxels;//
	float *dens;
	double **vec;
	int **xyz;
	float widthx,widthy,widthz;//len/m
	unsigned int Nact;//中间数据
	double dmax2,dsum,std,ave;
	double cent[3];
	double std_norm_ave;
} MRC;

int chkcmdline(int, char **,CMD *);
bool readmrc(MRC *,char *);
bool fastVEC(MRC *,MRC *);
void SetUpVoxSize(MRC *,MRC *,double,double);
void ShowVec(MRC *);
//	void ShowMrc(MRC *)

	
