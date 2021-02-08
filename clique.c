/*============================================================================
| (c) Copyright Arthur L. Corcoran, 1992, 1993.  All rights reserved.
|
| Genetic Algorithm Test Program 
============================================================================*/
#include "ga.h"
GA_run();
/* Global Variables*/
char graph[500][500];
int nnodes, nedges;
int clique=0;
//int largo;

/* function prototypes */
int obj_fun(Chrom_Ptr);
int read_instance();

int obj_fun();    /*--- Forward declaration ---*/

/*----------------------------------------------------------------------------
| main()
----------------------------------------------------------------------------*/
int main(int argc, char **argv) 
{
   if(argc > 0) {
        printf("%s\n", argv[1]);
        if(argc == 1){
            argv[1]="c1.txt";
        }
 
	GA_Info_Ptr ga_info;
	int i;

	/*--- Initialize the genetic algorithm ---*/
	ga_info = GA_config(argv[1],obj_fun);

   
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	read_instance("san200_0.7_1.clq.txt");
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   
	ga_info->chrom_len = nnodes;
	//largo = ga_info->chrom_len

	/*--- Run the GA ---*/
	GA_run(ga_info);

	printf("\nBest chrom:  ");
	for(i=0;i<ga_info->chrom_len;i++){
	  printf("%d",(int)ga_info->best->gene[i]);
	  clique +=(int)ga_info->best->gene[i];
	}
	printf("   (fitness: %g)\n\n",ga_info->best->fitness);
	printf("------------------------------------------------------------- TAMAÑO CLIQUE: %d\n",clique);
	}
}

/*----------------------------------------------------------------------------
| obj_fun() - user specified objective function
----------------------------------------------------------------------------*/
int obj_fun(Chrom_Ptr chrom) 
{
  int i,j; 
  double ones = 0.0;
  double edges = 0.0;
  
  for(i = 0; i < nnodes; i++){        
    if(chrom->gene[i] == 1){
      ones++; //cuento numero de unos o nodos
      for(j = i+1; j < nnodes; j++){
	if(chrom->gene[j] == 1){
	  if(graph[i][j] == 1){
	    edges++;} //cuento numero de aristas 
	  }
	}
      }
    }

  
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  chrom->fitness = ((ones*(ones-1)/2)/edges)*abs((ones*(ones-1)/2)-edges)-ones;
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //chrom->fitness = ones*abs((ones*(ones-1)/2)-edges)-ones;
  //chrom->fitness = ones-abs((ones*(ones-1)/2)-edges);
  //chrom->fitness = ones-((ones*(ones-1)/2)/edges)*abs((ones*(ones-1)/2)-edges);

  return 0;
  
}




// read DIMACS format
int read_instance(char *filename)
{
  char dummy1;
  char dummy2[100];
  int dummy3;
  int n1,n2;
  FILE *inputf;
  int i,j;
  

  nnodes = 0;
  nedges = 0;
  
  if( (inputf=fopen(filename,"rt")) == NULL )
    {
      printf( "Cannot open file %s\n",filename);
      exit(-1);
    }

  // lee la cabecera
  fscanf(inputf,"%c %s %d %d\n",&dummy1,dummy2,&nnodes,&nedges);
  
  printf("Opening %s (%d nodes, %d edges)\n",filename,nnodes,nedges);
   
    for(i=0;i<nnodes;i++){
        for(j=0;j<nnodes;j++){
            graph[i][j] = 0;
        }
    }
    

  // salta la lista de nodos
    for(i=0;i<nnodes;i++){
        fscanf(inputf,"%c  %d %d\n",&dummy1,&dummy3,&dummy3);
    }

  // lee los edges
  for(i=0;i<nedges;i++)
    {
      fscanf(inputf,"%c %d %d\n",&dummy1,&n1,&n2);
      graph[n1-1][n2-1] = 1;  // ojo que los vectores en C empiezan desde 0
      graph[n2-1][n1-1] = 1;
    }
  
    
  fclose(inputf);
  
}
