/***********************************************************/
/*                                                         */
/*               Executor of SP codes                      */
/*               (for double bucket Dijkstra)              */
/*                                                         */
/***********************************************************/


#include <stdio.h>
#include <stdlib.h>

#define YES 1
#define NO  0
#define STATISTIC YES

/* statistical variables */
long n_scans = 0;
long n_impr = 0;

/* definitions of types: node & arc */

#include "types_db.h"

/* parser for getting extended DIMACS format input and transforming the
   data to the internal representation */

#include "parser_db.c"

/* function 'timer()' for measuring processor time */

//#include "timer.c"

/* function for constructing shortest path tree */

#include "dikbd.c"


//main ()
int dikbd_run(char *spfile, int *distance, int *path)

{

float t;
arc *arp, *ta;
node *ndp, *source, *k;
long n, m, nmin, i; 
char name[21];
long mlen;
double sum_d = 0;
int cc;
FILE *fpout;
int first = 1;
node *nd;

 parse( &n, &m, &ndp, &arp, &source, &nmin, name, &mlen, spfile );


/*
printf ( "%s\nn= %ld, m= %ld, nmin= %ld, source = %ld, maxlen= %d\n",
        name,
        n,m,nmin, (source-ndp)+nmin, mlen );

printf ("\nordered arcs:\n");
for ( k = ndp; k< ndp + n; k++ )
  { i = (k-ndp)+nmin;
    for ( ta=k -> first; ta != (k+1)-> first; ta++ )
      printf ( " %2ld %2ld %4ld\n",
               i, ((ta->head)-ndp)+nmin, ta->len
             );

  }
*/

//t = timer();

cc = dikbd ( n, ndp, source, mlen );

//t = timer() - t;
t = 0;

i = 1;
for ( k= ndp; k< ndp + n; k++ ) {
    if ( k -> parent != (node*) NULL ) {
        sum_d += (double) (k -> dist);
		/*
		if (first) {
			first = 0;
			nd = k;
			while (nd->parent != (node*) NULL) {
				printf("i, status: %d %d %p\n",i,nd->status,nd);
				nd = nd->parent;
				i++;
			}
			printf("# of nodes on the path: %d\n",i);
		}
		*/
//   printf("To node: %d distance: %d\n",i, k->dist);
//   i++;
	}
}
printf ("\nDijkstra with double buckets -> problem name: %s\n\n\
Nodes: %ld    Arcs: %ld  cc: %d\n\
Number of scans: %ld\n Number of improvements: %ld\n\
Sum of distances: %8.0f\n\n", name, n, m, cc, n_scans, n_impr, sum_d);

//printf("problem name: %s Nodes: %ld Arcs: %ld  cc: %d scans: %ld improvements: %ld Sum: %8.0f\n\n", name, n, m, cc, n_scans, n_impr, sum_d); 

#define nd(ptr) (int)(ptr-ndp+nmin)

fpout = fopen("dikbd.out","w");
for ( k=ndp; k< ndp+n; k++ ) {
	fprintf (fpout," %d %d %d\n", nd(k), nd(k->parent), k->dist);
	path[nd(k)-1] = nd(k->parent) - 1;
	distance[nd(k)] = k->dist;
}
fclose(fpout);
return 0;
}
