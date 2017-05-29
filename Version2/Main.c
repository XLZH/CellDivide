/* Version 2
 * Date: 2017-05-28
*/

#include "Share.h"

void Usage( int argc )
{
    char *usage = "\n Usage: celldiv <ref.fa> <outfile> <depth>\n";

    if ( argc != 4 ) {
        fprintf(stderr,"%s", usage); exit(-1);
    }
}


int main( int argc, char **argv )
{
    time_t start, end; clock_t s, e;
    tree_t *tree = Inital();

    Usage(argc); time(&start);
    srand((unsigned)time(NULL)); /* random the seed */
    s = clock(); printf("Start dividing ...\n");
    for ( int i=0; i < GENERATION-1; i++ ) Divide(tree, i);
    printf("Time consum: %.2f(s)\n", (double)(clock()-s)/CLOCKS_PER_SEC);

    printf("Start to load the reference ...\n");
    uint8_t *fa = LoadRef(argv[1]); 

    s = clock(); printf("Start the pileup ...\n");
    pile_t *P = PileInit(fa, atoi(argv[3])); Pileup(P, tree, atoi(argv[3]));
    printf("Time consume: %.2f(s)\n", (double)(clock()-s)/CLOCKS_PER_SEC);

    printf("Write the snv info out ...\n");
    SnvWrite(argv[2], P);

    time(&end); printf("Total time consume is %.2f\n", difftime(end,start));
}

