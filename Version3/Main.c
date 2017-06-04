/* Main function
   Modified: 2017-06-03
   Fix:(BUG)92: depth = argv[3]
*/

#include "Share.h"

void Usage( void )
{
    char *usage =
        "\nUsage: celldiv <comand> [options]\n"
        "\n"
        "Commands:\n"
        "   normal              normal division\n"
        "       <ref>           reference for the cell type[.fa]\n"
        "       <out>           outfile name\n"
        "       <depth>         equal to sequencing depth\n"
        "   maline              dividing with mutation accumulation line\n"
        "       <ref>           reference for the cell type[.fa]\n"
        "       <out>           outfile name\n"
        "       <prolife>       proliferation num\n"
        "       <depth>         equal to sequencing depth\n\n";

    fprintf(stderr, "%s", usage); exit(-1);
}

int Divide_Normal( int argc, char **argv )
{
    mut_t *snv; time_t s, e;

    time(&s); srand((unsigned)time(NULL)); /* random the seed */
    tree_t *tree = Init(snv, 0);
    for ( int i=0; i < GENERATION-1; i++ ) {
        printf("\rStart dividing the %d generation ...", i+1); fflush(stdout);
        Divide(tree, i);
    } printf("\n");

    uint8_t *fa = LoadRef(argv[0]); /* Load the reference */

    printf("Start the pileup ...\n");
    pile_t *P = PileInit(fa, atoi(argv[2])); 
    Pileup(P, tree, atoi(argv[2]));

    SnvWrite(argv[1], P); time(&e);
    printf("Total time consume is %.2f\n\n", difftime(e,s));

    return 0;
}

int Divide_Maline( int argc, char **argv )
{
    mut_t snv; tree_t *tree;
    int pronum, *clist;
    time_t s, e;

    snv.altnum = 0; pronum = atoi(argv[2]);
    snv.pos = (uint32_t *)malloc(GENERATION * sizeof(uint32_t));
    snv.type = (uint8_t *)malloc(GENERATION * sizeof(uint8_t));
    if ( !snv.pos || !snv.type ) {
        fprintf(stderr, \
            "[Err::%s::%d] Failed to allocate memory!\n", __func__, __LINE__);
        return (-1);
    } time(&s);
    for ( int p=0; p < pronum; p++ ) {
        srand((unsigned)time(NULL)); /* random the seed */
        tree = Init(&snv, 1);
        printf("\rStart proliferaing for the %d times ...", p+1); fflush(stdout);
        for ( int i=0; i < GENERATION-1; i++ ) Divide(tree, i);
        clist = (int *)malloc(tree->typenum * sizeof(int));
        if ( !clist ) {
            fprintf(stderr, \
                "[Err::%s::%d] Failed to allocate memory!\n", __func__, __LINE__);
            return (-1);
        }
        clist[0] = tree->root[0].cellnum[GENERATION-1];
        for ( int i=1; i < tree->typenum; i++ )
            clist[i] = clist[i-1] + tree->root[i].cellnum[GENERATION-1];

        int pos = rand() % ( 1 << (GENERATION-1) ); 
        int index = LocIndex(clist, tree->typenum, pos);
        mut_t *tmp = tree->root[index].snv;

        snv.altnum = tmp->altnum;
        memcpy(snv.pos, tmp->pos, tmp->altnum * sizeof(uint32_t));
        memcpy(snv.type, tmp->type, tmp->altnum * sizeof(uint8_t));
        if ( p != pronum -1 ) Destory(tree); 
        free(clist);
    } printf("\n");

    uint8_t *fa = LoadRef(argv[0]); /* Load the reference */
    printf("Start the pileup ...\n");
    pile_t *P = PileInit(fa, atoi(argv[3])); Pileup(P, tree, atoi(argv[3]));
    SnvWrite(argv[1], P); time(&e);
    printf("Total time consume is %.2f\n\n", difftime(e,s));

    return 0;
}


int main( int argc, char **argv )
{
    if ( argc < 2 ) Usage();

    if ( strcmp(argv[1], "normal") == 0 ) {
        if ( argc != 5 )  Usage();
        else
           return Divide_Normal(argc-2, argv+2);
    }
    else if ( strcmp(argv[1], "maline") == 0 ) {
        if ( argc != 6 ) Usage();
        else
            return Divide_Maline(argc-2, argv+2);
    }
}
