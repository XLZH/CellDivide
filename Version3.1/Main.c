/* Main function
   Modified: 2017-06-03
   Fix:(BUG)92: depth = argv[3]
   Fix:(BUG)71&72: GENERATION => GSIZE
*/

#include "Share.h"


static int Divide_Normal( arg_t *Arg )
{
    mut_t *snv; time_t s, e;

    time(&s); srand((unsigned)time(NULL)); /* random the seed */
    tree_t *tree = Init(snv, 0);
    for ( int i=0; i < Arg->generation; i++ ) {
        printf("\rStart dividing the %d generation ...", i+1); fflush(stdout);
        Divide(tree, i, Arg);
    } printf("\n");

    uint8_t *fa = LoadRef(Arg); /* Load the reference */

    printf("Start the pileup ...\n");
    pile_t *P = PileInit(fa, Arg); 
    Pileup(P, tree, Arg);

    SnvWrite(P, Arg); time(&e);
    printf("Total time consume is %.2f\n\n", difftime(e,s));

    return 0;
}

static int Divide_Maline( arg_t *Arg )
{
    mut_t snv; tree_t *tree;
    int pronum, *clist;
    time_t s, e;

    snv.altnum = 0;
    snv.pos = (uint32_t *)malloc(Arg->gsize * sizeof(uint32_t));
    snv.type = (uint8_t *)malloc(Arg->gsize * sizeof(uint8_t));
    if ( !snv.pos || !snv.type ) {
        fprintf(stderr, \
            "[Err::%s::%d] Failed to allocate memory!\n", __func__, __LINE__);
        return (-1);
    } time(&s);

    srand((unsigned)time(NULL)); /* random the seed */
    for ( int p=0; p < Arg->transnum; p++ ) {
        tree = Init(&snv, 1);
        printf("\rStart proliferaing for the %d times ...", p+1); fflush(stdout);
        for ( int i=0; i < Arg->generation; i++ ) Divide(tree, i, Arg);
        clist = (int *)malloc(tree->typenum * sizeof(int));
        if ( !clist ) {
            fprintf(stderr, \
                "[Err::%s::%d] Failed to allocate memory!\n", __func__, __LINE__);
            return (-1);
        }
        clist[0] = tree->root[0].cellnum[Arg->generation];
        for ( int i=1; i < tree->typenum; i++ )
            clist[i] = clist[i-1] + tree->root[i].cellnum[Arg->generation];

        uint64_t pos = RAND64() % (1 << Arg->generation); 
        int index = LocIndex(clist, tree->typenum, pos);
        mut_t *tmp = tree->root[index].snv;

        snv.altnum = tmp->altnum;
        memcpy(snv.pos, tmp->pos, tmp->altnum * sizeof(uint32_t));
        memcpy(snv.type, tmp->type, tmp->altnum * sizeof(uint8_t));
        if ( p != Arg->transnum -1 ) Destory(tree); 
        free(clist);
    } printf("\n");

    uint8_t *fa = LoadRef(Arg); /* Load the reference */
    printf("Start the pileup ...\n");
    pile_t *P = PileInit(fa, Arg); Pileup(P, tree, Arg);
    SnvWrite(P, Arg); time(&e);
    printf("Total time consume is %.2f\n\n", difftime(e,s));

    return 0;
}


int main( int argc, char **argv )
{
    if ( argc < 2 ) Usage(0);

    if ( strcmp(argv[1], "normal") == 0 ) {
        if ( argc != 12 )  Usage(NORMAL);
        else {
            arg_t *Arg = ParseOpt(argc-1, argv+1);
            if ( !Arg->help ) {
                fprintf(stderr, "GenomeSize: %d\nGeneration:%d\nMutRate:%lu\n",
                            Arg->gsize, Arg->generation, Arg->rate);
                return Divide_Normal(Arg);
            }
            else Usage(NORMAL);
        }
    }
    else if ( strcmp(argv[1], "maline") == 0 ) {
        if ( argc != 14 ) Usage(MALINE);
        else {
            arg_t *Arg = ParseOpt(argc-1, argv+1);
            if ( !Arg->help ) {
                fprintf(stderr, "GenomeSize: %d\nGeneration:%d\nMutRate:%lu\n",
                            Arg->gsize, Arg->generation, Arg->rate);
                return Divide_Maline(Arg);
            }
            else Usage(MALINE);
        }
    }
    else Usage(0);
}
