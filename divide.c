/* Version: 3
 * New: remove the "search" function
 * date: 2017-5-13
 */

#include <stdio.h>
#include <time.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>

#define GENERATION 30          // the generation that cultured the cell
#define GSIZE 4746218          // the genome size of Ecoli
#define FREQLEVEL 1000000000   // the mutation frequence level for whole genome
#define SIZE 100000            // default size for memory realocated each time

typedef struct __mut_t {
    uint32_t    altnum;  // snv num for this cell type
    uint32_t    *pos;    // mutation position
    uint8_t     *type;   // mutation type [ A|T|C|G ]
} mut_t;

typedef struct __cell_t {
    uint32_t    cellnum[GENERATION]; // cellnum for each generation
    mut_t       *snv;                // mutation information
} cell_t;

typedef struct __tree_t {
    uint32_t    typenum;    // mutation type count
    cell_t      *root;      // pointer to the mutation type array
} tree_t;

typedef struct __pile_t {
    uint8_t     base[4];    // { refbase, altbase, [A|T|C|G], [A|T|C|G] }
    uint32_t    num[4];     // count for the base
} pile_t;


void InsertAlt( mut_t *new, const mut_t *f, uint32_t pos, uint8_t type )
{
    memcpy(new->pos, f->pos, f->altnum * sizeof(uint32_t));
    memcpy(new->type, f->type, f->altnum * sizeof(uint8_t));

    /* test whther the pos has existed in f.pos */
    for ( int i=0; i < f->altnum; i++ ) {
        if ( pos == new->pos[i] ) {
            new->type[i] = type; new->altnum = f->altnum; 
            return ;
        }
    }
    /* insert the pos and type in order */
    for ( int i=f->altnum -1; i >= -1; i-- ) { 
        if ( i != -1 && pos < new->pos[i] ) {
            new->pos[i+1] = new->pos[i]; new->type[i+1] = new->type[i];
        }
        else {
            new->pos[i+1] = pos; new->type[i+1] = type; new->altnum = f->altnum +1; 
            return ;
        }
    }
}

mut_t *Mutator( const mut_t *f )
{
    uint8_t type, mtype[4] = {65, 67, 71, 84};
    mut_t *snv; uint32_t pos = 0;

    pos = rand() % FREQLEVEL; type = mtype[rand() % 4];
    if ( pos <= GSIZE ) {
        snv = (mut_t *)malloc(sizeof(mut_t));
        snv->pos = (uint32_t *)malloc((f->altnum +1) * sizeof(uint32_t));
        snv->type = (uint8_t *)malloc((f->altnum +1) * sizeof(uint8_t));
        if ( !snv || !(snv->pos) || !(snv->type) ) {
            fprintf(stderr, \
                    "[Err::%s::%d] Failed to allocate memory\n", __func__, __LINE__);
            return NULL;
        }
        InsertAlt(snv, f, pos, type); return snv;
    }
    else return NULL;
}

void CreatCell( tree_t *tree, mut_t *new, int g )
{
    cell_t *start, *current;

    if ( !(tree->typenum % SIZE) ) {
        int newsize = tree->typenum + SIZE;
        start = (cell_t *)realloc(tree->root, newsize*sizeof(cell_t));
        if ( !start ) {
            fprintf(stderr, \
                    "[Err::%s::%d] Failed to realloc memory!\n", __func__, __LINE__);
            return ;
        }
        tree->root = start;
    }
    current = &tree->root[tree->typenum];
    memset(current, 0, sizeof(cell_t)); 
    current->snv = new; current->cellnum[g] = 1; tree->typenum++;

    return ;
}

void Divide( tree_t *tree, int g )
{
    cell_t *c, *pre; mut_t *tem; 
    int typenum = tree->typenum;

    for ( int i=0; i < typenum; i++ ) /* range from different mutation type */
    {
        c = &(tree->root[i]);
        for ( int j=0; j < c->cellnum[g]; j++ ) {   
            /* range form all the sametype cells */
            for ( int k=0; k < 2; k++ ) { /* the cell should divided by twice */
                tem = Mutator(c->snv);
                if ( tem ) { /* mutation occored */
                    CreatCell(tree, tem, g+1); c = &(tree->root[i]);
                }
                else  /* cell divided normally */
                    c->cellnum[g+1]++;
            }
        }
    }
}

tree_t *Inital( void )
{
    tree_t *tree; cell_t *father; mut_t *germ;

    tree = (tree_t *)malloc(sizeof(tree_t));
    tree->root = (cell_t *)malloc(SIZE * sizeof(cell_t));
    tree->root->snv = (mut_t *)malloc(sizeof(mut_t));
    if ( (!tree) || (!tree->root) || (!tree->root->snv) ) {
        fprintf(stderr, \
                "[Err::%s::%d] Failed to allocate memory!\n", __func__, __LINE__);
        return NULL;
    }
    father = tree->root; germ = father->snv;

    tree->typenum = 1; father->cellnum[0] = 1; germ->altnum = 0;
    germ->pos = NULL; germ->type = NULL;

    return tree;
}

uint8_t *LoadRef( char *filename )
{
    uint8_t *fa;
    FILE *fp = fopen(filename, "r");

    fa = (uint8_t *)malloc((GSIZE +8) * sizeof(uint8_t));
    if ( !fa ) {
        fprintf(stderr, \
                "[Err::%s::%d] Failed to allocate memory!\n", __func__, __LINE__);
        return NULL;
    }
    fgets(fa, (GSIZE +1), fp); fgets(fa, (GSIZE +1), fp); fclose(fp);

    return fa;
}

int LocIndex( tree_t *tree, int num )
{
    int curcell = 0;

    for ( int i=0; i < tree->typenum; i++ ) {
        curcell += tree->root[i].cellnum[GENERATION-1];
        if ( num < curcell )
            return i;
    }
    fprintf(stderr, "[Err::%s::%d] Never to be here!\n", __func__, __LINE__);
    return (-1);
}

int *RandSelect( tree_t *tree, int depth )
{
    int *count;
    int index, cellnum = 1 << (GENERATION-1);

    srand((unsigned)time(NULL));
    count = (int *)calloc(tree->typenum, sizeof(int));
    if ( !count ) {
        fprintf(stderr, \
                "[Err::%s::%d] Failed to allocate memory!\n", __func__, __LINE__);
        return NULL;
    }
    for ( int i=0; i < depth; i++ ) {
        index = LocIndex(tree, rand()%cellnum); 
        count[index]++;
    }

    return count;
}

void BaseSort( pile_t *pile )
{
    uint32_t num, i, j, *n = pile->num;
    uint8_t base, *b = pile->base;

    for ( i=1; i < 4; i++ ) {
        num = n[i]; base = b[i];
        for ( j=i; j > 0 && n[j-1] < num; j-- ) {
            n[j] = n[j-1]; b[j] = b[j-1];
        }
        n[j] = num; b[j] = base;
    }
}

void Pileup( pile_t *P, tree_t *tree, int *select )
{
    mut_t *snv;

    for ( int i=1; i < tree->typenum; i++ ) 
    {
        if ( select[i] ) {
            snv = tree->root[i].snv;
            for ( int j=0; j < snv->altnum; j++ ) 
            {
                uint32_t pos = snv->pos[j];
                for ( int k=0; k < 4; k++ ) 
                {
                    if ( snv->type[j] == P[pos].base[k] ) {
                        P[pos].num[k] += select[i]; P[pos].num[0] -= select[i];
                        break; 
                    } 
                }
            }
        } 
    } 
    for ( int i=0; i < GSIZE; i++ ) 
        BaseSort(&P[i]); 
}

pile_t *PileInit( uint8_t *ref, int depth )
{
    pile_t *P; uint8_t mtype[4] = {65, 67, 71, 84};

    P = (pile_t *)calloc(GSIZE, sizeof(pile_t));
    if ( !P ) {
        fprintf(stderr, \
                "[Err::%s::%d] Failed to allocate memory!\n", __func__, __LINE__);
        return NULL;
    }
    for ( int i=0; i < GSIZE; i++ ) {
        memcpy(P[i].base, mtype, 4*sizeof(uint8_t)); 
        P[i].num[0] = depth; P[i].base[0] = ref[i];

        for ( int j=1; j < 4; j++ ) {
            if ( P[i].base[j] == ref[i] ) P[i].base[j] = 65;
        }
    }
    return P;
}

void SnvWrite( char *snvname, pile_t *P )
{
    FILE *fp = fopen(snvname, "w");

    fprintf(fp, "Chrom\tPos\tRef\tAlt\tRefnum\tAltnum\n");
    for ( int i=0; i < GSIZE; i++ ) {
        if ( P[i].num[1] )
            fprintf(fp, "Ecoli\t%d\t%c\t%c\t%u\t%u\n", i+1, P[i].base[0], P[i].base[1], P[i].num[0], P[i].num[1]);
    }

    fclose(fp);
}

void Usage( int argc )
{
    char *usage = "\n Usage: celldiv <ref.fa> <outfile> <depth>\n";

    if ( argc != 4 ) {
        fprintf(stderr,"%s", usage); exit(-1);
    }
}


int main( int argc, char **argv )
{
    time_t start, end;
    tree_t *tree = Inital();

    Usage(argc); time(&start);
    srand((unsigned)time(NULL)); /* random the seed */
    for ( int i=0; i < GENERATION-1; i++ ) {
        printf("Start to divide the %d generation ...\n", i+2);
        Divide(tree, i);
    }
   
    printf("\nStart to load the reference ...\n");
    uint8_t *fa = LoadRef(argv[1]); 

    printf("Start to select the cell from the pool randomly ...\n");
    int *select = RandSelect(tree, atoi(argv[3])); 
    pile_t *P = PileInit(fa, atoi(argv[3])); Pileup(P, tree, select);

    printf("Write the snv info out ...\n");
    SnvWrite(argv[2], P);

    time(&end); printf("Time consume is %.2f\n", difftime(end,start));
}

