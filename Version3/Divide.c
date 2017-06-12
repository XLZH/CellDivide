/* Divide function
   Modified: 2017-6-9
*/
#include "Share.h"

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
    mut_t *snv; uint64_t pos = 0;

    pos = RAND64() % FREQLEVEL; type = mtype[rand() % 4];
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

tree_t *Init( mut_t *S, int flag )
{
    tree_t *tree; cell_t *father; mut_t *germ;

    tree = (tree_t *)malloc(sizeof(tree_t));
    tree->root = (cell_t *)malloc(SIZE * sizeof(cell_t));
    memset(tree->root, 0, sizeof(cell_t));
    tree->root->snv = (mut_t *)malloc(sizeof(mut_t));
    if ( (!tree) || (!tree->root) || (!tree->root->snv) ) {
        fprintf(stderr, \
            "[Err::%s::%d] Failed to allocate memory!\n", __func__, __LINE__);
        return NULL;
    }
    father = tree->root; germ = father->snv;

    tree->typenum = 1; father->cellnum[0] = 1;
    if ( flag == 0 ) {
        germ->altnum = 0; germ->pos = NULL; germ->type = NULL;
    }
    else {
        if ( S->altnum ) {
            germ->pos = (uint32_t *)malloc((S->altnum) * sizeof(uint32_t));
            germ->type = (uint8_t *)malloc((S->altnum) * sizeof(uint8_t));
            if ( !(germ->pos) || !(germ->type) ) {
                fprintf(stderr, \
                    "[Err::%s::%d] Failed to allocate memory!\n", __func__, __LINE__);
                return NULL;
            }
            memcpy(germ->pos, S->pos, S->altnum * sizeof(uint32_t));
            memcpy(germ->type, S->type, S->altnum * sizeof(uint8_t));
        }
        else { germ->pos = NULL; germ->type = NULL; }
        germ->altnum = S->altnum;
    }
    return tree;
}

void Destory( tree_t *tree )
{
    mut_t *snv;

    for ( int i=0; i < tree->typenum; i++ ) {
        snv = tree->root[i].snv;
        free(snv->pos); free(snv->type); free(snv);
    }
    free(tree->root); tree->root = NULL; free(tree);
}


