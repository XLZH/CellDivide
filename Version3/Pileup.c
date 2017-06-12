/* This is the source file of Pileup

   Modified: 2015-06-03
   Fix:(BUG)BaseSort: insert sort from index 2.
*/
#include "Share.h"

#define DEALSNV( _P, _snv, _start, _end )                               \
    for (int _i=0; _i < _snv->altnum; _i++) {                           \
        uint32_t _pos = _snv->pos[_i];                                  \
        if ( _pos >= _start && _pos < _end) {                           \
            for ( int _j=0; _j < 4; _j++ ) {                            \
                if ( _snv->type[_i] == _P[_pos].base[_j] ) {            \
                    _P[_pos].num[_j]++; _P[_pos].num[0]--; break; } } } \
    }

uint8_t *LoadRef( char *filename )
{
    uint8_t *fa, *cur, ch;

    FILE *fp = fopen(filename, "r");
    if ( !fp ) fprintf(stderr, \
            "[Err::%s::%d] Failed to open %s!\n", __func__, __LINE__, filename);
    fa = cur = (uint8_t *)malloc(GSIZE * sizeof(uint8_t));
    if ( !fa ) 
        fprintf(stderr, \
            "[Err::%s::%d] Failed to allocate memory!\n", __func__, __LINE__);
    while ( fgetc(fp) != '\n' ) ; // remove the fasta lines that starts with '>'
    while ( (ch = fgetc(fp)) != 0xff ) {
        if ( ch != '\r' && ch != '\n' ) *cur++ = ch;
    }
    return fa;
}

void BaseSort( pile_t *pile )
{
    uint32_t num, i, j, *n = pile->num;
    uint8_t base, *b = pile->base;

    for ( i=2; i < 4; i++ ) { // insert sort algorithm
        num = n[i]; base = b[i];
        for ( j=i; j > 1 && n[j-1] < num; j-- ) {
            n[j] = n[j-1]; b[j] = b[j-1];
        }
        n[j] = num; b[j] = base;
    }
}

int LocIndex( int *clist, int num, uint64_t pos )
{
    if ( pos <= clist[0] ) return 0;

    int low =1, mid, high =num -1;
    while ( low <= high ) {
        mid = (low + high) / 2;
        if ( pos < clist[mid] )
            high = mid -1;
        else if ( pos > clist[mid] )
            low = mid + 1;
        else return mid;
    }
    return low;
}

void Pileup( pile_t *P, tree_t *tree, int depth )
{
    mut_t *snv; 
    int *clist; uint64_t cellnum = 1 << (GENERATION-1);

    clist = (int *)malloc(tree->typenum * sizeof(int));
    if ( !clist )
        fprintf(stderr, \
            "[Err::%s::%d] Failed to allocate memory!\n", __func__, __LINE__);
    clist[0] = tree->root[0].cellnum[GENERATION-1];
	for ( int i=1; i < tree->typenum; i++ )
        clist[i] = clist[i-1] + tree->root[i].cellnum[GENERATION-1];

    srand((unsigned)time(NULL));
    for ( int b=0; b < GSIZE; b+=FRAGSIZE) {
        int start=b, end, index; uint64_t pos;
        end = b+FRAGSIZE > GSIZE ? GSIZE : b+FRAGSIZE;
        for ( int i=0; i < depth; i++ ) { /* random select the cell from pool */
            pos = RAND64() % cellnum; 
            index = LocIndex(clist, tree->typenum, pos);
            snv = tree->root[index].snv; DEALSNV(P, snv, start, end);
        }
    } free(clist);
    for ( int i=0; i < GSIZE; i++ )  BaseSort(&P[i]);
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
    float altfreq;

    fprintf(fp, "Chrom\tPos\tRef\tAlt\tRefnum\tAltnum\n");
    for ( int i=0; i < GSIZE; i++ ) {
        if ( P[i].num[1] ) {
            altfreq = (float)(P[i].num[1]) / (P[i].num[0] + P[i].num[1]) * 100;
            fprintf(fp, "Ecoli\t%d\t%c\t%c\t%u\t%u\t%.2f\n", \
                i+1, P[i].base[0], P[i].base[1], P[i].num[0], P[i].num[1], altfreq);
        }
    }
    fclose(fp);
}

