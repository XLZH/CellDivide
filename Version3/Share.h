/* Version: 3
*  Modified: 2017-6-9
*  Add: RAND64(): which range (0, 2^64 -1)
*/

#ifndef Divide_H
#define Divide_H

#include <stdio.h>
#include <time.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>

#define GENERATION 29          // the generation that cultured the cell
#define GSIZE 4639675          // the genome size of Ecoli
#define FREQLEVEL 4000000000   // the mutation frequence level for whole genome
#define SIZE 5000000           // default size for memory realocated each time
#define FRAGSIZE 500           // default fragment size

#define RAND64() ( (((uint64_t)rand() << 0) & 0x00000000FFFFFFFF) | \
                   (((uint64_t)rand() << 32) & 0xFFFFFFFF00000000) )

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

/* Divide function */
void InsertAlt( mut_t *, const mut_t *, uint32_t, uint8_t );
mut_t *Mutator( const mut_t * );
void CreatCell( tree_t *, mut_t *, int );
void Divide( tree_t *, int );
tree_t *Init( mut_t *, int );
void Destory( tree_t * );

/* Pileup function */
uint8_t *LoadRef( char * );
pile_t *PileInit( uint8_t *, int );
void Pileup( pile_t *, tree_t *, int );
int LocIndex( int *, int, uint64_t );
void BaseSort( pile_t * );
void SnvWrite( char *, pile_t * );

#endif
