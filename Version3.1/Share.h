/* Version: 3
*  date: 2017-6-22
*/

#ifndef Divide_H
#define Divide_H

#include <stdio.h>
#include <linux/limits.h>
#include <getopt.h>
#include <time.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>

#define GENERATION 64          // the maximum generation allowed
#define SIZE 5000000           // default size for memory realocated each time
#define FRAGSIZE 500           // default fragment size for pileup

#define RAND64() ( (((uint64_t)rand() << 0) & 0x00000000FFFFFFFF) | \
                   (((uint64_t)rand() << 32) & 0xFFFFFFFF00000000) )

/* two different division model 
 * Normal: The original cell is dividing for x generation
 * Maline: The original cell is dividing for x generation and then
 *         select one cell from the pool randomly and transfer for
 *         y times
*/
enum MODEL { NORMAL = 1, MALINE = 2 };


/*! @typedef mut_t
 @abstract structure for the mutation infomation
 @field altnum       snv num for this mutation type
 @field pos          pointer to the snv position array
 @field type         pointer to the snv type array [ A|T|C|G ]
*/
typedef struct __mut_t {
    uint32_t    altnum;
    uint32_t    *pos;
    uint8_t     *type;
} mut_t;


/*! @typedef cell_t
 @abstract structure for the mutation type
 @field cellnum      cellnum for each generation
 @field snv          pointer to the mutation infomation
*/
typedef struct __cell_t {
    uint32_t    cellnum[GENERATION];
    mut_t       *snv;
} cell_t;


/*! @typedef tree_t
 @abstract structure for the mutation tree
 @field typenum      total mutation type count
 @field root         pointer to the mutation type array
*/
typedef struct __tree_t {
    uint32_t    typenum;
    cell_t      *root;
} tree_t;


/*! @typedef pile_t
 @abstract structure for the pileup
 @field base         ordered with [ref, alt, other1, other2]
 @field num          coverage for the base
*/
typedef struct __pile_t {
    uint8_t     base[4];
    uint32_t    num[4];
} pile_t;


/*! @typedef arg_t
 @abstract structure for the args.
 @field ref          reference path for this cell type
 @field snvfile      snv filename for output
 @field depth        equal to sequencing depth
 @field generation   generation for each transfer
 @field transnum     equal to transfer times
 @field genomesize   genome size for the reference
 @field help         help status
 @field freq         default mutantion rate(scientific notation)
*/
typedef struct __arg_t {
    uint8_t ref[PATH_MAX];
    uint8_t snvfile[PATH_MAX];
    uint32_t depth, generation;
    uint32_t transnum, gsize, help;
    uint64_t rate;
} arg_t;

/* Parse the parmeters */
arg_t *ParseOpt( int, char ** );

/* Divide function */
void Usage( int );
void Divide( tree_t *, int, arg_t * );
tree_t *Init( mut_t *, int );
void Destory( tree_t * );

/* Pileup function */
uint8_t *LoadRef( arg_t * );
pile_t *PileInit( uint8_t *, arg_t * );
void Pileup( pile_t *, tree_t *, arg_t * );
int LocIndex( int *, int, uint64_t );
void SnvWrite( pile_t *, arg_t * );

#endif
