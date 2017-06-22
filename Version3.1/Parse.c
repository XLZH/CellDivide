#include "Share.h"


void Usage( int opt )
{
    char *usage_c = 
        "\nUsage: celldiv <comand> [options]\n"
        "\n"
        "Comands:\n"
        "    normal             Normal division\n"
        "                       Cell dividing in a normal way,which divided by\n"
        "                       specific generation for only one transfer\n"
        "    maline             Mutation accumulation line\n"
        "                       Cell dividing in a maline, which transfer for specific\n"
        "                       times and specific generation for each transfer\n\n";
    char *usage_n =
        "\nUsage: celldiv normal [options]\n"
        "\n"
        "Options:\n"
        "       -h|--help             print help infomation\n"
        "       -f|--reference        reference for the cell type[.fa]\n"
        "       -o|--output           output name for snv file\n"
        "       -d|--depth            [integer] equal to sequencing depth\n"
        "       -r|--mutrate          [scientific] mutation rate for cell division[eg: 2.5E-9]\n"
        "       -g|--generation       [integer] generation for each transfer\n\n";
    char *usage_m =
        "\nUsage: celldiv maline [options]\n"
        "\n"
        "Options:\n"
        "       -h|--help             print help infomation\n"
        "       -f|--reference        reference for the cell type[.fa]\n"
        "       -o|--output           output name for snv file\n"
        "       -d|--depth            [integer] equal to sequencing depth\n"
        "       -r|--mutrate          [scientific] mutation rate for cell division[eg: 2.5E-9]\n"
        "       -g|--generation       [integer] generation for each transfer\n"
        "       -t|--transfernum      [integer] transfer times\n\n";
    switch ( opt ) {
        case NORMAL:
            fprintf(stderr, "%s", usage_n); exit(-1); break;
        case MALINE:
            fprintf(stderr, "%s", usage_m); exit(-1); break;
        default:
            fprintf(stderr, "%s", usage_c); exit(-1); break;
    }
}

static const struct option long_options[] = 
{
    { "help", no_argument, NULL, 'h' },
    { "reference", required_argument, NULL, 'f' },
    { "output", required_argument, NULL, 'o' },
    { "depth", required_argument, NULL, 'd' },
    { "mutrate", required_argument, NULL, 'r' },
    { "generation", required_argument, NULL, 'g' },
    { "transfernum", optional_argument, NULL, 't' },
    { NULL, 0, NULL, 0 }
};


static uint64_t Parsefreq( const char *freq )
{
    char fstr[64];

    for ( int i=0; i < 64 && *freq; i++ ) {
        fstr[i] = ( *freq == '-' ) ? '+' : *freq;
        freq++;
    }
    return (uint64_t)(atof(fstr));
}

static uint32_t GetGsize( char *filename )
{
    uint32_t gsize =0;
    uint8_t ch; FILE *fp;

    fp = fopen(filename, "r");
    if ( !fp ) fprintf(stderr, \
            "[Err::%s::%d] Failed to open %s!\n", __func__, __LINE__, filename);
    while ( fgetc(fp) != '\n' ) ; //remove the fasta lines that starts with '>'
    while ( (ch = fgetc(fp)) != 0xff ) {
        if ( ch != '\r' && ch != '\n' ) gsize++;
    }
    fclose(fp); return gsize;
}

arg_t *ParseOpt( int argc, char **argv )
{
    int opt =0, opterr =0;

    arg_t *Arg = (arg_t *)calloc(1, sizeof(arg_t));
    if ( !Arg ) {
        fprintf(stderr, \
            "[Err::%s::%d] Failed to allocate memory!\n", __func__, __LINE__);
        return NULL;
    }
    while ( (opt = getopt_long(argc, argv, "f:o:d:r:g:t:h", long_options, NULL)) != -1 )
    {
        switch (opt) {
            case 'h': Arg->help = 1; break;
            case 'f': strcpy(Arg->ref, optarg); break;
            case 'o': strcpy(Arg->snvfile, optarg); break;
            case 'd': Arg->depth = atoi(optarg); break;
            case 'g': Arg->generation = atoi(optarg); break;
            case 'r': Arg->rate = Parsefreq(optarg); break;
            case 't': Arg->transnum = atoi(optarg); break;
            case '?': fprintf(stderr, "[Err::%s::%d] \
                                      Option error occour!.\n", __func__, __LINE__);
                      Arg->help = 1;
        }
    } Arg->gsize = GetGsize(Arg->ref);
    if ( Arg->generation > 60 ) {
        fprintf(stderr, "[Err::%s::%d] \
            Generation should be smaller than 60!\n", __func__, __LINE__);
        Arg->help = 1;
    }
    return Arg;
}

