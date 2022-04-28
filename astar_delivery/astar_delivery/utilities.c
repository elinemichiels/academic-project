#include "utilities.h"

void readNodeLine(char *buffer, node *nodes, int i){
    char *field;
    int j;
    field = strsep(&buffer, "|");
    if (strcmp(field, "node")){
        printf("No Node\n");
        exit(1);
    }
    else{
        field = strsep(&buffer, "|");
        sscanf(field, "%lu", &(nodes[i].id));
        for (j = 3; j<=9 ; j++){
            field = strsep(&buffer, "|");
        }
        field = strsep(&buffer, "|");
        sscanf(field, "%lf", &(nodes[i].lat));
        field = strsep(&buffer, "|");
        sscanf(field, "%lf", &(nodes[i].lon));
    }  
}

void printNodesList(node *nodes, int length){
    int i, j;
    for (i=0; i<length; i++){
        printf("= Node %d =\n", i+1);
        printf("ID: %lu\n", nodes[i].id);
        printf("Latitude: %lf\n", nodes[i].lat);
        printf("Longitude: %lf\n", nodes[i].lon);
        printf("Number of successors: %u\n", nodes[i].nsucc);
        printf("IDs of successors: ");
        for (j=0; j<nodes[i].nsucc; j++)
            printf("%lu, ", nodes[nodes[i].successors[j]].id);
        printf("\n");
    }
}

unsigned nodesearch(node *nodes, unsigned long id, unsigned long length){
    unsigned first = 0;
    unsigned last = length - 1;
    unsigned middle = (first+last)/2;

    while (first <= last) {
        if (nodes[middle].id < id)
            first = middle + 1;
        else if (nodes[middle].id == id) {
            return(middle);
        }
        else
            last = middle - 1;
        middle = (first + last)/2;
        }

    return(NOTFOUND);
}

int searchInSuccessors(node *nodes, unsigned prev_pos, unsigned long next_id){
    int k;
    unsigned *succs;
    succs = nodes[prev_pos].successors;
    for (k=0; k<nodes[prev_pos].nsucc; k++){
        if (nodes[succs[k]].id == next_id){
            return (1);
        }
    }
    return(0);
}

void linkNodes(node *nodes, unsigned prev_pos, unsigned next_pos){
    nodes[prev_pos].successors = (unsigned *) realloc(nodes[prev_pos].successors, (nodes[prev_pos].nsucc+1)*sizeof(unsigned));
    if (nodes[prev_pos].successors == NULL)
        ExitError("Reallocation failed while linking nodes", 14);
    nodes[prev_pos].successors[nodes[prev_pos].nsucc] = next_pos;
    nodes[prev_pos].nsucc++;
}

void computeValences(node *nodes, unsigned long length, int max_valence){
    int i;
    unsigned *valences;
    unsigned too_large = 0U;
    if ((valences = (unsigned*) malloc((max_valence+1)*sizeof(unsigned))) == NULL)
        ExitError("Allocation of memory while computing valences failed", 15);
    memset(valences, 0, (max_valence+1)*sizeof(unsigned));

    for (i=0; i<length; i++){
        if (nodes[i].nsucc <= max_valence)
            valences[nodes[i].nsucc]++;
        else
            too_large++;
    }

    printf("\n=== Valences of nodes found ===\n");
    for (i=0; i<=max_valence; i++){
        printf("Nodes with valence %d:\t%u\n", i, valences[i]);
    }
    printf("Nodes with valence higher than the max of %d: %u\n\n", max_valence, too_large);
}

void ExitError(const char *miss, int errcode) {
    fprintf (stderr, "\n\nERROR: %s.\n\nStopping...\n\n", miss); exit(errcode);
}

void writeBinary(node *nodes, unsigned long num_nodes, const char *filename){
    FILE *fin;
    char name[257];
    int i;
    printf("\nWriting binary file\n");

    /* Computing the total number of successors */
    unsigned long ntotnsucc=0UL;
    for(i=0; i < num_nodes; i++) ntotnsucc += nodes[i].nsucc;

    /* Setting the name of the binary file */
    strcpy(name, filename); strcpy(strrchr(name, '.'), ".bin");

    if ((fin = fopen (name, "wb")) == NULL)
        ExitError("the output binary data file cannot be opened", 31);

    /* Global data 􀀀􀀀􀀀 header */
    if( fwrite(&num_nodes, sizeof(unsigned long), 1, fin) +
        fwrite(&ntotnsucc, sizeof(unsigned long), 1, fin) != 2 )
            ExitError("when initializing the output binary data file", 32);

    /* Writing all nodes */
    if( fwrite(nodes, sizeof(node), num_nodes, fin) != num_nodes )
        ExitError("when writing nodes to the output binary data file", 32);

    /* Writing sucessors in blocks */
    for(i=0; i < num_nodes; i++) if(nodes[i].nsucc){
        if( fwrite(nodes[i].successors, sizeof(unsigned), nodes[i].nsucc, fin) != nodes[i].nsucc )
            ExitError("when writing edges to the output binary data file", 32);
    }
    fclose(fin);
}

node *readBinary(const char* filename, unsigned long *num_nodes_arg){
    FILE *fin;
    unsigned long ntotnsucc;
    unsigned long num_nodes;
    node *nodes;
    unsigned *allsuccessors;
    int i;

    if ((fin = fopen(filename, "r")) == NULL)
        ExitError("the data file does not exist or cannot be opened", 41);

    printf("Reading binary file %s\n", filename);

    /* Global data - header */
    if(fread(&num_nodes, sizeof(unsigned long), 1, fin) +
        fread(&ntotnsucc, sizeof(unsigned long), 1, fin) != 2)
            ExitError("when reading the header of the binary data file", 42);
    *num_nodes_arg = num_nodes;

    /* getting memory for all data */
    if((nodes = (node *) malloc(num_nodes*sizeof(node))) == NULL)
        ExitError("when allocating memory for the nodes vector", 43);
    if((allsuccessors = (unsigned *) malloc(ntotnsucc*sizeof(unsigned))) == NULL)
        ExitError("when allocating memory for the edges vector", 44);

    /* Reading all data from file */
    if(fread(nodes, sizeof(node), num_nodes, fin) != num_nodes)
        ExitError("when reading nodes from the binary data file", 45);
    if(fread(allsuccessors, sizeof(unsigned), ntotnsucc, fin) != ntotnsucc)
        ExitError("when reading sucessors from the binary data file", 46);

    fclose(fin);

    /* Setting pointers to successors */
    for(i=0; i < num_nodes; i++)
        if(nodes[i].nsucc){
            nodes[i].successors = allsuccessors; allsuccessors += nodes[i].nsucc;
        }
    
    return(nodes);
}