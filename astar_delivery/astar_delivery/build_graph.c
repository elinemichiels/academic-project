/*
Main function for reading the .csv-file and building the graph out of it.
Works for the spain.csv and cataluna.csv files.
*/

#include "utilities.h"

int main(int argc, char *argv[]){

    if (argc != 2) ExitError("Please specify the path to the .csv-file to read as an argument", 10);

    FILE *fp;
    size_t buffer_size = 0;
    size_t num_chars;
    char *buffer = NULL, *field = NULL;
    node *nodes;
    unsigned i;
    unsigned long num_nodes;
    int j;

    clock_t start_build = clock();

    /* Open the file */
    printf("Opening file %s\n\n", argv[1]);
    fp = fopen(argv[1], "r");
    if (!fp) ExitError("Unable to open file", 11);

    /* Count number of nodes and allocate memory for the nodes vector */
    for (i=0; i<3; i++)
        getline(&buffer, &buffer_size, fp); // Discard the first three lines

    num_nodes = 0UL;
    while (!feof(fp)){
        if (num_nodes%1000000 == 0)
            printf("Reading line %lu...\n", num_nodes);
        num_chars = getline(&buffer, &buffer_size, fp);
        if (strncmp(buffer, "node", strlen("node")))
            break;
        else
            num_nodes++;
    }

    printf("Found %lu lines.\n\n", num_nodes);

    if ((nodes = (node*) malloc(num_nodes*sizeof(node))) == NULL)
        ExitError("Unable to allocate memory for the nodes", 12);

    /* Read the information line by line and save it in the struct */
    rewind(fp);
    for (i=0; i<3; i++)
        getline(&buffer, &buffer_size, fp); // Discard the first three lines

    i = 0U;
    while(!feof(fp) && (i<num_nodes)){
        if (i%1000000 == 0)
            printf("Extracting info of node %d...\n", i);
        num_chars = getline(&buffer, &buffer_size, fp);
        readNodeLine(buffer, nodes, i);
        nodes[i].successors = NULL;
        nodes[i].nsucc = 0U;
        i++;
    }
    printf("Got info of all nodes. \n\n");

    /* Read the information of the ways */
    int oneway = 0;
    int count = 1;
    int found;
    unsigned long prev_id, next_id;
    unsigned prev_pos, next_pos;
    num_chars = getline(&buffer, &buffer_size, fp);
    while(!strncmp(buffer, "way", 3)){
        if (count%100000 == 0)
            printf("Reading way %d\n", count);
        for (j=1; j<=8; j++)
            field = strsep(&buffer, "|");
        if (!strcmp(field, "oneway"))
            oneway = 1;

        field = strsep(&buffer, "|");
        field = strsep(&buffer, "|");
        prev_id = 0UL;
        prev_pos = NOTFOUND;
        //Traverse the fields and link the nodes if they exist.
        while (field != NULL){
            if (!strcmp(field, ""))
                ExitError("Empty field while reading way", 13);
            sscanf(field, "%lu", &next_id);       
            field = strsep(&buffer, "|");
            next_pos = nodesearch(nodes, next_id, num_nodes);
            if (next_pos == NOTFOUND)
                continue;
            if (prev_pos != NOTFOUND){
                found = searchInSuccessors(nodes, prev_pos, next_id);
                if (!found)          
                    linkNodes(nodes, prev_pos, next_pos);
                // Do it the other way around if not oneway
                if (!oneway){
                    found = searchInSuccessors(nodes, next_pos, prev_id);
                    if (!found)             
                        linkNodes(nodes, next_pos, prev_pos);
                }
            }
            prev_pos = next_pos;
            prev_id = next_id;
        }
        oneway = 0;
        num_chars = getline(&buffer, &buffer_size, fp);
        count++;
    }
    computeValences(nodes, num_nodes, DEBUG_MAX_VALENCE);
    writeBinary(nodes, num_nodes, argv[1]);

    clock_t end_build = clock();
    printf("Done in %lf seconds.\n", ((double) end_build-start_build)/CLOCKS_PER_SEC);
    fclose(fp);
}