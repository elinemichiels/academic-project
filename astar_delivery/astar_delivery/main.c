#include "astar.h"

int main(int argc, char *argv[]){
    
    if (argc != 4) ExitError("Please specify the path to the binary file , start ID and end ID to read as arguments", 10);

    unsigned long num_nodes;
    clock_t start_binary = clock();
    node *nodes = readBinary(argv[1], &num_nodes);
    clock_t end_binary = clock();
    printf("Done in %lf seconds.\n\n", ((double) end_binary-start_binary)/CLOCKS_PER_SEC);

    unsigned long start;
    unsigned long end;
    sscanf(argv[2], "%lu", &start);
    sscanf(argv[3], "%lu", &end);

    char filename[257];
    char basename[50];
    sprintf(basename, "_%lu_%lu.csv", start, end);
    strcpy(filename, "results");
    strcat(filename, strrchr(argv[1], '/'));
    strcpy(strrchr(filename, '.'), basename);

    computeValences(nodes, num_nodes, DEBUG_MAX_VALENCE);

    unsigned start_index = nodesearch(nodes, start, num_nodes);
    if (start_index == NOTFOUND)
        ExitError("No start node with that ID", 91);
    unsigned end_index = nodesearch(nodes, end, num_nodes);
    if (end_index == NOTFOUND)
        ExitError("No end node with that ID", 92);

    printf("Initilazing AStar...\n");
    AStarData *PathData = astar_init(num_nodes);

    printf("Finding way from ID %lu to ID %lu...\n", start, end);
    clock_t start_astar = clock();
    astar(nodes, start_index, end_index, PathData);
    clock_t end_astar = clock();

    printf("Path found with distance %lfm in %.3lf seconds\n", PathData[end_index].g, ((double) end_astar-start_astar)/CLOCKS_PER_SEC);
    make_path(nodes, start_index, end_index, PathData, filename);
    printf("csv-file with the coordinates saved to %s\n", filename);
    
    return(0);
}