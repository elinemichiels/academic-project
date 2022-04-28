#include "astar.h"

const double R = 6371e3; //Radius of the earth.
double haversine(node *node1, node *node2){
    double lat1 = node1->lat;
    double lon1 = node1->lon;
    double lat2 = node2->lat;
    double lon2 = node2->lon;

    double phi1 = lat1*M_PI/180;
    double phi2 = lat2*M_PI/180;
    double delta_phi = (lat2-lat1)*M_PI/180;
    double delta_lambda = (lon2-lon1)*M_PI/180;

    double a = pow(sin(delta_phi/2), 2) + (cos(phi1)*cos(phi2)*pow(sin(delta_lambda/2), 2));
    double c = 2 * atan2(sqrt(a), sqrt(1-a));

    return (R*c);
}

double equirec(node *node1, node *node2){
    double lat1 = node1->lat;
    double lon1 = node1->lon;
    double lat2 = node2->lat;
    double lon2 = node2->lon;

    double phi1 = lat1*M_PI/180;
    double phi2 = lat2*M_PI/180;
    double delta_lambda = (lon2-lon1)*M_PI/180;

    double x = delta_lambda * cos((phi1+phi2)/2);
    double y = phi2-phi1;
    return ((sqrt(x*x + y*y))*R);
}

void astar(node *nodes, unsigned source, unsigned goal, AStarData *PathData){
    PriorityQueue Open = NULL;

    PathData[source].g = 0.0;
    PathData[source].parent = UINT_MAX;
    PathData[source].f = haversine(nodes+source, nodes+goal);
    add_with_priority(source, &Open, PathData);

    while(!IsEmpty(Open)){
        unsigned i;
        unsigned min_index = extract_min(&Open, PathData);
        if (min_index == goal)
            return;
        if (PathData[min_index].isOpen == 0)
            continue;
        for(i=0; i<nodes[min_index].nsucc; i++){
            unsigned adj = nodes[min_index].successors[i];
            double g_adj = PathData[min_index].g + equirec(nodes+min_index, nodes+adj);
            if (g_adj < PathData[adj].g){
                PathData[adj].parent = min_index;
                double h_adj;
                if (PathData[adj].g < DBL_MAX)
                    h_adj = PathData[adj].f - PathData[adj].g;
                else
                    h_adj = haversine(nodes+adj, nodes+goal);
                PathData[adj].f = g_adj + h_adj;
                PathData[adj].g = g_adj;
                if (!PathData[adj].isOpen)
                    add_with_priority(adj, &Open, PathData);
                else
                    requeue_with_priority(adj, &Open, PathData);
            }
        }
        PathData[min_index].isOpen = 0;        
    }
}

void make_path(node *nodes, unsigned start_index, unsigned end_index, AStarData *PathData, const char *filename){
    unsigned index;
    int path_count = 0;

    //Trace the path backwards and count how many nodes we need
    for (index=end_index; index != UINT_MAX; index = PathData[index].parent)
        path_count++;

    unsigned *path = (unsigned *) malloc(path_count*sizeof(unsigned));
    if (path==NULL)
        ExitError("Allocation of path index vector failed", 71);
    int j = 0;
    
    //Save the indices of the nodes in the path into an array
    for (index=end_index; index != UINT_MAX; index = PathData[index].parent){
        path[j] = index;
        j++;
    }

    FILE *fout = fopen(filename, "w");
    fprintf(fout, "lon,lat\n");
    node current;

    //Write the latitudes and longitudes in a file.
    for (j=path_count-1; j>=0; j--){
        current = nodes[path[j]];
        fprintf(fout, "%lf,%lf\n", current.lon, current.lat);
    }

}

AStarData *astar_init(unsigned long num_nodes){
    //Allocate memory for the AStarData vector
    AStarData *PathData = (AStarData *) malloc(num_nodes*sizeof(AStarData));
    if (PathData == NULL)
        ExitError("Initialization of the AStarData failed", 51);
    unsigned i;

    //Initialize the f- and g-values for each node to DBL_MAX, which represents infinity,
    //and set the isOpen-flag to zero.
    for (i=0; i<num_nodes; i++){
        PathData[i].f = DBL_MAX;
        PathData[i].g = DBL_MAX;
        PathData[i].isOpen = 0;
    }
    return PathData;
}