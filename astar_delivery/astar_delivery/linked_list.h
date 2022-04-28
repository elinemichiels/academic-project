/*
Defines the datatype and the functions needed to implement the priority queue
as a linked list
*/

#include <stdlib.h>
#include "utilities.h"

//Datatype for one element in the linked list
typedef struct QueueElementstructure {
    unsigned index; //Index of the node represented in the nodes vector.
    struct QueueElementstructure *next; //Pointer to the next node in the list.
    } QueueElement;
typedef QueueElement * PriorityQueue;

//Struct to save the data needed for each node to compute the A*-algorithm.
typedef struct{
    double f;
    double g;
    unsigned parent;
    char isOpen;
} AStarData;

//Returns true if the PriorityQueue is empty.
int IsEmpty( PriorityQueue Pq );

//Returns the index of the minimal element in the PriorityQueue and deletes it from the list.
unsigned extract_min(PriorityQueue *Pq, AStarData *PathData);

//Adds the index of a node to the linked list at the correct position according to its f-value.
void add_with_priority(unsigned index, PriorityQueue *Pq, AStarData *PathData);

//Changes the position of an element in the linked list if it's priority is updated during the algorithm.
void requeue_with_priority(unsigned index, PriorityQueue *Pq, AStarData *PathData);