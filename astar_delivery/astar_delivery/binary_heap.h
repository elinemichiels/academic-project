/*
Defines the dataype and functions needed to implement the
priority queue as a binary heap.
*/

#include <stdlib.h>
#include "utilities.h"

//Datatype of one element in the binary heap
typedef struct QueueElementstructure {
    unsigned index; //Index of the represented node in the nodes vector
    struct QueueElementstructure *left; //Left child
    unsigned left_size; //Total size of the left subtree
    struct QueueElementstructure *right; //Right child
    unsigned right_size; //Total size of the right subtree
    struct QueueElementstructure *parent;
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

//Returns the index of the element with the minimal priority then reorders the tree
//to be in the correct shape again.
unsigned extract_min(PriorityQueue *Pq, AStarData *PathData);

//Adds an element to the tree and moves its index to the correct position.
void add_with_priority(unsigned index, PriorityQueue *Pq, AStarData *PathData);

//Handles the update of the priority of an element.
void requeue_with_priority(unsigned index, PriorityQueue *Pq, AStarData *PathData);

//Restores the tree so that after insertion or a change in priority, the element at q is
//in the correct position.
void restore(QueueElement *q, AStarData *PathData);