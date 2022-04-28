/*
Defines utility functions needed for handling the nodes and building the graph.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

//We use an integer larger than the total number of nodes to indicate
//that a node has not been found or is not initalized in the linking process.
#define NOTFOUND 25000000
#define DEBUG_MAX_VALENCE 9
#define DEBUG_MAX_PRINT 5

typedef struct {
    unsigned long id; // Node identification
    double lat, lon; // Node position
    unsigned short nsucc; // Number of node successors; i. e. length of successors
    unsigned *successors; // Pointer to the field of the succesor's indices in the nodes vector
} node;

// Extract the info from the line in buffer and save it to the struct at position i in nodes.
void readNodeLine(char *buffer, node *nodes, int i);

// Print the info of the first length structs.
void printNodesList(node *nodes, int length);

// Search the position of the node with id id in the array of nodes.
// This is done using the binary search algorithm.
// If it is not found, NOTFOUND is returned.
unsigned nodesearch(node *first, unsigned long id, unsigned long length);

// Search, if the node with next_id is already in the succesors of the node at prev_pos in nodes.
// If yes, 1 is returned, otherwise 0.
int searchInSuccessors(node *first, unsigned prev_pos, unsigned long next_id);

// Adds the position next_pos to the list of successors of the node at prev_pos and increments
// nodes[prev_pos].nsucc by 1.
void linkNodes(node *nodes, unsigned prev_pos, unsigned next_pos);

// Computes and displays the distribution of valences of the nodes.
void computeValences(node *nodes, unsigned long length, int max_valence);

// Exit the program with code errcode, display the message miss and quit the program.
void ExitError(const char *miss, int errcode);

//Write the nodes vector of length num_nodes to a binary file at filename.
void writeBinary(node *nodes, unsigned long num_nodes, const char *filename);

//Read the nodes vector located at filename. Returns a pointer to the start of a nodes vector
//and saves the number of nodes in *num_nodes_arg
node *readBinary(const char* filename, unsigned long *num_nodes_arg);