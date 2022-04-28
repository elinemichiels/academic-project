===== A* Delivery - Eline Michiels and Jannik Jacobsen =====

We give all commands assuming the working directory of the terminal is the astar_delivery subfolder.

== FILE STRUCTURE ==
We split our functions over multiple files to make the code more easy to understand.

astar.c:	Main loop of the A* algorithm and distance functions.
binary_heap.c:	Functions and definitions to implement the priority queue as a binary heap.
build_graph.c:	Program to read the csv-file, build the graph and save it as a binary file.
linked_list.c:	Functions and definitions to implement the priority queue as a linked list.
main.c:		Main program, calls all functions to execute the A* algorithm after building the graph.
utilities.c:	Auxiliary definitions like the node structure or functions to save and read the binary file.


- Building the graph -
The code for building the graph is in build_graph.c. It uses some definitions from utilities.c, so it can be compiled using
gcc -Ofast build_graph.c utilities.c -o build_graph

The csv-file to be used must be specified as an command line argument.
We kept it in the ./data subfolder, so the command to execute the program reads
./build_graph data/spain.csv
The binary file will be placed in the same directory as the csv-file.

- Executing the algorithm -
The program in main.c calls all necessary functions to execute the algorithm.
The function for the A* algorithm is defined in astar.c.
The data structure for the priority queue is controlled in astar.h:
One of line 5 and line 6 is to be commented.
It also needs the auxiliary definitions in utilities.c, so the command for building the program is
gcc -Ofast main.c astar.c binary_heap.c utilities.c -o main
if the binary heap is to be used.
For the linked list, replace binary_heap.c by linked_list.c and change the command in astar.h.

It needs 3 command line arguments: The path to the binary file, the ID of the start and the ID of the end node.
This results in:
./main data/spain.bin 240949599 195977239

This will generate a csv-file with the coordinates of the nodes in the path found, which will be saved in the ./results subfolder.
With an R installation, the path can be visualized using the make_map.R script.

We provide the csv-file and map of the path as well as the executables of all programs built under Ubuntu.