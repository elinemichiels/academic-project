#include "binary_heap.h"

int IsEmpty(PriorityQueue Pq){
    return ( Pq == NULL );
}

unsigned extract_min(PriorityQueue *Pq, AStarData *PathData){
    PriorityQueue q = *Pq;
    unsigned index = q->index;
    if (q->left == NULL){
        (*Pq) = NULL;
        return index;
    }
    PriorityQueue prev = q;

    //Find the element of the heap to be put on top. Has to be the one such that the tree
    //is still balanced if we remove it from the bottom.
    while (q->left!=NULL){
        prev = q;
        if (q->left_size == q->right_size){
            q->right_size--;
            q = q->right;
        }
        else{
            q->left_size--;
            q = q->left;
        }
    }

    //Remove the element and copy its index to the root of the tree
    (*Pq)->index = q->index;
    if (prev->left==q){
        free(q);
        prev->left = NULL;
    }
    else{
        free(q);
        prev->right = NULL;
    }

    //Swap the indices of the heap elements until the correct order is restored.
    q = (*Pq);
    unsigned temp;
    PriorityQueue min_element = q;
    while(q->left!=NULL){
        if (PathData[q->left->index].f < PathData[q->index].f)
            min_element = q->left;
        if (q->right!=NULL){
            if (PathData[q->right->index].f < PathData[min_element->index].f)
                min_element = q->right;
        }
        if (min_element==q)
            break;
        else{
            temp = q->index;
            q->index = min_element->index;
            min_element->index = temp;
            q = min_element;
        }
    }
    return index;
}

void add_with_priority(unsigned index, PriorityQueue *Open, AStarData *PathData){
    QueueElement *aux = (QueueElement *) malloc(sizeof(QueueElement));
    if(aux == NULL)
        ExitError("Allocation of a new Queue element failed", 61);

    //Initilaize the data of the new element
    aux->index = index;
    aux->left = NULL;
    aux->left_size = 0;
    aux->right = NULL;
    aux->right_size = 0;
    aux->parent = NULL;
    double costv = PathData[index].f;
    PathData[index].isOpen = 1;

    if (*Open==NULL){
        *Open = aux;
        return;
    }
    QueueElement *q = *Open;
    QueueElement *next = q;

    //Find the parent where to insert the new element such that the tree is balanced
    //after the insertion.
    while(next != NULL){
        q = next;
        if (q->left_size==q->right_size){
            next = q->left;
            q->left_size++;
        }
        else{
            next = q->right;
            q->right_size++;
        }
    }
    if (q->left==NULL)
        q->left = aux;
    else
        q->right = aux;
    aux->parent = q;

    //Call the restore function to swap the index of the new element to the correct place.
    restore(aux, PathData);
    return;
}

void requeue_with_priority(unsigned index, PriorityQueue *Open, AStarData *PathData){
    if((*Open)->index == index)
        return;
    add_with_priority(index, Open, PathData);
    return;
}

void restore(QueueElement *q, AStarData *PathData){
    unsigned temp;
    //Swap the indices of child and parent until the correct order is reached.
    while((q!=NULL) && (q->parent!=NULL) && (PathData[q->index].f < PathData[q->parent->index].f)){
        temp = q->index;
        q->index = q->parent->index;
        q->parent->index = temp;
        q = q->parent;
    }
    return;
}

