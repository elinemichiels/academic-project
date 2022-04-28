#include "linked_list.h"

int IsEmpty(PriorityQueue Pq){
    return ( Pq == NULL );
}

unsigned extract_min(PriorityQueue *Pq, AStarData *PathData){
    PriorityQueue first = *Pq;
    unsigned index = first->index;
    *Pq = (*Pq)->next;
    free(first);
    return index;
}


void add_with_priority(unsigned index, PriorityQueue *Open, AStarData *PathData){
    QueueElement *aux = (QueueElement *) malloc(sizeof(QueueElement));
    if(aux == NULL)
        ExitError("Allocation of a new Queue element failed", 61);

    aux->index = index;
    double costv = PathData[index].f;
    PathData[index].isOpen = 1;

    //Special case if the new element has to be inserted at the start of the list
    if (*Open == NULL || !(costv > PathData[(*Open)->index].f)) {
        aux->next = *Open;
        *Open = aux;
        return;
    }
    register QueueElement * q;

    //Skip to the position where the element has to be inserted.
    for(q = *Open; q->next && PathData[q->next->index].f < costv; q = q->next );
    aux->next = q->next; q->next = aux;
    return;
}

void requeue_with_priority(unsigned index, PriorityQueue *Open, AStarData *PathData){
    if((*Open)->index == index)
        return;

    double costv = PathData[index].f;
    //Special case if the new element has to be moved to the start of the list.
    if(!(costv > PathData[(*Open)->index].f)){
        register QueueElement *prepv;
        //Skip to the element before the one that has the considered index.
        for(prepv = *Open; prepv->next->index != index; prepv = prepv->next);
        QueueElement *swap = *Open;
        *Open=prepv->next; prepv->next=prepv->next->next; (*Open)->next=swap;
        return;
    }

    register QueueElement *q, *prepv;
    //Skip to the element before the one that has the considered index.
    for(q = *Open; PathData[q->next->index].f < costv; q = q->next );
    if(q->next->index == index)
        return;
    //Skip to the position before the one where the element has to be moved to.
    for(prepv = q->next; prepv->next->index != index; prepv = prepv->next);
    QueueElement *pv = prepv->next;
    prepv->next = pv->next; pv->next = q->next; q->next = pv;
    return;
}