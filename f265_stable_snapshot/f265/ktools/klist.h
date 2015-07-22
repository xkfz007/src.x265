#ifndef __KLIST_H__
#define __KLIST_H__

#include "ktools.h"

// Node of the doubly-linked list.
typedef struct klist_node
{
    // Pointer to the previous node.
    struct klist_node *next;

    // Pointer to the next node.
    struct klist_node *prev;

    // Pointer to the data contained in this node.
    void *data;
} klist_node;

// Doubly-linked list class.
typedef struct klist
{
    // Number of elements that are present in the list.
    size_t size;

    // The root node, i.e the first and last element of the doubly-linked
    // list. If the list is empty, the 'next' and 'prev' pointers of the
    // root node point to the root node itself.
    klist_node root;
} klist;

#ifdef __cplusplus
extern "C" {
#endif
void klist_clear(klist *self);
void klist_assign_klist(klist *self, klist *other);
void klist_cleanse(klist *self, kcleaner cleaner);
klist_node* klist_add_after(klist *self, klist_node *node, void *data);
klist_node* klist_add_before(klist *self, klist_node *node, void *data);
void* klist_remove(klist *self, klist_node *node);
void* klist_iter_next(klist *self, klist_node **iter_handle);
    
static inline void klist_init(klist *self)
{
    self->root.prev = self->root.next = &self->root;
    self->size = 0;
}

/* This method returns true if the list is empty. */
static inline int klist_is_empty(klist *self)
{
    return (self->root.next == &self->root);
}

/* This method adds a new node at the beginning of the list. The new node is
 * returned.
 */
static inline klist_node* klist_prepend(klist *self, void *data)
{
    return klist_add_after(self, &self->root, data);
}

/* Same as above, but the new node is added at the end of the list. */
static inline klist_node* klist_append(klist *self, void *data)
{
    return klist_add_before(self, &self->root, data);
}

/* This method returns the first node, or NULL if there is none. */
static inline klist_node* klist_get_first_node(klist *self)
{
    return klist_is_empty(self) ? NULL : self->root.next;
}

/* This method returns the last node, or NULL if there is none. */
static inline klist_node* klist_get_last_node(klist *self)
{
    return klist_is_empty(self) ? NULL : self->root.prev;
}

/* This method removes the first node and returns the data pointer of
 * that node. If there is no such node, NULL is returned.
 */
static inline void* klist_remove_first_node(klist *self)
{
    klist_node *node = klist_get_first_node(self);
    if (node == NULL) return NULL;
    return klist_remove(self, node);
}

/* Complement of remove_first_node(). */
static inline void* klist_remove_last_node(klist *self)
{
    klist_node *node = klist_get_last_node(self);
    if (node == NULL) return NULL;
    return klist_remove(self, node);
}

/* This method returns an iterator node. Use it to iterate over all the
 * nodes in the list. The returned pointer is NULL if the list is empty.
 */
static inline klist_node* klist_iter_start(klist *self)
{
    if (! self->size) return NULL;
    return self->root.next;
}

#ifdef __cplusplus
};
#endif

#endif

