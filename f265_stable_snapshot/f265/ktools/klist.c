#include "klist.h"

/* This method deletes all the elements of the list. */
void klist_clear(klist *self)
{
    klist_node *node = self->root.next;

    while (node != &self->root)
    {
        klist_node *next = node->next;
        kfree(node);
        node = next;
    }

    self->root.prev = self->root.next = &self->root;
    self->size = 0;
}

void klist_assign_klist(klist *self, klist *other)
{
    klist_clear(self);

    self->size = other->size;
    klist_node *this_node = &self->root;
    klist_node *other_node = other->root.next;

    while (other_node != &other->root)
    {
        klist_node *new_node = (klist_node*)kmalloc(sizeof(klist_node));
        this_node->next = new_node;
        new_node->prev = this_node;
        new_node->data = other_node->data;
        this_node = new_node;
        other_node = other_node->next;
    }

    // Chain the last element.
    this_node->next = &self->root;
    self->root.prev = this_node;
}

/* This function destroys the objects contained in the list with the function
 * specified and clears the list.
 */
void klist_cleanse(klist *self, kcleaner cleaner)
{
    klist_node *node = self->root.next;

    while (node != &self->root)
    {
        klist_node *next = node->next;
        cleaner(node->data);
        kfree(node);
        node = next;
    }

    self->root.prev = self->root.next = &self->root;
    self->size = 0;
}

/* This method adds a new node after the node specified. The new node is
 * returned.
 */
klist_node* klist_add_after(klist *self, klist_node *node, void *data)
{
    klist_node *new_node = (klist_node*)kmalloc(sizeof(klist_node));
    new_node->next = node->next;
    new_node->prev = node;
    new_node->data = data;
    node->next->prev = new_node;
    node->next = new_node;
    self->size++;
    return new_node;
}

/* Same as above, but adds the new node before the pointed node. */
klist_node* klist_add_before(klist *self, klist_node *node, void *data)
{
    klist_node *new_node = (klist_node*)kmalloc(sizeof(klist_node));
    new_node->next = node;
    new_node->prev = node->prev;
    new_node->data = data;
    node->prev->next = new_node;
    node->prev = new_node;
    self->size++;
    return new_node;
}

/* This method removes the node specified from the list and returns the
 * associated data pointer.
 */
void* klist_remove(klist *self, klist_node *node)
{
    assert(node != &self->root);
    assert(node != NULL);
    node->next->prev = node->prev;
    node->prev->next = node->next;
    self->size--;
    void *data = node->data;
    kfree(node);
    return data;
}

/* This method makes the iterator pointer points to the next node in the list.
 * For convenience, the method also returns the value pointer of the current
 * node. Be careful not to iterate past the end of the list (don't call this
 * method if *iter_handle is NULL).
 */
void* klist_iter_next(klist *self, klist_node **iter_handle)
{
    assert(*iter_handle != NULL);
    void *value = (*iter_handle)->data;
    *iter_handle = (*iter_handle)->next;
    if (*iter_handle == &self->root) *iter_handle = NULL;
    return value;
}

