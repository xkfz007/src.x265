#ifndef __KRB_TREE_H__
#define __KRB_TREE_H__

#include "ktools.h"
#include "kutil.h"

/* Struct krb_tree is a red-black tree implementation. Such trees are always
 * balanced, and most operations over them are O(lg(n)) in the worst case. The
 * tree maintains a set of ordered key-value pairs. You can ask to get a node by
 * its key or its index (from 0 to size - 1). You can iterate on all the nodes
 * of the tree in O(n) time. To allow the tree to sort the keys, you must supply
 * a comparison function. This function receives two keys as its parameters and
 * returns the order of the keys. The function should returns -1 if the first
 * key comes before the second key, 1 if the first key comes after the second
 * key, and 0 if the keys are equal. By default, the code compares the keys by
 * address.
 *
 * Implementation note:
 * An effort was made to minimize the memory usage and initialization time of
 * empty trees. An empty tree contains only two pointers, the NULL root pointer
 * and the pointer to a function that compares keys by address. The nil node
 * required for the unification of the code is created when the root is created,
 * and destroyed when the root is destroyed. Furthermore, the code reserves the
 * nil node's right pointer to unify iterations: iter_start() uses the nil node
 * as the first node returned to the user (or NULL if the root does not exist).
 */

/* Red-black tree node. */
struct krb_node {
    
    /* Pointer to the parent, left and right nodes. */
    struct krb_node *parent;
    struct krb_node *left;
    struct krb_node *right;

    /* Key contained in this node. You can modify this directly, but if you do
     * make sure the order of the keys in the tree is still respected.
     */
    void *key;

    /* The value corresponding to the key. */
    void *value;

    /* Size of this node's subtree. */
    int size;

    /* Color of this node. */
    char color;
};

/* Red-black tree. */
typedef struct krb_tree {

    /* The pointer to the root node, or NULL if there is no root. Note that
     * root_node->parent is the nil node.
     */
    struct krb_node *root_node;

    /* The comparison function. */
    int (*cmp_func) (void *, void *);

} krb_tree;

#ifdef __cplusplus
extern "C" {
#endif
struct krb_node * krb_tree_get_node(krb_tree *self, void *key);
struct krb_node * krb_tree_get_node_by_index(krb_tree *self, int index);
int krb_tree_get_node_index(krb_tree *self, struct krb_node *node);
struct krb_node * krb_tree_get_successor(krb_tree *self, struct krb_node *node);
struct krb_node * krb_tree_iter_start(krb_tree *self);
void * krb_tree_iter_next(krb_tree *self, struct krb_node **iter_handle);
struct krb_node* krb_tree_add(krb_tree *self, void *key, void *value);
struct krb_node* krb_tree_add_fast(krb_tree *self, void *key, void *value);
void * krb_tree_remove(krb_tree *self, void *key);
void * krb_tree_remove_by_index(krb_tree *self, int index);
void krb_tree_remove_node(krb_tree *self, struct krb_node *node);
void krb_tree_reset(krb_tree *self);
void krb_tree_check_consistency(krb_tree *self);

static inline struct krb_node * krb_node_new() {
    return (struct krb_node *) kmalloc(sizeof(struct krb_node));
}

static inline void krb_node_destroy(struct krb_node *self) {
    if (self) {
        kfree(self);
    }
}

/* This function sets the fields of a node to the specified values, except the
 * size which is set to 1.
 */
static inline void krb_node_set(struct krb_node *self,
                                void *key,
                                void *value,
                                struct krb_node *parent,
                                struct krb_node *left,
                                struct krb_node *right,
                                int color) {
    self->key = key;
    self->value = value;
    self->parent = parent;
    self->left = left;
    self->right = right;
    self->color = color;
    self->size = 1;
}

/* This function initializes the tree. By default keys are compared by integers. */
static inline void krb_tree_init(krb_tree *self) {
    self->root_node = NULL;
    self->cmp_func = kutil_int32_cmp;
}

/* This function initializes the tree with the comparison function specified. */
static inline void krb_tree_init_func(krb_tree *self, int (*cmp_func) (void *, void *)) {
    self->root_node = NULL;
    self->cmp_func = cmp_func;
}

static inline void krb_tree_clean(krb_tree *self) {
    krb_tree_reset(self);
}

static inline krb_tree * krb_tree_new() {
    krb_tree *self = (krb_tree *) kmalloc(sizeof(krb_tree));
    krb_tree_init(self);
    return self;
}

static inline void krb_tree_destroy(krb_tree *self) {
    if (self) {
        krb_tree_clean(self);
        kfree(self);
    }
}

/* This function sets the compare function used to compare keys. */
static inline void krb_tree_set_func(krb_tree *self, int (*cmp_func) (void *, void *)) {
    self->cmp_func = cmp_func;
}

/* This method returns the tree size. */
static inline int krb_tree_size(krb_tree *self) {
    if (self->root_node == NULL) return 0;
    return self->root_node->size;
}

/* This method returns true if the key specified exists in the tree. */
static inline int krb_tree_exist(krb_tree *self, void *key) {
    return (krb_tree_get_node(self, key) != NULL);
}

/* This method returns the value corresponding to the key specified. The method
 * assumes that the key exists.
 */
static inline void * krb_tree_get_fast(krb_tree *self, void *key) {
    struct krb_node *node = krb_tree_get_node(self, key);
    assert(node != NULL);
    return node->value;
}

/* This method returns the value corresponding to the key specified if it
 * exists. Otherwise NULL is returned.
 */
static inline void * krb_tree_get(krb_tree *self, void *key) {
    struct krb_node *node = krb_tree_get_node(self, key);
    if (node == NULL) return NULL;
    return node->value;
}

/* This method returns the value corresponding to the index specified (from 0 to
 * size - 1). The method assumes that the index is valid.
 */
static inline void * krb_tree_get_by_index(krb_tree *self, int index) {
    return krb_tree_get_node_by_index(self, index)->value;
}
#ifdef __cplusplus
};
#endif

#endif
