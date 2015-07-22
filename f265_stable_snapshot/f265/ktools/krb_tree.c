#include "krb_tree.h"

/* This method rotates the tree left at the node specified. */
static void krb_tree_left_rotate(krb_tree *self, struct krb_node *node) {

    /* Swap subtree and link y's parent. */
    struct krb_node *y = node->right;
    node->right = y->left;

    /* y->left might be nil. Thus, we might be changing the nil node's parent
     * here. The insert and remove code is aware of this.
     */
    y->left->parent = node;
    y->parent = node->parent;

    /* Node is the root. */
    if (node == self->root_node) {
        self->root_node = y;
    }

    /* Node is the left node of its parent. */
    else if (node == node->parent->left) {
        node->parent->left = y;
    }

    /* Node is the right node of its parent. */
    else {
        node->parent->right = y;
    }

    /* Link node and y. */
    y->left = node;
    node->parent = y;

    /* Update sizes. */
    y->size = node->size;
    node->size = node->left->size + node->right->size + 1;
}

/* This method rotates the tree right at the node specified. */
static void krb_tree_right_rotate(krb_tree *self, struct krb_node *node) {

    /* Swap subtree and link y's parent. */
    struct krb_node *y = node->left;
    node->left = y->right;

    y->right->parent = node;
    y->parent = node->parent;

    /* Node is the root. */
    if (node == self->root_node) {
        self->root_node = y;
    }

    /* Node is the right node of its parent. */
    else if (node == node->parent->right) {
        node->parent->right = y;
    }

    /* Node is the left node of its parent. */
    else {
        node->parent->left = y;
    }

    /* Link node and y. */
    y->right = node;
    node->parent = y;

    /* Update sizes. */
    y->size = node->size;
    node->size = node->left->size + node->right->size + 1;
}

/* Helper method for krb_tree_reset(). */
static void krb_tree_reset_helper(struct krb_node *node, struct krb_node *nil) {
    
    /* Destroy left and right subtrees. */
    if (node->left != nil) krb_tree_reset_helper(node->left, nil);
    if (node->right != nil) krb_tree_reset_helper(node->right, nil);

    /* Delete this node. */
    krb_node_destroy(node);
}

/* Helper method for krb_tree_check_consistency(). */
static int krb_tree_check_consistency_helper(krb_tree *self, struct krb_node *node, struct krb_node *nil) {
    int left_black, right_black;
    
    /* Nil nodes are black. */
    if (node == nil) return 1;

    /* Node is red. Its children must be black. */
    if (node->color) {
        assert(node->left->color == 0 && node->right->color == 0);
    }

    /* Check chaining. */
    assert(node->left == nil || node->left->parent == node);
    assert(node->right == nil || node->right->parent == node);

    /* Check keys consistency. */
    assert(node->left == nil || self->cmp_func(node->key, node->left->key) > 0);
    assert(node->right == nil || self->cmp_func(node->key, node->right->key) < 0);

    /* Check size. */
    assert(node->size == node->left->size + node->right->size + 1);

    /* Check black-height. Line 'right_black += 0;' is used to shut up gcc. */
    left_black = krb_tree_check_consistency_helper(self, node->left, nil);
    right_black = krb_tree_check_consistency_helper(self, node->right, nil);
    right_black += 0;
    assert(left_black == right_black);

    /* Return black-height of this subtree. */
    if (node->color) return left_black;
    return (left_black + 1);
}

/* This method returns the node corresponding to the key specified, or NULL if
 * the key cannot be found.
 */
struct krb_node * krb_tree_get_node(krb_tree *self, void *key) {
    struct krb_node *node, *nil; 

    /* There is no root. */
    if (self->root_node == NULL) return NULL;

    /* Search for the key, starting at the root. */
    node = self->root_node;
    nil = self->root_node->parent;

    while (1) {
        int order = self->cmp_func(key, node->key);

        /* Go left. */
        if (order < 0) {
            if (node->left == nil) return NULL;
            node = node->left;
        }

        /* Go right. */
        else if (order > 0) {
            if (node->right == nil) return NULL;
            node = node->right;
        }

        /* We found it. */
        else return node;
    }
}

/* This method returns the node corresponding to the index specified (from 0 to
 * size - 1). The method assumes that the index is valid.
 */
struct krb_node * krb_tree_get_node_by_index(krb_tree *self, int index) {
    struct krb_node *node = self->root_node;

    assert(index >= 0 && index <= krb_tree_size(self) - 1);

    while (1) {

        /* We found it. */
        if (index == node->left->size) return node;

        /* The node we look for is at the left. */
        else if (index < node->left->size) {
            node = node->left;
        }

        /* The node we look for is at the right. */
        else {
            index -= node->left->size + 1;
            node = node->right;
        }
    }
}

/* This method returns the specified node's index. */
int krb_tree_get_node_index(krb_tree *self, struct krb_node *node) {
    int sum = 0;
    struct krb_node *nil = self->root_node->parent;

    while (1) {

        /* We reached the top. */
        if (node->parent == nil) {
            return sum;
        }

        /* We're the right child of our parent. Increment sum. */
        else if (node == node->parent->right) {
            sum += node->parent->left->size + 1;
        }

        /* Work our way toward the top. */
        node = node->parent;
    }
}

/* This method returns the successor of the node specified. The method assumes
 * that such a successor exists.
 */
struct krb_node * krb_tree_get_successor(krb_tree *self, struct krb_node *node) {
    struct krb_node *nil = self->root_node->parent;

    /* The node has a right subtree, get the minimum node in it. */
    if (node->right != nil) {
        node = node->right;
        while (node->left != nil) node = node->left;
        return node;
    }

    /* Get the first parent for which the current node is a left node. */
    assert(node != self->root_node);

    while (node == node->parent->right) {
        node = node->parent;
        assert(node != self->root_node);
    }

    return node->parent;
}

/* This method returns an iterator pointer. Use it to iterate over all the nodes
 * of the tree.
 */
struct krb_node * krb_tree_iter_start(krb_tree *self) {

    /* We use a trick here: we make nil->right points to the root. That way,
     * when the user calls krb_tree_get_successor(), he gets the first node of
     * the tree. This code is safe even if it is executed concurrently
     * (read-only access, of course).
     */
    if (self->root_node == NULL) return NULL;
    self->root_node->parent->right = self->root_node;
    return self->root_node->parent;
}

/* This method makes the iterator pointer points to the next node in the tree.
 * For convenience, the method also returns the value pointer of the node. It is
 * safe to call this method even if some of the keys in the tree are invalid,
 * e.g. if you freed the pointers to the objects used as the keys. Be careful
 * not to iterate past the end of the tree.
 */
void * krb_tree_iter_next(krb_tree *self, struct krb_node **iter_handle) {
    assert(*iter_handle != NULL);
    *iter_handle = krb_tree_get_successor(self, *iter_handle);
    return (*iter_handle)->value;
}

/* This method adds a key-value pair to the tree. If the pair already exists,
 * the method replaces it. The node is returned.
 */
struct krb_node* krb_tree_add(krb_tree *self, void *key, void *value) {

    /* Get the node, if any. */
    struct krb_node *node = krb_tree_get_node(self, key);

    /* The pair already exists. */
    if (node != NULL) {
        node->key = key;
        node->value = value;
    }

    /* Add the pair. */
    else {
        node = krb_tree_add_fast(self, key, value);
    }
    
    return node;
}

/* This method adds a key-value pair in the tree. The method assumes that the
 * key is not already in the tree. The node is returned.
 */
struct krb_node* krb_tree_add_fast(krb_tree *self, void *key, void *value) {	
    struct krb_node *ret_node, *node, *nil;

    /* There is no root. Create it (with nil node) and return. */
    if (self->root_node == NULL) {
        self->root_node = krb_node_new();
        nil = krb_node_new();

        /* nil is black and has size 0. The other fields are best left
         * uninitialized (so valgrind can trap incorrect accesses).
         */
        nil->color = 0;
        nil->size = 0;

        /* The root is black, and its parent, left and right pointers point to
         * the nil node.
         */
        krb_node_set(self->root_node, key, value, nil, nil, nil, 0);

        return self->root_node;
    }

    /* Find position to insert, starting at the root. */
    node = self->root_node;
    nil = self->root_node->parent;

    while (1) {

        /* Increment the size of the current node. */
        node->size++;

        /* Compare the key to add with the key of the current node. */
        int order = self->cmp_func(key, node->key);

        /* Go left. */
        if (order < 0) {
            if (node->left == nil) {
                node->left = krb_node_new();
                krb_node_set(node->left, key, value, node, nil, nil, 1);
                node = node->left;
                break;
            }

            node = node->left;
        }

        /* Go right. */
        else {
            assert(order != 0);

            if (node->right == nil) {
                node->right = krb_node_new();
                krb_node_set(node->right, key, value, node, nil, nil, 1);
                node = node->right;
                break;
            }

            node = node->right;
        }
    }
    
    /* Remember the node to return. */
    ret_node = node;

    /* Restore red-black properties. */
    while (node->parent->color) {

        /* The parent node of this node is the left node of its parent. */
        if (node->parent == node->parent->parent->left) {

            /* Get the uncle of this node. */
            struct krb_node *uncle = node->parent->parent->right;

            /* Case 1: uncle is red. Color the uncle's parent red and keep
             * going.
             */
            if (uncle->color) {
                node->parent->color = 0;
                uncle->color = 0;

                node = node->parent->parent;
                node->color = 1;
            }

            /* Case 2 and 3: uncle is black. Fix and stop. */
            else {

                /* Case 2. Node is the right child of its parent. We'll perform
                 * a swap to fall in case 3.
                 */
                if (node == node->parent->right) {
                    node = node->parent;
                    krb_tree_left_rotate(self, node);
                }

                /* Case 3. Node is the left child of its parent. Color correctly
                 * (except root perhaps).
                 */
                node->parent->color = 0;
                node->parent->parent->color = 1;
                krb_tree_right_rotate(self, node->parent->parent);
                break;
            }
        }

        /* The parent node of this node is the right node of its parent. This
         * code is symmetric with above.
         */
        else {
            /* Get the uncle of this node. */
            struct krb_node *uncle = node->parent->parent->left;

            /* Case 1: uncle is red. Color the uncle's parent red and keep
             * going.
             */
            if (uncle->color) {
                node->parent->color = 0;
                uncle->color = 0;

                node = node->parent->parent;
                node->color = 1;
            }

            /* Case 2 and 3: uncle is black. Fix and stop. */
            else {

                /* Case 2. Node is the left child of its parent. We'll perform a
                 * swap to fall in case 3.
                 */
                if (node == node->parent->left) {
                    node = node->parent;
                    krb_tree_right_rotate(self, node);
                }

                /* Case 3. Node is the right child of its parent. Color
                 * correctly (except root perhaps).
                 */
                node->parent->color = 0;
                node->parent->parent->color = 1;
                krb_tree_left_rotate(self, node->parent->parent);
                break;
            }
        }
    }

    /* We have to color the root black in case it was colored red. */
    self->root_node->color = 0;
    
    return ret_node;
}

/* This method removes the node corresponding to the key specified, if any. If
 * the node exists, its value is returned, otherwise NULL is returned.
 */
void * krb_tree_remove(krb_tree *self, void *key) {
    void *value;

    /* Get the node. */
    struct krb_node *node = krb_tree_get_node(self, key);

    /* Not found. */
    if (node == NULL) return NULL;

    /* Get the data. */
    value = node->value;

    /* Remove the node. */
    krb_tree_remove_node(self, node);

    return value;
}

/* This method removes the node corresponding to the index specified (from 0 to
 * size - 1). The method assumes that the index is valid. The value of the node
 * is returned.
 */
void * krb_tree_remove_by_index(krb_tree *self, int index) {

    /* Get the node and the data. */
    struct krb_node *node = krb_tree_get_node_by_index(self, index);
    void *value = node->value;

    /* Remove the node. */
    krb_tree_remove_node(self, node);

    return value;
}

/* This method removes a node from the tree. */
void krb_tree_remove_node(krb_tree *self, struct krb_node *node) {

    /* Get the nil node and get some temporary pointers. */
    struct krb_node *nil = self->root_node->parent;
    struct krb_node *y, *x;

    /* Get the node that we will really remove to fill the gap left by the key
     * removal. If this node has one or zero child, we use this node. Otherwise,
     * we use this node's successor.
     */
    if (node->left == nil || node->right == nil) {
        y = node;
    }

    else {
        /* The successor is guaranteed to have at most one child. */
        y = krb_tree_get_successor(self, node);
    }

    /* Decrement the size of the nodes from y->parent up to the root. */
    x = y->parent;

    while (x != nil) {
        x->size--;
        x= x->parent;
    }

    /* Get the child of the node we're about to remove (or nil). */
    if (y->left != nil) {
        x = y->left;
    }

    else {
        x = y->right;
    }

    /* Make the child (or nil) node points to the parent of the node we're about
     * to remove.
     */
    x->parent = y->parent;

    /* We're removing the root. */
    if (y->parent == nil) {

        /* If there's no child, we delete the nil and root nodes and clear the
         * root node pointer.
         */
        if (x == nil) {
            krb_node_destroy(nil);
            krb_node_destroy(self->root_node);
            self->root_node = NULL;
            return;
        }

        /* Otherwise, we make the child the root. */
        else {
            self->root_node = x;
            assert(self->root_node->parent == nil);
        }
    }

    /* The node we're removing is the left child of its parent. */
    else if (y == y->parent->left) {
        y->parent->left = x;
    }

    /* The node we're removing is the right child of its parent. */
    else {
        y->parent->right = x;
    }

    /* Copy the key and value if the node we're actually removing is not the
     * node we wanted to remove initially.
     */
    if (y != node) {
        node->key = y->key;
        node->value = y->value;
    }

    /* If the node we removed is black, the red-black properties are no longer
     * respected. Fix them.
     */
    if (y->color == 0) {

        /* This loop finishes any time we do not enter directly in case 2. */
        while (x != self->root_node && x->color == 0) {

            /* X is the left child of its parent. */
            if (x == x->parent->left) {

                /* Get X's parent and uncle. */
                struct krb_node *parent = x->parent;
                struct krb_node *uncle = parent->right;

                /* Case 1: the uncle is red. We'll exit the loop through cases
                 * 2, 3 or 4.
                 */
                if (uncle->color) {
                    uncle->color = 0;
                    parent->color = 1;
                    krb_tree_left_rotate(self, parent);
                    uncle = parent->right;
                }

                /* Case 2: Both the children of the uncle are black. */
                if (uncle->left->color == 0 && uncle->right->color == 0) {	
                    uncle->color = 1;
                    x = parent;
                }

                /* Case 3 and 4. */
                else {

                    /* Case 3: the right child of the uncle is black. */
                    if (uncle->right->color == 0) {
                        uncle->left->color = 0;
                        uncle->color = 1;
                        krb_tree_right_rotate(self, uncle);
                        uncle = parent->right;
                    }

                    /* Case 4: the right child of the uncle is red. */
                    uncle->color = parent->color;
                    parent->color = 0;
                    uncle->right->color = 0;
                    krb_tree_left_rotate(self, parent);
                    x = self->root_node;
                }
            }

            /* X is the right child of its parent. */
            else {
                /* Get X's parent and uncle. */
                struct krb_node *parent = x->parent;
                struct krb_node *uncle = parent->left;

                /* Case 1: the uncle is red. We'll exit the loop through cases
                 * 2, 3 or 4.
                 */
                if (uncle->color) {
                    uncle->color = 0;
                    parent->color = 1;
                    krb_tree_right_rotate(self, parent);
                    uncle = parent->left;
                }

                /* Case 2: Both the children of the uncle are black. */
                if (uncle->right->color == 0 && uncle->left->color == 0) {
                    uncle->color = 1;
                    x = parent;
                }

                /* Case 3 and 4. */
                else {

                    /* Case 3: the left child of the uncle is black. */
                    if (uncle->left->color == 0) {
                        uncle->right->color = 0;
                        uncle->color = 1;
                        krb_tree_left_rotate(self, uncle);
                        uncle = parent->left;
                    }

                    /* Case 4: the left child of the uncle is red. */
                    uncle->color = parent->color;
                    parent->color = 0;
                    uncle->left->color = 0;
                    krb_tree_right_rotate(self, parent);
                    x = self->root_node;
                }
            }
        }

        /* Color X black. */
        x->color = 0;
    }

    /* Delete the node we actually removed. */
    assert(y != nil && self->root_node->parent == nil && self->root_node != y);
    krb_node_destroy(y);
}

/* This method destroys all the nodes in the tree. */
void krb_tree_reset(krb_tree *self) {
    struct krb_node *nil;

    if (self->root_node == NULL) return;

    /* Get the nil node. */
    nil = self->root_node->parent;

    /* Clear the tree. */
    krb_tree_reset_helper(self->root_node, nil);

    /* Destroy the nil node. */
    krb_node_destroy(nil);

    /* Clear the root pointer. */
    self->root_node = NULL;
}

/* This method verifies the internal consistency of the tree. Call this is you
 * suspect the tree is corrupted.
 */
void krb_tree_check_consistency(krb_tree *self) {
    struct krb_node *nil;
    
    if (self->root_node == NULL) return;
    
    /* Get nil node. Verify color and size. */
    nil = self->root_node->parent;
    assert(nil->color == 0);
    assert(nil->size == 0);

    /* Check the tree. */
    krb_tree_check_consistency_helper(self, self->root_node, nil);
}

