#ifndef __KHASH_H__
#define __KHASH_H__

#include "ktools.h"

/* A cell in the hash table: a key and its associated value. */
typedef struct khash_cell {
    void *key;
    void *value;
} khash_cell;

/* The hash itself. */
typedef struct khash {
    
    /* The hashing function used by this hash. This function takes a key object
     * as its argument and returns an integer often unique for that object. By
     * default, we use the value of the pointer as the integer.
     */
    uint32_t (*key_func) (void *);

    /* The comparison function used by this hash. This function takes two key
     * objects as its arguments and returns true if the objects are the same. By
     * default, we compare the values of the pointers.
     */
    int (*eq_func) (void *, void *);

    /* The table containing the hash cells. */
    struct khash_cell *cell_array;

    /* Size of the table. */
    int alloc_size;

    /* Number of cells used in the table. */
    int size;

    /* Maximum number of cells that can be used before the hash is expanded. */
    int used_limit;

    /* Next index in the prime table to grow the hash. */
    int next_prime_index;

} khash;

#ifdef __cplusplus
extern "C" {
#endif
khash * khash_new();
void khash_destroy(khash *self);
void khash_init(khash *self);
void khash_init_func(khash *self, uint32_t (*key_func) (void *), int (*eq_func) (void *, void *));
void khash_set_func(khash *self, uint32_t (*key_func) (void *), int (*eq_func) (void *, void *));
void khash_clean(khash *self);
void khash_grow(khash *self);
int khash_exist(khash *self, void *key);
int khash_add(khash *self, void *key, void *value);
int khash_remove(khash *self, void *key);
void* khash_get(khash *self, void *key);
void* khash_get_key(khash *self, void *key);
void khash_reset(khash *self);
void khash_cleanse(khash *self, kcleaner key_cleaner, kcleaner value_cleaner);
void khash_iter_next_item(khash *self, int *index, void **key_handle, void **value_handle);
void * khash_iter_next_key(khash *self, int *index);
void * khash_iter_next_value(khash *self, int *index);
#ifdef __cplusplus
};
#endif

#endif

