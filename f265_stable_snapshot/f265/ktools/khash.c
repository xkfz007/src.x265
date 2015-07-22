#include "khash.h"
#include "kutil.h"

/* Prime table used to grow the hash. */
static int khash_prime_table[] = {
    23, 53, 97, 193, 389, 769, 1543, 3079, 6151, 12289, 24593, 49157, 98317,
    196613, 393241, 786433, 1572869, 3145739, 6291469, 12582917, 25165843, 50331653,
    100663319, 201326611, 402653189, 805306457, 1610612741, 4294967291u
};

/* The proportion of the cells that must be used before we decide to grow the 
 * hash table.
 */
#define KHASH_FILL_THRESHOLD 0.7

/* This function assumes that pos1 comes before pos2. It returns the distance between the two.
 * It is used internally.
 * Arguments:
 * Position 1.
 * Position 2.
 * Size of the hash (for modulus operation).
 */
static inline int khash_dist(int pos1, int pos2, int base) {
    return (pos1 <= pos2) ? pos2 - pos1 : base - (pos1 - pos2);
}

khash * khash_new() {
    khash *self = (khash *) kmalloc(sizeof(khash));
    khash_init(self);
    return self;

}

void khash_destroy(khash *self) {
    if (self)
        khash_clean(self);
    kfree(self);
}

/* This function initializes the hash. By default keys are hashed and compared
 * by integers.
 */
void khash_init(khash *self) {
    khash_init_func(self, kutil_pointer_hc, kutil_int32_eq);
}

void khash_init_func(khash *self, uint32_t (*key_func) (void *), int (*eq_func) (void *, void *)) {
    self->key_func = key_func;
    self->eq_func = eq_func;
    self->size = 0;
    self->used_limit = (int) (11 * KHASH_FILL_THRESHOLD);
    self->next_prime_index = 0;
    self->alloc_size = 11;
    self->cell_array = (khash_cell *) kcalloc(self->alloc_size * sizeof(khash_cell));
}

/* This function sets the key hash and compare functions used to hash objects. */
void khash_set_func(khash *self, uint32_t (*key_func) (void *), int (*eq_func) (void *, void *)) {
    self->key_func = key_func;
    self->eq_func = eq_func;
}

/* This function frees the content of the hash. */
void khash_clean(khash *self) {
    if (self == NULL)
    	return;
    
    kfree(self->cell_array);
}

/* This function increases the size of the hash. */
void khash_grow(khash *self) {
    int index;
    int new_alloc_size;
    khash_cell *new_cell_array;
    
    /* Get the new size. */
    new_alloc_size = khash_prime_table[self->next_prime_index];
    self->next_prime_index++;
    self->used_limit = (int) (new_alloc_size * KHASH_FILL_THRESHOLD);

    /* Allocate the new table. */
    new_cell_array = (khash_cell *) kcalloc(new_alloc_size * sizeof(khash_cell));
    
    /* Copy the elements. */
    for (index = 0; index < self->alloc_size; index++) {
        void *key = self->cell_array[index].key;
        
        if (key != NULL) {
            /* Put the key at the right place. */
            int i = self->key_func(key) % new_alloc_size;
            
            while (1) {
                /* OK, found an empty slot. */
                if (new_cell_array[i].key == NULL) {
                    /* Set the key / value pair. */
                    new_cell_array[i].key = key;
                    new_cell_array[i].value = self->cell_array[index].value;
                    break;
    	    	}

                i = (i + 1) % new_alloc_size;
    	    }
    	}
    }
    
    /* Free the old table. */
    kfree(self->cell_array);
    
    /* Assign the new table and the new size. */
    self->alloc_size = new_alloc_size;
    self->cell_array = new_cell_array;
}    

/* This function returns the position corresponding to the key in the hash, or -1
 * if it is not there.
 * Arguments:
 * Key to locate.
 */
static int khash_locate_key(khash *self, void *key) {
    int index;
    assert(key != NULL);
    
    /* Go where the key should be close. */
    index = self->key_func(key) % self->alloc_size;
        
    while (1) {
        /* Empty slot (the key is not there). */
        if (self->cell_array[index].key == NULL) {
            return -1;
        }
            
        /* Same slot. Must compare key values. */
        else if (self->eq_func(self->cell_array[index].key, key))
            return index;
        
        /* Not the same key. Advance to the next position, possibly looping back
	 * to the beginning.
	 */
        index = (index + 1) % self->alloc_size;
    }
}

/* This function returns true if the key is in the hash.
 * Arguments:
 * Key to look for.
 */
int khash_exist(khash *self, void *key) {
    return (khash_locate_key(self, key) != -1);
}

/* This function adds a key / value pair in the hash. If the key is already
 * present, it will be replaced. The key cannot be NULL.
 * Arguments:
 * Key to add.
 * Value to add.
 */
int khash_add(khash *self, void *key, void *value) {
    int index;
    assert(key != NULL);
    
    /* Grow hash if it is too small. */
    if (self->size >= self->used_limit)
        khash_grow(self);

    /* Go where the key should be close. */
    index = self->key_func(key) % self->alloc_size;

    while (1) {
    
        /* Empty slot. It's a new key. */
        if (self->cell_array[index].key == NULL) {
	
            /* Increment use count. */
            self->size++;
            
            /* Set the key / value pair. */
            self->cell_array[index].key = key;
            self->cell_array[index].value = value;
            return 0;
    	}
        
        /* Not the same key. Advance to the next position, possibly looping back
	 * to the beginning.
	 */
        index = (index + 1) % self->alloc_size;
    }
}

/* This function removes the key / value pair from the hash (if any).
 * Arguments:
 * Key to remove.
 */
int khash_remove(khash *self, void *key) {
    int index = khash_locate_key(self, key);
    int gap_position = index;
    int scanned_pos = index;
    
    /* Key is not present in the hash. */
    if(index == -1)
        return -1;
    
    /* We must ensure that other keys remain contiguous when we remove a key.
     * The situation where we need to move a key up is when the position of the
     * key given by key_func() is farther than the actual position of the key in
     * the hash, counting  from the position of the gap. Picture:
     *
     * The key wants to be here. (Distance between gap and pos wanted is far.)
     * A gap is here.            (So we move the key here.)
     * The key is here.          (Distance between gap and key is close.)
     *
     * In this situation, we don't move:
     * The gap is here.        
     * The key wants to be here. (Distance between gap and pos wanted is short.)
     * The key is here.          (Distance between gap and key is far.)
     *
     * If the gap position matches the wanted pos, we must move the key to fill
     * the gap.
     *
     * So here's the plan:
     * First we locate the key to remove. Removing it causes a gap. We start
     * scanning the keys coming next. If we meet a NULL, we're done. If we meet
     * a key, we check where it wants to be. If it wants to be before the gap,
     * we move it there. Then the gap is now at the position of the key we
     * moved, and we continue at the next position. Otherwise, we just continue 
     * with the next position.
     */
    while (1) {
    	int wanted_pos_dist;
	int key_dist;
	
        /* Scan the next position. */
        scanned_pos = (scanned_pos + 1) % self->alloc_size;

        /* We're done. Just set the gap to NULL. */
        if (self->cell_array[scanned_pos].key == NULL) {
            self->cell_array[gap_position].key = NULL;
            break;
    	}

        /* Calculate the distances. */
        wanted_pos_dist = khash_dist(gap_position,
	    	    	    	     self->key_func(self->cell_array[scanned_pos].key) % self->alloc_size,
				     self->alloc_size);
        key_dist = khash_dist(gap_position, scanned_pos, self->alloc_size);    

        /* Situations where we must move key (and value). */
        if (wanted_pos_dist > key_dist || wanted_pos_dist == 0) {
            self->cell_array[gap_position].key = self->cell_array[scanned_pos].key;
            self->cell_array[gap_position].value = self->cell_array[scanned_pos].value;
            gap_position = scanned_pos;
    	}
    }

    /* Decrement the usage count. */
    self->size--;
    return 0;
}

/* This function returns the value corresponding to the key, or NULL if the key
 * is not in the hash.
 */
void* khash_get(khash *self, void *key) {
    int index = khash_locate_key(self, key);
    if (index == -1) return NULL;
    return self->cell_array[index].value;
}

/* This function returns the key in the hash corresponding to the key, or NULL
 * if it is not in the hash.
 */
void* khash_get_key(khash *self, void *key) {
    int index = khash_locate_key(self, key);
    if (index == -1) return NULL;
    return self->cell_array[index].key;
}

/* This function clears all entries in the hash. */
void khash_reset(khash *self) {
    self->size = 0;
    self->used_limit = (int) (11 * KHASH_FILL_THRESHOLD);
    self->next_prime_index = 0;
    self->alloc_size = 11;
    self->cell_array = (khash_cell *) krealloc(self->cell_array,
                                               self->alloc_size * sizeof(khash_cell));
    memset(self->cell_array, 0, self->alloc_size * sizeof(khash_cell));
}

/* This function destroys the objects contained in the hash with the functions
 * specified and clears the hash. Either function may be NULL.
 */
void khash_cleanse(khash *self, kcleaner key_cleaner, kcleaner value_cleaner)
{
    int i;
    for (i = 0; i < self->alloc_size; i++)
    {
        if (self->cell_array[i].key)
        {
            if (key_cleaner) key_cleaner(self->cell_array[i].key);
            if (value_cleaner) value_cleaner(self->cell_array[i].value);
        }
    }
    
    khash_reset(self);
}

/* This function sets the key / value pair of the next element in the
 * enumeration. It is safe to call this function even if some of the keys in
 * the hash are invalid, e.g. if you freed the pointers to the objects used as
 * the keys. Be careful not to iterate past the end of the hash.
 * Arguments:
 * Pointer to iterator index, which should be initialized to -1 prior to the 
 *   first call.
 * Pointer to the location where you wish the key to be set; can be NULL.
 * Pointer to the location where you wish the value to be set; can be NULL.
 */
void khash_iter_next_item(khash *self, int *index, void **key_handle, void **value_handle) {
    
    for ((*index)++; *index < self->alloc_size; (*index)++) {
        if (self->cell_array[*index].key != NULL) {
            if (key_handle != NULL)
                *key_handle = self->cell_array[*index].key;
            
            if (value_handle != NULL)
                *value_handle = self->cell_array[*index].value;
    	    
            return;
        }
    }
    
    assert(0);
}

/* Same as above, except that it returns only the next key. */
void * khash_iter_next_key(khash *self, int *index) {

    for ((*index)++; *index < self->alloc_size; (*index)++)
        if (self->cell_array[*index].key != NULL)
            return self->cell_array[*index].key;
    
    assert(0);
    return NULL;
}

/* Same as above, except that it returns only the next value. */
void * khash_iter_next_value(khash *self, int *index) {
    
    for ((*index)++; *index < self->alloc_size; (*index)++)
        if (self->cell_array[*index].key != NULL)
            return self->cell_array[*index].value;
    
    assert(0);
    return NULL;
}

