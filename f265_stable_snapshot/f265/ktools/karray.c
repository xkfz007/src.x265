#include "karray.h"
#include "kutil.h"

karray * karray_new() {
    karray *self = (karray *) kmalloc(sizeof(karray));
    karray_init(self);
    return self;
}

void karray_destroy(karray *self) {
    karray_clean(self);
    kfree(self);
}

void karray_init(karray *self) {
    self->alloc_size = 0;
    self->size = 0;
    self->data = NULL;
}

void karray_init_karray(karray *self, karray *init_array) {
    self->alloc_size = init_array->size;
    self->size = init_array->size;
    self->data = (void**)kmalloc(self->size * sizeof(void *));
    memcpy(self->data, init_array->data, self->size * sizeof(void *));
}

void karray_reset(karray *self) {
    self->size = 0;
}

void karray_clean(karray *self) {
    if (self == NULL) return;
    kfree(self->data);
}

/* This function increases the size of the array so that it may contain at least
 * 'min_len' elements.
 */
void karray_grow(karray *self, size_t min_len) {
    if (min_len > self->alloc_size) {
    
        /* Compute the snapped size for a given requested size. By snapping to powers
         * of 2 like this, repeated reallocations are avoided.
         */
        if (min_len < 4) {
            self->alloc_size = 4;
        }
        
        else {
	    self->alloc_size = kutil_next_power_of_2(min_len);
        }
        
        assert(self->alloc_size >= min_len);    
        self->data = (void**)krealloc(self->data, self->alloc_size * sizeof(void *));
    }
}

/* This function adds an element at the end of the array. */
void karray_push(karray *self, void *elem) {
    karray_set(self, self->size, elem);
}

/* This function remove the element at the end of the array and returns it. */
void *karray_pop(karray *self) {
    self->size -= 1;
    return self->data[self->size];
}

/* This function sets an element at the specified position in the array. */
void karray_set(karray *self, int pos, void *elem) {
    assert(pos >= 0);

    /* Make sure the array is big enough. */
    karray_grow(self, pos + 1);

    /* Increase the array size to contain at least this position. */
    if (self->size < pos + 1) {
    	self->size = pos + 1;
    }

    self->data[pos] = elem;
}

/* This function assigns a karray to this array. */
void karray_assign_karray(karray *self, karray *assign_array) {
    self->size = assign_array->size;
    karray_grow(self, self->size);
    memcpy(self->data, assign_array->data, self->size * sizeof(void *));    
}

/* This function appends a karray to this array. */
void karray_append_karray(karray *self, karray *append_array) {
    karray_grow(self, self->size + append_array->size);
    memcpy(self->data + self->size, append_array->data, append_array->size * sizeof(void *));
    self->size += append_array->size;
}

/* This function destroys the objects contained in the array with the function
 * specified and clears the array.
 */
void karray_cleanse(karray *self, kcleaner cleaner)
{
    size_t i;
    for (i = 0; i < self->size; i++) cleaner(self->data[i]);
    self->size = 0;
}

