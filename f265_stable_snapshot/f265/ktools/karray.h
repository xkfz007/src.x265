#ifndef __KARRAY_H__
#define __KARRAY_H__

#include "ktools.h"

typedef struct karray {
    
    /* The allocated array size (in terms of elements). */
    size_t alloc_size;
    
    /* The number of elements in the array. */ 
    size_t size;
    
    /* The element array. */
    void **data;
} karray;

#ifdef __cplusplus
extern "C" {
#endif
karray * karray_new();
void karray_destroy(karray *self);
void karray_init(karray *self);
void karray_init_karray(karray *self, karray *init_array);
void karray_reset(karray *self);
void karray_clean(karray *self);
void karray_grow(karray *self, size_t min_len);
void karray_push(karray *self, void *elem);
void *karray_pop(karray *self);
void karray_set(karray *self, int pos, void *elem);
void karray_assign_karray(karray *self, karray *assign_array);
void karray_append_karray(karray *self, karray *append_array);
void karray_cleanse(karray *self, kcleaner cleaner);
#ifdef __cplusplus
};
#endif

#endif

