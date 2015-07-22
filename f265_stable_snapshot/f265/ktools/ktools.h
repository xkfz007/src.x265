#ifndef __KTOOLS_H__
#define __KTOOLS_H__

#include <stddef.h>
#include <stdint.h>
#include <stdarg.h>
#include <assert.h>

// Cleaner function interface.
typedef void (*kcleaner)(void *);

// malloc() and free() wrappers. The malloc functions provided must not return
// if the memory allocation fails.
typedef void* (*kmalloc_func)(size_t size);
typedef void* (*kcalloc_func)(size_t size);
typedef void* (*krealloc_func)(void *ptr, size_t size);
typedef void (*kfree_func)(void *ptr);
extern kmalloc_func kmalloc;
extern kcalloc_func kcalloc;
extern krealloc_func krealloc;
extern kfree_func kfree;

#ifdef __cplusplus
extern "C" {
#endif

void ktools_init(kmalloc_func m, kcalloc_func c, krealloc_func r, kfree_func f);

#ifdef __cplusplus
};
#endif

#endif

