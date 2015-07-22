#include "ktools.h"

kmalloc_func kmalloc;
kcalloc_func kcalloc;
krealloc_func krealloc;
kfree_func kfree;

void ktools_init(kmalloc_func m, kcalloc_func c, krealloc_func r, kfree_func f)
{
    kmalloc = m;
    kcalloc = c;
    krealloc = r;
    kfree = f;
}

