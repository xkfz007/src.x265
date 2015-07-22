#ifndef __KSTR_H__
#define __KSTR_H__

#include "ktools.h"

typedef struct kstr
{
    /* The allocated buffer size. */
    int mlen;
    
    /* The string length, not including the final '0'. */
    int slen;
    
    /* The character buffer, always terminated by a '0'.
     * Note that there may be other '0' in the string.
     */
    char *data;
} kstr;

#ifdef __cplusplus
extern "C" {
#endif
kstr* kstr_new();
kstr* kstr_new_kstr(kstr *str);
kstr* kstr_new_cstr(char *str);
void kstr_destroy(kstr *str);
void kstr_init(kstr *self);
void kstr_init_cstr(kstr *self, const char *init_str);
void kstr_init_kstr(kstr *self, kstr *init_str);
void kstr_init_buf(kstr *self, const void *buf, int buf_len);
void kstr_init_sfv(kstr *self, const char *format, va_list args);
void kstr_init_sf(kstr *self, const char *format, ...);
void kstr_clean(kstr *self);
void kstr_grow(kstr *self, int min_slen);
void kstr_reset(kstr *self);
void kstr_shrink(kstr *self, int max_size);
void kstr_assign_cstr(kstr *self, const char *assign_str);
void kstr_assign_kstr(kstr *self, kstr *assign_str);
void kstr_assign_buf(kstr *self, const void *buf, int buf_len);
void kstr_append_char(kstr *self, char c);
void kstr_append_cstr(kstr *self, const char *append_str);
void kstr_append_kstr(kstr *self, kstr *append_str);
void kstr_append_buf(kstr *self, const void *buf, int buf_len);
void kstr_append_sf(kstr *self, const char *format, ...);
void kstr_append_sfv(kstr *self, const char *format, va_list arg);
void kstr_sf(kstr *self, const char *format, ...);
void kstr_sfv(kstr *self, const char *format, va_list arg);
void kstr_tolower(kstr *str);
void kstr_mid(kstr *self, kstr *mid_str, int begin_pos, int size);
void kstr_replace(kstr *self, char *from, char *to);
int kstr_equal_cstr(kstr *first, const char *second);
int kstr_equal_kstr(kstr *first, kstr *second);
#ifdef __cplusplus
};
#endif

#endif 

