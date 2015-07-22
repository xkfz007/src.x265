#include <stdio.h>
#include <ctype.h>
#include "kstr.h"
#include "kutil.h"

kstr* kstr_new() {
    kstr *self = (kstr *) kmalloc(sizeof(kstr));
    kstr_init(self);
    return self;
}

kstr* kstr_new_kstr(kstr *str) {
    kstr *self = (kstr *) kmalloc(sizeof(kstr));
    kstr_init_kstr(self, str);
    return self;
}

kstr* kstr_new_cstr(char *str) {
    kstr *self = (kstr *) kmalloc(sizeof(kstr));
    kstr_init_cstr(self, str);
    return self;
}

void kstr_destroy(kstr *str) {
    kstr_clean(str);
    kfree(str);
}

void kstr_init(kstr *self) {
    self->slen = 0;
    self->mlen = 8;
    self->data = (char *) kmalloc(self->mlen);
    self->data[0] = 0;
}

void kstr_init_cstr(kstr *self, const char *init_str) {
    if (init_str == NULL) {
        init_str = "";
    }
    
    kstr_init_buf(self, init_str, strlen(init_str));
}

void kstr_init_kstr(kstr *self, kstr *init_str) {
    kstr_init_buf(self, init_str->data, init_str->slen);
}

void kstr_init_buf(kstr *self, const void *buf, int buf_len) {
    self->slen = buf_len;
    self->mlen = buf_len + 1;
    self->data = (char *) kmalloc(self->mlen);
    memcpy(self->data, buf, buf_len);
    self->data[buf_len] = 0;
}

void kstr_init_sfv(kstr *self, const char *format, va_list args) {
    kstr_init(self);
    kstr_sfv(self, format, args);
}

void kstr_init_sf(kstr *self, const char *format, ...) {
    va_list args;
    va_start(args, format);
    kstr_init_sfv(self, format, args);
    va_end(args);
}


void kstr_clean(kstr *self) {
    
    if (self == NULL)
    	return;

    kfree(self->data);
}

/* This function increases the size of the memory containing the string so that
 * it may contain at least 'min_slen' characters (not counting the terminating
 * '0').
 */
void kstr_grow(kstr *self, int min_slen) {
    assert(min_slen >= 0);
    
    if (min_slen >= self->mlen) {
    
        /* Compute the snapped size for a given requested size. By snapping to powers
         * of 2 like this, repeated reallocations are avoided.
         */
        if (min_slen < 8) {
            self->mlen = 8;
        }
        
        else {
	    self->mlen = kutil_next_power_of_2(min_slen);
        }
        
        assert(self->mlen > min_slen);    
        self->data = (char *) krealloc(self->data, self->mlen);
    }
}

void kstr_reset(kstr *self) {
    self->slen = 0;
    self->data[0] = 0;
}

/* This function ensures that the string specified does not get too large. If
 * the internal memory associated to the string is bigger than the threshold
 * specified, the memory associated to the string is released and a new, small
 * buffer is allocated for the string. In all cases, the string is cleared.
 */
void kstr_shrink(kstr *self, int max_size) {
    if (self->slen > max_size) {
    	kstr_clean(self);
	kstr_init(self);
    }
    
    kstr_reset(self);
}

void kstr_assign_cstr(kstr *self, const char *assign_str) {
    if (assign_str == NULL) {
        assign_str = "";
    }
    
    kstr_assign_buf(self, assign_str, strlen(assign_str)); 
}

void kstr_assign_kstr(kstr *self, kstr *assign_str) {
    kstr_assign_buf(self, assign_str->data, assign_str->slen);   
}

void kstr_assign_buf(kstr *self, const void *buf, int buf_len) {
    kstr_grow(self, buf_len);
    memcpy(self->data, buf, buf_len);
    self->data[buf_len] = 0;
    self->slen = buf_len;
}

void kstr_append_char(kstr *self, char c) {
    self->slen++;
    kstr_grow(self, self->slen);
    self->data[self->slen - 1] = c;
    self->data[self->slen] = 0;
}

void kstr_append_cstr(kstr *self, const char *append_str) {
    kstr_append_buf(self, append_str, strlen(append_str));
}

void kstr_append_kstr(kstr *self, kstr *append_str) {
    kstr_append_buf(self, append_str->data, append_str->slen);
}

void kstr_append_buf(kstr *self, const void *buf, int buf_len) {
    kstr_grow(self, self->slen + buf_len);
    memcpy(self->data + self->slen, buf, buf_len);
    self->slen += buf_len;
    self->data[self->slen] = 0;
}

void kstr_append_sf(kstr *self, const char *format, ...) {
    va_list arg;
    va_start(arg, format);
    kstr_append_sfv(self, format, arg);
    va_end(arg);
}

void kstr_append_sfv(kstr *self, const char *format, va_list arg) {

    /* vsnprintf() is defined in a standard. Some UNIX systems implement it
     * correctly.
     */
    #ifndef __WINDOWS__
    
    /* Determine the size of the resulting string. */
    int print_size;
    int pos;
    va_list arg2;
    
    va_copy(arg2, arg);
    print_size = vsnprintf(NULL, 0, format, arg2);
    va_end(arg2);
    
    pos = self->slen;
    self->slen += print_size;
    kstr_grow(self, self->slen);
    
    /* Do the sprintf(). */
    print_size = vsnprintf(self->data + pos, self->slen - pos + 1, format, arg);
    /* Windows doesn't support it correctly, though. */
    #else
    while (1) {
        /* Use _vsnprintf() with its ugly semantics. */
        va_list arg2;
        int print_size;
        
        va_copy(arg2, arg);
        print_size = _vsnprintf(self->data + self->slen, self->mlen - self->slen, format, arg2);
        va_end(arg2);
        
        if (print_size == -1 || print_size == self->mlen - self->slen) {
            kstr_grow(self, self->mlen * 2);
        }
        
        else {
            self->slen += print_size;
            break;
        }
    }
    #endif
}

void kstr_sf(kstr *self, const char *format, ...) {
    va_list arg;
    va_start(arg, format);
    kstr_sfv(self, format, arg);
    va_end(arg);
}

void kstr_sfv(kstr *self, const char *format, va_list arg) {
    kstr_reset(self);
    kstr_append_sfv(self, format, arg);
}

/* This function puts all the characters of a string in lowercase. */
void kstr_tolower(kstr *str) {
    int i;
    for (i = 0 ; i < str->slen ; i++) 
        str->data[i] = tolower(str->data[i]);
}

/* This function extracts a substring from the string and places it in 'mid_str'.
 * Arguments:
 * Source string.
 * String that will contain the substring.
 * Beginning of the substring in this string.
 * Size of the substring.
 */
void kstr_mid(kstr *self, kstr *mid_str, int begin_pos, int size) {
    assert(begin_pos + size <= self->slen);
    kstr_grow(mid_str, size);
    memcpy(mid_str->data, self->data + begin_pos, size);
    mid_str->data[size] = 0;
    mid_str->slen = size;
}

/* Replace all occurences of the string 'from' with the string 'to' in the
 * string specified.
 */
void kstr_replace(kstr *self, char *from, char *to) {
    kstr tmp;
    int i = 0;
    int l = strlen(from);
    
    assert(l);
    kstr_init(&tmp);
    
    while (i < self->slen) {
        if (! strncmp(self->data + i, from, l)) {
            kstr_append_cstr(&tmp, to);
            i += l;
        }
        
        else {
            kstr_append_char(&tmp, self->data[i]);
            i++;
        }
    }
    
    kstr_assign_kstr(self, &tmp);
    kstr_clean(&tmp);
}

int kstr_equal_cstr(kstr *first, const char *second) {
    return (strcmp(first->data, second) == 0);
}

int kstr_equal_kstr(kstr *first, kstr *second) {
    return (strcmp(first->data, second->data) == 0);
}

