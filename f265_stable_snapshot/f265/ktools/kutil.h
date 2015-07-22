#ifndef __KUTIL_H__
#define __KUTIL_H__

#include <string.h>
#include "ktools.h"

#ifdef __cplusplus
extern "C" {
#endif
uint32_t kutil_pointer_hc(void *k);
int kutil_pointer_eq(void *k1, void *k2);
int kutil_pointer_cmp(void *k1, void *k2);
uint32_t kutil_int8_hc(void *k);
int kutil_int8_eq(void *k1, void *k2);
int kutil_int8_cmp(void *k1, void *k2);
int kutil_uint8_cmp(void *k1, void *k2);
uint32_t kutil_int16_hc(void *k);
int kutil_int16_eq(void *k1, void *k2);
int kutil_int16_cmp(void *k1, void *k2);
int kutil_uint16_cmp(void *k1, void *k2);
uint32_t kutil_int32_hc(void *k);
int kutil_int32_eq(void *k1, void *k2);
int kutil_int32_cmp(void *k1, void *k2);
int kutil_uint32_cmp(void *k1, void *k2);
uint32_t kutil_int64_hc(void *k);
int kutil_int64_eq(void *k1, void *k2);
int kutil_int64_cmp(void *k1, void *k2);
int kutil_uint64_cmp(void *k1, void *k2);
uint32_t kutil_cstr_hc(void *k);
int kutil_cstr_eq(void *k1, void *k2);
int kutil_cstr_cmp(void *k1, void *k2);
uint32_t kutil_string_hc(void *k);
int kutil_string_eq(void *k1, void *k2);
int kutil_string_cmp(void *k1, void *k2);
size_t kutil_next_power_of_2(size_t val);
#ifdef __cplusplus
};
#endif

#endif

