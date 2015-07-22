#include "kutil.h"
#include "kstr.h"

/* List of hashing, equality and ordering functions frequently used. */
uint32_t kutil_pointer_hc(void *k) { return (long)k * 3; }
int kutil_pointer_eq(void *k1, void *k2) { return k1 == k2; }
int kutil_pointer_cmp(void *k1, void *k2) { if (k1 < k2) return -1; if (k1 > k2) return 1; return 0; }

uint32_t kutil_int8_hc(void *k) { return (int) (*(uint8_t *) k) * 3; }
int kutil_int8_eq(void *k1, void *k2) { return *(int8_t *) k1 == *(int8_t *) k2; }
int kutil_int8_cmp(void *k1, void *k2)
{
	int8_t c1 = *(int8_t *) k1;
	int8_t c2 = *(int8_t *) k2;
	if (c1 < c2) return -1;
	if (c1 > c2) return 1;
	return 0;
}
int kutil_uint8_cmp(void *k1, void *k2)
{
	uint8_t c1 = *(uint8_t *) k1;
	uint8_t c2 = *(uint8_t *) k2;
	if (c1 < c2) return -1;
	if (c1 > c2) return 1;
	return 0;
}

uint32_t kutil_int16_hc(void *k) { return (int) (*(uint16_t *) k) * 3; }
int kutil_int16_eq(void *k1, void *k2) { return *(int16_t *) k1 == *(int16_t *) k2; }
int kutil_int16_cmp(void *k1, void *k2)
{
	int16_t c1 = *(int16_t *) k1;
	int16_t c2 = *(int16_t *) k2;
	if (c1 < c2) return -1;
	if (c1 > c2) return 1;
	return 0;
}
int kutil_uint16_cmp(void *k1, void *k2)
{
	uint16_t c1 = *(uint16_t *) k1;
	uint16_t c2 = *(uint16_t *) k2;
	if (c1 < c2) return -1;
	if (c1 > c2) return 1;
	return 0;
}

uint32_t kutil_int32_hc(void *k) { return (int) (*(int *) k) * 3; }
int kutil_int32_eq(void *k1, void *k2) { return *(int32_t *) k1 == *(int32_t *) k2; }
int kutil_int32_cmp(void *k1, void *k2)
{
	int32_t c1 = *(int32_t *) k1;
	int32_t c2 = *(int32_t *) k2;
	if (c1 < c2) return -1;
	if (c1 > c2) return 1;
	return 0;
}
int kutil_uint32_cmp(void *k1, void *k2)
{
	int c1 = *(int *) k1;
	int c2 = *(int *) k2;
	if (c1 < c2) return -1;
	if (c1 > c2) return 1;
	return 0;
}

uint32_t kutil_int64_hc(void *k) { return (int) (*(uint64_t *) k) * 3; }
int kutil_int64_eq(void *k1, void *k2) { return *(int64_t *) k1 == *(int64_t *) k2; }
int kutil_int64_cmp(void *k1, void *k2)
{
	int64_t c1 = *(int64_t *) k1;
	int64_t c2 = *(int64_t *) k2;
	if (c1 < c2) return -1;
	if (c1 > c2) return 1;
	return 0;
}
int kutil_uint64_cmp(void *k1, void *k2)
{
	uint64_t c1 = *(uint64_t *) k1;
	uint64_t c2 = *(uint64_t *) k2;
	if (c1 < c2) return -1;
	if (c1 > c2) return 1;
	return 0;
}

uint32_t kutil_cstr_hc(void *k)
{
        // Code taken from http://www.cse.yorku.ca/~oz/hash.html.
	char *c = (char *) k; 
	int h = 0;
	while (*c) h = *c++ + (h << 6) + (h << 16) - h;
	return h;
}

int kutil_cstr_eq(void *k1, void *k2)
{
	char *c1 = (char *) k1;
	char *c2 = (char *) k2;
	return (! strcmp(c1, c2));
}

int kutil_cstr_cmp(void *k1, void *k2)
{
	char *c1 = (char *) k1;
	char *c2 = (char *) k2;
	return strcmp(c1, c2);
}

uint32_t kstr_hc(void *k)
{
	kstr *c = (kstr *) k; 
	return kutil_cstr_hc(c->data);
}

int kstr_eq(void *k1, void *k2)
{
	kstr *c1 = (kstr *) k1;
	kstr *c2 = (kstr *) k2;
	return !strcmp(c1->data, c2->data);
}

int kstr_cmp(void *k1, void *k2)
{
	kstr *c1 = (kstr *) k1;
	kstr *c2 = (kstr *) k2;
	return strcmp(c1->data, c2->data);
}

/* This function returns the power of 2 just over the value specified (up to 32
 * bits).
 */
size_t kutil_next_power_of_2(size_t val)
{
    /* Fill every bit at the right of the leftmost bit set:
     * 0001 0000 0000
     * 0001 1000 0000
     * 0001 1110 0000
     * 0001 1111 1110
     * 0001 1111 1111
     *
     * Then add one:
     * 0010 0000 0000
     *
     * So only the bit to the left of the leftmost bit set is set after those 
     * operations.
     */
    val |= (val >>  1);
    val |= (val >>  2);
    val |= (val >>  4);
    val |= (val >>  8);
    val |= (val >> 16);
    val += 1;

    return val;
}

