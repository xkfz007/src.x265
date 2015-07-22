#!/usr/bin/python

# This script auto-generates the quarterpel interpolation code.


#################
### UTILITIES ###
#################

import os, sys, re, string, time, random, subprocess, errno, readline, getopt, stat, tempfile, ConfigParser
from subprocess import *

# Write the content of the file specified.
def write_file(path, data):
    f = open(path, "wb")
    f.write(data)
    f.close()

# Ordered dictionary.
class odict(dict):
    def __init__(self, data=None):
        self._keys = []
        dict.__init__(self)

        if data:
            # we were provided a regular dict
            if isinstance(data, dict):
                self.append_from_dict(data)

            # we were provided a tuple list
            elif type(data) == list:
                self.append_from_plist(data)

            # we were provided invalid input
            else:
                raise Exception("expected a dict or a tuple list")

    def append_from_dict(self, dict):
        map(self.__setitem__, dict.keys(), dict.values())

    def append_from_plist(self, plist):
        for pair in plist:
            if len(pair) != 2:
                raise Exception("invalid pairs list")
        for (k, v) in plist:
            self.__setitem__(k, v)

    def __delitem__(self, key):
        if not key in self._keys:
            raise KeyError, key
        dict.__delitem__(self, key)
        self._keys.remove(key)

    def __setitem__(self, key, item):
        dict.__setitem__(self, key, item)
        if key not in self._keys:
            self._keys.append(key)

    def clear(self):
        dict.clear(self)
        self._keys = []

    def copy(self):
        return odict(self.plist())

    def items(self):
        return zip(self._keys, self.values())

    def keys(self):
        return list(self._keys) # return a copy of the list

    def values(self):
        return map(self.get, self._keys)

    def plist(self):
        p = []
        for k, v in self.items():
            p.append( (k, v) )
        return p

    def __str__(self):
        s = "{"
        l = len(self._keys)
        for k, v in self.items():
            l -= 1
            strkey = str(k)
            if isinstance(k, basestring): strkey = "'"+strkey+"'"
            strval = str(v)
            if isinstance(v, basestring): strval = "'"+strval+"'"
            s += strkey + ":" + strval
            if l > 0: s += ", "
        s += "}"
        return s


######################
### IMPLEMENTATION ###
######################

def get_one_interpol_c(comp, dst, frac):
    func_name = "venc_interpol_%s_qpel_%s_%s_c" % (comp, dst, frac)
    dst_type = "f265_pix" if dst == "pix" else "int16_t"
    frac_and = "3" if comp == "luma" else "7"
    frac_shift = "2" if comp == "luma" else "3"
    frac_str = "frac&%s" % (frac_and) if frac == "h" else "frac>>%s" % (frac_shift)
    factors = "f265_if_%s[%s]" % (comp, frac_str)
    if dst == "pix":
        shift = "20 - bd" if frac == "d" else "6"
    else:
        shift = "6" if frac == "d" else "bd - 8"

    s = ""
    s += "void %s(%s *dst, int dst_stride, f265_pix *src, int src_stride, int frac,\n" % (func_name, dst_type)
    for i in range(6+len(func_name)): s += " "
    s += "int packed_dims, uint8_t *spill)"
    func_decl = s
    proto = "%s;\n" % (func_decl)

    s = func_decl + "\n"
    s += "{\n"
    s += "    #ifdef F265_LBD\n"
    s += "    int bd = 8, width = (packed_dims>>8)&0xff, height = packed_dims&0xff;\n"
    s += "    #else\n"
    s += "    int bd = (packed_dims>>16)&0xff, width = (packed_dims>>8)&0xff, height = packed_dims&0xff;\n"
    s += "    #endif\n"
    s += "    int shift = %s;\n" % (shift)
    if dst == "pix": s += "    int bias = 1<<(shift-1), max_val = (1<<bd)-1;\n";
    s += "    const int8_t *factors = %s;\n" % (factors)
    s += "\n"

    if frac == "d":
        if comp == "luma":
            row_above = 3
            row_below = 4
        else:
            row_above = 1
            row_below = 2
        call_name = "venc_interpol_%s_qpel_s16_h_c" % (comp)
        s += "    // Interpolate %d rows above the block and %d rows below.\n" % (row_above, row_below)
        s += "    int16_t *h_buf = (int16_t*)spill;\n"
        s += "    %s(h_buf, width, src - %d*src_stride, src_stride, frac,\n" % (call_name, row_above)
        for i in range(5+len(call_name)): s += " "
        s += "((bd-8)<<16)|(width<<8)|(height+%d), NULL);\n" % (row_above+row_below)
        s += "    h_buf += %d*width;\n" % (row_above)
        s += "\n"

    space_8 = "        "
    space_12 = "            "
    if frac == "d": s += "    for (int y = 0; y < height; y++, dst += dst_stride, h_buf += width)\n"
    else: s += "    for (int y = 0; y < height; y++, dst += dst_stride, src += src_stride)\n"
    s += space_8 + "for (int x = 0; x < width; x++)\n"

    filter_src = "h_buf" if frac == "d" else "src"
    if frac == "h": filter_stride = "1"
    elif frac == "v": filter_stride = "src_stride"
    else: filter_stride = "width"
    taps = 8 if comp == "luma" else 4
    filter_str = "FILTER_%dTAPS(%s+x, %s, factors)" % (taps, filter_src, filter_stride)
    if dst == "pix":
        s += space_8 + "{\n"
        s += space_12 + "int pix = (%s + bias)>>shift;\n" % (filter_str)
        s += space_12 + "dst[x] = F265_CLAMP(pix, 0, max_val);\n"
        s += space_8 + "}\n"
    else:
        s += space_12 + "dst[x] = %s>>shift;\n" % (filter_str)

    s += "}\n\n"
    define = s

    return (proto, define)

def get_interpol_c():
    proto = ""
    define = ""
    for comp in ["luma", "chroma"]:
        for frac in ["h", "v", "d"]:
            for dst in ["pix", "s16"]:
                ret = get_one_interpol_c(comp, dst, frac)
                proto += ret[0]
                define += ret[1]

    return "%s\n%s" % (proto, define)

def main():
    print(get_interpol_c())

main()

