#!/usr/bin/python

# This script auto-generates the assembly glue code. See declare_all() below.
#
# THIS SCRIPT MUST BE RUN FROM THE SNIPPETS DIRECTORY. IT WRITES FILES IN PLACE.


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

# Function declaration properties. Properties are the declare_func arguments.
class Function:
    def __init__(self):
        pass

# Dictionary containing the declarations. Keys are function names, in order.
declare_dict = odict()

# Declare the function called "name" having the return type and the arguments
# specified, optionally for both low and high bit depths, and optionally in
# array form for the indices specified. An index with value "X" has no function
# mapped in that slot. For each supported architecture, the arguments "arch",
# "arch_lbd" and "arch_hbd" specifies the slots (as applicable) that are
# currently implemented in assembly. If "arch" is specified, it shadows
# "arch_lbd" and "arch_hbd". If single_c=true, the same C function is mapped to
# every slot.
def declare_func(name, ret="void", args="", bd=0, indices=None, avx2=None, avx2_lbd=None, avx2_hbd=None, single_c=True):
    f = Function()
    f.name = name
    f.ret = ret
    f.args = args
    f.bd = bd
    f.indices = indices
    f.avx2 = avx2
    f.avx2_lbd = avx2_lbd
    f.avx2_hbd = avx2_hbd
    f.single_c = single_c
    declare_dict[name] = f

### Declare all functions. ###
def declare_all():
    df = declare_func

    amp_indices =        ["2", "4", "8", "16", "32", "64", "6", "12", "24", "48"]
    luma_amp_indices_x = ["X", "4", "8", "16", "32", "64", "X", "12", "24", "48"]
    luma_amp_indices =   ["4", "8", "16", "32", "64", "12", "24", "48"]
    luma_qpel_indices = []
    luma_qpel_indices_avx2 = []
    for index in luma_amp_indices_x:
        for frac in [ "h", "v", "d"]:
            luma_qpel_indices.append("X" if index == "X" else "%s_%s" % (index, frac))
            if index != "X" and int(index) % 8 == 0:
                luma_qpel_indices_avx2.append("%s_%s" % (index, frac))

    # Declarations go here.
    df("dct", bd=1, single_c=False,
       args="int16_t *dst, f265_pix *src, int src_stride, f265_pix *pred, int pred_stride, uint8_t *spill",
       indices=["4", "8", "16", "32", "dst"], avx2_lbd=["4", "8", "16", "32", "dst"])
    df("idct", bd=1, single_c=False,
       args="f265_pix *dst, int dst_stride, f265_pix *pred, int pred_stride, int16_t *coeffs, uint8_t *spill",
       indices=["4", "8", "16", "32", "dst"], avx2_lbd=["4", "8", "16", "32", "dst"])

    df("quant", bd=1,
       ret = "int", args="int16_t *dst, int16_t *src, int bs, int mult, int add, int shift",
       indices=["4", "8", "16", "32"], avx2_lbd=["4", "8", "16", "32"])

    df("dequant", bd=1,
       ret = "void", args="int16_t *dst, int16_t *src, int bs, int mult, int add, int shift",
       indices=["4", "8", "16", "32"], avx2_lbd=["4", "8", "16", "32"])

    df("fsad", bd=1,
       ret = "int", args="f265_pix *src, int src_stride, f265_pix *ref, int ref_stride, int packed_dims",
       indices=amp_indices,
       avx2_lbd=luma_amp_indices)

    df("sad3", bd=1,
       args="int *costs, f265_pix *src, int src_stride, f265_pix **refs, int ref_stride, int packed_dims",
       indices=luma_amp_indices_x,
       avx2_lbd=1)

    df("sad4", bd=1,
       args="int *costs, f265_pix *src, int src_stride, f265_pix **refs, int ref_stride, int packed_dims",
       indices=luma_amp_indices_x,
       avx2_lbd=1)

    df("avg_pix", bd=1,
       args="f265_pix *dst, f265_pix *src0, int src0_stride, f265_pix *src1, int src1_stride, int packed_dims",
       indices=luma_amp_indices_x,
       avx2_lbd=["4", "8", "16", "32", "64", "12", "24", "48"])

    df("interpol_luma_qpel_pix", bd=1,
       args="f265_pix *dst, int dst_stride, f265_pix *src, int src_stride, int frac, int packed_dims, uint8_t *spill",
       indices=luma_qpel_indices,
       avx2_lbd=luma_qpel_indices_avx2)


### AVX2 SAD special code. ###
def avx2_sad_special_code():

    # Handle sizes 12, 24, 48 by calling two assembly functions.
    base_code =\
        "int f265_lbd_fsad_SZ0_avx2(f265_pix *src, int src_stride, f265_pix *ref, int ref_stride, int packed_dims)\n"\
        "{\n"\
        "    return f265_lbd_fsad_SZ1_avx2(src, src_stride, ref, ref_stride, packed_dims) +\n"\
        "           f265_lbd_fsad_SZ2_avx2(src+SZ1, src_stride, ref+SZ1, ref_stride, packed_dims);\n"\
        "}\n"\
        "\n"\
        "void f265_lbd_sad3_SZ0_avx2(int *costs, f265_pix *src, int src_stride, f265_pix **refs, int ref_stride, int packed_dims)\n"\
        "{\n"\
        "    f265_pix *refs2[3] = { refs[0]+SZ1, refs[1]+SZ1, refs[2]+SZ1 };\n"\
        "    int costs2[4];\n"\
        "    f265_lbd_sad3_SZ1_avx2(costs, src, src_stride, refs, ref_stride, packed_dims);\n"\
        "    f265_lbd_sad3_SZ2_avx2(costs2, src+SZ1, src_stride, refs2, ref_stride, packed_dims);\n"\
        "    for (int i = 0; i < 3; i++) costs[i] += costs2[i];\n"\
        "}\n"\
        "\n"\
        "void f265_lbd_sad4_SZ0_avx2(int *costs, f265_pix *src, int src_stride, f265_pix **refs, int ref_stride, int packed_dims)\n"\
        "{\n"\
        "    f265_pix *refs2[4] = { refs[0]+SZ1, refs[1]+SZ1, refs[2]+SZ1, refs[3]+SZ1 };\n"\
        "    int costs2[4];\n"\
        "    f265_lbd_sad4_SZ1_avx2(costs, src, src_stride, refs, ref_stride, packed_dims);\n"\
        "    f265_lbd_sad4_SZ2_avx2(costs2, src+SZ1, src_stride, refs2, ref_stride, packed_dims);\n"\
        "    for (int i = 0; i < 4; i++) costs[i] += costs2[i];\n"\
        "}\n"\
        "\n"

    s = ""
    for sz in [12, 24, 48]:
        sz1 = (sz*2)/3
        sz2 = sz/3
        s += base_code.replace("f265_pix", "uint8_t").replace("SZ0", str(sz)).\
             replace("SZ1", str(sz1)).replace("SZ2", str(sz2))
    return s

### AVX2 luma qpel interpolation special code. ###
def avx2_interpol_luma_qpel_special_code():
    # Call the 8x8 function multiple times.
    # We only care about the height in assembly, so packed_dims doesn't need to be updated.
    base_code =\
        "void f265_lbd_interpol_luma_qpel_pix_SZ_FRAC_avx2(f265_pix *dst, int dst_stride, f265_pix *src, int src_stride, int frac,\n"\
        "                                               int packed_dims, uint8_t *spill)\n"\
        "{\n"\
        "    int nb_blocks = ((packed_dims>>8)&0xff)>>3;\n"\
        "    for (int i = 0; i < nb_blocks; i++)\n"\
        "        f265_lbd_interpol_luma_qpel_pix_8_FRAC_avx2(dst + 8*i, dst_stride, src + 8*i, src_stride,\n"\
        "                                                 frac, packed_dims, spill);\n"\
        "}\n"\
        "\n"
    s = ""
    for sz in [16, 32, 64, 24, 48]:
        for frac in [ "h", "v", "d"]:
            s += base_code.replace("f265_pix", "uint8_t").replace("SZ", str(sz)).replace("FRAC", frac)
    return s


### Special code injection. ###
def get_special_code():
    s = ""
    s += "#ifdef F265_HAVE_ASM\n"
    s += avx2_sad_special_code()
    s += avx2_interpol_luma_qpel_special_code()
    s += "#endif\n"
    return s

# Generate the output text for the C/header file.
def get_output():

    # Program prefix.
    prog = "f265"

    # List of assembly architectures supported.
    asm_arch_list = [ "avx2" ]

    # Full list of architectures.
    arch_list = ["c"] + asm_arch_list

    # Prototype text.
    proto_text = ""

    # Typedef text.
    typedef_text = ""

    # Global variable text. Group by bit depth to help caching.
    global_var_text_lbd = ""
    global_var_text_hbd = ""
    extern_var_text = ""

    # Arch-specific assignment text.
    assign_text = {}
    for arch in arch_list:
        assign_text[arch] = ""

    # Pass every function.
    for f in declare_dict.values():

        # Base function name.
        base_func_name = "%s_%s" % (prog, f.name)

        # Iterate on the bit depths, if any.
        bd_list = ["lbd", "hbd"] if f.bd else [None]
        for bd in bd_list:

            # Adjust the function name for the bit depth.
            bd_func_name = base_func_name
            if bd != None: bd_func_name = "%s_%s_%s" % (prog, bd, f.name)

            # Do the substitutions for the bit depth in the arguments and the
            # return type.
            pix_type = "uint8_t" if bd == "lbd" else "int16_t"
            func_ret_str = f.ret.replace("f265_pix", pix_type)
            func_args_str = f.args.replace("f265_pix", pix_type)

            # Declare the typedef.
            typedef_text += "typedef %s(*%s_func)(%s);\n" % (func_ret_str, bd_func_name, func_args_str)

            # Declare the global variable, along with documentation.
            var_indice_str = "[%d]" % (len(f.indices)) if f.indices else ""
            global_str = "%s_func %s%s;\n" % (bd_func_name, bd_func_name, var_indice_str)
            if bd != "hbd": global_var_text_lbd += global_str
            else: global_var_text_hbd += global_str
            if f.indices != None: extern_var_text += "// Indices: %s.\n" % (", ".join(f.indices))
            extern_var_text += "extern " + global_str + "\n";

            # Iterate on the indices, if any.
            index_list = f.indices if f.indices != None else [None]
            for index_pos in range(len(index_list)):
                index = index_list[index_pos]

                # Adjust the function name for the index.
                index_func_name = bd_func_name
                if index != None: index_func_name += "_" + index

                # Iterate on the architectures.
                for arch in arch_list:

                    # Adjust the function name for the architecture.
                    arch_func_name = index_func_name + "_" + arch
                    if f.single_c and arch == "c":
                        arch_func_name = bd_func_name + "_c"

                    # Check whether the architecture supports this function.

                    # Skipped slot.
                    if index == "X":
                        support_flag = 0

                    # C always supports the function.
                    elif arch == "c":
                        support_flag = 1

                    # Handle assembly.
                    else:
                        # Get the relevant fields.
                        bdi_field = getattr(f, arch)
                        bd_field = bdi_field if bd == None else getattr(f, "%s_%s" % (arch, bd))

                        # Do the shadowing.
                        field = bd_field if bdi_field == None else bdi_field

                        # Explicitly supported.
                        support_flag = field == 1 or type(field) is list and index in field

                    # Declare the prototype.
                    if (arch == "c" and f.single_c and index_pos == 0) or\
                       (support_flag and (arch != "c" or not f.single_c)):
                        # Kludge for the interpolation functions.
                        if arch_func_name.find("interpol") != -1 and arch == "c":
                            for frac in [ "h", "v", "d"]:
                                proto_text += "%s %s_%s_c(%s);\n" % (func_ret_str, bd_func_name, frac, func_args_str);
                        # Normal declaration.
                        else:
                            proto_text += "%s %s(%s);\n" % (func_ret_str, arch_func_name, func_args_str);

                    # Not supported, skip.
                    if not support_flag: continue

                    # Do the assignments.
                    assign_tabs = "    "
                    if arch != "c": assign_tabs += "    "
                    assign_index_str = "[%d]" % (index_pos) if f.indices else ""
                    assign_val = arch_func_name
                    # Kludge for the interpolation functions.
                    if arch_func_name.find("interpol") != -1 and arch == "c":
                        assign_val = "%s_%s_c" % (arch_func_name[:-2], index[-1])
                    assign_text[arch] += "%s%s%s = %s;\n" % (assign_tabs, bd_func_name, assign_index_str, assign_val)

        proto_text += "\n"


    # Perform the final assemblage.
    s = ""
    s += "// This file was auto-generated by snippets/asm.py.\n"
    s += "// It handles the linkage of the assembly functions.\n"
    top_msg = s

    s = top_msg
    s += "\n"
    s += "// Prototypes.\n"
    s += proto_text
    s += "// Special code.\n"
    s += get_special_code() + "\n"
    s += "// Globals.\n"
    s += global_var_text_lbd + "\n"
    s += global_var_text_hbd + "\n"
    s += "// Linkage at runtime.\n"
    s += "static void f265_link_asm(int avx2_flag)\n"
    s += "{\n"
    s += assign_text["c"] + "\n"
    s += "    #ifdef F265_HAVE_ASM\n"
    s += "    if (avx2_flag)\n"
    s += "    {\n"
    s += assign_text["avx2"]
    s += "    }\n"
    s += "    #endif\n"
    s += "}\n\n"
    c_content = s

    s = top_msg
    s += "\n"
    s += "// Typedefs.\n"
    s += typedef_text + "\n"
    s += "// Globals.\n\n"
    s += extern_var_text + "\n"
    h_content = s

    return (c_content, h_content)

def main():
    declare_all()
    (c_content, h_content) = get_output()

    # Sanity check.
    if not os.path.exists("../f265/asm.c"):
        print "**asm.c not found, not touching filesystem**\n"
        print c_content
        print h_content
        return

    write_file("../f265/asm.c", c_content)
    write_file("../f265/asm.h", h_content)

main()

