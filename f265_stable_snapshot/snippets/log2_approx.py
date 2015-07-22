#!/usr/bin/python

# Code to generate some fast math function approximations.


#################
### UTILITIES ###
#################

import os, sys, re, string, time, random, subprocess, errno, readline, getopt, stat, tempfile, ConfigParser
from subprocess import *

# This class creates an object having the specified attributes.
# Example: person = Namespace(name="Mickey", age=18)
class Namespace(object):
    def __init__(self, **kwds): self.__dict__ = kwds

# This class creates an object in which it is possible to add fields
# dynamically. Example: store = PropStore(); store.foo = "bar"
class PropStore(object):
    
    def __setattr__(self, name, value):
        self.__dict__[name] = value
    
    def __getattr__(self, name):
        if not self.__dict__.has_key(name): raise AttributeError, name
        return self.__dict__[name]
    
    def __setitem__(self, name, value):
        self.__dict__[name] = value
    
    def __getitem__(self, name):
        if not self.__dict__.has_key(name): raise KeyError, name
        return self.__dict__[name]
    
    def __delitem__(self, name):
        del self.__dict__[name]
    
    def has_key(self, name):
        return self.__dict__.has_key(name)

    def keys(self):
        return self.__dict__.keys()

    def values(self):
        return self.__dict__.values()

# This function returns true if the string specified is alphanumeric (a-z, 0-9,
# '_'). This is a workaround for Python's unexpected implementation.
def isalpha(s):
    return s.replace('_', '').isalnum()

# This function adds spaces to the string specified until the total length of
# the string is at least 'min'.
def fill_string(s, min):
    while len(s) < min: s += ' '
    return s
    
# This function converts a string to hexadecimal. Function taken from the Python
# Cookbook.
def str_to_hex(s):
    lst = []
    for ch in s:
        hv = hex(ord(ch)).replace('0x', '')
        if len(hv) == 1:
            hv = '0'+hv
        lst.append(hv)
    
    return reduce(lambda x,y:x+y, lst)

# This function converts an hexadecimal number to a string.
def hex_to_str(s):
    return s and chr(string.atoi(s[:2], base=16)) + hex_to_str(s[2:]) or ''

# This function generates a random string, suitable for a username or password.
def gen_random(nb):
    generator = random.SystemRandom()
    s = ""
    
    for i in range(nb):
        s += generator.choice(['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 
                               'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z',
                               '0', '1', '2', '3', '4', '5', '6', '7', '8', '9'])
    return s

# This function checks if the exception specified corresponds to EINTR or
# EAGAIN. If not, the exception is raised.
def check_interrupted_ex(e):
    if e.args[0] != errno.EINTR and e.args[0] != errno.EAGAIN: raise e

# This function is a wrapper around select.select() to give it sane semantics.
def select_wrapper(rlist, wlist, xlist, timeout):
    try: return select.select(rlist, wlist, xlist, timeout)
    except select.error, e:
        check_interrupted_ex(e)
        return ([], [], [])

# Read the content of the file specified.
def read_file(path):
    f = open(path)
    data = f.read()
    f.close()
    return data

# Write the content of the file specified.
def write_file(path, data):
    f = open(path, "wb")
    f.write(data)
    f.close()

# This function writes the content of a file atomically by using a temporary
# file. The permissions of the destination file are preserved by default,
# otherwise the mode 644 is used.
def write_file_atom(path, data, preserve_flag=1):
    
    # Note that mkstemp creates the file with the mode 600.
    (unix_fileno, tmp_path) = tempfile.mkstemp(dir=os.path.dirname(path))
    tmp_file = os.fdopen(unix_fileno, "wb")
    tmp_file.write(data)
    tmp_file.close()
    
    try:
        if os.path.isfile(path):
            dest_stat = os.stat(path)
            os.chmod(tmp_path, stat.S_IMODE(dest_stat.st_mode))
            os.chown(tmp_path, dest_stat.st_uid, dest_stat.st_gid)
        else:
            os.chmod(tmp_path, 0644)
    except: pass
    
    os.rename(tmp_path, path)

# Move a file atomically and preserve the permissions and ownership information
# of the destination file, if it exists, if requested (disabled by default) and
# if possible.
def move_file(src, dest, preserve_flag=0):
    if preserve_flag and os.path.isfile(dest):
        try:
            dest_stat = os.stat(dest)
            os.chmod(src, stat.S_IMODE(dest_stat.st_mode))
            os.chown(src, dest_stat.st_uid, dest_stat.st_gid)
        except: pass
    os.rename(src, dest)

# Create the directory specified if it does not exist.
def create_dir(path):
    if not os.path.isdir(path): get_cmd_output(["mkdir", "-p", path])

# Delete the directory specified recursively.
def delete_dir(path):
    if os.path.isdir(path): get_cmd_output(["rm", "-rf", path])

# Append a slash to the path specified if the path is not "" and it does not end
# with a slash.
def append_trailing_slash(path):
    if path != "" and not path.endswith("/"): return path + "/"
    return path

# Remove any trailing slash from the path specified unless the path is '/'.
def strip_trailing_slash(path):
    if path != "/" and path.endswith("/"): return path[:-1]
    return path

# This function reads a ConfigParser object representing an INI file.
def read_ini_file(path):
    f = open(path, "rb")
    parser = ConfigParser.ConfigParser()
    parser.readfp(f)
    f.close()
    return parser

# This function writes the content of a ConfigParser object into an INI file
# atomically. The permissions of the destination file are preserved by default.
def write_ini_file(path, parser, preserve_flag=1):
    (unix_fileno, tmp_path) = tempfile.mkstemp()
    tmp_file = os.fdopen(unix_fileno, "wb")
    parser.write(tmp_file)
    tmp_file.close()
    move_file(tmp_path, path, preserve_flag)

# This function filters the specified file with the specified list of (pattern,
# replacement) pairs. If a pattern matches, the current line is replaced by the
# associated replacement string. If the replacement string is 'None', the
# current line is discarded. If no pattern matches the current line, the line is
# added back as-is.
def filter_generic_config_file(path, pair_list):
    data = ""
    f = open(path, "rb")
    
    for line in f.readlines():
        line = line.rstrip("\n")
        matched_flag = 0
        
        for pair in pair_list:
            regex = re.compile(pair[0])
            if regex.search(line):
                matched_flag = 1
                if pair[1] != None: data += regex.sub(pair[1], line) + "\n"
                break
         
        if not matched_flag:
            data += line + "\n"
    
    f.close()
    write_file_atom(path, data)

# This function returns an escaped and quoted shell argument.
def escape_shell_arg(arg):
    if type(arg) ==  None or arg == "":
        return "\'\'"

    arg = str(arg)
    arg = arg.replace("\'", "\'\"\'\"\'")

    return "\'" + arg + "\'"

# Helper function for get_cmd_output() and show_cmd_output(). Join a list of
# strings with whitespaces.
def cmd_output_join_with_whitespace(l):
    res = ""
    for e in l:
        if res != "": res += " "
        res += e
    return res

# Helper function for get_cmd_output() and show_cmd_output(). Return a list of
# strings or a single string depending on the value of 'shell_flag'.
def cmd_output_adjust_arg_list(arg_list, shell_flag):
    args = arg_list
    if shell_flag and type(arg_list) == list: args = cmd_output_join_with_whitespace(arg_list)
    elif not shell_flag and type(arg_list) == str: args = arg_list.split()
    return args

# This function executes the command specified and returns the standard output
# of the command.
# 
# By default, the function does not execute the command in the context of a
# shell. This can be overriden by setting 'shell_flag' to true.
# 
# The first argument expected by the function is a list of strings or a single
# string. If a single string is provided and the command is not executed in the
# context of a shell, the string is split at whitespaces to obtain the list of
# arguments. Conversely, if a list of strings is provided and the command is
# executed in the context of a shell, the list of strings is joined with
# whitespaces.
# 
# The 'err_behavior' value controls the behavior of the function when an error
# occurs. If 'err_behavior' is set to 'brief', an exception is thrown containing
# the error string generated by the command. If 'err_behavior' is set to 'full',
# an exception is thrown containing both the text of the command and the error
# string generated by the command. If 'err_behavior' is set to 'ignore', no
# exception is thrown. By default, 'err_behavior' is set to 'brief'.
# 
# If 'input_str' is non-null, it is written to the standard input of the
# process.
def get_cmd_output(arg_list, err_behavior="brief", shell_flag=0, input_str=None):
    args = cmd_output_adjust_arg_list(arg_list, shell_flag)
    
    try:
        stdin = None
        if input_str != None: stdin=PIPE
        proc = Popen(args=args, stdin=stdin, stdout=PIPE, stderr=PIPE, shell=shell_flag)
        (out_text, err_text) = proc.communicate(input_str)
        
        # An error occurred.
        if proc.returncode != 0 and err_behavior != "ignore":
            
            # Strip the surrounding whitespaces and the trailing '.' of both
            # streams.
            out_text = out_text.strip().rstrip('.')
            err_text = err_text.strip().rstrip('.')
            
            # If err_text is empty, use out_text if not empty, otherwise use
            # 'unknown error'.
            if len(err_text): msg = err_text
            elif len(out_text): msg = out_text
            else: msg = 'unknown error'
            
            raise Exception(msg)
        
        return out_text
        
    except Exception, e:
        
        # We're ignoring errors and it seems the command could not be executed.
        # Return an empty string.
        if err_behavior == "ignore": return ""
        
        err_msg = str(e)
        if err_behavior == "full":
            if type(args) == list: cmd_text = cmd_output_join_with_whitespace(args)
            else: cmd_text = args
            err_msg = "command '%s' failed: %s" % (cmd_text, err_msg)
        raise Exception(err_msg)

# This function is similar to get_cmd_output(), with the difference that the
# output of the command (stdout, stderr) is not redirected. The function throws
# an exception if the command fails if requested.
def show_cmd_output(arg_list, ignore_error=0, shell_flag=0, input_str=None):
    args = cmd_output_adjust_arg_list(arg_list, shell_flag)
    
    if type(arg_list) == str: cmd_name = arg_list.split()[0]
    else: cmd_name = arg_list[0]
    
    try:
        # Flush stdout and stderr since Python is buffering those streams and
        # the Popen() call is going to write to those streams directly,
        # resulting in out-of-order output when the streams are not redirected
        # to a terminal.
        sys.stdout.flush()
        sys.stderr.flush()
        
        stdin = None
        if input_str != None: stdin=PIPE
        proc = Popen(args=args, stdin=stdin, shell=shell_flag)
        proc.communicate(input_str)
        if proc.returncode != 0 and not ignore_error: raise Exception("command " + cmd_name + " failed")
        
    except Exception, e:
        if not ignore_error: raise Exception("command " + cmd_name + " failed")
    
# This class setups command completion in readline.
class readline_completer:
    def __init__(self, words):
        self.words = words
        self.prefix = None
        
    def complete(self, prefix, index):
        if prefix != self.prefix:
            self.matching_words = [ w for w in self.words if w.startswith(prefix) ]
            self.prefix = prefix
        try:
            return self.matching_words[index]
        except IndexError:
            return None

# This function prompts the user for a confirmation (y/n). It returns true if
# the confirmation was given. Note: I wrote this on a friday evening.
def get_confirm(prompt):
    try:
        while 1:
            res = raw_input(prompt + " ")
            res = string.lower(res)
            
            if (res == "yes" or res == "aye" or res == "sure" or res == "of course" or\
                res == "go ahead" or res == "why not" or res == "yeah" or res == "y"): return 1
            if (res == "no" or res == "nay" or res == "nah" or res == "never" or res == "n"): return 0
            
            print "Please answer with 'y' or 'n'.\n"
            
    except Exception:
        print ""
        raise KeyboardInterrupt

# This function prompts the user for a string. It returns the string entered,
# which can be "". The string is stripped of its surrounding whitespaces.
def prompt_string(prompt):
    try: return raw_input(prompt + " ").strip()
    except Exception:
        print ""
        raise KeyboardInterrupt

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

# Original topological sort code written by Ofer Faigon (www.bitformation.com)
# and used with permission.
def topological_sort(items, partial_order):
    """Perform topological sort.
       items is a list of items to be sorted.
       partial_order is a list of pairs. If pair (a,b) is in it, it means
       that item a should appear before item b.
       Returns a list of the items in one of the possible orders, or None
       if partial_order contains a loop.
    """

    def add_node(graph, node):
        """Add a node to the graph if not already exists."""
        if not graph.has_key(node):
            graph[node] = [0] # 0 = number of arcs coming into this node.

    def add_arc(graph, fromnode, tonode):
        """Add an arc to a graph. Can create multiple arcs.
           The end nodes must already exist."""
        graph[fromnode].append(tonode)
        # Update the count of incoming arcs in tonode.
        graph[tonode][0] = graph[tonode][0] + 1

    # step 1 - create a directed graph with an arc a->b for each input
    # pair (a,b).
    # The graph is represented by a dictionary. The dictionary contains
    # a pair item:list for each node in the graph. /item/ is the value
    # of the node. /list/'s 1st item is the count of incoming arcs, and
    # the rest are the destinations of the outgoing arcs. For example:
    #           {'a':[0,'b','c'], 'b':[1], 'c':[1]}
    # represents the graph:   c <-- a --> b
    # The graph may contain loops and multiple arcs.
    # Note that our representation does not contain reference loops to
    # cause GC problems even when the represented graph contains loops,
    # because we keep the node names rather than references to the nodes.
    graph = {}
    for v in items:
        add_node(graph, v)
    for a,b in partial_order:
        add_arc(graph, a, b)

    # Step 2 - find all roots (nodes with zero incoming arcs).
    roots = [node for (node,nodeinfo) in graph.items() if nodeinfo[0] == 0]

    # step 3 - repeatedly emit a root and remove it from the graph. Removing
    # a node may convert some of the node's direct children into roots.
    # Whenever that happens, we append the new roots to the list of
    # current roots.
    sorted = []
    while len(roots) != 0:
        # If len(roots) is always 1 when we get here, it means that
        # the input describes a complete ordering and there is only
        # one possible output.
        # When len(roots) > 1, we can choose any root to send to the
        # output; this freedom represents the multiple complete orderings
        # that satisfy the input restrictions. We arbitrarily take one of
        # the roots using pop(). Note that for the algorithm to be efficient,
        # this operation must be done in O(1) time.
        root = roots.pop()
        sorted.append(root)
        for child in graph[root][1:]:
            graph[child][0] = graph[child][0] - 1
            if graph[child][0] == 0:
                roots.append(child)
        del graph[root]
    if len(graph.items()) != 0:
        # There is a loop in the input.
        return None
    return sorted


######################
### IMPLEMENTATION ###
######################

import math, scipy, scipy.linalg

# Evaluate the approximating polynomial at a position.
def eval_poly(x, coeffs):
    sum = 0.0
    for i in range(len(coeffs)): sum += coeffs[i] * x**i;
    return sum

# Evaluate the real function at a position.
def eval_real(x):
    return 2.0**x

def eval_real_log(x):
    return math.log(x, 2)

# Evaluate the approximation error at a position.
def eval_poly_diff(x, coeffs):
    return eval_poly(x, coeffs) - eval_real(x)

# Find the extrema positions and the error.
def find_poly_diff_extrema(coeffs, start_x, end_x):
    def update_max_error(max_error, x, y):
        y_abs = abs(y)
        if y_abs > max_error[1]:
            max_error[0] = x
            max_error[1] = y_abs
    
    # Brute force evaluation steps over the interval.
    steps = 10000
    delta = (end_x - start_x) / float(steps)
    
    # Maximum error (X,Y).
    max_error = [0, 0.0];
    
    # We don't know whether the error is going up or down initially. We detect
    # this by evaluating if the difference increases between start and start +
    # delta.
    start_y = eval_poly_diff(start_x, coeffs)
    start_delta_y = eval_poly_diff(start_x + delta, coeffs)
    up_flag = start_delta_y > start_y
    
    # Process the first extrema.
    update_max_error(max_error, start_x, start_y)
    
    # Find the middle extrema.
    extrema_list = []
    last_y = start_y
    for i in range(1, steps):
        
        # Evaluate the difference at this step.
        cur_x = start_x + i*delta
        cur_y = eval_poly_diff(cur_x, coeffs)
        diff_y = cur_y - last_y
        update_max_error(max_error, cur_x, cur_y)
        
        # Detect extrema.
        pos_flag = diff_y > 0
        extrema_flag = (pos_flag != up_flag)
        if extrema_flag:
            extrema_list.append(cur_x)
            up_flag = not up_flag
        
        last_y = cur_y
    
    # Process the last extrema.
    end_y = eval_poly_diff(end_x, coeffs)
    update_max_error(max_error, end_x, end_y)
    
    return (extrema_list, max_error[1])

# Solve the polynomial equation system.
def solve_poly(pos_list):
    row_list = []
    y_list = []
    for i in range(len(pos_list)):
        x = pos_list[i]
        row = []
        
        # Polynomial coefficients.
        for j in range(len(pos_list) - 1): row.append(x**j)
        
        # Error term.
        row.append((-1)**i)
        
        row_list.append(row)
        y_list.append(eval_real(x))
    
    row_array = scipy.array(row_list)
    y_array = scipy.array(y_list)
    solution = scipy.linalg.solve(row_array, y_array)
    out = [float(x) for x in solution]
    return out

# Gnuplot command for printing the error function.
def gnuplot_string(coeffs):
    s = "plot [0:1] ("
    for i in range(len(coeffs)):
        v = float(coeffs[i])
        if i: s += " + "
        s += "%f*x**%d" % (v, i)
    s += ") - 2**x"
    return s

def gnuplot_string_log(coeffs):
    s = "plot [1:2] ("
    for i in range(len(coeffs)):
        v = float(coeffs[i])
        if i: s += " + "
        s += "%f*x**%d" % (v, i)
    s += ") - (log(x)/log(2))"
    return s

# Find the minimax solution given the provided starting positions with the Remez
# algorithm. The order of the polynomial is derived implicitly from the number
# of positions.
def find_minimax(init_pos_list):
    
    # Best observed error so far.
    best_error = 100000.0
    
    # Interval of interest.
    start_x = init_pos_list[0]
    end_x = init_pos_list[-1]
    
    # Current polynomial positions and coefficients.
    cur_pos_list = init_pos_list
    coeffs = []
    
    # Iterate until the solution stabilizes.
    for i in range(0, 10):
        
        # Get the solution. This includes the coefficients and the wanted error.
        sol = solve_poly(cur_pos_list)
        coeffs = sol[0:-1]
        wanted_error = sol[-1]
        
        # Find the extrema positions of the polynomial and the actual error.
        (extrema, actual_error) = find_poly_diff_extrema(coeffs, start_x, end_x)
        
        # Verify convergence.
        print""
        if actual_error < best_error:
            print("Converging.")
            best_error = actual_error
        elif actual_error > best_error: print("Non-convergence warning.")
        else: print("Stable.")
        print("Iter %d wanted %f actual %f best %f." % (i, wanted_error, actual_error, best_error))
        print("Coefficients: " + str(coeffs))
        print("Extrema: " + str(extrema))
            
        # Update the positions.
        cur_pos_list = [ start_x ] + extrema + [ end_x ]
        
        # Failure, we're getting more/less extrema than expected.
        if len(cur_pos_list) != len(init_pos_list):
            print("Extrema number mismatch, bailing out.")
            break
    
    print("")
    print("Final coefficients: " + str(coeffs))
    print("Absolute error: %f." % (best_error))
    print("Gnuplot: " + gnuplot_string(coeffs))

def main():
        #find_minimax([0.0, 0.5, 1.0])
        find_minimax([0.0, 0.25, 0.75, 1.0])
        #find_minimax([0.0, 0.2, 0.4, 0.8, 1.0])

def main_log():
        #find_minimax([1.0, 1.5, 2.0])
        find_minimax([1.0, 1.25, 1.75, 2.0])
        #find_minimax([1.0, 1.2, 1.4, 1.8, 2.0])
     
main()

