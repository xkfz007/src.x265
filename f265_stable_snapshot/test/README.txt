This document describes the f265test and f265bench utilities.


* f265test.

f265test is a script that performs benchmarks and regression tests. It is
heavily tailored for UNIX-like systems. It is possible to use it on Windows in
a MinGW/MSYS environment. However, some functionalities may not work. We did
manage to use f265test on Windows 8.1 with raw YUV input sequences. See the
notes at the end of this document for the required steps.

Example:

$ ./f265test tests.ini
Processing tests.ini ...
  Run test foo, video crash 640x480 @ QP 30 (881 kb/s):
    f265     secs 0.76  bytes 24500  psnr 50.72 56.81 56.21  ssim [...]
    x264     secs 0.36  bytes 11451  psnr 54.42 59.06 58.61  ssim [...]
    HM       secs 0.46  bytes 24500  psnr 50.72 56.81 56.21  ssim [...]

In a nutshell, f265test parses every file specified on the command line and
executes the tests contained within in declaration order. Each file may define
individual tests, test batches, test results and video definitions.

An individual test executes one or several encoders with a specific video,
bitrate, encoder parameters, execution environment and compilation settings.
Several test modes are supported:
- Execute the encoders and report the run time, output size and quality.
- Compare the bitstream of two encoders for bit-exactness.
- Compare the bitstream with previously recorded results for bit-exactness.
- Run a unit test in an encoder.

A test batch executes individual tests for several videos and bitrates. The
batch modes are as follow:
- Run the tests with different QPs.
- Run the tests with different bitrates (can be scaled for each video).

For the f265 encoder, f265test can verify if the assembly code produces the
same output as its C equivalent, verify if the YUV reconstruction matches the
decoded bitstream, and verify if Valgrind finds errors during the execution.

F265test supports the following use cases:
- Fast regression tests before commit.
- Long regression tests during the night build.
- Video/bitrate coverage tests during development.


** Processing.

The processing is purely data-driven. The configuration and the test data are
stored in INI files (a bunch of sections containing key-value pairs).

The INI files are processed in the order specified on the command line, after
the processing of the configuration files. Whatever data is defined during the
processing of the current file is kept around for use by the next files. New
declarations clobber old declarations. Hence you can mix-and-match the standard
test files with your own custom test files.

The name of an INI file section determines how the section is interpreted:
- config:   local script configuration.
- aliases:  set up file aliases, e.g. full => f265/regression/full.ini.
- video_X:  describe video "X".
- test_X:   define individual test "X".
- result_X: results of test "X".
- batch:    define a test batch.

There can be at most one batch section per INI file. If a batch section exists,
all the individual tests in the file become part of the batch, overriding their
mode, video, QP, etc. with the batch settings.


** First-time setup.

From this directory (i.e. f265/test):
$ cp data/f265testrc ~/.f265testrc
$ vim ~/.f265testrc                             # Edit and set up as needed.
$ cp data/f265_enc.sh data/hm_enc.sh $HOME      # Or make your own scripts.
$ ./f265test -v data/test_syntax.ini            # Run a quick test. This
                                                # clobbers the van_cfg.h files.


** Usage notes.

f265test can be told to ignore database results and print the actual results in
INI format after a file has been executed. Thus, updating a file can be done by
deleting the results from the file (located at the end) and running the script:
$ f265test -iq file1 >> file1

f265test feeds input files directly to the encoder if it accepts them,
otherwise GStreamer is used to write data in a named pipe to fake a raw YUV
input file. At this time only H.264/MP4 can be faked in this way.


** f265test on Windows.

To get a minimum amount of functionality on your Windows plateform, here are
four tips to help you on your way. We did manage to get f265test to work on
Windows 8.1. However, YUV inputs are required when using the HM reference
software.
1. Favor Windows-style paths over Unix-style paths.
2. Converting shell scripts to batch files.
3. Installing pyreadline.
4. Installing libdl
5. Calling f265test.

1. Favor Windows-style paths over Unix-style paths.

f265test is written in Python. To minimize problems, we strongly encourage you
to use Windows-style paths (e.g. C:\path\to\resource) instead of Unix-style
paths (e.g. /path/to/resource). While the Unix-style paths should work in the
MinGW/MSYS environment, we have encountered situations where resources could
not be found, but where correctly located using the Windows-style paths. This
is especially true when using absolute paths.

2. Converting shell scripts to batch files.

Shell scripts will work in a MinGW/MSYS environment. However, when called from
Python, shell scripts are not recognized as executables on a Windows plateform.
The simplest workaround is to convert the 'hm_enc.sh' and 'f265_enc.sh' shell
scripts to batch files. Here are two short working examples. Simply adapt the
paths for your system.

hm_enc.bat

set curr_dir=%cd%
chdir /D C:\path\to\hm\folder
make
chdir /D %curr_dir%

f265_enc.bat

set curr_dir=%cd%
chdir /D C:\path\to\f265\folder
if "%1" == "-c" (
    scons asm=0
) else (
    scons asm=1
)
chdir /D %curr_dir%

In the f265testrc file, you will need to replace the ".sh" instances by ".bat".

3. Installing pyreadline

f265test imports readline, a useful GNU package not found on Windows.
Fortunately, a Python implementation of readline called pyreadline can be
installed to replicate readline's functionalities. The package can be found at
https://pypi.python.org/pypi/pyreadline.

We installed pyreadline 2.0. We are unaware if the older 1.7.1 version works.
Given that we are running Python 2.7.6 32-bit, we installed the win32 version
using the MS Windows installer.

4. Installing libdl

MinGW/MSYS does not come with libdl, which is required to compile the HM using
the makefile found in the build/linux folder. You will need to build and
install the dlfcn-win32 library found at https://code.google.com/p/dlfcn-win32.

From the MSYS shell:
$ tar -xjf dlfcn-win32-r19.tar.bz2
$ cd dlfcn-win32-r19
$ ./configure
$ make
$ make install

5. Calling f265test

Additionally, given that the script was written on a Linux plateform, you may
encounter newline character problems (a.k.a. bad interpreter). A safer bet is
to use the following syntax (instead of the one used in the example above):

$ python f265test -v data/test_syntax.ini


* f265bench.

f265bench is a layer on top of f265test to execute tests and display the
results in a web page. This is a work in progress. The documentation is incomplete.


** Quick setup.

Create a benchmark directory, e.g. "bench".
Copy test/data/bench.ini in the benchmark directory.
Edit bench.ini and add your tests.
Make sure the path to f265test is set correctly in bench.ini.
Execute f265bench on the benchmark directory
$ f265bench -v bench


** Design document (to integrate properly).

Goals.

1) Fit right into the code-compile-test routine.
2) Allow us to share benchmark reports easily.

We are busy people; every second spent dealing with a clunky tool is wasted.

Use cases.

Compile the codec and launch a quick test to determine the impact on quality,
check whether the result is bit-exact with the prior run, and check if there
are reconstruction or Valgrind errors.

Compare the current results with those of a prior run, for quality and
performance. Scrap the results if they aren't good or keep them as reference.

Manage the result set by copying or deleting files. There is one file per
result. The files are kept in filesystem directories, in a shallow hierarchy.
You want some results gone, you delete their files or directories. There is no
database, no metadata, and no housekeeping.

Manage the test setup by editing a single INI file. The INI file declares the
videos, the number of frames and the QP range. There is one test per INI file
section. The section name identifies the test name and the encoder used. A
section field sets the encoder parameters. The user modifies the test setup by
copying text around and commenting out stuff.

f265bench executes the tests in order. We use the existing f265test script for
the execution. For each video, f265bench writes a "batch" definition (iteration
over the QPs) for consumption by f265test. The script output is written in a
text file. That file is the test result for the video. It contains, for each
QP, the md5sum and the size of the output file, the average PSNR/SSIM
statistics and the execution time.

It's possible to use all the cores of the machine to execute the tests, if
there are enough test videos. If requested, f265bench will run several
instances of f265test in parallel, up to one instance per test video. In "fair"
mode, the parallelism is limited to the videos from the same test. In "unfair"
mode, there is no such limitation. There is no support to parallelize f265test
itself, i.e. the QP range for a video is processed serially.

The number of frames can be adjusted per video. This can be exploited to keep
all cores busy by setting fewer frames for larger videos. This is useful to
obtain as much data as possible in a short time frame, during development.

The test results are cached by default. If a result file exists, the test is
not executed again. This behavior can be overridden by an option in the test
section. In a typical development scenario, we want to establish a large test
bed with different setups and encoders, then focus on optimizing one test in
particular. In that case, we can configure that particular test to auto-clobber
and keep on testing until the code starts working. To preserve results, either
the test name can be changed or the directory containing the test results can
be renamed.

The encoder binaries known to f265test are used by default. It's also possible
to specify another path to use a precompiled binary in a test.

f265bench collects MD5 sums at every step to help the user detect if the output
of the encoder has changed. For each test, and for each video, f265bench
outputs the MD5 sum of the MD5 sums appearing in a result file (there is one
MD5 sum per QP in a result file). For each test, f265bench outputs the MD5 sum
of the MD5 sum of the test videos. Finally, f265bench outputs the master MD5
sum, which is the MD5 sum of every test MD5 sum. The user notes the master MD5
sum in a file, and checks that it still matches while doing changes that
shouldn't impact quality such as assembly programming. Typically it's
sufficient to remember the first few characters, so it's practical solution.

All binaries are copied in a work directory before executing any of the tests.
Hence the user can launch the tests and still be able to compile and run the
code. The work directory used by f265bench can be named with the process ID
(PID) to make it unique. If so desired, the user can launch a test batch,
change the code, compile, change the test name in the config file and run
another batch of tests without waiting for the first batch to complete. This is
living dangerously, but it works if you're careful.

The main output of f265bench is a self-contained web page that summarizes the
test result. By default all test results are compared together, for each video.
f265bench doesn't care where the results come from. It merely scans the dataset
directory and picks up result files from there. Hence if the user deletes all
tests from the config file, f265bench will compare whatever results are cached
in the dataset directory. The layout of the dataset directory is
dataset/test_name/video_name.result.

It's the user's responsibility to ensure that the results are comparable. The
main thing to avoid is changing the number of frames encoded for a video.
Normally, for the same video, all results are comparable. f265bench will record
and check the number of frames in the result file for mismatches.

The reports outputted by f265bench can be tweaked by adding report sections in
the config file. Each report section is essentially a bunch of regular
expressions to match result files.

f265bench can be instructed to open a new window in Firefox when the report is
ready, so you can launch a test and go back to work until it pops up. The web
page contains HTML anchors so that the user can view exactly the desired
result.

Sharing results is important. It must be easy to send/view the results for both
the sender and the recipients, otherwise neither will bother. Zipping and
unzipping an archive takes too much effort. The easiest way I can imagine is to
copy the report directory in the local web server directory and send out the
resulting URL to the recipients, then leave it there gathering dust.

The web page will be organized by report sections. Each section has three
graphs per video. The main graph contains the usual RD curves. It plots the
file size (not the bitrate, not the QP, just the file size) against the PSNR or
SSIM metric. That graph becomes an unreadable mess when it contains more than 4
lines, so the two other graphs are histograms that scale for more results.

The histograms are obtained by selecting a region of interest (ROI) on the
curves and splitting it in intervals. The ROI is the part of the graph where
there is data for each line. For example,
 -------------  
----------
  ---------
  ^      ^
  | ROI  |


The histograms display the SSIM and encoding time at those intervals by
interpolating the stats on either side of these points. The histograms are
empty if the ROI is empty because the RD curves do not overlap for any common
point. There is one important degenerate case that we handle. If at most one RD
curve contains a single point (e.g. for a test at QP 30), then the histograms
contain a single line for that point. That is useful to do a quicky-and-dirty
test at one QP after having plotted the rest of the curve once. The histograms
contain numeric values so there is no guesswork.

