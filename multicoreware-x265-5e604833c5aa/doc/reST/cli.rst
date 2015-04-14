*********************
Command Line Options
*********************

.. _string-options-ref:

Note that unless an option is listed as **CLI ONLY** the option is also
supported by x265_param_parse(). The CLI uses getopt to parse the
command line options so the short or long versions may be used and the
long options may be truncated to the shortest unambiguous abbreviation.
Users of the API must pass x265_param_parse() the full option name.

Preset and tune have special implications. The API user must call
x265_param_default_preset() with the preset and tune parameters they
wish to use, prior to calling x265_param_parse() to set any additional
fields. The CLI does this for the user implicitly, so all CLI options
are applied after the user's preset and tune choices, regardless of the
order of the arguments on the command line.

If there is an extra command line argument (not an option or an option
value) the CLI will treat it as the input filename.  This effectively
makes the :option:`--input` specifier optional for the input file. If
there are two extra arguments, the second is treated as the output
bitstream filename, making :option:`--output` also optional if the input
filename was implied. This makes :command:`x265 in.y4m out.hevc` a valid
command line. If there are more than two extra arguments, the CLI will
consider this an error and abort.

Generally, when an option expects a string value from a list of strings
the user may specify the integer ordinal of the value they desire. ie:
:option:`--log-level` 3 is equivalent to :option:`--log-level` debug.

Standalone Executable Options
=============================

.. option:: --help, -h

	Display help text

	**CLI ONLY**

.. option:: --version, -V

	Display version details

	**CLI ONLY**

.. option:: --asm <integer:false:string>, --no-asm

	x265 will use all detected CPU SIMD architectures by default. You can
	disable all assembly by using :option:`--no-asm` or you can specify
	a comma separated list of SIMD architectures to use, matching these
	strings: MMX2, SSE, SSE2, SSE3, SSSE3, SSE4, SSE4.1, SSE4.2, AVX, XOP, FMA4, AVX2, FMA3

	Some higher architectures imply lower ones being present, this is
	handled implicitly.

	One may also directly supply the CPU capability bitmap as an integer.

.. option:: --threads <integer>

	Number of threads to allocate for the worker thread pool  This pool
	is used for WPP and for distributed analysis and motion search:
	:option:`--wpp` :option:`--pmode` and :option:`--pme` respectively.

	If :option:`--threads`=1 is specified, then no thread pool is
	created. When no thread pool is created, all the thread pool
	features are implicitly disabled. If all the pool features are
	disabled by the user, then the pool is implicitly disabled.

	Default 0, one thread is allocated per detected hardware thread
	(logical CPU cores)

.. option:: --pmode, --no-pmode

	Parallel mode decision, or distributed mode analysis. When enabled
	the encoder will distribute the analysis work of each CU (merge,
	inter, intra) across multiple worker threads. Only recommended if
	x265 is not already saturating the CPU cores. In RD levels 3 and 4
	it will be most effective if --rect was enabled. At RD levels 5 and
	6 there is generally always enough work to distribute to warrant the
	overhead, assuming your CPUs are not already saturated.
	
	--pmode will increase utilization without reducing compression
	efficiency. In fact, since the modes are all measured in parallel it
	makes certain early-outs impractical and thus you usually get
	slightly better compression when it is enabled (at the expense of
	not skipping improbable modes).

	This feature is implicitly disabled when no thread pool is present.

	Default disabled

.. option:: --pme, --no-pme

	Parallel motion estimation. When enabled the encoder will distribute
	motion estimation across multiple worker threads when more than two
	references require motion searches for a given CU. Only recommended
	if x265 is not already saturating CPU cores. :option:`--pmode` is
	much more effective than this option, since the amount of work it
	distributes is substantially higher. With --pme it is not unusual
	for the overhead of distributing the work to outweigh the
	parallelism benefits.
	
	This feature is implicitly disabled when no thread pool is present.

	--pme will increase utilization on many core systems with no effect
	on the output bitstream.
	
	Default disabled

.. option:: --preset, -p <integer|string>

	Sets parameters to preselected values, trading off compression efficiency against 
	encoding speed. These parameters are applied before all other input parameters are 
	applied, and so you can override any parameters that these values control.

	0. ultrafast
	1. superfast
	2. veryfast
	3. faster
	4. fast
	5. medium **(default)**
	6. slow
	7. slower
	8. veryslow
	9. placebo

.. option:: --tune, -t <string>

	Tune the settings for a particular type of source or situation. The changes will
	be applied after :option:`--preset` but before all other parameters. Default none

	**Values:** psnr, ssim, zero-latency, fast-decode.

.. option:: --frame-threads, -F <integer>

	Number of concurrently encoded frames. Using a single frame thread
	gives a slight improvement in compression, since the entire reference
	frames are always available for motion compensation, but it has
	severe performance implications. Default is an autodetected count
	based on the number of CPU cores and whether WPP is enabled or not.

	Over-allocation of frame threads will not improve performance, it
	will generally just increase memory use.

.. option:: --log-level <integer|string>

	Logging level. Debug level enables per-frame QP, metric, and bitrate
	logging. If a CSV file is being generated, debug level makes the log
	be per-frame rather than per-encode. Full level enables hash and
	weight logging. -1 disables all logging, except certain fatal
	errors, and can be specified by the string "none".

	0. error
	1. warning
	2. info **(default)**
	3. debug
	4. full

.. option:: --csv <filename>

	Writes encoding results to a comma separated value log file. Creates
	the file if it doesnt already exist, else adds one line per run.  if
	:option:`--log-level` is debug or above, it writes one line per
	frame. Default none

.. option:: --cu-stats, --no-cu-stats

	Records statistics on how each CU was coded (split depths and other
	mode decisions) and reports those statistics at the end of the
	encode. Default disabled

.. option:: --output, -o <filename>

	Bitstream output file name. If there are two extra CLI options, the
	first is implicitly the input filename and the second is the output
	filename, making the :option:`--output` option optional.

	The output file will always contain a raw HEVC bitstream, the CLI
	does not support any container file formats.

	**CLI ONLY**

.. option:: --no-progress

	Disable CLI periodic progress reports

	**CLI ONLY**

Quality reporting metrics
=========================

.. option:: --ssim, --no-ssim

	Calculate and report Structural Similarity values. It is
	recommended to use :option:`--tune` ssim if you are measuring ssim,
	else the results should not be used for comparison purposes.
	Default disabled

.. option:: --psnr, --no-psnr

	Calculate and report Peak Signal to Noise Ratio.  It is recommended
	to use :option:`--tune` psnr if you are measuring PSNR, else the
	results should not be used for comparison purposes.  Default
	disabled

Input Options
=============

.. option:: --input <filename>

	Input filename, only raw YUV or Y4M supported. Use single dash for
	stdin. This option name will be implied for the first "extra"
	command line argument.

	**CLI ONLY**

.. option:: --y4m

	Parse input stream as YUV4MPEG2 regardless of file extension,
	primarily intended for use with stdin (ie: :option:`--input` -
	:option:`--y4m`).  This option is implied if the input filename has
	a ".y4m" extension

	**CLI ONLY**

.. option:: --input-depth <integer>

	YUV only: Bit-depth of input file or stream

	**Values:** any value between 8 and 16. Default is internal depth.

	**CLI ONLY**

.. option:: --dither

	Enable high quality downscaling. Dithering is based on the diffusion
	of errors from one row of pixels to the next row of pixels in a
	picture. Only applicable when the input bit depth is larger than
	8bits and internal bit depth is 8bits. Default disabled

	**CLI ONLY**

.. option:: --nr <integer>

	Noise reduction - an adaptive deadzone applied after DCT
	(subtracting from DCT coefficients), before quantization, on inter
	blocks. It does no pixel-level filtering, doesn't cross DCT block
	boundaries, has no overlap, doesn't affect intra blocks. The higher
	the strength value parameter, the more aggressively it will reduce
	noise.

	Enabling noise reduction will make outputs diverge between different
	numbers of frame threads. Outputs will be deterministic but the
	outputs of -F2 will no longer match the outputs of -F3, etc.

	**Values:** any value in range of 100 to 1000. Default disabled.

.. option:: --input-res <wxh>

	YUV only: Source picture size [w x h]

	**CLI ONLY**

.. option:: --input-csp <integer|string>

	YUV only: Source color space. Only i420, i422, and i444 are
	supported at this time. The internal color space is always the
	same as the source color space (libx265 does not support any color
	space conversions).

	0. i400
	1. i420 **(default)**
	2. i422
	3. i444
	4. nv12
	5. nv16

.. option:: --fps <integer|float|numerator/denominator>

	YUV only: Source frame rate

	**Range of values:** positive int or float, or num/denom

.. option:: --interlaceMode <false|tff|bff>, --no-interlaceMode

	**EXPERIMENTAL** Specify interlace type of source pictures. 
	
	0. progressive pictures **(default)**
	1. top field first 
	2. bottom field first

	HEVC encodes interlaced content as fields. Fields must be provided to
	the encoder in the correct temporal order. The source dimensions
	must be field dimensions and the FPS must be in units of fields per
	second. The decoder must re-combine the fields in their correct
	orientation for display.

.. option:: --seek <integer>

	Number of frames to skip at start of input file. Default 0

	**CLI ONLY**

.. option:: --frames, -f <integer>

	Number of frames to be encoded. Default 0 (all)

	**CLI ONLY**

.. option:: --qpfile <filename>

	Specify a text file which contains frametypes and QPs for some or
	all frames. The format of each line is:

	framenumber frametype QP

	Frametype can be one of [I,i,P,B,b]. **B** is a referenced B frame,
	**b** is an unreferenced B frame.  **I** is a keyframe (random
	access point) while **i** is a I frame that is not a keyframe
	(references are not broken).

	Specifying QP (integer) is optional, and if specified they are
	clamped within the encoder to qpmin/qpmax.

.. option:: --scaling-list <filename>

	Quantization scaling lists. HEVC supports 6 quantization scaling
	lists to be defined; one each for Y, Cb, Cr for intra prediction and
	one each for inter prediction.

	x265 does not use scaling lists by default, but this can also be
	made explicit by :option:`--scaling-list` *off*.

	HEVC specifies a default set of scaling lists which may be enabled
	without requiring them to be signaled in the SPS. Those scaling
	lists can be enabled via :option:`--scaling-list` *default*.
    
	All other strings indicate a filename containing custom scaling
	lists in the HM format. The encode will abort if the file is not
	parsed correctly. Custom lists must be signaled in the SPS

.. option:: --lambda-file <filename>

	Specify a text file containing values for x265_lambda_tab and
	x265_lambda2_tab. Each table requires MAX_MAX_QP+1 (70) float
	values.
	
	The text file syntax is simple. Comma is considered to be
	white-space. All white-space is ignored. Lines must be less than 2k
	bytes in length. Content following hash (#) characters are ignored.
	The values read from the file are logged at :option:`--log-level`
	debug.

	Note that the lambda tables are process-global and so the new values
	affect all encoders running in the same process. 
	
	Lambda values affect encoder mode decisions, the lower the lambda
	the more bits it will try to spend on signaling information (motion
	vectors and splits) and less on residual. This feature is intended
	for experimentation.

Profile, Level, Tier
====================

.. option:: --profile <string>

	Enforce the requirements of the specified profile, ensuring the
	output stream will be decodable by a decoder which supports that
	profile.  May abort the encode if the specified profile is
	impossible to be supported by the compile options chosen for the
	encoder (a high bit depth encoder will be unable to output
	bitstreams compliant with Main or Mainstillpicture).

	API users must use x265_param_apply_profile() after configuring
	their param structure. Any changes made to the param structure after
	this call might make the encode non-compliant.

	**Values:** main, main10, mainstillpicture, main422-8, main422-10, main444-8, main444-10

	**CLI ONLY**

.. option:: --level-idc <integer|float>

	Minimum decoder requirement level. Defaults to 0, which implies
	auto-detection by the encoder. If specified, the encoder will
	attempt to bring the encode specifications within that specified
	level. If the encoder is unable to reach the level it issues a
	warning and aborts the encode. If the requested requirement level is
	higher than the actual level, the actual requirement level is
	signaled.

	Beware, specifying a decoder level will force the encoder to enable
	VBV for constant rate factor encodes, which may introduce
	non-determinism.

	The value is specified as a float or as an integer with the level
	times 10, for example level **5.1** is specified as "5.1" or "51",
	and level **5.0** is specified as "5.0" or "50".

	Annex A levels: 1, 2, 2.1, 3, 3.1, 4, 4.1, 5, 5.1, 5.2, 6, 6.1, 6.2

.. option:: --high-tier, --no-high-tier

	If :option:`--level-idc` has been specified, the option adds the
	intention to support the High tier of that level. If your specified
	level does not support a High tier, a warning is issued and this
	modifier flag is ignored.

.. note::
	:option:`--profile`, :option:`--level-idc`, and
	:option:`--high-tier` are only intended for use when you are
	targeting a particular decoder (or decoders) with fixed resource
	limitations and must constrain the bitstream within those limits.
	Specifying a profile or level may lower the encode quality
	parameters to meet those requirements but it will never raise
	them.

Quad-Tree analysis
==================

.. option:: --wpp, --no-wpp

	Enable Wavefront Parallel Processing. The encoder may begin encoding
	a row as soon as the row above it is at least two CTUs ahead in the
	encode process. This gives a 3-5x gain in parallelism for about 1%
	overhead in compression efficiency. Default: Enabled

.. option:: --ctu, -s <64|32|16>

	Maximum CU size (width and height). The larger the maximum CU size,
	the more efficiently x265 can encode flat areas of the picture,
	giving large reductions in bitrate. However this comes at a loss of
	parallelism with fewer rows of CUs that can be encoded in parallel,
	and less frame parallelism as well. Because of this the faster
	presets use a CU size of 32. Default: 64

.. option:: --tu-intra-depth <1..4>

	The transform unit (residual) quad-tree begins with the same depth
	as the coding unit quad-tree, but the encoder may decide to further
	split the transform unit tree if it improves compression efficiency.
	This setting limits the number of extra recursion depth which can be
	attempted for intra coded units. Default: 1, which means the
	residual quad-tree is always at the same depth as the coded unit
	quad-tree
	
	Note that when the CU intra prediction is NxN (only possible with
	8x8 CUs), a TU split is implied, and thus the residual quad-tree
	begins at 4x4 and cannot split any futhrer.

.. option:: --tu-inter-depth <1..4>

	The transform unit (residual) quad-tree begins with the same depth
	as the coding unit quad-tree, but the encoder may decide to further
	split the transform unit tree if it improves compression efficiency.
	This setting limits the number of extra recursion depth which can be
	attempted for inter coded units. Default: 1. which means the
	residual quad-tree is always at the same depth as the coded unit
	quad-tree unless the CU was coded with rectangular or AMP
	partitions, in which case a TU split is implied and thus the
	residual quad-tree begins one layer below the CU quad-tree.

Temporal / motion search options
================================

.. option:: --me <integer|string>

	Motion search method. Generally, the higher the number the harder
	the ME method will try to find an optimal match. Diamond search is
	the simplest. Hexagon search is a little better. Uneven
	Multi-Hexegon is an adaption of the search method used by x264 for
	slower presets. Star is a three step search adapted from the HM
	encoder: a star-pattern search followed by an optional radix scan
	followed by an optional star-search refinement. Full is an
	exhaustive search; an order of magnitude slower than all other
	searches but not much better than umh or star.

	0. dia
	1. hex **(default)**
	2. umh
	3. star
	4. full

.. option:: --subme, -m <0..7>

	Amount of subpel refinement to perform. The higher the number the
	more subpel iterations and steps are performed. Default 2

	+----+------------+-----------+------------+-----------+-----------+
	| -m | HPEL iters | HPEL dirs | QPEL iters | QPEL dirs | HPEL SATD |
	+====+============+===========+============+===========+===========+
	|  0 | 1          | 4         | 0          | 4         | false     |
	+----+------------+-----------+------------+-----------+-----------+
	|  1 | 1          | 4         | 1          | 4         | false     |
	+----+------------+-----------+------------+-----------+-----------+
	|  2 | 1          | 4         | 1          | 4         | true      |
	+----+------------+-----------+------------+-----------+-----------+
	|  3 | 2          | 4         | 1          | 4         | true      |
	+----+------------+-----------+------------+-----------+-----------+
	|  4 | 2          | 4         | 2          | 4         | true      |
	+----+------------+-----------+------------+-----------+-----------+
	|  5 | 1          | 8         | 1          | 8         | true      |
	+----+------------+-----------+------------+-----------+-----------+
	|  6 | 2          | 8         | 1          | 8         | true      |
	+----+------------+-----------+------------+-----------+-----------+
	|  7 | 2          | 8         | 2          | 8         | true      |
	+----+------------+-----------+------------+-----------+-----------+

.. option:: --merange <integer>

	Motion search range. Default 57

	The default is derived from the default CTU size (64) minus the luma
	interpolation half-length (4) minus maximum subpel distance (2)
	minus one extra pixel just in case the hex search method is used. If
	the search range were any larger than this, another CTU row of
	latency would be required for reference frames.

	**Range of values:** an integer from 0 to 32768

.. option:: --max-merge <1..5>

	Maximum number of neighbor (spatial and temporal) candidate blocks
	that the encoder may consider for merging motion predictions. If a
	merge candidate results in no residual, it is immediately selected
	as a "skip".  Otherwise the merge candidates are tested as part of
	motion estimation when searching for the least cost inter option.
	The max candidate number is encoded in the SPS and determines the
	bit cost of signaling merge CUs. Default 2

.. option:: --temporal-mvp, --no-temporal-mvp

	Enable temporal motion vector predictors in P and B slices.
	This enables the use of the motion vector from the collocated block
	in the previous frame to be used as a predictor. Default is enabled

Spatial/intra options
=====================

.. option:: --rdpenalty <0..2>

	When set to 1, transform units of size 32x32 are given a 4x bit cost
	penalty compared to smaller transform units, in intra coded CUs in P
	or B slices.

	When set to 2, transform units of size 32x32 are not even attempted,
	unless otherwise required by the maximum recursion depth.  For this
	option to be effective with 32x32 intra CUs,
	:option:`--tu-intra-depth` must be at least 2.  For it to be
	effective with 64x64 intra CUs, :option:`--tu-intra-depth` must be
	at least 3.

	Note that in HEVC an intra transform unit (a block of the residual
	quad-tree) is also a prediction unit, meaning that the intra
	prediction signal is generated for each TU block, the residual
	subtracted and then coded. The coding unit simply provides the
	prediction modes that will be used when predicting all of the
	transform units within the CU. This means that when you prevent
	32x32 intra transform units, you are preventing 32x32 intra
	predictions.

	Default 0, disabled.

	**Values:** 0:disabled 1:4x cost penalty 2:force splits

.. option:: --b-intra, --no-b-intra

	Enables the evaluation of intra modes in B slices. Default disabled.

.. option:: --tskip, --no-tskip

	Enable evaluation of transform skip (bypass DCT but still use
	quantization) coding for 4x4 TU coded blocks.

	Only effective at RD levels 3 and above, which perform RDO mode
	decisions. Default disabled

.. option:: --tskip-fast, --no-tskip-fast

	Only evaluate transform skip for NxN intra predictions (4x4 blocks).
	Only applicable if transform skip is enabled. For chroma, only
	evaluate if luma used tskip. Inter block tskip analysis is
	unmodified. Default disabled

.. option:: --strong-intra-smoothing, --no-strong-intra-smoothing

	Enable strong intra smoothing for 32x32 intra blocks. Default enabled

.. option:: --constrained-intra, --no-constrained-intra

	Constrained intra prediction. When generating intra predictions for
	blocks in inter slices, only intra-coded reference pixels are used.
	Inter-coded reference pixels are replaced with intra-coded neighbor
	pixels or default values. The general idea is to block the
	propagation of reference errors that may have resulted from lossy
	signals. Default disabled

Mode decision / Analysis
========================

.. option:: --rect, --no-rect

	Enable analysis of rectangular motion partitions Nx2N and 2NxN
	(50/50 splits, two directions). Default disabled

.. option:: --amp, --no-amp

	Enable analysis of asymmetric motion partitions (75/25 splits, four
	directions). At RD levels 0 through 4, AMP partitions are only
	considered at CU sizes 32x32 and below. At RD levels 5 and 6, it
	will only consider AMP partitions as merge candidates (no motion
	search) at 64x64, and as merge or inter candidates below 64x64.

	The AMP partitions which are searched are derived from the current
	best inter partition. If Nx2N (vertical rectangular) is the best
	current prediction, then left and right asymmetrical splits will be
	evaluated. If 2NxN (horizontal rectangular) is the best current
	prediction, then top and bottom asymmetrical splits will be
	evaluated, If 2Nx2N is the best prediction, and the block is not a
	merge/skip, then all four AMP partitions are evaluated.

	This setting has no effect if rectangular partitions are disabled.
	Default disabled

.. option:: --early-skip, --no-early-skip

	Measure full CU size (2Nx2N) merge candidates first; if no residual
	is found the analysis is short circuited. Default disabled

.. option:: --fast-cbf, --no-fast-cbf

	Short circuit analysis if a prediction is found that does not set
	the coded block flag (aka: no residual was encoded).  It prevents
	the encoder from perhaps finding other predictions that also have no
	residual but require less signaling bits or have less distortion.
	Only applicable for RD levels 5 and 6. Default disabled

.. option:: --fast-intra, --no-fast-intra

	Perform an initial scan of every fifth intra angular mode, then
	check modes +/- 2 distance from the best mode, then +/- 1 distance
	from the best mode, effectively performing a gradient descent. When
	enabled 10 modes in total are checked. When disabled all 33 angular
	modes are checked.  Only applicable for :option:`--rd` levels 3 and
	below (medium preset and faster).

.. option:: --weightp, -w, --no-weightp

	Enable weighted prediction in P slices. This enables weighting
	analysis in the lookahead, which influences slice decisions, and
	enables weighting analysis in the main encoder which allows P
	reference samples to have a weight function applied to them prior to
	using them for motion compensation.  In video which has lighting
	changes, it can give a large improvement in compression efficiency.
	Default is enabled

.. option:: --weightb, --no-weightb

	Enable weighted prediction in B slices. Default disabled

.. option:: --rd <0..6>

	Level of RDO in mode decision. The higher the value, the more
	exhaustive the analysis and the more rate distortion optimization is
	used. The lower the value the faster the encode, the higher the
	value the smaller the bitstream (in general). Default 3

	Note that this table aims for accuracy, but is not necessarily our
	final target behavior for each mode.

	+-------+---------------------------------------------------------------+
	| Level | Description                                                   |
	+=======+===============================================================+
	| 0     | sa8d mode and split decisions, intra w/ source pixels         |
	+-------+---------------------------------------------------------------+
	| 1     | recon generated (better intra), RDO merge/skip selection      |
	+-------+---------------------------------------------------------------+
	| 2     | RDO splits and merge/skip selection                           |
	+-------+---------------------------------------------------------------+
	| 3     | RDO mode and split decisions                                  |
	+-------+---------------------------------------------------------------+
	| 4     | Adds RDO Quant                                                |
	+-------+---------------------------------------------------------------+
	| 5     | Adds RDO prediction decisions                                 |
	+-------+---------------------------------------------------------------+
	| 6     | Currently same as 5                                           |
	+-------+---------------------------------------------------------------+

	**Range of values:** 0: least .. 6: full RDO analysis

.. option:: --cu-lossless, --no-cu-lossless

	For each CU, evaluate lossless (transform and quant bypass) encode
	of the best non-lossless mode option as a potential rate distortion
	optimization. If the global option :option:`--lossless` has been
	specified, all CUs will be encoded as lossless unconditionally
	regardless of whether this option was enabled. Default disabled.

	Only effective at RD levels 3 and above, which perform RDO mode
	decisions.

.. option:: --signhide, --no-signhide

	Hide sign bit of one coeff per TU (rdo). The last sign is implied.
	This requires analyzing all the coefficients to determine if a sign
	must be toggled, and then to determine which one can be toggled with
	the least amount of distortion. Default enabled
 
Psycho-visual options
=====================

Left to its own devices, the encoder will make mode decisions based on a
simple rate distortion formula, trading distortion for bitrate. This is
generally effective except for the manner in which this distortion is
measured. It tends to favor blurred reconstructed blocks over blocks
which have wrong motion. The human eye generally prefers the wrong
motion over the blur and thus x265 offers psycho-visual adjustments to
the rate distortion algorithm.

:option:`--psy-rd` will add an extra cost to reconstructed blocks which
do not match the visual energy of the source block. The higher the
strength of :option:`--psy-rd` the more strongly it will favor similar
energy over blur and the more aggressively it will ignore rate
distortion. If it is too high, it will introduce visal artifacts and
increase bitrate enough for rate control to increase quantization
globally, reducing overall quality. psy-rd will tend to reduce the use
of blurred prediction modes, like DC and planar intra and bi-directional
inter prediction.

:option:`--psy-rdoq` will adjust the distortion cost used in
rate-distortion optimized quantization (RDO quant), enabled in
:option:`--rd` 4 and above, favoring the preservation of energy in the
reconstructed image.  :option:`--psy-rdoq` prevents RDOQ from blurring
all of the encoding options which psy-rd has to chose from.  At low
strength levels, psy-rdoq will influence the quantization level
decisions, favoring higher AC energy in the reconstructed image. As
psy-rdoq strength is increased, more non-zero coefficient levels are
added and fewer coefficients are zeroed by RDOQ's rate distortion
analysis. High levels of psy-rdoq can double the bitrate which can have
a drastic effect on rate control, forcing higher overall QP, and can
cause ringing artifacts. psy-rdoq is less accurate than psy-rd, it is
biasing towards energy in general while psy-rd biases towards the energy
of the source image. But very large psy-rdoq values can sometimes be
beneficial, preserving film grain for instance.

As a general rule, when both psycho-visual features are disabled, the
encoder will tend to blur blocks in areas of difficult motion. Turning
on small amounts of psy-rd and psy-rdoq will improve the perceived
visual quality. Increasing psycho-visual strength further will improve
quality and begin introducing artifacts and increase bitrate, which may
force rate control to increase global QP. Finding the optimal
psycho-visual parameters for a given video requires experimentation. Our
recommended defaults (1.0 for both) are generally on the low end of the
spectrum. And generally the lower the bitrate, the lower the optimal
psycho-visual settings.

.. option:: --psy-rd <float>

	Influence rate distortion optimizated mode decision to preserve the
	energy of the source image in the encoded image at the expense of
	compression efficiency. It only has effect on presets which use
	RDO-based mode decisions (:option:`--rd` 3 and above).  1.0 is a
	typical value. Default disabled.  Experimental

	**Range of values:** 0 .. 2.0

.. option:: --psy-rdoq <float>

	Influence rate distortion optimized quantization by favoring higher
	energy in the reconstructed image. This generally improves perceived
	visual quality at the cost of lower quality metric scores.  It only
	has effect on slower presets which use RDO Quantization
	(:option:`--rd` 4, 5 and 6). 1.0 is a typical value. Default
	disabled. High values can be beneficial in preserving high-frequency
	detail like film grain. Experimental

	**Range of values:** 0 .. 50.0


Slice decision options
======================

.. option:: --open-gop, --no-open-gop

	Enable open GOP, allow I-slices to be non-IDR. Default enabled

.. option:: --keyint, -I <integer>

	Max intra period in frames. A special case of infinite-gop (single
	keyframe at the beginning of the stream) can be triggered with
	argument -1. Use 1 to force all-intra. Default 250

.. option:: --min-keyint, -i <integer>

	Minimum GOP size. Scenecuts closer together than this are coded as I
	or P, not IDR. Minimum keyint is clamped to be at least half of
	:option:`--keyint`. If you wish to force regular keyframe intervals
	and disable adaptive I frame placement, you must use
	:option:`--no-scenecut`.

	**Range of values:** >=0 (0: auto)

.. option:: --scenecut <integer>, --no-scenecut

	How aggressively I-frames need to be inserted. The higher the
	threshold value, the more aggressive the I-frame placement.
	:option:`--scenecut` 0 or :option:`--no-scenecut` disables adaptive
	I frame placement. Default 40

.. option:: --rc-lookahead <integer>

	Number of frames for slice-type decision lookahead (a key
	determining factor for encoder latency). The longer the lookahead
	buffer the more accurate scenecut decisions will be, and the more
	effective cuTree will be at improving adaptive quant. Having a
	lookahead larger than the max keyframe interval is not helpful.
	Default 20

	**Range of values:** Between the maximum consecutive bframe count (:option:`--bframes`) and 250

.. option:: --b-adapt <integer>

	Adaptive B frame scheduling. Default 2

	**Values:** 0:none; 1:fast; 2:full(trellis)

.. option:: --bframes, -b <0..16>

	Maximum number of consecutive b-frames. Use :option:`--bframes` 0 to
	force all P/I low-latency encodes. Default 4. This parameter has a
	quadratic effect on the amount of memory allocated and the amount of
	work performed by the full trellis version of :option:`--b-adapt`
	lookahead.

.. option:: --bframe-bias <integer>

	Bias towards B frames in slicetype decision. The higher the bias the
	more likely x265 is to use B frames. Can be any value between -90
	and 100 and is clipped to that range. Default 0

.. option:: --b-pyramid, --no-b-pyramid

	Use B-frames as references, when possible. Default enabled

.. option:: --ref <1..16>

	Max number of L0 references to be allowed. This number has a linear
	multiplier effect on the amount of work performed in motion search,
	but will generally have a beneficial affect on compression and
	distortion. Default 3

Quality, rate control and rate distortion options
=================================================

.. option:: --bitrate <integer>

	Enables single-pass ABR rate control. Specify the target bitrate in
	kbps. Default is 0 (CRF)

	**Range of values:** An integer greater than 0

.. option:: --crf <0..51.0>

	Quality-controlled variable bitrate. CRF is the default rate control
	method; it does not try to reach any particular bitrate target,
	instead it tries to achieve a given uniform quality and the size of
	the bitstream is determined by the complexity of the source video.
	The higher the rate factor the higher the quantization and the lower
	the quality. Default rate factor is 28.0.

.. option:: --crf-max <0..51.0>

	Specify an upper limit to the rate factor which may be assigned to
	any given frame (ensuring a max QP).  This is dangerous when CRF is
	used in combination with VBV as it may result in buffer underruns.
	Default disabled
        
.. option:: --crf-min <0..51.0>

	Specify an lower limit to the rate factor which may be assigned to
	any given frame (ensuring a min QP).  This is dangerous when CRF is
	used in combination with VBV as it may result in buffer underruns.
	Default disabled

.. option:: --vbv-bufsize <integer>

	Specify the size of the VBV buffer (kbits). Enables VBV in ABR
	mode.  In CRF mode, :option:`--vbv-maxrate` must also be specified.
	Default 0 (vbv disabled)

.. option:: --vbv-maxrate <integer>

	Maximum local bitrate (kbits/sec). Will be used only if vbv-bufsize
	is also non-zero. Both vbv-bufsize and vbv-maxrate are required to
	enable VBV in CRF mode. Default 0 (disabled)

.. option:: --vbv-init <float>

	Initial buffer occupancy. The portion of the decode buffer which
	must be full before the decoder will begin decoding.  Determines
	absolute maximum frame size. May be specified as a fractional value
	between 0 and 1, or in kbits. In other words these two option pairs
	are equivalent::

	:option:`--vbv-bufsize` 1000 :option:`--vbv-init` 900
	:option:`--vbv-bufsize` 1000 :option:`--vbv-init` 0.9

	Default 0.9

	**Range of values:** fractional: 0 - 1.0, or kbits: 2 .. bufsize

.. option:: --qp, -q <integer>

	Specify base quantization parameter for Constant QP rate control.
	Using this option enables Constant QP rate control. The specified QP
	is assigned to P slices. I and B slices are given QPs relative to P
	slices using param->rc.ipFactor and param->rc.pbFactor unless QP 0
	is specified, in which case QP 0 is used for all slice types.  Note
	that QP 0 does not cause lossless encoding, it only disables
	quantization. Default disabled (CRF)

	**Range of values:** an integer from 0 to 51

.. option:: --ipratio <float>

	QP ratio factor between I and P slices. This ratio is used in all of
	the rate control modes. Some :option:`--tune` options may change the
	default value. It is not typically manually specified. Default 1.4

.. option:: --pbratio <float>

	QP ratio factor between P and B slices. This ratio is used in all of
	the rate control modes. Some :option:`--tune` options may change the
	default value. It is not typically manually specified. Default 1.3

.. option:: --lossless, --no-lossless

	Enables true lossless coding by bypassing scaling, transform,
	quantization and in-loop filter processes. This is used for
	ultra-high bitrates with zero loss of quality. Reconstructed output
	pictures are bit-exact to the input pictures. Lossless encodes
	implicitly have no rate control, all rate control options are
	ignored. Slower presets will generally achieve better compression
	efficiency (and generate smaller bitstreams). Default disabled.

.. option:: --aq-mode <0|1|2>

	Adaptive Quantization operating mode. Raise or lower per-block
	quantization based on complexity analysis of the source image. The
	more complex the block, the more quantization is used. This offsets
	the tendency of the encoder to spend too many bits on complex areas
	and not enough in flat areas.

	0. disabled
	1. AQ enabled
	2. AQ enabled with auto-variance **(default)**

.. option:: --aq-strength <float>

	Adjust the strength of the adaptive quantization offsets. Setting
	:option:`--aq-strength` to 0 disables AQ. Default 1.0.

	**Range of values:** 0.0 to 3.0

.. option:: --cutree, --no-cutree

	Enable the use of lookahead's lowres motion vector fields to
	determine the amount of reuse of each block to tune adaptive
	quantization factors. CU blocks which are heavily reused as motion
	reference for later frames are given a lower QP (more bits) while CU
	blocks which are quickly changed and are not referenced are given
	less bits. This tends to improve detail in the backgrounds of video
	with less detail in areas of high motion. Default enabled

.. option:: --cbqpoffs <integer>

	Offset of Cb chroma QP from the luma QP selected by rate control.
	This is a general way to spend more or less bits on the chroma
	channel.  Default 0

	**Range of values:** -12 to 12

.. option:: --crqpoffs <integer>

	Offset of Cr chroma QP from the luma QP selected by rate control.
	This is a general way to spend more or less bits on the chroma
	channel.  Default 0

	**Range of values:**  -12 to 12

.. option:: --pass <integer>

	Enable multipass rate control mode. Input is encoded multiple times,
	storing the encoded information of each pass in a stats file from which
	the consecutive pass tunes the qp of each frame to improve the quality
	of the output. Default disabled

	1. First pass, creates stats file
	2. Last pass, does not overwrite stats file
	3. Nth pass, overwrites stats file

	**Range of values:** 1 to 3

.. option:: --slow-firstpass, --no-slow-firstpass

	Enable a slow and more detailed first pass encode in Multipass rate
	control mode.  Speed of the first pass encode is slightly lesser and
	quality midly improved when compared to the default settings in a
	multipass encode. Default disabled (turbo mode enabled)

	When **turbo** first pass is not disabled, these options are
	set on the first pass to improve performance:
	
	* :option:`--fast-intra`
	* :option:`--no-rect`
	* :option:`--no-amp`
	* :option:`--early-skip`
	* :option:`--ref` = 1
	* :option:`--max-merge` = 1
	* :option:`--me` = DIA
	* :option:`--subme` = MIN(2, :option:`--subme`)
	* :option:`--rd` = MIN(2, :option:`--rd`)

.. option:: --analysis-mode <string|int>

	Specify whether analysis information of each frame is output by encoder
	or input for reuse. By reading the analysis data writen by an
	earlier encode of the same sequence, substantial redundant work may
	be avoided.

	The following data may be stored and reused:
	I frames   - split decisions and luma intra directions of all CUs.
	P/B frames - motion vectors are dumped at each depth for all CUs.

	**Values:** off(0), save(1): dump analysis data, load(2): read analysis data

.. option:: --analysis-file <filename>

	Specify a filename for analysis data (see :option:`--analysis-mode`)
	If no filename is specified, x265_analysis.dat is used.

Loop filters
============

.. option:: --lft, --no-lft

	Toggle deblocking loop filter, default enabled

.. option:: --sao, --no-sao

	Toggle Sample Adaptive Offset loop filter, default enabled

.. option:: --sao-non-deblock, --no-sao-non-deblock

	Specify how to handle depencency between SAO and deblocking filter.
	When enabled, non-deblocked pixels are used for SAO analysis. When
	disabled, SAO analysis skips the right/bottom boundary areas.
	Default disabled

VUI (Video Usability Information) options
=========================================

x265 emits a VUI with only the timing info by default. If the SAR is
specified (or read from a Y4M header) it is also included.  All other
VUI fields must be manually specified.

.. option:: --sar <integer|w:h>

	Sample Aspect Ratio, the ratio of width to height of an individual
	sample (pixel). The user may supply the width and height explicitly
	or specify an integer from the predefined list of aspect ratios
	defined in the HEVC specification.  Default undefined (not signaled)

	1. 1:1 (square)
	2. 12:11
	3. 10:11
	4. 16:11
	5. 40:33
	6. 24:11
	7. 20:11
	8. 32:11
	9. 80:33
	10. 18:11
	11. 15:11
	12. 64:33
	13. 160:99
	14. 4:3
	15. 3:2
	16. 2:1

.. option:: --crop-rect <left,top,right,bottom>

	Define the (overscan) region of the image that does not contain
	information because it was added to achieve certain resolution or
	aspect ratio. The decoder may be directed to crop away this region
	before displaying the images via the :option:`--overscan` option.
	Default undefined (not signaled)

.. option:: --overscan <show|crop>

	Specify whether it is appropriate for the decoder to display or crop
	the overscan area. Default unspecified (not signaled)

.. option:: --videoformat <integer|string>

	Specify the source format of the original analog video prior to
	digitizing and encoding. Default undefined (not signaled)

	0. component
	1. pal
	2. ntsc
	3. secam
	4. mac
	5. undefined

.. option:: --range <full|limited>

	Specify output range of black level and range of luma and chroma
	signals. Default undefined (not signaled)

.. option:: --colorprim <integer|string>

	Specify color primitive to use when converting to RGB. Default
	undefined (not signaled)

	1. bt709
	2. undef
	3. **reserved**
	4. bt470m
	5. bt470bg
	6. smpte170m
	7. smpte240m
	8. film
	9. bt2020

.. option:: --transfer <integer|string>

	Specify transfer characteristics. Default undefined (not signaled)

	1. bt709
	2. undef
	3. **reserved**
	4. bt470m
	5. bt470bg
	6. smpte170m
	7. smpte240m
	8. linear
	9. log100
	10. log316
	11. iec61966-2-4
	12. bt1361e
	13. iec61966-2-1
	14. bt2020-10
	15. bt2020-12

.. option:: --colormatrix <integer|string>

	Specify color matrix setting i.e set the matrix coefficients used in
	deriving the luma and chroma. Default undefined (not signaled)

	0. GBR
	1. bt709
	2. undef 
	3. **reserved**
	4. fcc
	5. bt470bg
	6. smpte170m
	7. smpte240m
	8. YCgCo
	9. bt2020nc
	10. bt2020c

.. option:: --chromalocs <0..5>

	Specify chroma sample location for 4:2:0 inputs. Consult the HEVC
	specification for a description of these values. Default undefined
	(not signaled)

Bitstream options
=================

.. option:: --repeat-headers, --no-repeat-headers

	If enabled, x265 will emit VPS, SPS, and PPS headers with every
	keyframe. This is intended for use when you do not have a container
	to keep the stream headers for you and you want keyframes to be
	random access points. Default disabled

.. option:: --info, --no-info

	Emit an informational SEI with the stream headers which describes
	the encoder version, build info, and encode parameters. This is very
	helpful for debugging purposes but encoding version numbers and
	build info could make your bitstreams diverge and interfere with
	regression testing. Default enabled

.. option:: --hrd, --no-hrd

	Enable the signalling of HRD parameters to the decoder. The HRD
	parameters are carried by the Buffering Period SEI messages and
	Picture Timing SEI messages providing timing information to the
	decoder. Default disabled

.. option:: --aud, --no-aud

	Emit an access unit delimiter NAL at the start of each slice access
	unit. If option:`--repeat-headers` is not enabled (indicating the
	user will be writing headers manually at the start of the stream)
	the very first AUD will be skipped since it cannot be placed at the
	start of the access unit, where it belongs. Default disabled

.. option:: --hash <integer>

	Emit decoded picture hash SEI, so the decoder may validate the
	reconstructed pictures and detect data loss. Also useful as a
	debug feature to validate the encoder state. Default None

	1. MD5
	2. CRC
	3. Checksum

Debugging options
=================

.. option:: --recon, -r <filename>

	Output file containing reconstructed images in display order. If the
	file extension is ".y4m" the file will contain a YUV4MPEG2 stream
	header and frame headers. Otherwise it will be a raw YUV file in the
	encoder's internal bit depth.

	**CLI ONLY**

.. option:: --recon-depth <integer>

	Bit-depth of output file. This value defaults to the internal bit
	depth and currently cannot to be modified.

	**CLI ONLY**

.. vim: noet
