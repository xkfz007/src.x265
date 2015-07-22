// This file is automatically written by compilation scripts.

// Force the use of intra vertical 16x16. Soon-to-be obsolete and badly-named
// (backward compatibility), use only for regression.
//#define VAN_DUMP_INTRA

// Print the bitstream as it is encoded.
//#define VAN_TRACE_SYNTAX

// Load the CTB analysis.
//#define VAN_LOAD_CTB_ANALYSIS

// Verify if the analysis matches the encoding (full RDO case only).
//#define VAN_VERIFY_RDO_ANALYSIS

// Debug the deblocking filter with HM.
//#define VAN_DUMP_DEBLOCK_FILTER

// Test inter PMV, inter merge candidate and intra MPM.
// Requires VAN_LOAD_CTB_ANALYSIS to be defined.
//#define VAN_VALIDATE_MODE_PRED

// Trace the frame encoding progression.
//#define VAN_TRACE_FRAME_ENCODE

// Force the use of Full Search over TZ search.
//#define VAN_HM_RDO_USE_FS

// Trace the motion vector costs (Full Search).
//#define VAN_TRACE_FS_MV_COST

