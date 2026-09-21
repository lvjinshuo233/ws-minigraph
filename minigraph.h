#ifndef MINIGRAPH2_H
#define MINIGRAPH2_H

#include <stdint.h>
#include "gfa.h"
#include "bloom_filter/bloom_filter.hpp"
#include "collect-I.h"
#define MG_VERSION "0.21-r606"
#define MM_I_HPC          0x1
#define MG_M_SPLICE       0x10
#define MG_M_SR           0x20
#define MG_M_FRAG_MODE    0x40
#define MG_M_FRAG_MERGE   0x80
#define MG_M_FOR_ONLY     0x100
#define MG_M_REV_ONLY     0x200
#define MG_M_HEAP_SORT    0x400
#define MG_M_VERTEX_COOR  0x800
#define MG_M_ALL_CHAINS   0x1000
#define MG_M_PRINT_2ND    0x2000
#define MG_M_CAL_COV      0x4000
#define MG_M_RMQ          0x8000
#define MG_M_COPY_COMMENT 0x10000
#define MG_M_INDEPEND_SEG 0x20000
#define MG_M_NO_QUAL      0x40000
#define MG_M_2_IO_THREADS 0x80000
#define MG_M_SHOW_UNMAP   0x100000
#define MG_M_NO_COMP_PATH 0x200000
#define MG_M_NO_DIAG      0x400000
#define MG_M_WRITE_LCHAIN 0x800000
#define MG_M_WRITE_MZ     0x1000000
#define MG_M_SKIP_GCHECK  0x2000000
#define MG_M_CIGAR        0x4000000

#define MG_G_NONE         0
#define MG_G_GGSIMPLE     1

#define MG_G_NO_QOVLP     0x1
#define MG_G_CAL_COV      0x2
#define MG_G_NO_INV       0x4
#define MG_G_CALL         0x8

#define SUPER_WINDOW_SIZE 1000000
extern int buf_weight, cur_weights;

extern long int n_weights, n_weights_centro, n_sync_centro;
extern long int n_syncmers, n_w_k, window, n_strobes;
extern double index_start_time,	index_end_time, map_start_time, map_end_time;
extern FILE *f_freq;
extern FILE *f_seq;
extern FILE *f_check;
extern FILE *f_samplingInfo;
extern FILE *f_rep;
extern long int sync_distance, dist_temp, dist_count;
extern long int lines, bases;
typedef struct { uint64_t x, y; } mg128_t;
typedef struct { size_t n, m; mg128_t *a; } mg128_v;
typedef struct { int32_t n, m; uint32_t *a; } mg32_v;
typedef struct { int32_t n, m; uint64_t *a; } mg64_v;

typedef struct {
	int w, k, s, t, syncmer, min_rep_len, max_rep_len;
	int bucket_bits;
} mg_idxopt_t;

typedef struct{
	long int n_kmers, n_bases;
	int n_windows;
	int n_sketch;
	int n_types_kmers;
	char *seq;
	long int kmer_aft_filtered;
	long int kmer;
}SuperWindow;

extern SuperWindow SW;

typedef struct {
	uint64_t flag;
	int64_t mini_batch_size;
	int seed;
	int max_qlen;
	int pe_ori;
	int occ_max1, occ_max1_cap;
	float occ_max1_frac;
	int bw, bw_long;
	int rmq_size_cap;
	int rmq_rescue_size;
	float rmq_rescue_ratio;
	int max_gap_pre, max_gap, max_gap_ref, max_frag_len;
	float div;
	float chn_pen_gap, chn_pen_skip;
	int max_lc_skip, max_lc_iter, max_gc_skip;
	int min_lc_cnt, min_lc_score;
	int min_gc_cnt, min_gc_score;
	int gdp_max_ed, lc_max_trim, lc_max_occ;
	float mask_level;
	int sub_diff;
	int best_n;
	float pri_ratio;
	int ref_bonus;
	int64_t cap_kalloc;
	int min_cov_mapq, min_cov_blen;
	const char *kmer_freq_filename; //file name containing k-mer frequencies
} mg_mapopt_t;

typedef struct {
	uint64_t flag;
	int algo;
	int min_mapq;
	int min_map_len, min_depth_len;
	int min_var_len, match_pen;
	// parameters specific to ggsimple/ggs
	int ggs_shrink_pen;
	int ggs_min_end_cnt;
	float ggs_min_end_frac;
	// scoring for SW check
	float ggs_max_iden, ggs_min_inv_iden;
} mg_ggopt_t;

typedef struct {
	int syncmer;
	const gfa_t *g;
	gfa_edseq_t *es;
	int32_t b, w, k, flag, n_seg, s, t, min_rep_len, max_rep_len;
	struct mg_idx_bucket_s *B; // index (hidden)
	bloom_filter* downFilter;  // bloom filter for down-weighted kmers
	void *strobe_idx;          // 新增：strobemer 索引指针（strobe_hash_t*）
	mg128_v all_syncs;   // 新增：保存所有参考 syncmer 的有序数组
} mg_idx_t;

typedef struct {
	int32_t off, cnt:31, inner_pre:1;
	uint32_t v;
	int32_t rs, re, qs, qe;
	int32_t score, dist_pre;
	uint32_t hash_pre;
} mg_lchain_t;

typedef struct {
	int32_t off, cnt;
	uint32_t v;
	int32_t score;
	int32_t ed;
} mg_llchain_t;

typedef struct {
	int32_t n_cigar, mlen, blen, aplen, ss, ee; // ss: start on the start vertex; ee: end on the end vertex
	uint64_t cigar[];
} mg_cigar_t;

typedef struct {
	int32_t len, n_off, *off;
	char *ds;
} mg_ds_t;

typedef struct {
	int32_t id, parent;
	int32_t off, cnt;
	int32_t n_anchor, score;
	int32_t qs, qe;
	int32_t plen, ps, pe;
	int32_t blen, mlen;
	float div;
	uint32_t hash;
	int32_t subsc, n_sub;
	uint32_t mapq:8, flt:1, dummy:23;
	mg_cigar_t *p;
	mg_ds_t ds;
} mg_gchain_t;

typedef struct {
	void *km;
	int32_t n_gc, n_lc, n_a, rep_len;
	mg_gchain_t *gc;
	mg_llchain_t *lc;
	mg128_t *a; // minimizer positions; see comments above mg_update_anchors() for details
	int32_t *filt_pos;   // 查询上被过滤掉的 minimizer 的位置（0‑based）
    int32_t n_filt;      // 被过滤的数量
	    int32_t n_mini;      // 总 minimizer 数
    int32_t *all_mini_pos;   // 新增：每个minimizer的查询起始位置（0‑based）
    int32_t *all_mini_info;  // 新增：bit0=过滤标志（1=过滤），高位=出现次数（occ）
    int64_t *all_ref_pos;   // 修改类型：存储 cr[0] (seg_id<<32 | (pos<<1) | strand)

} mg_gchains_t;

typedef struct mg_tbuf_s mg_tbuf_t;

extern int mg_verbose, mg_dbg_flag;
extern double mg_realtime0;

#ifdef __cplusplus
extern "C" {
#endif

// options
int mg_opt_set(const char *preset, mg_idxopt_t *io, mg_mapopt_t *mo, mg_ggopt_t *go);
int mg_opt_check(const mg_idxopt_t *io, const mg_mapopt_t *mo, const mg_ggopt_t *go);
void mg_opt_update(const mg_idx_t *gi, mg_mapopt_t *mo, mg_ggopt_t *go);

// index operations
mg_idx_t *mg_index(gfa_t *g, const mg_idxopt_t *io, int n_threads, mg_mapopt_t *mo, const char* kmer_freq_filename, int flag, collect_rep_t *collect_rep); // combine mg_index_core() and mg_opt_update()
void mg_idx_destroy(mg_idx_t *gi);

// mapping
mg_tbuf_t *mg_tbuf_init(void);
void mg_tbuf_destroy(mg_tbuf_t *b);
mg_gchains_t *mg_map(const mg_idx_t *gi, int qlen, const char *seq, mg_tbuf_t *b, const mg_mapopt_t *opt, const char *qname);
void mg_map_frag(const mg_idx_t *gi, int n_segs, const int *qlens, const char **seqs, mg_gchains_t **gcs, mg_tbuf_t *b, const mg_mapopt_t *opt, const char *qname, collect_rep_t *collect_rep);

// high-level mapping APIs
int mg_map_files(gfa_t *g, int n_fn, const char **fn, const mg_idxopt_t *ipt, const mg_mapopt_t *opt0, int n_threads, const char* kmer_freq_filename);

// graph generation
int mg_ggen(gfa_t *g, int32_t n_fn, const char **fn, const mg_idxopt_t *ipt, const mg_mapopt_t *opt0, const mg_ggopt_t *go, int n_threads);

#ifdef __cplusplus
}
#endif

#endif
