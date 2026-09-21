#ifndef WS_STROBE_H
#define WS_STROBE_H

#include <stdint.h>
#include "mgpriv.h"      // 需要 mg_idx_t, mg128_v, mg128_t 等
#include "collect-I.h"   // 需要 seg_rep_t, rep_interval_t

#ifdef __cplusplus
extern "C" {
#endif

typedef struct strobe_hash_s strobe_hash_t;

/* ---------- 新增：比对阶段使用的 strobemer 结构 ---------- */
typedef struct {
    uint64_t h1, h2;
    int32_t s1_start, s2_start;
    int32_t dist;
    int32_t freq1, freq2;
    int32_t s1_qspan, s2_qspan;   // q_span 信息
} read_strobe_t;

/* 查询 read 上的 strobemer，生成 read_strobe_t 数组 */
int collect_ws_strobes(const mg_idx_t *gi, const mg128_v *mv, seg_rep_t *reads_rep,
                       read_strobe_t **out_strobes);

/* 释放 read_rep 结构（包括链表） */
void free_reads_rep(seg_rep_t *reads_rep, int n_segs);

/* 将 read 上的 strobemer 与参考索引匹配，生成锚点并添加到锚点数组 */
void add_interval_anchors(void *km, mg128_t **pa, int64_t *pn_a,
                          const mg_idx_t *gi,
                          const mg128_v *mv,
                          const seg_rep_t *reads_rep,
                          const void *strobe_idx,
                          const char *query_seq, int qlen,
                          int max_occ);   /* 新增参数 */

/* 原有接口 */
strobe_hash_t *ws_strobe_build(const void *gi, const void *collect_rep, const mg128_v *all_syncs);
int ws_strobe_get(const strobe_hash_t *h, uint64_t h1, uint64_t h2, uint64_t *val);
int ws_strobe_get_all(const strobe_hash_t *h, uint64_t h1, uint64_t h2, uint64_t **vals, int *n);
void ws_strobe_destroy(strobe_hash_t *h);
void remove_anchors_in_ref_rep(void *km, mg128_t **pa, int64_t *pn_a,
                               const mg_idx_t *gi, const collect_rep_t *collect_rep,
                               int64_t *removed_count);

strobe_hash_t *build_ref_interval_index(const void *gi, const void *collect_rep);

void print_seg_high_intervals(const collect_rep_t *collect_rep, const char *sname);
void print_seg_indexed_syncmers(const mg_idx_t *gi, const char *sname);

void second_sampling_and_anchoring(void *km, mg128_t **pa, int64_t *pn_a,
                                    const mg_idx_t *gi,
                                    const mg128_v *mv,
                                    const seg_rep_t *reads_rep,
                                    const void *strobe_idx,
                                    const char *query_seq, int qlen,
                                    int max_occ);
#ifdef __cplusplus
}
#endif

#endif