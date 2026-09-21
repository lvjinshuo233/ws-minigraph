#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdio.h>
#include "ws-strobe.h"
#include "mgpriv.h"
#include "collect-I.h"
#include "khashl.h"
#include <cassert>
extern volatile long long total_gap_on_read;
extern volatile long long n_high_itv_seen;
extern volatile long long n_high_itv_processed;
extern volatile long long total_high_itv_len;
extern volatile long long n_kmer_brigde, total_bk_freq;
extern volatile long long n_itv_no_low_region;
extern volatile long long n_itv_no_bridge_cand;
extern volatile long long n_bridge_t0;
extern volatile long long n_cand_empty;
extern volatile long long n_p_loop_total;
extern volatile long long n_low_region_no_syncmer;
extern volatile long long n_low_region_all_zero;
extern volatile long long total_n_tried;
extern volatile long long n_bridge_all_failed;
/* 缓存状态（每个 read 一份） */
typedef struct {
    int      active;
    int32_t  rid;           /* 缓存的参考 segment id */
    int32_t  ref_lo;        /* 0-based，缓存窗口起点 */
    int32_t  ref_hi;        /* 0-based，缓存窗口终点（exclusive） */
    int32_t  read_hi;       /* 0-based，缓存对应的 read 结束位置 */
} cache_state_t;

/* 在 ws-strobe.c 顶部加上新增统计的 extern */
extern volatile long long n_cache_started;
extern volatile long long n_cache_continued;
extern volatile long long n_cache_broken;
/* 比对阶段使用的辅助结构体，保持内部可见 */
typedef struct {
    uint64_t hash;
    int32_t start;
    int32_t freq;
    int32_t qspan;
} read_sync_info_t;
/* ---- 桥接 k-mer 候选（用于按频次迭代尝试） ---- */
typedef struct {
    uint64_t hash;
    int32_t  q_end;
    int32_t  q_span;
    int32_t  freq;
} bridge_cand_t;

/* 用于 print_seg_indexed_syncmers 的辅助结构体和比较器 */
typedef struct {
    int32_t start_1b;
    int32_t freq;
} sync_pos_t;

static int cmp_sync_pos(const void *a, const void *b) {
    const sync_pos_t *sa = (const sync_pos_t*)a;
    const sync_pos_t *sb = (const sync_pos_t*)b;
    return (sa->start_1b > sb->start_1b) - (sa->start_1b < sb->start_1b);
}
static int cmp_int32(const void *a, const void *b) {
    int32_t x = *(const int32_t*)a;
    int32_t y = *(const int32_t*)b;
    return (x > y) - (x < y);
}
static int cmp_bridge_cand(const void *a, const void *b) {
    const bridge_cand_t *ca = (const bridge_cand_t*)a;
    const bridge_cand_t *cb = (const bridge_cand_t*)b;
    return ca->freq - cb->freq;    /* 频次升序 */
}
static int cmp_read_sync(const void *a, const void *b) {
    const read_sync_info_t *sa = (const read_sync_info_t*)a;
    const read_sync_info_t *sb = (const read_sync_info_t*)b;
    return sa->start - sb->start;
}
/* 重新定义 mg_idx_bucket_t，与 index.c 中的定义保持一致 */
typedef struct mg_idx_bucket_s {
    mg128_v a;
    int32_t n;
    uint64_t *p;
    void *h;
} mg_idx_bucket_t;

/* 重新生成 idxhash_t 类型 */
#define idx_hash(a) ((a)>>1)
#define idx_eq(a, b) ((a)>>1 == (b)>>1)
KHASHL_MAP_INIT(KH_LOCAL, idxhash_t, mg_hidx, uint64_t, uint64_t, idx_hash, idx_eq)

/* strobemer 值：动态数组，保存所有匹配位置 */
typedef struct {
    uint64_t *vals;
    int n, m;
} strobe_val_t;

/* strobemer 键 */
typedef struct {
    uint64_t h1, h2;
} strobe_key_t;

static inline uint64_t strobe_hash(strobe_key_t key) {
    return key.h1 ^ (key.h2 + 0x9e3779b9 + (key.h1 << 6) + (key.h1 >> 2));
}
static inline int strobe_eq(strobe_key_t a, strobe_key_t b) {
    return a.h1 == b.h1 && a.h2 == b.h2;
}

KHASHL_MAP_INIT(KH_LOCAL, strobehash_t, strobe_h, strobe_key_t, strobe_val_t,
                strobe_hash, strobe_eq)

typedef struct {
    int32_t *start;   // 1-based 起始
    int32_t *end;     // 1-based 结束（exclusive）

    int n;
} rep_interval_arr_t;

struct strobe_hash_s {
    strobehash_t *table;
    rep_interval_arr_t *segs;   // 每个 segment 一个
    int n_seg;
};

/* 编码宏 */
/**
 * strobe_val 位置信息编码 
 * 低 32 位：参考 s2 的末位碱基位置（0-based）

中间 16 位：s1 与 s2 的距离 dist（s2_start - s1_start）

高 16 位：参考段 ID rid
 */
#define STROBE_VAL(end_pos, dist, rid) \
    (((uint64_t)(rid) << 48) | ((uint64_t)(dist) << 32) | (uint32_t)(end_pos))
#define STROBE_RID(val)   ((uint32_t)((val) >> 48))
#define STROBE_DIST(val)  ((uint16_t)((val) >> 32))
#define STROBE_END(val)   ((uint32_t)(val))

/* 提取的 syncmer 信息 */
typedef struct {
    uint64_t hash;
    uint32_t rid;
    int32_t start;
    int32_t freq;
} sync_info_t;

static sync_info_t *syncs = NULL;
static size_t n_syncs = 0, m_syncs = 0;

static void syncs_push(sync_info_t s) {
    if (n_syncs == m_syncs) {
        m_syncs = m_syncs ? m_syncs << 1 : 1024;
        syncs = (sync_info_t*)realloc(syncs, m_syncs * sizeof(sync_info_t));
    }
    syncs[n_syncs++] = s;
}

static int cmp_sync(const void *a, const void *b) {
    const sync_info_t *sa = (const sync_info_t*)a;
    const sync_info_t *sb = (const sync_info_t*)b;
    uint64_t ka = ((uint64_t)sa->rid << 32) | (uint32_t)sa->start;
    uint64_t kb = ((uint64_t)sb->rid << 32) | (uint32_t)sb->start;
    return (ka > kb) - (ka < kb);
}

static void extract_syncmers(const mg_idx_t *gi) {
    n_syncs = 0;
    int k = gi->k;
    for (uint32_t b = 0; b < (1U << gi->b); ++b) {
        mg_idx_bucket_t *bucket = &gi->B[b];
        idxhash_t *h = (idxhash_t*)bucket->h;
        if (!h) continue;
        for (khint_t it = 0; it < kh_end(h); ++it) {
            if (!kh_exist(h, it)) continue;
            uint64_t key = kh_key(h, it);
            uint64_t val = kh_val(h, it);
            uint64_t hash = ((key >> 1) << gi->b) | b;
            if (key & 1) {
                uint64_t y = val;
                uint32_t rid = (uint32_t)(y >> 32);
                uint32_t last_base = (uint32_t)((y >> 1) & 0x7fffffff);
                int32_t start = (int32_t)last_base - k + 1;
                syncs_push((sync_info_t){hash, rid, start, 1});
            } else {
                uint32_t count = (uint32_t)(val & 0xffffffff);
                uint32_t offset = (uint32_t)(val >> 32);
                const uint64_t *positions = &bucket->p[offset];
                for (uint32_t j = 0; j < count; ++j) {
                    uint64_t y = positions[j];
                    uint32_t rid = (uint32_t)(y >> 32);
                    uint32_t last_base = (uint32_t)((y >> 1) & 0x7fffffff);
                    int32_t start = (int32_t)last_base - k + 1;
                    syncs_push((sync_info_t){hash, rid, start, (int32_t)count});
                }
            }
        }
    }
    qsort(syncs, n_syncs, sizeof(sync_info_t), cmp_sync);
}

/* 插入 strobemer 值，支持多值 */
static void strobe_insert(strobehash_t *h, strobe_key_t key, uint64_t val) {
    int absent;
    khint_t it = strobe_h_put(h, key, &absent);
    if (absent) {
        kh_val(h, it).vals = (uint64_t*)malloc(4 * sizeof(uint64_t));
        kh_val(h, it).n = 0;
        kh_val(h, it).m = 4;
    }
    strobe_val_t *sv = &kh_val(h, it);
    if (sv->n == sv->m) {
        sv->m <<= 1;
        sv->vals = (uint64_t*)realloc(sv->vals, sv->m * sizeof(uint64_t));
    }
    sv->vals[sv->n++] = val;
}

strobe_hash_t *ws_strobe_build(const void *gi_ptr, const void *collect_rep_ptr, const mg128_v *all_syncs) {
    const mg_idx_t *gi = (const mg_idx_t*)gi_ptr;
    const collect_rep_t *collect_rep = (const collect_rep_t*)collect_rep_ptr;
    if (!gi || !collect_rep || !all_syncs || all_syncs->n == 0) return NULL;

    strobe_hash_t *sh = (strobe_hash_t*)calloc(1, sizeof(strobe_hash_t));
    if (!sh) return NULL;
    sh->table = strobe_h_init();

    int k = gi->k;
    size_t total_intervals = 0;
    const int max_dist_limit = (gi->max_rep_len > 0) ? (2 * gi->max_rep_len) : 0xFFFF;

    size_t cursor = 0;   // 全局游标，指向 all_syncs

    for (uint32_t rid = 0; rid < (uint32_t)collect_rep->len; ++rid) {
        const seg_rep_t *seg = &collect_rep->segs[rid];
        const rep_interval_t *itv = seg->high_first;
        const rep_interval_t *prev_itv = NULL;

        // 统计区间总数
        for (const rep_interval_t *tmp = itv; tmp; tmp = tmp->next)
            total_intervals++;

        // 跳过 rid 小于当前的 syncmer
        while (cursor < all_syncs->n && (uint32_t)(all_syncs->a[cursor].y >> 32) < rid)
            cursor++;

        size_t seg_start = cursor;   // 当前 segment 的 syncmer 起始索引

        while (itv) {
            // ---------- 定位整个高频区间内的 syncmer 范围 [whole_left, whole_right) ----------
            size_t whole_left = seg_start;
            while (whole_left < all_syncs->n &&
                   (uint32_t)(all_syncs->a[whole_left].y >> 32) == rid &&
                   ((int32_t)((all_syncs->a[whole_left].y >> 1) & 0x7fffffff) - k + 1) + 1 < itv->start)
                whole_left++;
            size_t whole_right = whole_left;
            while (whole_right < all_syncs->n &&
                   (uint32_t)(all_syncs->a[whole_right].y >> 32) == rid &&
                   ((int32_t)((all_syncs->a[whole_right].y >> 1) & 0x7fffffff) - k + 1) + 1 < itv->end)
                whole_right++;

            size_t n_total = whole_right - whole_left;
            // 修改：只要有一个候选 syncmer 就可以选择 s2（不再要求至少 3 个）
            if (n_total < 1) {
                prev_itv = itv;
                itv = itv->next;
                continue;
            }

            int32_t interval_len = itv->end - itv->start;

            if (gi->max_rep_len <= 0 || interval_len <= gi->max_rep_len) {
                // ================= 不分段：整个区间作为一个整体 =================
                // 提取候选数组
                read_sync_info_t *cand = (read_sync_info_t*)malloc(n_total * sizeof(read_sync_info_t));
                for (size_t i = 0; i < n_total; ++i) {
                    mg128_t *p = &all_syncs->a[whole_left + i];
                    cand[i].hash = p->x >> 8;
                    cand[i].start = (int32_t)((p->y >> 1) & 0x7fffffff) - k + 1;
                    int freq;
                    mg_idx_get(gi, cand[i].hash, &freq);
                    cand[i].freq = freq;
                    cand[i].qspan = p->x & 0xff;
                }

                // 修改：直接选择频次最低的 syncmer 作为 s2，不再限制中间三分之一
                const read_sync_info_t *best_s2 = &cand[0];
                for (int i = 1; i < n_total; ++i)
                    if (cand[i].freq < best_s2->freq) best_s2 = &cand[i];

                // 选择 s1：来自前一个低频区
                if (prev_itv) {
                    int region_lo = prev_itv->end - 1;   // 0-based inclusive
                    int region_hi = itv->start - 1;      // 0-based exclusive

                    size_t r_left = seg_start;
                    while (r_left < whole_left &&
                           (uint32_t)(all_syncs->a[r_left].y >> 32) == rid &&
                           ((int32_t)((all_syncs->a[r_left].y >> 1) & 0x7fffffff) - k + 1) < region_lo)
                        r_left++;
                    size_t r_right = r_left;
                    while (r_right < whole_left &&
                           (uint32_t)(all_syncs->a[r_right].y >> 32) == rid &&
                           ((int32_t)((all_syncs->a[r_right].y >> 1) & 0x7fffffff) - k + 1) < region_hi)
                        r_right++;

                    if (r_right > r_left) {
                        read_sync_info_t *s1_cand = (read_sync_info_t*)malloc((r_right - r_left) * sizeof(read_sync_info_t));
                        for (size_t i = r_left; i < r_right; ++i) {
                            mg128_t *p = &all_syncs->a[i];
                            s1_cand[i - r_left].hash = p->x >> 8;
                            s1_cand[i - r_left].start = (int32_t)((p->y >> 1) & 0x7fffffff) - k + 1;
                            int freq;
                            mg_idx_get(gi, s1_cand[i - r_left].hash, &freq);
                            s1_cand[i - r_left].freq = freq;
                            s1_cand[i - r_left].qspan = p->x & 0xff;
                        }
                        const read_sync_info_t *best_s1 = &s1_cand[0];
                        for (size_t i = 1; i < r_right - r_left; ++i)
                            if (s1_cand[i].freq < best_s1->freq) best_s1 = &s1_cand[i];

                        int dist = best_s2->start - best_s1->start;
                        if (dist >= 0 && dist <= max_dist_limit) {
                            uint32_t end_pos = best_s2->start + k - 1;
                            uint64_t encoded = STROBE_VAL(end_pos, dist, rid);
                            strobe_key_t key = {best_s1->hash, best_s2->hash};
                            strobe_insert(sh->table, key, encoded);
                        }
                        free(s1_cand);
                    }
                }
                free(cand);
            } else {
                // ================= 分段处理 =================
                int n_segments = (interval_len + gi->max_rep_len - 1) / gi->max_rep_len;
                // 用于保存上一个子段选出的 s2 的结构体副本
                read_sync_info_t prev_sel_s2_copy;
                int has_prev_s2 = 0;

                for (int seg_idx = 0; seg_idx < n_segments; ++seg_idx) {
                    int32_t sub_start = itv->start + seg_idx * gi->max_rep_len;
                    int32_t sub_end   = sub_start + gi->max_rep_len;
                    if (sub_end > itv->end) sub_end = itv->end;

                    // 定位子段内的 syncmer 范围 [sub_left, sub_right)
                    size_t sub_left = whole_left;
                    while (sub_left < whole_right &&
                           ((int32_t)((all_syncs->a[sub_left].y >> 1) & 0x7fffffff) - k + 1) + 1 < sub_start)
                        sub_left++;
                    size_t sub_right = sub_left;
                    while (sub_right < whole_right &&
                           ((int32_t)((all_syncs->a[sub_right].y >> 1) & 0x7fffffff) - k + 1) + 1 < sub_end)
                        sub_right++;

                    size_t n_sub = sub_right - sub_left;
                    // 修改：子段内不足 1 个 syncmer 才跳过（即只要有 1 个就选它）
                    if (n_sub < 1) continue;

                    // 提取子段候选信息
                    read_sync_info_t *sub_cand = (read_sync_info_t*)malloc(n_sub * sizeof(read_sync_info_t));
                    for (size_t i = 0; i < n_sub; ++i) {
                        mg128_t *p = &all_syncs->a[sub_left + i];
                        sub_cand[i].hash = p->x >> 8;
                        sub_cand[i].start = (int32_t)((p->y >> 1) & 0x7fffffff) - k + 1;
                        int freq;
                        mg_idx_get(gi, sub_cand[i].hash, &freq);
                        sub_cand[i].freq = freq;
                        sub_cand[i].qspan = p->x & 0xff;
                    }

                    // 修改：直接选择频次最低的 syncmer 作为 s2，不再限制中间三分之一
                    const read_sync_info_t *best_s2 = &sub_cand[0];
                    for (size_t i = 1; i < n_sub; ++i)
                        if (sub_cand[i].freq < best_s2->freq) best_s2 = &sub_cand[i];

                    // 确定 s1（可能来自前一个低频区或上一个子段的 s2）
                    read_sync_info_t best_s1_copy;
                    int has_s1 = 0;

                    if (seg_idx == 0) {
                        // 第一个子段：s1 来自前一个低频区（可能截断）
                        if (prev_itv) {
                            int region_lo = prev_itv->end - 1;
                            int region_hi = itv->start - 1;
                            if (region_hi - region_lo > gi->max_rep_len)
                                region_lo = region_hi - gi->max_rep_len;

                            size_t r_left = seg_start;
                            while (r_left < sub_left &&
                                   (uint32_t)(all_syncs->a[r_left].y >> 32) == rid &&
                                   ((int32_t)((all_syncs->a[r_left].y >> 1) & 0x7fffffff) - k + 1) < region_lo)
                                r_left++;
                            size_t r_right = r_left;
                            while (r_right < sub_left &&
                                   (uint32_t)(all_syncs->a[r_right].y >> 32) == rid &&
                                   ((int32_t)((all_syncs->a[r_right].y >> 1) & 0x7fffffff) - k + 1) < region_hi)
                                r_right++;

                            if (r_right > r_left) {
                                read_sync_info_t *s1_cand = (read_sync_info_t*)malloc((r_right - r_left) * sizeof(read_sync_info_t));
                                for (size_t i = r_left; i < r_right; ++i) {
                                    mg128_t *p = &all_syncs->a[i];
                                    s1_cand[i - r_left].hash = p->x >> 8;
                                    s1_cand[i - r_left].start = (int32_t)((p->y >> 1) & 0x7fffffff) - k + 1;
                                    int freq;
                                    mg_idx_get(gi, s1_cand[i - r_left].hash, &freq);
                                    s1_cand[i - r_left].freq = freq;
                                    s1_cand[i - r_left].qspan = p->x & 0xff;
                                }
                                const read_sync_info_t *best_s1 = &s1_cand[0];
                                for (size_t i = 1; i < r_right - r_left; ++i)
                                    if (s1_cand[i].freq < best_s1->freq) best_s1 = &s1_cand[i];
                                best_s1_copy = *best_s1;
                                has_s1 = 1;
                                free(s1_cand);
                            }
                        }
                    } else {
                        // 后续子段：s1 是上一个子段选出的 s2
                        if (has_prev_s2) {
                            best_s1_copy = prev_sel_s2_copy;
                            has_s1 = 1;
                        }
                    }

                    // 生成 strobemer
                    if (has_s1) {
                        int dist = best_s2->start - best_s1_copy.start;
                        if (dist >= 0 && dist <= max_dist_limit) {
                            uint32_t end_pos = best_s2->start + k - 1;
                            uint64_t encoded = STROBE_VAL(end_pos, dist, rid);
                            strobe_key_t key = {best_s1_copy.hash, best_s2->hash};
                            strobe_insert(sh->table, key, encoded);
                        }
                    }

                    // 保存当前子段选出的 s2 作为下一个子段的 s1（拷贝副本）
                    prev_sel_s2_copy = *best_s2;
                    has_prev_s2 = 1;

                    free(sub_cand);
                }
            }

            prev_itv = itv;
            itv = itv->next;
        }

        // 更新游标到当前 segment 结束
        cursor = seg_start;
        while (cursor < all_syncs->n && (uint32_t)(all_syncs->a[cursor].y >> 32) <= rid)
            cursor++;
    }

    // 调试输出（保留原逻辑）
    int dbg_cnt = 0;
    size_t total_strobes = 0;
    for (khint_t i = 0; i < kh_end(sh->table); ++i) {
        if (!kh_exist(sh->table, i)) continue;
        strobe_key_t key = kh_key(sh->table, i);
        strobe_val_t *sv = &kh_val(sh->table, i);
        total_strobes += sv->n;
        if (dbg_cnt < 100) {
            uint64_t val = sv->vals[0];
            uint32_t rid = STROBE_RID(val);
            uint16_t dist = STROBE_DIST(val);
            uint32_t end_pos = STROBE_END(val);
            int32_t s2_start = (int32_t)(end_pos - k + 1);
            int32_t s1_start = s2_start - dist;
            /*fprintf(stderr,
                    "S1_hash=%016llx  S2_hash=%016llx  rid=%u  dist=%d  "
                    "s1_start=%d  s2_start=%d  end=%u  total_vals=%d\n",
                    (unsigned long long)key.h1,
                    (unsigned long long)key.h2,
                    rid, dist, s1_start, s2_start, end_pos, sv->n);*/
            dbg_cnt++;
        }
    }
    fprintf(stderr, "Total strobemers: %zu, Total keys: %zu, Total intervals: %zu\n",
            total_strobes, (size_t)kh_size(sh->table), total_intervals);
    n_syncmers += (long int)total_strobes;
    n_strobes += (long int)total_strobes;
        // 预计算每个 segment 的高频区间数组（用于后续快速二分查找）
    sh->n_seg = (int)collect_rep->len;
    sh->segs = (rep_interval_arr_t*)calloc(sh->n_seg, sizeof(rep_interval_arr_t));
    for (uint32_t rid = 0; rid < (uint32_t)collect_rep->len; ++rid) {
        int n = 0;
        for (const rep_interval_t *itv = collect_rep->segs[rid].high_first; itv; itv = itv->next) n++;
        if (n == 0) continue;
        sh->segs[rid].start = (int32_t*)malloc(n * sizeof(int32_t));
        sh->segs[rid].end = (int32_t*)malloc(n * sizeof(int32_t));
        sh->segs[rid].n = n;
        int i = 0;
        for (const rep_interval_t *itv = collect_rep->segs[rid].high_first; itv; itv = itv->next) {
            sh->segs[rid].start[i] = itv->start;
            sh->segs[rid].end[i] = itv->end;
            i++;
        }
    }
    return sh;
}

int ws_strobe_get(const strobe_hash_t *h, uint64_t h1, uint64_t h2, uint64_t *val) {
    if (!h || !h->table) return 0;
    strobe_key_t key = {h1, h2};
    khint_t it = strobe_h_get(h->table, key);
    if (it == kh_end(h->table)) return 0;
    strobe_val_t *sv = &kh_val(h->table, it);
    if (sv->n == 0) return 0;
    *val = sv->vals[0];
    return 1;
}

int ws_strobe_get_all(const strobe_hash_t *h, uint64_t h1, uint64_t h2, uint64_t **vals, int *n) {
    if (!h || !h->table) return 0;
    strobe_key_t key = {h1, h2};
    khint_t it = strobe_h_get(h->table, key);
    if (it == kh_end(h->table)) return 0;
    strobe_val_t *sv = &kh_val(h->table, it);
    if (sv->n == 0) return 0;
    *vals = sv->vals;
    *n = sv->n;
    return 1;
}

void ws_strobe_destroy(strobe_hash_t *h) {
    if (!h) return;
    if (h->table) {
        for (khint_t i = 0; i < kh_end(h->table); ++i) {
            if (!kh_exist(h->table, i)) continue;
            strobe_val_t *sv = &kh_val(h->table, i);
            free(sv->vals);
        }
        strobe_h_destroy(h->table);
    }
    if (h->segs) {
    for (int i = 0; i < h->n_seg; ++i) {
        free(h->segs[i].start);
        free(h->segs[i].end);
    }
    free(h->segs);
}
    free(h);
}
static read_sync_info_t* extract_read_syncmers(const mg_idx_t *gi, const mg128_v *mv, int *n_out) {
    int n = mv->n;
    read_sync_info_t *arr = (read_sync_info_t*)malloc(n * sizeof(read_sync_info_t));
    int cnt = 0;
    for (int i = 0; i < n; ++i) {
        uint64_t hash = mv->a[i].x >> 8;  // 去掉低8位k-mer长度
        int t;
        mg_idx_get(gi, hash, &t);
        if (t <= 0) continue;
        int32_t end_pos = (int32_t)(mv->a[i].y) >> 1;
        int32_t start = end_pos - gi->k + 1;
        if (start < 0) continue;
        arr[cnt].hash = hash;
        arr[cnt].start = start;
        arr[cnt].freq = t;
        arr[cnt].qspan = mv->a[i].x & 0xff;
        cnt++;
    }
    qsort(arr, cnt, sizeof(read_sync_info_t), cmp_read_sync);
    *n_out = cnt;
    return arr;
}

int collect_ws_strobes(const mg_idx_t *gi, const mg128_v *mv, seg_rep_t *reads_rep,
                       read_strobe_t **out_strobes) {
    if (mv->n == 0) {
        *out_strobes = NULL;
        return 0;
    }

    int total_strobes = 0;
    int cap = 0;
    read_strobe_t *strobes = NULL;
    const int max_dist_limit = (gi->max_rep_len > 0) ? (2 * gi->max_rep_len) : 0xFFFF;

    rep_interval_t *prev_itv = NULL;
    size_t cursor = 0;   // 指向 mv 数组的游标

    for (rep_interval_t *itv = reads_rep->high_first; itv; itv = itv->next) {
        int32_t interval_len = itv->end - itv->start;

        // 定位整个高频区间内的 syncmer 范围 [whole_left, whole_right)
        size_t whole_left = cursor;
        while (whole_left < mv->n && 
               ((int32_t)((mv->a[whole_left].y >> 1) & 0x7fffffff) - gi->k + 1) + 1 < itv->start)
            whole_left++;
        size_t whole_right = whole_left;
        while (whole_right < mv->n &&
               ((int32_t)((mv->a[whole_right].y >> 1) & 0x7fffffff) - gi->k + 1) + 1 < itv->end)
            whole_right++;

        size_t n_total = whole_right - whole_left;
        if (n_total < 3) {
            prev_itv = itv;
            cursor = whole_right;   // 更新游标
            continue;
        }

        if (gi->max_rep_len <= 0 || interval_len <= gi->max_rep_len) {
            // ---------- 不分段 ----------
            // 提取候选数组（直接使用 mv 解码，避免额外分配大数组？但为了简化，保留局部数组）
            read_sync_info_t *cand = (read_sync_info_t*)malloc(n_total * sizeof(read_sync_info_t));
            for (size_t i = 0; i < n_total; ++i) {
                mg128_t *p = &mv->a[whole_left + i];
                cand[i].hash = p->x >> 8;
                cand[i].start = (int32_t)((p->y >> 1) & 0x7fffffff) - gi->k + 1;
                int freq;
                mg_idx_get(gi, cand[i].hash, &freq);
                cand[i].freq = freq;
                cand[i].qspan = p->x & 0xff;
            }

            int w_min = (n_total + 2) / 3;
            int w_max = (2 * n_total + 2) / 3;
            const read_sync_info_t *best_s2 = &cand[w_min - 1];
            for (int i = w_min - 1; i < w_max; ++i)
                if (cand[i].freq < best_s2->freq) best_s2 = &cand[i];

            if (prev_itv) {
                int region_lo = prev_itv->end - 1;
                int region_hi = itv->start - 1;

                size_t r_left = cursor;  // 注意：cursor 在上一个区间后已经更新
                while (r_left < whole_left &&
                       ((int32_t)((mv->a[r_left].y >> 1) & 0x7fffffff) - gi->k + 1) < region_lo)
                    r_left++;
                size_t r_right = r_left;
                while (r_right < whole_left &&
                       ((int32_t)((mv->a[r_right].y >> 1) & 0x7fffffff) - gi->k + 1) < region_hi)
                    r_right++;

                if (r_right > r_left) {
                    read_sync_info_t *s1_cand = (read_sync_info_t*)malloc((r_right - r_left) * sizeof(read_sync_info_t));
                    for (size_t i = r_left; i < r_right; ++i) {
                        mg128_t *p = &mv->a[i];
                        s1_cand[i - r_left].hash = p->x >> 8;
                        s1_cand[i - r_left].start = (int32_t)((p->y >> 1) & 0x7fffffff) - gi->k + 1;
                        int freq;
                        mg_idx_get(gi, s1_cand[i - r_left].hash, &freq);
                        s1_cand[i - r_left].freq = freq;
                        s1_cand[i - r_left].qspan = p->x & 0xff;
                    }
                    const read_sync_info_t *best_s1 = &s1_cand[0];
                    for (size_t i = 1; i < r_right - r_left; ++i)
                        if (s1_cand[i].freq < best_s1->freq) best_s1 = &s1_cand[i];

                    int dist = best_s2->start - best_s1->start;
                    if (dist >= 0 && dist <= max_dist_limit) {
                        if (total_strobes == cap) {
                            cap = cap ? cap * 2 : 16;
                            strobes = (read_strobe_t*)realloc(strobes, cap * sizeof(read_strobe_t));
                        }
                        strobes[total_strobes].h1 = best_s1->hash;
                        strobes[total_strobes].h2 = best_s2->hash;
                        strobes[total_strobes].s1_start = best_s1->start;
                        strobes[total_strobes].s2_start = best_s2->start;
                        strobes[total_strobes].dist = dist;
                        strobes[total_strobes].freq1 = best_s1->freq;
                        strobes[total_strobes].freq2 = best_s2->freq;
                        strobes[total_strobes].s1_qspan = best_s1->qspan;
                        strobes[total_strobes].s2_qspan = best_s2->qspan;
                        total_strobes++;
                    }
                    free(s1_cand);
                }
            }
            free(cand);
        } else {
            // ---------- 分段处理 ----------
            int n_segments = (interval_len + gi->max_rep_len - 1) / gi->max_rep_len;
            read_sync_info_t prev_sel_s2_copy;
            int has_prev_s2 = 0;

            for (int seg_idx = 0; seg_idx < n_segments; ++seg_idx) {
                int32_t sub_start = itv->start + seg_idx * gi->max_rep_len;
                int32_t sub_end   = sub_start + gi->max_rep_len;
                if (sub_end > itv->end) sub_end = itv->end;

                size_t sub_left = whole_left;
                while (sub_left < whole_right &&
                       ((int32_t)((mv->a[sub_left].y >> 1) & 0x7fffffff) - gi->k + 1) + 1 < sub_start)
                    sub_left++;
                size_t sub_right = sub_left;
                while (sub_right < whole_right &&
                       ((int32_t)((mv->a[sub_right].y >> 1) & 0x7fffffff) - gi->k + 1) + 1 < sub_end)
                    sub_right++;

                size_t n_sub = sub_right - sub_left;
                if (n_sub < 2) continue;

                read_sync_info_t *sub_cand = (read_sync_info_t*)malloc(n_sub * sizeof(read_sync_info_t));
                for (size_t i = 0; i < n_sub; ++i) {
                    mg128_t *p = &mv->a[sub_left + i];
                    sub_cand[i].hash = p->x >> 8;
                    sub_cand[i].start = (int32_t)((p->y >> 1) & 0x7fffffff) - gi->k + 1;
                    int freq;
                    mg_idx_get(gi, sub_cand[i].hash, &freq);
                    sub_cand[i].freq = freq;
                    sub_cand[i].qspan = p->x & 0xff;
                }

                const read_sync_info_t *best_s2 = &sub_cand[0];
                for (size_t i = 1; i < n_sub; ++i)
                    if (sub_cand[i].freq < best_s2->freq) best_s2 = &sub_cand[i];

                read_sync_info_t best_s1_copy;
                int has_s1 = 0;
                if (seg_idx == 0) {
                    if (prev_itv) {
                        int region_lo = prev_itv->end - 1;
                        int region_hi = itv->start - 1;
                        if (region_hi - region_lo > gi->max_rep_len)
                            region_lo = region_hi - gi->max_rep_len;

                        size_t r_left = cursor;
                        while (r_left < sub_left &&
                               ((int32_t)((mv->a[r_left].y >> 1) & 0x7fffffff) - gi->k + 1) < region_lo)
                            r_left++;
                        size_t r_right = r_left;
                        while (r_right < sub_left &&
                               ((int32_t)((mv->a[r_right].y >> 1) & 0x7fffffff) - gi->k + 1) < region_hi)
                            r_right++;

                        if (r_right > r_left) {
                            read_sync_info_t *s1_cand = (read_sync_info_t*)malloc((r_right - r_left) * sizeof(read_sync_info_t));
                            for (size_t i = r_left; i < r_right; ++i) {
                                mg128_t *p = &mv->a[i];
                                s1_cand[i - r_left].hash = p->x >> 8;
                                s1_cand[i - r_left].start = (int32_t)((p->y >> 1) & 0x7fffffff) - gi->k + 1;
                                int freq;
                                mg_idx_get(gi, s1_cand[i - r_left].hash, &freq);
                                s1_cand[i - r_left].freq = freq;
                                s1_cand[i - r_left].qspan = p->x & 0xff;
                            }
                            const read_sync_info_t *best_s1 = &s1_cand[0];
                            for (size_t i = 1; i < r_right - r_left; ++i)
                                if (s1_cand[i].freq < best_s1->freq) best_s1 = &s1_cand[i];
                            best_s1_copy = *best_s1;
                            has_s1 = 1;
                            free(s1_cand);
                        }
                    }
                } else {
                    if (has_prev_s2) {
                        best_s1_copy = prev_sel_s2_copy;
                        has_s1 = 1;
                    }
                }

                if (has_s1) {
                    int dist = best_s2->start - best_s1_copy.start;
                    if (dist >= 0 && dist <= max_dist_limit) {
                        if (total_strobes == cap) {
                            cap = cap ? cap * 2 : 16;
                            strobes = (read_strobe_t*)realloc(strobes, cap * sizeof(read_strobe_t));
                        }
                        strobes[total_strobes].h1 = best_s1_copy.hash;
                        strobes[total_strobes].h2 = best_s2->hash;
                        strobes[total_strobes].s1_start = best_s1_copy.start;
                        strobes[total_strobes].s2_start = best_s2->start;
                        strobes[total_strobes].dist = dist;
                        strobes[total_strobes].freq1 = best_s1_copy.freq;
                        strobes[total_strobes].freq2 = best_s2->freq;
                        strobes[total_strobes].s1_qspan = best_s1_copy.qspan;
                        strobes[total_strobes].s2_qspan = best_s2->qspan;
                        total_strobes++;
                    }
                }

                prev_sel_s2_copy = *best_s2;
                has_prev_s2 = 1;
                free(sub_cand);
            }
        }

        prev_itv = itv;
        cursor = whole_right;   // 更新游标
    }

    *out_strobes = strobes;
    return total_strobes;
}
void free_reads_rep(seg_rep_t *reads_rep, int n_segs) {
    if (!reads_rep) return;
    for (int i = 0; i < n_segs; ++i) {
        seg_rep_t *seg = &reads_rep[i];
        rep_interval_t *itv = seg->high_reps_list->next;
        while (itv) {
            rep_interval_t *next = itv->next;
            free(itv);
            itv = next;
        }
        free(seg->high_reps_list);
        if (seg->low_reps_list) {
            itv = seg->low_reps_list->next;
            while (itv) {
                rep_interval_t *next = itv->next;
                free(itv);
                itv = next;
            }
            free(seg->low_reps_list);
        }
    }
    free(reads_rep);
}

void add_strobe_anchors(void *km, mg128_t **pa, int64_t *pn_a,
                        const read_strobe_t *strobes, int n_strobes,
                        const mg_idx_t *gi,
                        const char *query_seq, int qlen) {
    if (n_strobes == 0 || !gi->strobe_idx) return;

    int64_t old_n = *pn_a;
    int64_t new_anchors = 0;
    int *match_counts = (int*)malloc(n_strobes * sizeof(int));

    for (int i = 0; i < n_strobes; ++i) {
        uint64_t *vals = NULL;
        int n_vals = 0;
        if (ws_strobe_get_all((strobe_hash_t*)gi->strobe_idx, strobes[i].h1, strobes[i].h2, &vals, &n_vals)) {
            match_counts[i] = n_vals;
            new_anchors += n_vals;
        } else {
            match_counts[i] = 0;
        }
    }
    if (new_anchors == 0) {
        free(match_counts);
        return;
    }

    int64_t new_n = old_n + new_anchors;
    mg128_t *a = (mg128_t*)krealloc(km, *pa, new_n * sizeof(mg128_t));
    if (!a) {
        fprintf(stderr, "[E::%s] failed to realloc for strobe anchors\n", __func__);
        abort();
    }
    *pa = a;

    int k = gi->k;
    int64_t idx = old_n;
    for (int i = 0; i < n_strobes; ++i) {
        if (match_counts[i] == 0) continue;
        uint64_t *vals = NULL;
        int n_vals = 0;
        ws_strobe_get_all((strobe_hash_t*)gi->strobe_idx, strobes[i].h1, strobes[i].h2, &vals, &n_vals);
        for (int j = 0; j < n_vals; ++j) {
            uint64_t val = vals[j];
            uint32_t rid = (uint32_t)(val >> 48);
            uint32_t ref_end = (uint32_t)(val & 0xFFFFFFFF);
            int32_t q_end = strobes[i].s2_start + k - 1;
            int32_t q_span_orig = strobes[i].s2_qspan;

            a[idx].x = (uint64_t)rid << 33 | ref_end;
            a[idx].y = (uint64_t)q_span_orig << 32 | q_end;
            a[idx].y |= (uint64_t)1 << MG_SEED_OCC_SHIFT;

            // 向前扩展
            if (gi->g && rid < (uint32_t)gi->g->n_seg && query_seq && qlen > 0) {
                const char *ref_seq = gi->g->seg[rid].seq;
                int q_start = q_end - q_span_orig + 1;
                int ref_start = ref_end - q_span_orig + 1;
                int ext = 0;
                int max_ext = 255 - q_span_orig;
                while (ext < max_ext) {
                    int q_pos = q_start - ext - 1;
                    int r_pos = ref_start - ext - 1;
                    if (q_pos < 0 || r_pos < 0) break;
                    if (query_seq[q_pos] != ref_seq[r_pos]) break;
                    ext++;
                }
                if (ext > 0) {
                    q_span_orig += ext;
                    a[idx].y = (uint64_t)q_span_orig << 32 | q_end;
                    a[idx].y |= (uint64_t)1 << MG_SEED_OCC_SHIFT;
                }
            }
            idx++;
        }
    }
    *pn_a = new_n;
    free(match_counts);
}

void remove_anchors_in_ref_rep(void *km, mg128_t **pa, int64_t *pn_a,
                               const mg_idx_t *gi, const collect_rep_t *collect_rep,
                               int64_t *removed_count) {
    if (*pn_a == 0 || !gi || !collect_rep) {
        if (removed_count) *removed_count = 0;
        return;
    }

    mg128_t *a = *pa;
    int64_t old_n = *pn_a;
    int64_t write = 0;
    int64_t n_removed = 0;
    int k = gi->k;

    for (int64_t i = 0; i < old_n; ++i) {
        uint32_t rid = (uint32_t)(a[i].x >> 33);
        uint32_t r_end = (uint32_t)(a[i].x & 0x7fffffff);
        int32_t r_start = (int32_t)r_end - k + 1;
        int32_t r_start_1b = r_start + 1;

        int del = 0;
        if (rid < collect_rep->len) {
            const seg_rep_t *seg = &collect_rep->segs[rid];
            const rep_interval_t *itv = seg->high_first;
            while (itv) {
                if (r_start_1b >= itv->start && r_start_1b < itv->end) {
                    del = 1;
                    break;
                }
                itv = itv->next;
            }
        }
        if (del) {
            n_removed++;
        } else {
            if (write != i) a[write] = a[i];
            write++;
        }
    }
    *pn_a = write;
    if (removed_count) *removed_count = n_removed;
}

/* ===== 二分查找辅助 ===== */
static int find_first_interval_start_ge(const rep_interval_arr_t *arr, int32_t target) {
    int lo = 0, hi = arr->n;
    while (lo < hi) {
        int mid = (lo + hi) / 2;
        if (arr->start[mid] < target) lo = mid + 1;
        else hi = mid;
    }
    return lo;
}

static int find_last_interval_end_le(const rep_interval_arr_t *arr, int32_t target) {
    int lo = 0, hi = arr->n;
    while (lo < hi) {
        int mid = (lo + hi) / 2;
        if (arr->end[mid] <= target) lo = mid + 1;
        else hi = mid;
    }
    return lo - 1;  // -1 表示没有
}

/* 在 mv 中二分查找第一个 0-based 起始位置 >= target 的 syncmer 索引 */
static size_t find_first_syncmer_ge(const mg128_v *mv, int k, int32_t target) {
    size_t lo = 0, hi = mv->n;
    while (lo < hi) {
        size_t mid = (lo + hi) / 2;
        int32_t end_0b = (int32_t)((mv->a[mid].y >> 1) & 0x7fffffff);
        int32_t start_0b = end_0b - k + 1;
        if (start_0b < target) lo = mid + 1;
        else hi = mid;
    }
    return lo;
}

/* 在 all_syncs 中二分查找第一个 (sid >= target_sid 且该 sid 内 start_0b >= target_start_0b) 的条目 */
static size_t find_first_syncmer_in_ref(const mg128_v *all_syncs, int k,
                                        int32_t target_sid, int32_t target_start_0b)
{
    size_t lo = 0, hi = all_syncs->n;
    while (lo < hi) {
        size_t mid = (lo + hi) / 2;
        uint32_t sid = (uint32_t)(all_syncs->a[mid].y >> 32);
        int32_t end_0b = (int32_t)((all_syncs->a[mid].y >> 1) & 0x7fffffff);
        int32_t start_0b = end_0b - k + 1;
        if (sid < (uint32_t)target_sid ||
            (sid == (uint32_t)target_sid && start_0b < target_start_0b)) {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    return lo;
}

typedef struct {
    uint64_t hash;
    int32_t  end_0b;
} rf_kmer_t;

static int cmp_rf_kmer(const void *a, const void *b)
{
    const rf_kmer_t *pa = (const rf_kmer_t*)a;
    const rf_kmer_t *pb = (const rf_kmer_t*)b;
    if (pa->hash < pb->hash) return -1;
    if (pa->hash > pb->hash) return  1;
    return 0;
}

void add_interval_anchors(void *km, mg128_t **pa, int64_t *pn_a,
                          const mg_idx_t *gi,
                          const mg128_v *mv,
                          const seg_rep_t *reads_rep,
                          const void *strobe_idx,
                          const char *query_seq, int qlen,
                          int max_occ)
{
    if (!mv || mv->n == 0 || !reads_rep || !gi || !strobe_idx) return;
    if (gi->max_rep_len <= 0) return;

    const strobe_hash_t *sh = (const strobe_hash_t*)strobe_idx;
    if (!sh->segs || sh->n_seg == 0) return;

    const int k = gi->k;
    const int max_dist = 2 * gi->max_rep_len;
    const int min_low_len = gi->min_rep_len / 2;
    const int freq_threshold = 2 * max_occ;   /* 频次阈值 */

    int64_t old_n = *pn_a;
    mg128_t *new_a = NULL;
    int64_t added = 0;
    int64_t cap_added = 0;

    const rep_interval_t *prev_itv = NULL;
    for (const rep_interval_t *itv = reads_rep->high_first; itv; itv = itv->next) {
        int32_t rs = itv->start;   /* 1-based */
        int32_t re = itv->end;
        if (re - rs < 10) {
            prev_itv = itv;
            continue;
        }

        /* ===== 两侧尝试：前侧优先；若选出的 k-mer 频次 > 2*max_occ，换后侧 ===== */
        int found = 0;
        int use_pre_final = 0;
        uint64_t best_hash_final = 0;
        int32_t bridge_q_end_final  = 0;
        int32_t bridge_q_span_final = 0;
        int best_freq_final = 0;

        for (int attempt = 0; attempt < 2 && !found; ++attempt) {
            int use_pre = 0;
            int32_t reg_lo_1b = 0, reg_hi_1b = 0;

            if (attempt == 0) {
                /* ---- 前侧低频区 ---- */
                if (!prev_itv) continue;
                int32_t pre_lo = prev_itv->end;
                int32_t pre_hi = rs;
                if (pre_hi - pre_lo < min_low_len) continue;
                if (pre_hi - pre_lo > gi->max_rep_len) pre_lo = pre_hi - gi->max_rep_len;
                reg_lo_1b = pre_lo;
                reg_hi_1b = pre_hi;
                use_pre = 1;
            } else {
                /* ---- 后侧低频区 ---- */
                int32_t post_lo = re;
                int32_t post_hi = itv->next ? itv->next->start : (qlen + 1);
                if (post_hi - post_lo < min_low_len) continue;
                if (post_hi - post_lo > gi->max_rep_len) post_hi = post_lo + gi->max_rep_len;
                reg_lo_1b = post_lo;
                reg_hi_1b = post_hi;
                use_pre = 0;
            }

            int32_t reg_lo_0b = reg_lo_1b - 1;
            int32_t reg_hi_0b = reg_hi_1b - 1;

            /* ---- 收集该侧低频区内所有候选 syncmer ---- */
            size_t lo = find_first_syncmer_ge(mv, k, reg_lo_0b);

            bridge_cand_t *cands = NULL;
            int n_cand = 0, cap_cand = 0;

            for (size_t i = lo; i < mv->n; ++i) {
                int32_t end_0b   = (int32_t)((mv->a[i].y >> 1) & 0x7fffffff);
                int32_t start_0b = end_0b - k + 1;
                if (start_0b >= reg_hi_0b) break;

                uint64_t h = mv->a[i].x >> 8;
                int freq;
                mg_idx_get(gi, h, &freq);

                if (n_cand == cap_cand) {
                    cap_cand = cap_cand ? cap_cand * 2 : 8;
                    cands = (bridge_cand_t*)realloc(cands, cap_cand * sizeof(bridge_cand_t));
                }
                cands[n_cand].hash   = h;
                cands[n_cand].q_end  = end_0b;
                cands[n_cand].q_span = mv->a[i].x & 0xff;
                cands[n_cand].freq   = freq;
                n_cand++;
            }

            if (n_cand == 0) {
                free(cands);
                continue;   /* 该侧没有候选，试另一侧 */
            }

            qsort(cands, n_cand, sizeof(bridge_cand_t), cmp_bridge_cand);

            /* ---- 依次尝试候选，找第一个能在参考候选区内连上的 ---- */
            int chosen = -1;
            int n_tried = 0;
            int32_t chosen_gap = 0;

            for (int ci = 0; ci < n_cand; ++ci) {
                if (cands[ci].freq <= 0) continue;
                n_tried++;

                const uint64_t *cr;
                int t;
                cr = mg_idx_get(gi, cands[ci].hash, &t);
                if (t <= 0 || t > freq_threshold) continue;
                int32_t cand_q_end   = cands[ci].q_end;
                int32_t cand_q_span  = cands[ci].q_span;
                int32_t cand_q_start = cand_q_end - cand_q_span + 1;
                int32_t gap_on_read;
                if (use_pre) gap_on_read = rs - cand_q_end - 2;
                else         gap_on_read = cand_q_start - re + 1;
                if (gap_on_read < 0) gap_on_read = 0;

                int usable = 0;
                for (int p = 0; p < t && !usable; ++p) {
                    uint32_t rid_ref    = (uint32_t)(cr[p] >> 32);
                    uint32_t pos_ref_0b = (uint32_t)(cr[p] >> 1) & 0x7fffffff;
                    int32_t  pos_ref_1b = (int32_t)pos_ref_0b + 1;
                    if (rid_ref >= (uint32_t)sh->n_seg) continue;
                    const rep_interval_arr_t *arr = &sh->segs[rid_ref];
                    if (arr->n == 0) continue;

                    int cand_lo, cand_hi;
                    if (use_pre) {
                        cand_lo = find_first_interval_start_ge(arr, pos_ref_1b);
                        cand_hi = arr->n;
                        while (cand_hi > cand_lo &&
                               arr->start[cand_hi - 1] - pos_ref_1b > max_dist) cand_hi--;
                    } else {
                        int last = find_last_interval_end_le(arr, pos_ref_1b);
                        cand_lo = 0;
                        cand_hi = last + 1;
                        while (cand_lo < cand_hi &&
                               pos_ref_1b - arr->end[cand_lo] > max_dist) cand_lo++;
                    }
                    if (cand_hi > cand_lo) {
                        usable = 1;
                        chosen_gap = gap_on_read;
                    }
                }
                if (usable) { chosen = ci; break; }
            }

            if (chosen < 0) {
                free(cands);
                continue;   /* 该侧没有能连上的候选，试另一侧 */
            }

            /* ---- 频次检查：选出的 k-mer 若 > 2*max_occ，换另一侧 ---- */
            if (cands[chosen].freq > freq_threshold) {
                free(cands);
                continue;   /* 换另一侧重新尝试 */
            }

            /* 该侧成功 */
            found               = 1;
            use_pre_final       = use_pre;
            best_hash_final     = cands[chosen].hash;
            bridge_q_end_final  = cands[chosen].q_end;
            bridge_q_span_final = cands[chosen].q_span;
            best_freq_final     = cands[chosen].freq;
            free(cands);
        }

        if (!found) {
            /* 两侧都没有合格候选（无可连上的、或连上的频次都超阈值）：舍弃该区间 */
            prev_itv = itv;
            continue;
        }

        /* ---- 查参考位置 ---- */
        int use_pre = use_pre_final;
        uint64_t best_hash      = best_hash_final;
        int32_t bridge_q_end    = bridge_q_end_final;
        int32_t bridge_q_span   = bridge_q_span_final;
        int32_t bridge_q_start  = bridge_q_end - bridge_q_span + 1;
        int best_freq           = best_freq_final;
        (void)best_freq;  /* 若不需要统计，可去掉本行 */

        const uint64_t *cr;
        int t;
        cr = mg_idx_get(gi, best_hash, &t);
        if (t <= 0) {
            prev_itv = itv;
            continue;
        }

        /* ---- 预计算 I_r 内 syncmer 范围 ---- */
        size_t ir_lo = find_first_syncmer_ge(mv, k, rs - 1);
        size_t ir_hi = find_first_syncmer_ge(mv, k, re - 1);

        int bridge_dir = use_pre ? +1 : -1;
        int inner_dir  = -bridge_dir;

        /* ---- 对每个参考出现位置生成锚点 ---- */
        for (int p = 0; p < t; ++p) {
            uint32_t rid_ref    = (uint32_t)(cr[p] >> 32);
            uint32_t pos_ref_0b = (uint32_t)(cr[p] >> 1) & 0x7fffffff;
            int32_t  pos_ref_1b = (int32_t)pos_ref_0b + 1;

            if (rid_ref >= (uint32_t)sh->n_seg) continue;
            const rep_interval_arr_t *arr = &sh->segs[rid_ref];
            if (arr->n == 0) continue;
            const gfa_seg_t *ref_seg = &gi->g->seg[rid_ref];
            const char *ref_seq = ref_seg->seq;
            int ref_len = ref_seg->len;

            int cand_lo, cand_hi;
            if (use_pre) {
                cand_lo = find_first_interval_start_ge(arr, pos_ref_1b);
                cand_hi = arr->n;
                while (cand_hi > cand_lo &&
                       arr->start[cand_hi - 1] - pos_ref_1b > max_dist) cand_hi--;
            } else {
                int last = find_last_interval_end_le(arr, pos_ref_1b);
                cand_lo = 0;
                cand_hi = last + 1;
                while (cand_lo < cand_hi &&
                       pos_ref_1b - arr->end[cand_lo] > max_dist) cand_lo++;
            }
            if (cand_hi <= cand_lo) continue;

            /* ---- 桥接 k-mer 锚点（带扩展） ---- */
            {
                int32_t q_span_final = bridge_q_span;
                int32_t q_end_final  = bridge_q_end;
                int32_t r_end_final  = (int32_t)pos_ref_0b;

                if (bridge_dir > 0) {
                    int q_pos = bridge_q_end + 1;
                    int r_pos = (int)pos_ref_0b + 1;
                    int ext = 0, max_ext = 255 - q_span_final;
                    while (ext < max_ext && q_pos < qlen && r_pos < ref_len) {
                        if (query_seq[q_pos] != ref_seq[r_pos]) break;
                        q_pos++; r_pos++; ext++;
                    }
                    if (ext > 0) {
                        q_span_final += ext;
                        q_end_final   = bridge_q_end + ext;
                        r_end_final   = (int32_t)pos_ref_0b + ext;
                    }
                } else {
                    int q_pos = bridge_q_start - 1;
                    int r_pos = (int)pos_ref_0b - bridge_q_span;
                    int ext = 0, max_ext = 255 - q_span_final;
                    while (ext < max_ext && q_pos >= 0 && r_pos >= 0) {
                        if (query_seq[q_pos] != ref_seq[r_pos]) break;
                        q_pos--; r_pos--; ext++;
                    }
                    if (ext > 0) q_span_final += ext;
                }

                if (added == cap_added) {
                    cap_added = cap_added ? cap_added * 2 : 16;
                    new_a = (mg128_t*)realloc(new_a, cap_added * sizeof(mg128_t));
                }
                new_a[added].x  = (uint64_t)rid_ref << 33 | (uint32_t)r_end_final;
                new_a[added].y  = (uint64_t)q_span_final << 32 | (uint32_t)q_end_final;
                new_a[added].y |= (uint64_t)1 << MG_SEED_OCC_SHIFT;
                added++;
            }

            /* ---- I_r 内共享 k-mer 锚点（带扩展） ---- */
            for (int ci2 = cand_lo; ci2 < cand_hi; ++ci2) {
                int32_t is    = arr->start[ci2];
                int32_t ie    = arr->end[ci2];
                int32_t is_0b = is - 1;
                int32_t ie_0b = ie - 1;

                for (size_t j = ir_lo; j < ir_hi; ++j) {
                    uint64_t h      = mv->a[j].x >> 8;
                    int32_t q_span  = mv->a[j].x & 0xff;
                    int32_t q_end   = (int32_t)((mv->a[j].y >> 1) & 0x7fffffff);
                    int32_t q_start = q_end - q_span + 1;

                    const uint64_t *cr2;
                    int t2;
                    cr2 = mg_idx_get(gi, h, &t2);
                    if (t2 <= 0) continue;

                    for (int p2 = 0; p2 < t2; ++p2) {
                        uint32_t rid_r2    = (uint32_t)(cr2[p2] >> 32);
                        uint32_t pos_r2_0b = (uint32_t)(cr2[p2] >> 1) & 0x7fffffff;
                        if (rid_r2 != rid_ref) continue;
                        int32_t r2_start = (int32_t)pos_r2_0b - k + 1;
                        if (r2_start < is_0b || r2_start >= ie_0b) continue;

                        int32_t q_span_final = q_span;
                        int32_t q_end_final  = q_end;
                        int32_t r_end_final  = (int32_t)pos_r2_0b;

                        if (inner_dir > 0) {
                            int q_pos = q_end + 1;
                            int r_pos = (int)pos_r2_0b + 1;
                            int ext = 0, max_ext = 255 - q_span_final;
                            while (ext < max_ext && q_pos < qlen && r_pos < ref_len) {
                                if (query_seq[q_pos] != ref_seq[r_pos]) break;
                                q_pos++; r_pos++; ext++;
                            }
                            if (ext > 0) {
                                q_span_final += ext;
                                q_end_final   = q_end + ext;
                                r_end_final   = (int32_t)pos_r2_0b + ext;
                            }
                        } else {
                            int q_pos = q_start - 1;
                            int r_pos = (int)pos_r2_0b - q_span;
                            int ext = 0, max_ext = 255 - q_span_final;
                            while (ext < max_ext && q_pos >= 0 && r_pos >= 0) {
                                if (query_seq[q_pos] != ref_seq[r_pos]) break;
                                q_pos--; r_pos--; ext++;
                            }
                            if (ext > 0) q_span_final += ext;
                        }

                        if (added == cap_added) {
                            cap_added = cap_added ? cap_added * 2 : 16;
                            new_a = (mg128_t*)realloc(new_a, cap_added * sizeof(mg128_t));
                        }
                        new_a[added].x  = (uint64_t)rid_ref << 33 | (uint32_t)r_end_final;
                        new_a[added].y  = (uint64_t)q_span_final << 32 | (uint32_t)q_end_final;
                        new_a[added].y |= (uint64_t)1 << MG_SEED_OCC_SHIFT;
                        added++;
                    }
                }
            }
        }

        prev_itv = itv;
    }

    if (added > 0) {
        int64_t new_n = old_n + added;
        mg128_t *a = (mg128_t*)krealloc(km, *pa, new_n * sizeof(mg128_t));
        if (!a) {
            fprintf(stderr, "[E::%s] 锚点 realloc 失败\n", __func__);
            free(new_a);
            abort();
        }
        *pa = a;
        memcpy(&a[old_n], new_a, added * sizeof(mg128_t));
        *pn_a = new_n;
    }

    free(new_a);
}
strobe_hash_t *build_ref_interval_index(const void *gi_ptr, const void *collect_rep_ptr) {
    const collect_rep_t *collect_rep = (const collect_rep_t*)collect_rep_ptr;
    if (!collect_rep) return NULL;

    strobe_hash_t *sh = (strobe_hash_t*)calloc(1, sizeof(strobe_hash_t));
    if (!sh) return NULL;

    sh->n_seg = (int)collect_rep->len;
    sh->segs = (rep_interval_arr_t*)calloc(sh->n_seg, sizeof(rep_interval_arr_t));

    for (uint32_t rid = 0; rid < (uint32_t)collect_rep->len; ++rid) {
        int n = 0;
        for (const rep_interval_t *itv = collect_rep->segs[rid].high_first; itv; itv = itv->next) n++;
        if (n == 0) continue;
        sh->segs[rid].start = (int32_t*)malloc(n * sizeof(int32_t));
        sh->segs[rid].end   = (int32_t*)malloc(n * sizeof(int32_t));
        sh->segs[rid].n = n;
        int i = 0;
        for (const rep_interval_t *itv = collect_rep->segs[rid].high_first; itv; itv = itv->next) {
            sh->segs[rid].start[i] = itv->start;
            sh->segs[rid].end[i]   = itv->end;
            i++;
        }
    }
    return sh;
}

void print_seg_high_intervals(const collect_rep_t *collect_rep, const char *sname)
{
    if (!collect_rep || !sname) return;

    for (uint32_t i = 0; i < collect_rep->len; ++i) {
        const seg_rep_t *seg = &collect_rep->segs[i];
        if (!seg->sname || strcmp(seg->sname, sname) != 0) continue;

        fprintf(stderr,"=== %s: high-frequency intervals ===\n", seg->sname);
        int idx = 0;
        for (const rep_interval_t *itv = seg->high_first; itv; itv = itv->next) {
            fprintf(stderr, "[%d, %d)\tlen=%d\n",
                   itv->start, itv->end, itv->end - itv->start);
            idx++;
        }
        fprintf(stderr, "total: %d intervals\n", idx);
        return;
    }
    fprintf(stderr, "segment %s not found\n", sname);
}

void print_seg_indexed_syncmers(const mg_idx_t *gi, const char *sname)
{
    if (!gi || !gi->g || !sname) return;
    const gfa_t *g = gi->g;

    /* 1. 找 sname 对应的 rid */
    int32_t rid = -1;
    for (uint32_t i = 0; i < g->n_seg; ++i) {
        if (g->seg[i].name && strcmp(g->seg[i].name, sname) == 0) {
            rid = (int32_t)i;
            break;
        }
    }
    if (rid < 0) {
        fprintf(stderr, "segment %s not found in graph\n", sname);
        return;
    }

    const int k = gi->k;
    const char *seq = g->seg[rid].seq;
    const int32_t seq_len = g->seg[rid].len;

    /* 2. 遍历所有 bucket，收集属于 rid 的 syncmer 起始位置 + 频次 */
    int64_t cap = 1024, n = 0;
    sync_pos_t *positions = (sync_pos_t*)malloc(cap * sizeof(sync_pos_t));
    if (!positions) return;

    for (uint32_t b = 0; b < (1U << gi->b); ++b) {
        mg_idx_bucket_t *bucket = &gi->B[b];
        idxhash_t *h = (idxhash_t*)bucket->h;
        if (!h) continue;

        for (khint_t it = 0; it < kh_end(h); ++it) {
            if (!kh_exist(h, it)) continue;
            uint64_t key = kh_key(h, it);
            uint64_t val = kh_val(h, it);

            if (key & 1) {
                uint64_t y = val;
                if ((uint32_t)(y >> 32) != (uint32_t)rid) continue;
                uint32_t end_0b = (uint32_t)((y >> 1) & 0x7fffffff);
                int32_t start_1b = (int32_t)end_0b - k + 2;
                if (n == cap) {
                    cap <<= 1;
                    positions = (sync_pos_t*)realloc(positions, cap * sizeof(sync_pos_t));
                }
                positions[n].start_1b = start_1b;
                positions[n].freq     = 1;
                n++;
            } else {
                uint32_t count  = (uint32_t)(val & 0xffffffff);
                uint32_t offset = (uint32_t)(val >> 32);
                const uint64_t *p = &bucket->p[offset];
                for (uint32_t j = 0; j < count; ++j) {
                    uint64_t y = p[j];
                    if ((uint32_t)(y >> 32) != (uint32_t)rid) continue;
                    uint32_t end_0b = (uint32_t)((y >> 1) & 0x7fffffff);
                    int32_t start_1b = (int32_t)end_0b - k + 2;
                    if (n == cap) {
                        cap <<= 1;
                        positions = (sync_pos_t*)realloc(positions, cap * sizeof(sync_pos_t));
                    }
                    positions[n].start_1b = start_1b;
                    positions[n].freq     = (int32_t)count;
                    n++;
                }
            }
        }
    }

    /* 3. 排序 */
    qsort(positions, n, sizeof(sync_pos_t), cmp_sync_pos);

    /* 4. 输出：位置、频次、kmer 序列 */
    fprintf(stderr, "=== %s: indexed weighted syncmer start positions (1-based) ===\n", sname);
    for (int64_t i = 0; i < n; ++i) {
        int32_t s1b = positions[i].start_1b;
        int32_t start_0b = s1b - 1;
        if (start_0b >= 0 && start_0b + k <= seq_len) {
            fprintf(stderr, "%d\t%d\t%.*s\n",
                    s1b, positions[i].freq, k, seq + start_0b);
        } else {
            fprintf(stderr, "%d\t%d\t<OUT_OF_RANGE>\n",
                    s1b, positions[i].freq);
        }
    }
    fprintf(stderr, "total: %lld\n", (long long)n);

    free(positions);
}

void second_sampling_and_anchoring(void *km, mg128_t **pa, int64_t *pn_a,
                                    const mg_idx_t *gi,
                                    const mg128_v *mv,
                                    const seg_rep_t *reads_rep,
                                    const void *strobe_idx,
                                    const char *query_seq, int qlen,
                                    int max_occ)
{
    if (!mv || mv->n == 0 || !reads_rep || !gi || !strobe_idx) return;
    if (gi->max_rep_len <= 0) return;
    if (!gi->all_syncs.a || gi->all_syncs.n == 0) return;

    const strobe_hash_t *sh = (const strobe_hash_t*)strobe_idx;
    if (!sh->segs || sh->n_seg == 0) return;

    const int k = gi->k;
    const int max_dist    = 2 * gi->max_rep_len;
    const int min_low_len = gi->min_rep_len / 2;
    const int freq_threshold = 2 * max_occ;

    int64_t old_n = *pn_a;
    mg128_t *new_a = NULL;
    int64_t added = 0;
    int64_t cap_added = 0;

    const rep_interval_t *prev_itv = NULL;
    for (const rep_interval_t *itv = reads_rep->high_first; itv; itv = itv->next) {
        int32_t rs = itv->start;
        int32_t re = itv->end;
        if (re - rs < 10) {
            prev_itv = itv;
            continue;
        }

        /* ===== 选桥接 k-mer ===== */
        int found = 0;
        int use_pre_final = 0;
        uint64_t best_hash_final = 0;
        int32_t bridge_q_end_final  = 0;
        int32_t bridge_q_span_final = 0;

        for (int attempt = 0; attempt < 2 && !found; ++attempt) {
            int use_pre = 0;
            int32_t reg_lo_1b = 0, reg_hi_1b = 0;

            if (attempt == 0) {
                if (!prev_itv) continue;
                int32_t pre_lo = prev_itv->end;
                int32_t pre_hi = rs;
                if (pre_hi - pre_lo < min_low_len) continue;
                if (pre_hi - pre_lo > gi->max_rep_len) pre_lo = pre_hi - gi->max_rep_len;
                reg_lo_1b = pre_lo;
                reg_hi_1b = pre_hi;
                use_pre = 1;
            } else {
                int32_t post_lo = re;
                int32_t post_hi = itv->next ? itv->next->start : (qlen + 1);
                if (post_hi - post_lo < min_low_len) continue;
                if (post_hi - post_lo > gi->max_rep_len) post_hi = post_lo + gi->max_rep_len;
                reg_lo_1b = post_lo;
                reg_hi_1b = post_hi;
                use_pre = 0;
            }

            int32_t reg_lo_0b = reg_lo_1b - 1;
            int32_t reg_hi_0b = reg_hi_1b - 1;

            size_t lo = find_first_syncmer_ge(mv, k, reg_lo_0b);

            bridge_cand_t *cands = NULL;
            int n_cand = 0, cap_cand = 0;

            for (size_t i = lo; i < mv->n; ++i) {
                int32_t end_0b   = (int32_t)((mv->a[i].y >> 1) & 0x7fffffff);
                int32_t start_0b = end_0b - k + 1;
                if (start_0b >= reg_hi_0b) break;

                uint64_t h = mv->a[i].x >> 8;
                int freq;
                mg_idx_get(gi, h, &freq);

                if (n_cand == cap_cand) {
                    cap_cand = cap_cand ? cap_cand * 2 : 8;
                    cands = (bridge_cand_t*)realloc(cands, cap_cand * sizeof(bridge_cand_t));
                }
                cands[n_cand].hash   = h;
                cands[n_cand].q_end  = end_0b;
                cands[n_cand].q_span = mv->a[i].x & 0xff;
                cands[n_cand].freq   = freq;
                n_cand++;
            }

            if (n_cand == 0) {
                free(cands);
                continue;
            }

            qsort(cands, n_cand, sizeof(bridge_cand_t), cmp_bridge_cand);

            int chosen = -1;
            for (int ci = 0; ci < n_cand; ++ci) {
                if (cands[ci].freq <= 0) continue;

                const uint64_t *cr;
                int t;
                cr = mg_idx_get(gi, cands[ci].hash, &t);
                if (t <= 0 || t > freq_threshold) continue;

                int32_t cand_q_end   = cands[ci].q_end;
                int32_t cand_q_span  = cands[ci].q_span;
                int32_t cand_q_start = cand_q_end - cand_q_span + 1;
                int32_t gap_on_read;
                if (use_pre) gap_on_read = rs - cand_q_end - 2;
                else         gap_on_read = cand_q_start - re + 1;
                if (gap_on_read < 0) gap_on_read = 0;

                int usable = 0;
                for (int p = 0; p < t && !usable; ++p) {
                    uint32_t rid_ref    = (uint32_t)(cr[p] >> 32);
                    uint32_t pos_ref_0b = (uint32_t)(cr[p] >> 1) & 0x7fffffff;
                    int32_t  pos_ref_1b = (int32_t)pos_ref_0b + 1;
                    if (rid_ref >= (uint32_t)sh->n_seg) continue;
                    const rep_interval_arr_t *arr = &sh->segs[rid_ref];
                    if (arr->n == 0) continue;

                    int cand_lo, cand_hi;
                    if (use_pre) {
                        cand_lo = find_first_interval_start_ge(arr, pos_ref_1b);
                        cand_hi = arr->n;
                        while (cand_hi > cand_lo &&
                               arr->start[cand_hi - 1] - pos_ref_1b > max_dist) cand_hi--;
                    } else {
                        int last = find_last_interval_end_le(arr, pos_ref_1b);
                        cand_lo = 0;
                        cand_hi = last + 1;
                        while (cand_lo < cand_hi &&
                               pos_ref_1b - arr->end[cand_lo] > max_dist) cand_lo++;
                    }
                    if (cand_hi > cand_lo) usable = 1;
                }
                if (usable) { chosen = ci; break; }
            }

            if (chosen < 0) { free(cands); continue; }
            if (cands[chosen].freq > freq_threshold) { free(cands); continue; }

            found               = 1;
            use_pre_final       = use_pre;
            best_hash_final     = cands[chosen].hash;
            bridge_q_end_final  = cands[chosen].q_end;
            bridge_q_span_final = cands[chosen].q_span;
            free(cands);
        }

        if (!found) {
            prev_itv = itv;
            continue;
        }

        int use_pre = use_pre_final;
        uint64_t best_hash      = best_hash_final;
        int32_t bridge_q_end    = bridge_q_end_final;
        int32_t bridge_q_span   = bridge_q_span_final;
        int32_t bridge_q_start  = bridge_q_end - bridge_q_span + 1;

        const uint64_t *cr;
        int t;
        cr = mg_idx_get(gi, best_hash, &t);
        if (t <= 0) { prev_itv = itv; continue; }

        size_t ir_lo = find_first_syncmer_ge(mv, k, rs - 1);
        size_t ir_hi = find_first_syncmer_ge(mv, k, re - 1);

        int bridge_dir = use_pre ? +1 : -1;
        int inner_dir  = -bridge_dir;

        for (int p = 0; p < t; ++p) {
            uint32_t rid_ref    = (uint32_t)(cr[p] >> 32);
            uint32_t pos_ref_0b = (uint32_t)(cr[p] >> 1) & 0x7fffffff;
            int32_t  pos_ref_1b = (int32_t)pos_ref_0b + 1;

            if (rid_ref >= (uint32_t)sh->n_seg) continue;
            const rep_interval_arr_t *arr = &sh->segs[rid_ref];
            if (arr->n == 0) continue;
            const gfa_seg_t *ref_seg = &gi->g->seg[rid_ref];
            const char *ref_seq = ref_seg->seq;
            int ref_len = ref_seg->len;

            int cand_lo, cand_hi;
            if (use_pre) {
                cand_lo = find_first_interval_start_ge(arr, pos_ref_1b);
                cand_hi = arr->n;
                while (cand_hi > cand_lo &&
                       arr->start[cand_hi - 1] - pos_ref_1b > max_dist) cand_hi--;
            } else {
                int last = find_last_interval_end_le(arr, pos_ref_1b);
                cand_lo = 0;
                cand_hi = last + 1;
                while (cand_lo < cand_hi &&
                       pos_ref_1b - arr->end[cand_lo] > max_dist) cand_lo++;
            }
            if (cand_hi <= cand_lo) continue;

            /* ---- 桥接 k-mer 锚点 ---- */
            {
                int32_t q_span_final = bridge_q_span;
                int32_t q_end_final  = bridge_q_end;
                int32_t r_end_final  = (int32_t)pos_ref_0b;

                if (bridge_dir > 0) {
                    int q_pos = bridge_q_end + 1;
                    int r_pos = (int)pos_ref_0b + 1;
                    int ext = 0, max_ext = 255 - q_span_final;
                    while (ext < max_ext && q_pos < qlen && r_pos < ref_len) {
                        if (query_seq[q_pos] != ref_seq[r_pos]) break;
                        q_pos++; r_pos++; ext++;
                    }
                    if (ext > 0) {
                        q_span_final += ext;
                        q_end_final   = bridge_q_end + ext;
                        r_end_final   = (int32_t)pos_ref_0b + ext;
                    }
                } else {
                    int q_pos = bridge_q_start - 1;
                    int r_pos = (int)pos_ref_0b - bridge_q_span;
                    int ext = 0, max_ext = 255 - q_span_final;
                    while (ext < max_ext && q_pos >= 0 && r_pos >= 0) {
                        if (query_seq[q_pos] != ref_seq[r_pos]) break;
                        q_pos--; r_pos--; ext++;
                    }
                    if (ext > 0) q_span_final += ext;
                }

                /* ===== 新增：坐标校验，跳过无效锚点 ===== */
                if (q_span_final <= 0 || q_span_final > 255) goto skip_bridge;
                if (q_end_final < 0 || q_end_final >= qlen) goto skip_bridge;
                if (r_end_final < 0 || r_end_final >= ref_len) goto skip_bridge;
                if (q_end_final - q_span_final + 1 < 0) goto skip_bridge;
                if (r_end_final - q_span_final + 1 < 0) goto skip_bridge;

                if (added == cap_added) {
                    cap_added = cap_added ? cap_added * 2 : 16;
                    new_a = (mg128_t*)realloc(new_a, cap_added * sizeof(mg128_t));
                }
                new_a[added].x  = (uint64_t)rid_ref << 33 | (uint32_t)r_end_final;
                new_a[added].y  = (uint64_t)q_span_final << 32 | (uint32_t)q_end_final;
                new_a[added].y |= (uint64_t)1 << MG_SEED_OCC_SHIFT;
                added++;
            }
skip_bridge:;

            /* ===== 共享 k-mer 锚点：用 all_syncs 枚举 R_f 内的 k-mer ===== */
            for (int ci2 = cand_lo; ci2 < cand_hi; ++ci2) {
                int32_t is    = arr->start[ci2];
                int32_t ie    = arr->end[ci2];
                int32_t is_0b = is - 1;
                int32_t ie_0b = ie - 1;

                size_t rf_lo = find_first_syncmer_in_ref(&gi->all_syncs, k, (int32_t)rid_ref, is_0b);
                size_t rf_hi = rf_lo;
                while (rf_hi < gi->all_syncs.n) {
                    uint32_t sid = (uint32_t)(gi->all_syncs.a[rf_hi].y >> 32);
                    if (sid != rid_ref) break;
                    int32_t end_0b = (int32_t)((gi->all_syncs.a[rf_hi].y >> 1) & 0x7fffffff);
                    int32_t start_0b = end_0b - k + 1;
                    if (start_0b >= ie_0b) break;
                    rf_hi++;
                }
                if (rf_hi <= rf_lo) continue;

                int n_rf = (int)(rf_hi - rf_lo);
                rf_kmer_t *rf_kmers = (rf_kmer_t*)malloc(n_rf * sizeof(rf_kmer_t));
                if (!rf_kmers) continue;
                for (size_t jj = rf_lo; jj < rf_hi; ++jj) {
                    rf_kmers[jj - rf_lo].hash   = gi->all_syncs.a[jj].x >> 8;
                    rf_kmers[jj - rf_lo].end_0b = (int32_t)((gi->all_syncs.a[jj].y >> 1) & 0x7fffffff);
                }
                qsort(rf_kmers, n_rf, sizeof(rf_kmer_t), cmp_rf_kmer);

                for (size_t j = ir_lo; j < ir_hi; ++j) {
                    uint64_t h      = mv->a[j].x >> 8;
                    int32_t q_span  = mv->a[j].x & 0xff;
                    int32_t q_end   = (int32_t)((mv->a[j].y >> 1) & 0x7fffffff);
                    int32_t q_start = q_end - q_span + 1;

                    int lo2 = 0, hi2 = n_rf;
                    while (lo2 < hi2) {
                        int mid2 = (lo2 + hi2) / 2;
                        if (rf_kmers[mid2].hash < h) lo2 = mid2 + 1;
                        else hi2 = mid2;
                    }
                    int first = lo2;

                    lo2 = first; hi2 = n_rf;
                    while (lo2 < hi2) {
                        int mid2 = (lo2 + hi2) / 2;
                        if (rf_kmers[mid2].hash <= h) lo2 = mid2 + 1;
                        else hi2 = mid2;
                    }
                    int last = lo2;

                    for (int idx = first; idx < last; ++idx) {
                        int32_t pos_r2_0b = rf_kmers[idx].end_0b;

                        int32_t q_span_final = q_span;
                        int32_t q_end_final  = q_end;
                        int32_t r_end_final  = pos_r2_0b;

                        if (inner_dir > 0) {
                            int q_pos = q_end + 1;
                            int r_pos = (int)pos_r2_0b + 1;
                            int ext = 0, max_ext = 255 - q_span_final;
                            while (ext < max_ext && q_pos < qlen && r_pos < ref_len) {
                                if (query_seq[q_pos] != ref_seq[r_pos]) break;
                                q_pos++; r_pos++; ext++;
                            }
                            if (ext > 0) {
                                q_span_final += ext;
                                q_end_final   = q_end + ext;
                                r_end_final   = (int32_t)pos_r2_0b + ext;
                            }
                        } else {
                            int q_pos = q_start - 1;
                            int r_pos = (int)pos_r2_0b - q_span;
                            int ext = 0, max_ext = 255 - q_span_final;
                            while (ext < max_ext && q_pos >= 0 && r_pos >= 0) {
                                if (query_seq[q_pos] != ref_seq[r_pos]) break;
                                q_pos--; r_pos--; ext++;
                            }
                            if (ext > 0) q_span_final += ext;
                        }

                        /* ===== 新增：坐标校验，跳过无效锚点 ===== */
                        if (q_span_final <= 0 || q_span_final > 255) continue;
                        if (q_end_final < 0 || q_end_final >= qlen) continue;
                        if (r_end_final < 0 || r_end_final >= ref_len) continue;
                        if (q_end_final - q_span_final + 1 < 0) continue;
                        if (r_end_final - q_span_final + 1 < 0) continue;

                        if (added == cap_added) {
                            cap_added = cap_added ? cap_added * 2 : 16;
                            new_a = (mg128_t*)realloc(new_a, cap_added * sizeof(mg128_t));
                        }
                        new_a[added].x  = (uint64_t)rid_ref << 33 | (uint32_t)r_end_final;
                        new_a[added].y  = (uint64_t)q_span_final << 32 | (uint32_t)q_end_final;
                        new_a[added].y |= (uint64_t)1 << MG_SEED_OCC_SHIFT;
                        added++;
                    }
                }

                free(rf_kmers);
            }
        }

        prev_itv = itv;
    }

    if (added > 0) {
        int64_t new_n = old_n + added;
        mg128_t *a = (mg128_t*)krealloc(km, *pa, new_n * sizeof(mg128_t));
        if (!a) {
            fprintf(stderr, "[E::%s] 锚点 realloc 失败\n", __func__);
            free(new_a);
            abort();
        }
        *pa = a;
        memcpy(&a[old_n], new_a, added * sizeof(mg128_t));
        *pn_a = new_n;
    }

    free(new_a);
}