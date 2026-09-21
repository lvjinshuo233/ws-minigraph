#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "kalloc.h"
#include "mgpriv.h"
#include "khashl.h"
#include "sys.h"
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include "khashl.h"   
#include "ws-strobe.h" 
__thread const char *g_qname = NULL;
// 改为自定义哈希函数
extern volatile long long total_anchors;
extern volatile long long total_anchors_rem;
static pthread_mutex_t anchor_mutex = PTHREAD_MUTEX_INITIALIZER;
static kh_inline khint_t pos_int_hash(int32_t x) { return (khint_t)x; }
static kh_inline int pos_int_eq(int32_t a, int32_t b) { return a == b; }
KHASHL_MAP_INIT(static, pos2idx_t, pos2idx, int32_t, int32_t, pos_int_hash, pos_int_eq)
struct mg_tbuf_s {
	void *km;
	int frag_gap;
};

mg_tbuf_t *mg_tbuf_init(void)
{
	mg_tbuf_t *b;
	b = (mg_tbuf_t*)calloc(1, sizeof(mg_tbuf_t));
	if (!(mg_dbg_flag & MG_DBG_NO_KALLOC)) b->km = km_init();
	return b;
}

void mg_tbuf_destroy(mg_tbuf_t *b)
{
	if (b == 0) return;
	if (b->km) km_destroy(b->km);
	free(b);
}

void *mg_tbuf_get_km(mg_tbuf_t *b)
{
	return b->km;
}

static void collect_minimizers(void *km, const mg_mapopt_t *opt, const mg_idx_t *gi, int n_segs, const int *qlens, const char **seqs, mg128_v *mv, seg_rep_t *reads_rep)
{
	int i, n, sum = 0;
	mv->n = 0;
	
	for (i = n = 0; i < n_segs; ++i) {
		size_t j;
		open_sync_sketch(km, seqs[i], qlens[i], gi->k, gi->s, gi->t, gi->flag&MM_I_HPC, i, mv, gi, 3, reads_rep, gi->min_rep_len);
        //fix_seg_rep(reads_rep, gi->s, gi->min_rep_len);
        //show_seg_rep(reads_rep);
		for (j = n; j < mv->n; ++j)
			mv->a[j].y += sum << 1;
		sum += qlens[i], n = mv->n;
	}
}

#include "ksort.h"
#define heap_lt(a, b) ((a).x > (b).x)
KSORT_INIT(heap, mg128_t, heap_lt)

typedef struct {
	uint32_t n;
	uint32_t q_pos, q_span;
	uint32_t seg_id:31, is_tandem:1;
	const uint64_t *cr;
} mg_match_t;

static int cmp_int32(const void *a, const void *b) {
    return (*(int32_t*)a - *(int32_t*)b);
}
typedef struct {
    int32_t *read_starts;
    int n_read;
    int32_t **seg_pos;
    int *seg_n;
    const mg_idx_t *gi;
    int k;
    const char *method;
    uint64_t anchor_y;   // 新增：当前检查的锚点的 y 值（用于错误输出）
} match_check_ctx_t;
static void validate_kmer_match(const match_check_ctx_t *ctx,
                                int32_t qpos_end,
                                uint32_t seg_id,
                                uint32_t rpos_end)
{
    int32_t q_start = qpos_end - ctx->k + 1;
    int32_t r_start = rpos_end - ctx->k + 1;

    // 检查 read 端
    if (q_start >= 0) {
        int lo = 0, hi = ctx->n_read - 1, found = 0;
        while (lo <= hi) {
            int mid = (lo + hi) / 2;
            if (ctx->read_starts[mid] == q_start) { found = 1; break; }
            else if (ctx->read_starts[mid] < q_start) lo = mid + 1;
            else hi = mid - 1;
        }
        if (!found) {
            fprintf(stderr, "ERROR: read start position %d (end %d) not found in reads_index (method=%s, read=%s)\n",
                    q_start, qpos_end, ctx->method, g_qname);
            abort();
        }
    }

    // 检查 ref 端（懒加载 + 二分查找）
    if (seg_id >= (uint32_t)ctx->gi->g->n_seg) return;

    if (ctx->seg_pos[seg_id] == NULL) {
        char path[2048];
        snprintf(path, sizeof(path),
                 "/home/lvjinshuo/BIBM_paper/scripts/new_centromere_fraction/ref_index/%s_m/%s_%s_sample.txt",
                 ctx->method, ctx->method, ctx->gi->g->seg[seg_id].name);
        FILE *fp = fopen(path, "r");
        if (!fp) {
            fprintf(stderr, "ERROR: cannot open ref_index file %s\n", path);
            abort();
        }
        int32_t *tmp_arr = NULL;
        int n_tmp = 0;
        int32_t tmp;
        while (fscanf(fp, "%d", &tmp) == 1) {
            tmp_arr = (int32_t*)realloc(tmp_arr, (n_tmp + 1) * sizeof(int32_t));
            tmp_arr[n_tmp++] = tmp;
        }
        fclose(fp);
        if (n_tmp > 0) {
            qsort(tmp_arr, n_tmp, sizeof(int32_t), cmp_int32);
            int u = 1;
            for (int j = 1; j < n_tmp; ++j)
                if (tmp_arr[j] != tmp_arr[j - 1])
                    tmp_arr[u++] = tmp_arr[j];
            n_tmp = u;
        }
        int32_t **pp = (int32_t**)&ctx->seg_pos[seg_id];
        int *pn = (int*)&ctx->seg_n[seg_id];
        *pp = tmp_arr;
        *pn = n_tmp;
    }

    if (ctx->seg_n[seg_id] > 0) {
        int lo = 0, hi = ctx->seg_n[seg_id] - 1, found = 0;
        while (lo <= hi) {
            int mid = (lo + hi) / 2;
            if (ctx->seg_pos[seg_id][mid] == r_start) { found = 1; break; }
            else if (ctx->seg_pos[seg_id][mid] < r_start) lo = mid + 1;
            else hi = mid - 1;
        }
        if (!found) {
            fprintf(stderr, "ERROR: ref start position %d (end %u) not found in ref_index for segment %s (method=%s, read=%s)\n",
                    r_start, rpos_end, ctx->gi->g->seg[seg_id].name, ctx->method, g_qname);
            abort();
        }
    }
}
static mg_match_t *collect_matches(void *km, int *_n_m, int max_occ, const mg_idx_t *gi,
                                   const mg128_v *mv, int64_t *n_a, int *rep_len, int *n_mini_pos,
                                   int32_t **mini_pos, int32_t **filt_pos, int32_t *n_filt,
                                   int32_t **all_mini_info, int64_t **all_ref_pos)
{
    int rep_st = 0, rep_en = 0, n_m;
    size_t i;
    mg_match_t *m;
    *n_mini_pos = 0;
    KMALLOC(km, *mini_pos, mv->n);
    m = (mg_match_t*)kmalloc(km, mv->n * sizeof(mg_match_t));
    int32_t *fp = 0;
    int n_fp = 0;
    KMALLOC(km, fp, mv->n);
    KMALLOC(km, *all_mini_info, mv->n);
    KMALLOC(km, *all_ref_pos, mv->n);
    int64_t *aref = *all_ref_pos;
    int32_t *ainfo = *all_mini_info;

    for (i = 0, n_m = 0, *rep_len = 0, *n_a = 0; i < mv->n; ++i) {
        const uint64_t *cr;
        mg128_t *p = &mv->a[i];
        uint32_t q_pos = (uint32_t)p->y, q_span = p->x & 0xff;
        int t;
        cr = mg_idx_get(gi, p->x>>8, &t);
        int64_t ref_val = -1;
        if (t > 0) ref_val = (int64_t)cr[0];
        aref[i] = ref_val;
        if (t >= max_occ) {
            ainfo[i] = 1;
            fp[n_fp++] = q_pos >> 1;
            int en = (q_pos >> 1) + 1, st = en - q_span;
            if (st > rep_en) {
                *rep_len += rep_en - rep_st;
                rep_st = st, rep_en = en;
            } else rep_en = en;
        } else {
            ainfo[i] = (t << 1);
            mg_match_t *q = &m[n_m++];
            q->q_pos = q_pos, q->q_span = q_span, q->cr = cr, q->n = t, q->seg_id = p->y >> 32;
            q->is_tandem = 0;
            if (i > 0 && p->x>>8 == mv->a[i - 1].x>>8) q->is_tandem = 1;
            if (i < mv->n - 1 && p->x>>8 == mv->a[i + 1].x>>8) q->is_tandem = 1;
            *n_a += q->n;
            (*mini_pos)[(*n_mini_pos)++] = q_pos >> 1;
        }
    }
    *rep_len += rep_en - rep_st;
    *_n_m = n_m;
    *filt_pos = fp;
    *n_filt = n_fp;
    return m;
}
/* ---------- 修改后的 collect_matches ---------- */
/*static mg_match_t *collect_matches(void *km, int *_n_m, int max_occ, const mg_idx_t *gi,
                                   const mg128_v *mv, int64_t *n_a, int *rep_len, int *n_mini_pos,
                                   int32_t **mini_pos, int32_t **filt_pos, int32_t *n_filt,
                                   int32_t **all_mini_info, int64_t **all_ref_pos)
{
    int rep_st = 0, rep_en = 0, n_m;
    size_t i;
    mg_match_t *m;
    *n_mini_pos = 0;
    KMALLOC(km, *mini_pos, mv->n);
    m = (mg_match_t*)kmalloc(km, mv->n * sizeof(mg_match_t));
    int32_t *fp = 0;
    int n_fp = 0;
    KMALLOC(km, fp, mv->n);
    KMALLOC(km, *all_mini_info, mv->n);
    KMALLOC(km, *all_ref_pos, mv->n);
    int64_t *aref = *all_ref_pos;
    int32_t *ainfo = *all_mini_info;

    // ---------- 验证上下文（仅对 read S1_21 启用） ----------
    static int validated_count = 0;
    static match_check_ctx_t ctx;
    static int ctx_initialized = 0;
    if (g_qname && gi && strcmp(g_qname, "S1_21") == 0) {
        const char *method = "ws";    // 可根据版本改为 "s"
        if (!ctx_initialized) {
            ctx.method = method;
            ctx.gi = gi;
            ctx.k = gi->k;
            int n_seg = gi->g->n_seg;
            ctx.seg_pos = (int32_t**)calloc(n_seg, sizeof(int32_t*));
            ctx.seg_n = (int*)calloc(n_seg, sizeof(int));

            char path[2048];
            snprintf(path, sizeof(path),
                     "/home/lvjinshuo/BIBM_paper/scripts/new_centromere_fraction/reads_index/%s_m/%s_%s_sample.txt",
                     method, method, g_qname);
            FILE *fp_r = fopen(path, "r");
            if (fp_r) {
                int32_t tmp;
                while (fscanf(fp_r, "%d", &tmp) == 1) {
                    ctx.read_starts = (int32_t*)realloc(ctx.read_starts, (ctx.n_read + 1) * sizeof(int32_t));
                    ctx.read_starts[ctx.n_read++] = tmp;
                }
                fclose(fp_r);
                if (ctx.n_read > 0) {
                    qsort(ctx.read_starts, ctx.n_read, sizeof(int32_t), cmp_int32);
                    int u = 1;
                    for (size_t ii = 1; ii < ctx.n_read; ++ii)
                        if (ctx.read_starts[ii] != ctx.read_starts[ii-1])
                            ctx.read_starts[u++] = ctx.read_starts[ii];
                    ctx.n_read = u;
                }
            } else {
                fprintf(stderr, "WARNING: reads_index file %s not found, skip validation\n", path);
            }
            ctx_initialized = 1;
        }
    }

    // ---------- 原有匹配收集循环 ----------
    for (i = 0, n_m = 0, *rep_len = 0, *n_a = 0; i < mv->n; ++i) {
        const uint64_t *cr;
        mg128_t *p = &mv->a[i];
        uint32_t q_pos = (uint32_t)p->y, q_span = p->x & 0xff;
        int t;
        cr = mg_idx_get(gi, p->x>>8, &t);
        int64_t ref_val = -1;
        if (t > 0) ref_val = (int64_t)cr[0];
        aref[i] = ref_val;

        // 仅 S1_21 时验证
        if (g_qname && strcmp(g_qname, "S1_21") == 0 && ctx_initialized && ctx.read_starts) {
            int32_t qpos_end = q_pos >> 1;
            for (int k = 0; k < t; ++k) {
                uint32_t sid = (uint32_t)(cr[k] >> 32);
                uint32_t rp = (uint32_t)(cr[k] >> 1) & 0x7FFFFFFF;
                validate_kmer_match(&ctx, qpos_end, sid, rp);
                fprintf(stderr, "VALIDATE OK #%d: read=%s, read_start=%d (end=%d) -> seg=%s, ref_start=%u (end=%u)\n",
                        ++validated_count, g_qname,
                        qpos_end - gi->k + 1, qpos_end,
                        (sid < gi->g->n_seg) ? gi->g->seg[sid].name : "?",
                        rp - gi->k + 1, rp);
            }
        }

        // 原有匹配逻辑
        if (t >= max_occ) {
            ainfo[i] = 1;
            fp[n_fp++] = q_pos >> 1;
            int en = (q_pos >> 1) + 1, st = en - q_span;
            if (st > rep_en) {
                *rep_len += rep_en - rep_st;
                rep_st = st, rep_en = en;
            } else rep_en = en;
        } else {
            ainfo[i] = (t << 1);
            mg_match_t *q = &m[n_m++];
            q->q_pos = q_pos, q->q_span = q_span, q->cr = cr, q->n = t, q->seg_id = p->y >> 32;
            q->is_tandem = 0;
            if (i > 0 && p->x>>8 == mv->a[i - 1].x>>8) q->is_tandem = 1;
            if (i < mv->n - 1 && p->x>>8 == mv->a[i + 1].x>>8) q->is_tandem = 1;
            *n_a += q->n;
            (*mini_pos)[(*n_mini_pos)++] = q_pos >> 1;
        }
    }
    *rep_len += rep_en - rep_st;
    *_n_m = n_m;
    *filt_pos = fp;
    *n_filt = n_fp;
    return m;
}*/
//检查ref和read文件是否与工具存储一致的版本
/*static mg_match_t *collect_matches(void *km, int *_n_m, int max_occ, const mg_idx_t *gi, const mg128_v *mv, int64_t *n_a, int *rep_len, int *n_mini_pos, int32_t **mini_pos, int32_t **filt_pos, int32_t *n_filt,int32_t **all_mini_info, int64_t **all_ref_pos)
{
    int rep_st = 0, rep_en = 0, n_m;
    size_t i;
    mg_match_t *m;
    *n_mini_pos = 0;
    KMALLOC(km, *mini_pos, mv->n);
    m = (mg_match_t*)kmalloc(km, mv->n * sizeof(mg_match_t));
    int32_t *fp = 0;
    int n_fp = 0;
    KMALLOC(km, fp, mv->n);
    KMALLOC(km, *all_mini_info, mv->n);
    KMALLOC(km, *all_ref_pos, mv->n);
    int64_t *aref = *all_ref_pos;
    int32_t *ainfo = *all_mini_info;
// ========== 调试：检查 read 和 ref 的 k‑mer 采样位置是否与索引文件一致（纯 C） ==========
if (g_qname && gi && gi->g && gi->g->n_seg > 0) {
    const char *method = "ws";
    char path[2048];
    int32_t *mv_starts_arr = NULL;
    int n_mv_starts = 0;
    int32_t *file_starts_arr = NULL;
    int n_file_starts = 0;
    int i, j, cmp_ok;

    // ---- 1. 提取 mv 中的起始位置并去重排序 ----
    for (i = 0; i < (int)mv->n; ++i) {
        int32_t end_pos = (int32_t)(mv->a[i].y) >> 1;
        int32_t start_pos = end_pos - gi->k + 1;
        if (start_pos >= 0) {
            mv_starts_arr = (int32_t*)realloc(mv_starts_arr, (n_mv_starts + 1) * sizeof(int32_t));
            mv_starts_arr[n_mv_starts++] = start_pos;
        }
    }
    if (n_mv_starts > 0) {
        qsort(mv_starts_arr, n_mv_starts, sizeof(int32_t), cmp_int32);
        int unique = 1;
        for (i = 1; i < n_mv_starts; ++i) {
            if (mv_starts_arr[i] != mv_starts_arr[i - 1])
                mv_starts_arr[unique++] = mv_starts_arr[i];
        }
        n_mv_starts = unique;
    }

    // ---- 2. 读取 reads_index 文件 ----
    snprintf(path, sizeof(path),
             "/home/lvjinshuo/BIBM_paper/scripts/new_centromere_fraction/reads_index/%s_m/%s_%s_sample.txt",
             method, method, g_qname);
    FILE *fp = fopen(path, "r");
    if (!fp) {
        fprintf(stderr, "ERROR: cannot open reads_index file %s\n", path);
        abort();
    }
    int32_t tmp;
    while (fscanf(fp, "%d", &tmp) == 1) {
        file_starts_arr = (int32_t*)realloc(file_starts_arr, (n_file_starts + 1) * sizeof(int32_t));
        file_starts_arr[n_file_starts++] = tmp;
    }
    fclose(fp);
    if (n_file_starts > 0) {
        qsort(file_starts_arr, n_file_starts, sizeof(int32_t), cmp_int32);
        int unique = 1;
        for (i = 1; i < n_file_starts; ++i) {
            if (file_starts_arr[i] != file_starts_arr[i - 1])
                file_starts_arr[unique++] = file_starts_arr[i];
        }
        n_file_starts = unique;
    }

    // ---- 3. 比较两者 ----
    cmp_ok = (n_mv_starts == n_file_starts);
    if (cmp_ok) {
        for (i = 0; i < n_mv_starts; ++i) {
            if (mv_starts_arr[i] != file_starts_arr[i]) {
                cmp_ok = 0;
                break;
            }
        }
    }
    if (!cmp_ok) {
        fprintf(stderr, "MISMATCH for read %s: mv has %d unique starts, file has %d\n",
                g_qname, n_mv_starts, n_file_starts);
        i = j = 0;
        int diff_cnt = 0;
        while ((i < n_mv_starts || j < n_file_starts) && diff_cnt < 20) {
            if (i < n_mv_starts && j < n_file_starts && mv_starts_arr[i] == file_starts_arr[j]) {
                ++i; ++j;
            } else {
                if (i < n_mv_starts && (j == n_file_starts || mv_starts_arr[i] < file_starts_arr[j])) {
                    fprintf(stderr, "  MV only: %d\n", mv_starts_arr[i]); ++i;
                } else if (j < n_file_starts && (i == n_mv_starts || file_starts_arr[j] < mv_starts_arr[i])) {
                    fprintf(stderr, "  FILE only: %d\n", file_starts_arr[j]); ++j;
                }
                ++diff_cnt;
            }
        }
        abort();
    } else {
        fprintf(stderr, "Read %s: PASS (read kmer positions match)\n", g_qname);
    }
    free(mv_starts_arr);
    free(file_starts_arr);

    // ---- 4. 检查 ref 端 ----
    int n_seg = gi->g->n_seg;
    int32_t **seg_pos = (int32_t**)calloc(n_seg, sizeof(int32_t*));
    int *seg_n = (int*)calloc(n_seg, sizeof(int));
    for (i = 0; i < n_seg; ++i) {
        snprintf(path, sizeof(path),
                 "/home/lvjinshuo/BIBM_paper/scripts/new_centromere_fraction/ref_index/ws_m/ws_%s_sample.txt",
                 gi->g->seg[i].name);
        fp = fopen(path, "r");
        if (!fp) {
            fprintf(stderr, "WARNING: cannot open ref_index file %s, skipped\n", path);
            continue;
        }
        while (fscanf(fp, "%d", &tmp) == 1) {
            seg_pos[i] = (int32_t*)realloc(seg_pos[i], (seg_n[i] + 1) * sizeof(int32_t));
            seg_pos[i][seg_n[i]++] = tmp;
        }
        fclose(fp);
        if (seg_n[i] > 0) {
            qsort(seg_pos[i], seg_n[i], sizeof(int32_t), cmp_int32);
            int u = 1;
            for (j = 1; j < seg_n[i]; ++j) {
                if (seg_pos[i][j] != seg_pos[i][j - 1])
                    seg_pos[i][u++] = seg_pos[i][j];
            }
            seg_n[i] = u;
        }
    }

    int ref_ok = 1;
    for (i = 0; i < (int)mv->n && ref_ok; ++i) {
        int t;
        const uint64_t *cr = mg_idx_get(gi, mv->a[i].x >> 8, &t);
        if (t <= 0) continue;
        for (int k = 0; k < t; ++k) {
            uint32_t sid = (uint32_t)(cr[k] >> 32);
            uint32_t rp_end = (uint32_t)(cr[k] >> 1) & 0x7FFFFFFF;
            int32_t rp_start = rp_end - gi->k + 1;
            if (rp_start < 0) continue;
            if (sid >= (uint32_t)n_seg || !seg_pos[sid]) {
                fprintf(stderr, "REF MISMATCH: seg_id=%u not found\n", sid);
                ref_ok = 0;
                break;
            }
            int lo = 0, hi = seg_n[sid] - 1, found = 0;
            while (lo <= hi) {
                int mid = (lo + hi) / 2;
                if (seg_pos[sid][mid] == rp_start) { found = 1; break; }
                else if (seg_pos[sid][mid] < rp_start) lo = mid + 1;
                else hi = mid - 1;
            }
            if (!found) {
                fprintf(stderr, "REF MISMATCH: seg=%s, read_pos=%d, ref_start=%d not in ref_index file\n",
                        gi->g->seg[sid].name, (int32_t)(mv->a[i].y) >> 1, rp_start);
                ref_ok = 0;
                break;
            }
        }
    }
    if (ref_ok) {
        fprintf(stderr, "Read %s: PASS (ref kmer positions match)\n", g_qname);
    } else {
        abort();
    }

    for (i = 0; i < n_seg; ++i) free(seg_pos[i]);
    free(seg_pos);
    free(seg_n);

}
    for (i = 0, n_m = 0, *rep_len = 0, *n_a = 0; i < mv->n; ++i) {
        const uint64_t *cr;
        mg128_t *p = &mv->a[i];
        uint32_t q_pos = (uint32_t)p->y, q_span = p->x & 0xff;
        int t;
        cr = mg_idx_get(gi, p->x>>8, &t);
        int64_t ref_val = -1;
        if (t > 0) ref_val = (int64_t)cr[0];
        aref[i] = ref_val;
        if (t >= max_occ) {
            ainfo[i] = 1; 
            fp[n_fp++] = q_pos >> 1;   // 结束位置
            int en = (q_pos >> 1) + 1, st = en - q_span;
            if (st > rep_en) {
                *rep_len += rep_en - rep_st;
                rep_st = st, rep_en = en;
            } else rep_en = en;
        } else {
            ainfo[i] = (t << 1);
            mg_match_t *q = &m[n_m++];
            q->q_pos = q_pos, q->q_span = q_span, q->cr = cr, q->n = t, q->seg_id = p->y >> 32;
            q->is_tandem = 0;
            if (i > 0 && p->x>>8 == mv->a[i - 1].x>>8) q->is_tandem = 1;
            if (i < mv->n - 1 && p->x>>8 == mv->a[i + 1].x>>8) q->is_tandem = 1;
            *n_a += q->n;
            (*mini_pos)[(*n_mini_pos)++] = q_pos >> 1;   // 结束位置
        }
    }
    *rep_len += rep_en - rep_st;
    *_n_m = n_m;
    *filt_pos = fp;
    *n_filt = n_fp;
    return m;
}*/

/* 旧版 collect_matches 函数（已注释） */
/*static mg_match_t *collect_matches(void *km, int *_n_m, int max_occ, const mg_idx_t *gi, const mg128_v *mv, int64_t *n_a, int *rep_len, int *n_mini_pos, int32_t **mini_pos, int32_t **filt_pos, int32_t *n_filt)
{
    const int keep_step = 1024;   // 每1024bp保留一个高频kmer
    const int retained_occ = max_occ - 1;   

    int rep_en = 0, n_m;
    size_t i;
    mg_match_t *m;
    *n_mini_pos = 0;
    KMALLOC(km, *mini_pos, mv->n);
    m = (mg_match_t*)kmalloc(km, mv->n * sizeof(mg_match_t));
    int32_t *fp = 0;
    int n_fp = 0;
    KMALLOC(km, fp, mv->n);

    mg_match_t *rep_buf = 0;
    int rep_buf_n = 0, rep_buf_m = 0;
    int cur_rep_st = 0, cur_rep_en = 0;

    for (i = 0, n_m = 0, *rep_len = 0, *n_a = 0; i < mv->n; ++i) {
        const uint64_t *cr;
        mg128_t *p = &mv->a[i];
        uint32_t q_pos = (uint32_t)p->y, q_span = p->x & 0xff;
        int t;
        cr = mg_idx_get(gi, p->x>>8, &t);

        if (t >= max_occ) { // 高频kmer
            int en = (q_pos >> 1) + 1, st = en - q_span;

            if (rep_buf_n == 0) {
                cur_rep_st = st;
                cur_rep_en = en;
            } else {
                if (st > rep_en) { // 不连续，结束前一区域
                    int region_len = cur_rep_en - cur_rep_st;
                    int keep_n = (region_len + keep_step - 1) / keep_step;
                    if (keep_n > 0 && rep_buf_n > 0) {
                        for (int k = 0; k < keep_n; ++k) {
                            int idx = (int)((double)k / keep_n * rep_buf_n + 0.5);
                            if (idx >= rep_buf_n) idx = rep_buf_n - 1;
                            mg_match_t *q = &m[n_m++];
                            *q = rep_buf[idx];
                            q->n = retained_occ;   // 设置为 max_occ - 1
                            *n_a += q->n;
                            (*mini_pos)[(*n_mini_pos)++] = q->q_pos >> 1;
                        }
                    }
                    *rep_len += region_len;
                    rep_buf_n = 0;
                    cur_rep_st = st;
                    cur_rep_en = en;
                } else {
                    if (en > cur_rep_en) cur_rep_en = en;
                }
            }

            // 暂存当前高频kmer
            if (rep_buf_n == rep_buf_m) {
                int new_m = rep_buf_m == 0 ? 64 : rep_buf_m * 2;
                rep_buf = (mg_match_t*)krealloc(km, rep_buf, new_m * sizeof(mg_match_t));
                rep_buf_m = new_m;
            }
            mg_match_t *tmp_q = &rep_buf[rep_buf_n++];
            tmp_q->q_pos = q_pos;
            tmp_q->q_span = q_span;
            tmp_q->cr = cr;
            tmp_q->n = t;          // 暂存原始值，不影响
            tmp_q->seg_id = p->y >> 32;
            tmp_q->is_tandem = 0;
            if (i > 0 && p->x>>8 == mv->a[i - 1].x>>8) tmp_q->is_tandem = 1;
            if (i < mv->n - 1 && p->x>>8 == mv->a[i + 1].x>>8) tmp_q->is_tandem = 1;

            rep_en = en;
        } else { // 非高频kmer
            if (rep_buf_n > 0) {
                int region_len = cur_rep_en - cur_rep_st;
                int keep_n = (region_len + keep_step - 1) / keep_step;
                if (keep_n > 0) {
                    for (int k = 0; k < keep_n; ++k) {
                        int idx = (int)((double)k / keep_n * rep_buf_n + 0.5);
                        if (idx >= rep_buf_n) idx = rep_buf_n - 1;
                        mg_match_t *q = &m[n_m++];
                        *q = rep_buf[idx];
                        q->n = retained_occ;
                        *n_a += q->n;
                        (*mini_pos)[(*n_mini_pos)++] = q->q_pos >> 1;
                    }
                }
                for (int j = 0; j < rep_buf_n; ++j)
                    fp[n_fp++] = rep_buf[j].q_pos >> 1;
                *rep_len += region_len;
                rep_buf_n = 0;
            }

            mg_match_t *q = &m[n_m++];
            q->q_pos = q_pos, q->q_span = q_span, q->cr = cr, q->n = t, q->seg_id = p->y >> 32;
            q->is_tandem = 0;
            if (i > 0 && p->x>>8 == mv->a[i - 1].x>>8) q->is_tandem = 1;
            if (i < mv->n - 1 && p->x>>8 == mv->a[i + 1].x>>8) q->is_tandem = 1;
            *n_a += q->n;
            (*mini_pos)[(*n_mini_pos)++] = q_pos>>1;

            rep_en = 0;
        }
    }

    if (rep_buf_n > 0) {
        int region_len = cur_rep_en - cur_rep_st;
        int keep_n = (region_len + keep_step - 1) / keep_step;
        if (keep_n > 0) {
            for (int k = 0; k < keep_n; ++k) {
                int idx = (int)((double)k / keep_n * rep_buf_n + 0.5);
                if (idx >= rep_buf_n) idx = rep_buf_n - 1;
                mg_match_t *q = &m[n_m++];
                *q = rep_buf[idx];
                q->n = retained_occ;
                *n_a += q->n;
                (*mini_pos)[(*n_mini_pos)++] = q->q_pos >> 1;
            }
        }
        for (int j = 0; j < rep_buf_n; ++j)
            fp[n_fp++] = rep_buf[j].q_pos >> 1;
        *rep_len += region_len;
    }

    kfree(km, rep_buf);
    *_n_m = n_m;
    *filt_pos = fp;
    *n_filt = n_fp;
    return m;
}*/   /* <--- 注意这里添加了闭合注释的 */ 

static mg128_t *collect_seed_hits_heap(void *km, const mg_mapopt_t *opt, int max_occ, const mg_idx_t *gi, const char *qname, const mg128_v *mv, int qlen, int64_t *n_a, int *rep_len,
								  int *n_mini_pos, int32_t **mini_pos,  int32_t **filt_pos, int32_t *n_filt,int32_t **all_mini_info, int64_t **all_ref_pos)
{
	int i, n_m, heap_size = 0;
	int64_t n_for = 0, n_rev = 0;
	mg_match_t *m;
	mg128_t *a, *heap;

	m = collect_matches(km, &n_m, max_occ, gi, mv, n_a, rep_len, n_mini_pos, mini_pos, filt_pos, n_filt,all_mini_info, all_ref_pos);

	heap = (mg128_t*)kmalloc(km, n_m * sizeof(mg128_t));
	a = (mg128_t*)kmalloc(km, *n_a * sizeof(mg128_t));

	for (i = 0, heap_size = 0; i < n_m; ++i) {
		if (m[i].n > 0) {
			heap[heap_size].x = m[i].cr[0];
			heap[heap_size].y = (uint64_t)i<<32;
			++heap_size;
		}
	}
	ks_heapmake_heap(heap_size, heap);
	while (heap_size > 0) {
		mg_match_t *q = &m[heap->y>>32];
		mg128_t *p;
		uint64_t r = heap->x;
		int32_t rpos = (uint32_t)r >> 1;
		// TODO: skip anchor if MG_F_NO_DIAL
		if ((r&1) == (q->q_pos&1)) { // forward strand
			p = &a[n_for++];
			p->x = r>>32<<33 | rpos;
		} else { // reverse strand; TODO: more testing needed for this block
			p = &a[(*n_a) - (++n_rev)];
			p->x = r>>32<<33 | 1ULL<<32 | (gi->g->seg[r>>32].len - (rpos + 1 - q->q_span) - 1);
		}
		p->y = (uint64_t)q->q_span << 32 | q->q_pos >> 1;
		p->y |= (uint64_t)q->seg_id << MG_SEED_SEG_SHIFT;
		if (q->is_tandem) p->y |= MG_SEED_TANDEM;
		p->y |= (uint64_t)(q->n < 255? q->n : 255) << MG_SEED_OCC_SHIFT;
		// update the heap
		if ((uint32_t)heap->y < q->n - 1) {
			++heap[0].y;
			heap[0].x = m[heap[0].y>>32].cr[(uint32_t)heap[0].y];
		} else {
			heap[0] = heap[heap_size - 1];
			--heap_size;
		}
		ks_heapdown_heap(0, heap_size, heap);
	}
	kfree(km, m);
	kfree(km, heap);

	// reverse anchors on the reverse strand, as they are in the descending order
	if (*n_a > n_for + n_rev) {
		memmove(a + n_for, a + (*n_a) - n_rev, n_rev * sizeof(mg128_t));
		*n_a = n_for + n_rev;
	}
	return a;
}

static mg128_t *collect_seed_hits(void *km, const mg_mapopt_t *opt, int max_occ, const mg_idx_t *gi, const char *qname, const mg128_v *mv, int qlen, int64_t *n_a, int *rep_len,
								  int *n_mini_pos, int32_t **mini_pos, int32_t **filt_pos, int32_t *n_filt,int32_t **all_mini_info, int64_t **all_ref_pos)
{
	int i, n_m;
	mg_match_t *m;
	mg128_t *a;
	m = collect_matches(km, &n_m, max_occ, gi, mv, n_a, rep_len, n_mini_pos, mini_pos, filt_pos, n_filt,all_mini_info, all_ref_pos);
	a = (mg128_t*)kmalloc(km, *n_a * sizeof(mg128_t));
	for (i = 0, *n_a = 0; i < n_m; ++i) {
		mg_match_t *q = &m[i];
		const uint64_t *r = q->cr;
		uint32_t k;
		for (k = 0; k < q->n; ++k) {
			int32_t rpos = (uint32_t)r[k] >> 1;
			mg128_t *p;
			if (qname && (opt->flag & MG_M_NO_DIAG)) {
				const gfa_seg_t *s = &gi->g->seg[r[k]>>32];
				const char *gname = s->snid >= 0 && gi->g->sseq? gi->g->sseq[s->snid].name : s->name;
				int32_t g_pos;
				if (s->snid >= 0 && gi->g->sseq)
					gname = gi->g->sseq[s->snid].name, g_pos = s->soff + (uint32_t)r[k];
				else
					gname = s->name, g_pos = (uint32_t)r[k];
				if ((uint32_t)g_pos == q->q_pos && strcmp(qname, gname) == 0)
					continue;
			}
			p = &a[(*n_a)++];
			if ((r[k]&1) == (q->q_pos&1)) // forward strand
				p->x = r[k]>>32<<33 | rpos;
			else // reverse strand
				p->x = r[k]>>32<<33 | 1ULL<<32 | (gi->g->seg[r[k]>>32].len - (rpos + 1 - q->q_span) - 1);
			p->y = (uint64_t)q->q_span << 32 | q->q_pos >> 1;
			p->y |= (uint64_t)q->seg_id << MG_SEED_SEG_SHIFT;
			if (q->is_tandem) p->y |= MG_SEED_TANDEM;
			p->y |= (uint64_t)(q->n < 255? q->n : 255) << MG_SEED_OCC_SHIFT;
		}
	}
	kfree(km, m);
	radix_sort_128x(a, a + (*n_a));
	return a;
}

static void mm_fix_bad_ends(const mg128_t *a, int32_t lc_max_occ, int32_t lc_max_trim, int32_t *as, int32_t *cnt)
{
	int32_t i, k, as0 = *as, cnt0 = *cnt;
	for (i = as0 + cnt0 - 1, k = 0; k < lc_max_trim && k < cnt0; ++k, --i)
		if (a[i].y>>MG_SEED_OCC_SHIFT <= (uint32_t)lc_max_occ)
			break;
	*cnt -= k;
	for (i = as0, k = 0; k < *cnt && k < lc_max_trim; ++i, ++k)
		if (a[i].y>>MG_SEED_OCC_SHIFT <= (uint32_t)lc_max_occ)
			break;
	*as += k, *cnt -= k;
}

static void mm_fix_bad_ends_alt(const mg128_t *a, int32_t score, int bw, int min_match, int32_t *as, int32_t *cnt)
{
	int32_t i, l, m, as0 = *as, cnt0 = *cnt;
	if (cnt0 < 3) return;
	m = l = a[as0].y >> 32 & 0xff;
	for (i = as0 + 1; i < as0 + cnt0 - 1; ++i) {
		int32_t lq, lr, min, max;
		int32_t q_span = a[i].y >> 32 & 0xff;
		lr = (int32_t)a[i].x - (int32_t)a[i-1].x;
		lq = (int32_t)a[i].y - (int32_t)a[i-1].y;
		min = lr < lq? lr : lq;
		max = lr > lq? lr : lq;
		if (max - min > l >> 1) *as = i;
		l += min;
		m += min < q_span? min : q_span;
		if (l >= bw << 1 || (m >= min_match && m >= bw) || m >= score>>1) break;
	}
	*cnt = as0 + cnt0 - *as;
	m = l = a[as0 + cnt0 - 1].y >> 32 & 0xff;
	for (i = as0 + cnt0 - 2; i > *as; --i) {
		int32_t lq, lr, min, max;
		int32_t q_span = a[i+1].y >> 32 & 0xff;
		lr = (int32_t)a[i+1].x - (int32_t)a[i].x;
		lq = (int32_t)a[i+1].y - (int32_t)a[i].y;
		min = lr < lq? lr : lq;
		max = lr > lq? lr : lq;
		if (max - min > l >> 1) *cnt = i + 1 - *as;
		l += min;
		m += min < q_span? min : q_span;
		if (l >= bw << 1 || (m >= min_match && m >= bw) || m >= score>>1) break;
	}
}

static int *collect_long_gaps(void *km, int as1, int cnt1, mg128_t *a, int min_gap, int *n_)
{
	int i, n, *K;
	*n_ = 0;
	for (i = 1, n = 0; i < cnt1; ++i) {
		int gap = ((int32_t)a[as1 + i].y - a[as1 + i - 1].y) - ((int32_t)a[as1 + i].x - a[as1 + i - 1].x);
		if (gap < -min_gap || gap > min_gap) ++n;
	}
	if (n <= 1) return 0;
	K = (int*)kmalloc(km, n * sizeof(int));
	for (i = 1, n = 0; i < cnt1; ++i) {
		int gap = ((int32_t)a[as1 + i].y - a[as1 + i - 1].y) - ((int32_t)a[as1 + i].x - a[as1 + i - 1].x);
		if (gap < -min_gap || gap > min_gap)
			K[n++] = i;
	}
	*n_ = n;
	return K;
}

static void mm_filter_bad_seeds(void *km, int as1, int cnt1, mg128_t *a, int min_gap, int diff_thres, int max_ext_len, int max_ext_cnt)
{
	int max_st, max_en, n, i, k, max, *K;
	K = collect_long_gaps(km, as1, cnt1, a, min_gap, &n);
	if (K == 0) return;
	max = 0, max_st = max_en = -1;
	for (k = 0;; ++k) {
		int gap, l, n_ins = 0, n_del = 0, qs, rs, max_diff = 0, max_diff_l = -1;
		if (k == n || k >= max_en) {
			if (max_en > 0)
				for (i = K[max_st]; i < K[max_en]; ++i)
					a[as1 + i].y |= MG_SEED_IGNORE;
			max = 0, max_st = max_en = -1;
			if (k == n) break;
		}
		i = K[k];
		gap = ((int32_t)a[as1 + i].y - (int32_t)a[as1 + i - 1].y) - (int32_t)(a[as1 + i].x - a[as1 + i - 1].x);
		if (gap > 0) n_ins += gap;
		else n_del += -gap;
		qs = (int32_t)a[as1 + i - 1].y;
		rs = (int32_t)a[as1 + i - 1].x;
		for (l = k + 1; l < n && l <= k + max_ext_cnt; ++l) {
			int j = K[l], diff;
			if ((int32_t)a[as1 + j].y - qs > max_ext_len || (int32_t)a[as1 + j].x - rs > max_ext_len) break;
			gap = ((int32_t)a[as1 + j].y - (int32_t)a[as1 + j - 1].y) - (int32_t)(a[as1 + j].x - a[as1 + j - 1].x);
			if (gap > 0) n_ins += gap;
			else n_del += -gap;
			diff = n_ins + n_del - abs(n_ins - n_del);
			if (max_diff < diff)
				max_diff = diff, max_diff_l = l;
		}
		if (max_diff > diff_thres && max_diff > max)
			max = max_diff, max_st = k, max_en = max_diff_l;
	}
	kfree(km, K);
}

static void mm_filter_bad_seeds_alt(void *km, int as1, int cnt1, mg128_t *a, int min_gap, int max_ext)
{
	int n, k, *K;
	K = collect_long_gaps(km, as1, cnt1, a, min_gap, &n);
	if (K == 0) return;
	for (k = 0; k < n;) {
		int i = K[k], l;
		int gap1 = ((int32_t)a[as1 + i].y - (int32_t)a[as1 + i - 1].y) - ((int32_t)a[as1 + i].x - (int32_t)a[as1 + i - 1].x);
		int re1 = (int32_t)a[as1 + i].x;
		int qe1 = (int32_t)a[as1 + i].y;
		gap1 = gap1 > 0? gap1 : -gap1;
		for (l = k + 1; l < n; ++l) {
			int j = K[l], gap2, q_span_pre, rs2, qs2, m;
			if ((int32_t)a[as1 + j].y - qe1 > max_ext || (int32_t)a[as1 + j].x - re1 > max_ext) break;
			gap2 = ((int32_t)a[as1 + j].y - (int32_t)a[as1 + j - 1].y) - (int32_t)(a[as1 + j].x - a[as1 + j - 1].x);
			q_span_pre = a[as1 + j - 1].y >> 32 & 0xff;
			rs2 = (int32_t)a[as1 + j - 1].x + q_span_pre;
			qs2 = (int32_t)a[as1 + j - 1].y + q_span_pre;
			m = rs2 - re1 < qs2 - qe1? rs2 - re1 : qs2 - qe1;
			gap2 = gap2 > 0? gap2 : -gap2;
			if (m > gap1 + gap2) break;
			re1 = (int32_t)a[as1 + j].x;
			qe1 = (int32_t)a[as1 + j].y;
			gap1 = gap2;
		}
		if (l > k + 1) {
			int j, end = K[l - 1];
			for (j = K[k]; j < end; ++j)
				a[as1 + j].y |= MG_SEED_IGNORE;
			a[as1 + end].y |= MG_SEED_FIXED;
		}
		k = l;
	}
	kfree(km, K);
}

static double print_time(double t0, int stage, const char *qname)
{
	double t;
	t = realtime();
	fprintf(stderr, "Q%d\t%s\t%.3f\n", stage, qname, t - t0);
	return t;
}

void mg_map_frag(const mg_idx_t *gi, int n_segs, const int *qlens, const char **seqs, mg_gchains_t **gcs, mg_tbuf_t *b, const mg_mapopt_t *opt, const char *qname, collect_rep_t *collect_rep)
{
	int i, l, rep_len, qlen_sum, n_lc, n_gc, n_mini_pos;
	int max_chain_gap_qry, max_chain_gap_ref, is_splice = !!(opt->flag & MG_M_SPLICE), is_sr = !!(opt->flag & MG_M_SR);
	uint32_t hash;
	int64_t n_a;
	uint64_t *u;
	int32_t *mini_pos;
	mg128_t *a;
	mg128_v mv = {0,0,0};
	mg_lchain_t *lc;
	char *seq_cat;
	km_stat_t kmst;
	float tmp, chn_pen_gap, chn_pen_skip;
	double t = 0.0;
	int64_t *all_ref_pos = NULL;  // 新增局部变量

    g_qname = qname;
    //mkdir("/home/lvjinshuo/BIBM_paper/scripts/new_centromere_fraction/reads_index/ws_m", 0755);

	for (i = 0, qlen_sum = 0; i < n_segs; ++i)
		qlen_sum += qlens[i], gcs[i] = 0;

	if (qlen_sum == 0 || n_segs <= 0 || n_segs > MG_MAX_SEG) return;
	if (opt->max_qlen > 0 && qlen_sum > opt->max_qlen) return;

	hash  = qname? kh_hash_str(qname) : 0;
	hash ^= kh_hash_uint32(qlen_sum) + kh_hash_uint32(opt->seed);
	hash  = kh_hash_uint32(hash);
	/**
	 * Collect minimizers
	 */
    seg_rep_t *reads_rep = (seg_rep_t*)calloc(n_segs, sizeof(seg_rep_t));
    if(!(init_seg_rep(reads_rep, (char*)qname))) {
        fprintf(stderr, "Error: init_reads_rep failed\n");
        assert(0);
    }
	// 收集 minimizers 和 read strobemers
collect_minimizers(b->km, opt, gi, n_segs, qlens, seqs, &mv, reads_rep);
//read_strobe_t *read_strobes = NULL;
//int n_read_strobes = collect_ws_strobes(gi, &mv, reads_rep, &read_strobes);
//free_reads_rep(reads_rep, n_segs);

    // 可在此处使用 read_strobes，例如查询参考索引或输出调试
/*for (int i = 0; i < n_read_strobes; ++i) {
    if (i < 100) {
        fprintf(stderr, "ReadStrobe: h1=%016llx h2=%016llx s1_start=%d s2_start=%d dist=%d freq1=%d freq2=%d\n",
                (unsigned long long)read_strobes[i].h1,
                (unsigned long long)read_strobes[i].h2,
                read_strobes[i].s1_start, read_strobes[i].s2_start,
                read_strobes[i].dist, read_strobes[i].freq1, read_strobes[i].freq2);
        }
    }
    free(read_strobes);*/
	/**
	 * Collect anchors
	 */
	int32_t *all_mini_pos = NULL;
	int32_t *all_mini_info = NULL;
	int32_t n_mini_total = mv.n;   // 保存原始 minimizer 总数
	int32_t *filt_pos = 0;
	int32_t n_filt = 0;
	KMALLOC(b->km, all_mini_pos, n_mini_total);
	for (int k = 0; k < n_mini_total; ++k)
		all_mini_pos[k] = (int32_t)(mv.a[k].y) >> 1;   // 结束位置 (0‑based)
	if (opt->flag & MG_M_HEAP_SORT) a = collect_seed_hits_heap(b->km, opt, opt->occ_max1, gi, qname, &mv, qlen_sum, &n_a, &rep_len, &n_mini_pos, &mini_pos, &filt_pos, &n_filt,&all_mini_info, &all_ref_pos);
	else a = collect_seed_hits(b->km, opt, opt->occ_max1, gi, qname, &mv, qlen_sum, &n_a, &rep_len, &n_mini_pos, &mini_pos, &filt_pos, &n_filt,&all_mini_info, &all_ref_pos);
int64_t removed_anchors = 0;
//remove_anchors_in_ref_rep(b->km, &a, &n_a, gi, collect_rep, &removed_anchors);
//fprintf(stderr, "REMOVED_ANCHORS\t%s\t%lld\n", qname, (long long)removed_anchors);
//__sync_fetch_and_add(&total_anchors_rem, removed_anchors);
    // 添加 strobemer 锚点
    // 添加 strobemer 锚点

    	KMALLOC(b->km, seq_cat, qlen_sum);
	for (i = l = 0; i < n_segs; ++i) {
		strncpy(&seq_cat[l], seqs[i], qlens[i]);
		l += qlens[i];
	}
//add_strobe_anchors(b->km, &a, &n_a, read_strobes, n_read_strobes, gi, seq_cat, qlen_sum);
    /* add_interval_anchors(b->km, &a, &n_a, gi, &mv, reads_rep, gi->strobe_idx,
                           seq_cat, qlen_sum, opt->occ_max1); */

    /* 新版本：用参考侧 M 数组（gi->all_syncs）定位 R_f 内的 k-mer */
    second_sampling_and_anchoring(b->km, &a, &n_a, gi, &mv, reads_rep, gi->strobe_idx,
                                  seq_cat, qlen_sum, opt->occ_max1);
    free_reads_rep(reads_rep, n_segs);

// 统计匹配情况并输出
    int n_matched_strobes = 0;
    int total_matches = 0;
/*if (gi->strobe_idx) {
    for (int i = 0; i < n_read_strobes; ++i) {
        uint64_t *vals = NULL;
        int nv = 0;
        if (ws_strobe_get_all((strobe_hash_t*)gi->strobe_idx, read_strobes[i].h1, read_strobes[i].h2, &vals, &nv)) {
            n_matched_strobes++;
            total_matches += nv;
        }
    }
}*/
//fprintf(stderr, "NUM_ANCHORS\t%s\t%lld\n", qname, (long long)n_a);
//fprintf(stderr, "STROBE_MATCH\t%s\t%d\t%d\t%d\n", qname, n_read_strobes, n_matched_strobes, total_matches);

//free(read_strobes);

// 重新排序锚点并重建 mini_pos
// 排序
if (n_a > 0) {
    radix_sort_128x(a, a + n_a);

    // 去重：去除完全相同的锚点
    int64_t write = 0;
    for (int64_t i = 0; i < n_a; ++i) {
        if (i == 0 || a[i].x != a[write - 1].x || a[i].y != a[write - 1].y) {
            a[write++] = a[i];
        }
    }
    n_a = write;
    n_mini_pos = n_a; // 去重后重新生成 mini_pos

    // 重新生成 mini_pos
    if (mini_pos) kfree(b->km, mini_pos);
    KMALLOC(b->km, mini_pos, n_a);
    for (int64_t i = 0; i < n_a; ++i) {
        mini_pos[i] = (int32_t)(a[i].y) & 0x7fffffff;
    }
    qsort(mini_pos, n_a, sizeof(int32_t), cmp_int32);
}
    //__sync_fetch_and_add(&total_anchors, n_a);

	if (mg_dbg_flag & MG_DBG_SEED) {
		fprintf(stderr, "RS\t%d\n", rep_len);
		for (i = 0; i < n_a; ++i)
			fprintf(stderr, "SD\t%s\t%d\t%c\t%d\t%d\t%d\n", gi->g->seg[a[i].x>>33].name, (int32_t)a[i].x, "+-"[a[i].x>>32&1], (int32_t)a[i].y, (int32_t)(a[i].y>>32&0xff),
					i == 0? 0 : ((int32_t)a[i].y - (int32_t)a[i-1].y) - ((int32_t)a[i].x - (int32_t)a[i-1].x));
	}

	// set max chaining gap on the query and the reference sequence
	if (is_sr)
		max_chain_gap_qry = qlen_sum > opt->max_gap? qlen_sum : opt->max_gap;
	else max_chain_gap_qry = opt->max_gap;
	if (opt->max_gap_ref > 0) {
		max_chain_gap_ref = opt->max_gap_ref;
	} else if (opt->max_frag_len > 0) {
		max_chain_gap_ref = opt->max_frag_len - qlen_sum;
		if (max_chain_gap_ref < opt->max_gap) max_chain_gap_ref = opt->max_gap;
	} else max_chain_gap_ref = opt->max_gap;

	tmp = expf(-opt->div * gi->k);
	chn_pen_gap = opt->chn_pen_gap * tmp;
	chn_pen_skip = opt->chn_pen_skip * tmp;

	if (mg_dbg_flag & MG_DBG_QNAME) t = realtime();

	if (n_a == 0) {
		if (a) kfree(b->km, a);
		a = 0, n_lc = 0, u = 0;
		if (filt_pos) {
			kfree(b->km, filt_pos);
			filt_pos = NULL;
		}
		n_filt = 0;
	} else {
		if (opt->flag & MG_M_RMQ) {
			a = mg_lchain_rmq(opt->max_gap, opt->max_gap_pre, opt->bw, opt->max_lc_skip, opt->rmq_size_cap, opt->min_lc_cnt, opt->min_lc_score,
							  chn_pen_gap, chn_pen_skip, n_a, a, &n_lc, &u, b->km);
		} else {
			a = mg_lchain_dp(max_chain_gap_ref, max_chain_gap_qry, opt->bw, opt->max_lc_skip, opt->max_lc_iter, opt->min_lc_cnt, opt->min_lc_score,
							 chn_pen_gap, chn_pen_skip, is_splice, n_segs, n_a, a, &n_lc, &u, b->km);
		}
	}
	if (mg_dbg_flag & MG_DBG_QNAME) t = print_time(t, 1, qname);

	if (opt->bw_long > opt->bw && (opt->flag & (MG_M_SPLICE|MG_M_SR)) == 0 && n_segs == 1 && n_lc > 1) {
		int32_t st = (int32_t)a[0].y, en = (int32_t)a[(int32_t)u[0] - 1].y;
		if (qlen_sum - (en - st) > opt->rmq_rescue_size || qlen_sum - (en - st) > qlen_sum * opt->rmq_rescue_ratio) {
			int32_t i;
			for (i = 0, n_a = 0; i < n_lc; ++i) n_a += (int32_t)u[i];
			kfree(b->km, u);
			radix_sort_128x(a, a + n_a);
			a = mg_lchain_rmq(opt->max_gap, opt->max_gap_pre, opt->bw_long, opt->max_lc_skip, opt->rmq_size_cap, opt->min_lc_cnt, opt->min_lc_score,
							  chn_pen_gap, chn_pen_skip, n_a, a, &n_lc, &u, b->km);
		}
	}

	b->frag_gap = max_chain_gap_ref;
	kfree(b->km, mv.a);

	if (n_lc) {
		lc = mg_lchain_gen(b->km, hash, qlen_sum, n_lc, u, a);
		if (n_lc > 1) {
			int32_t n_lc_new = 0;
			for (i = 0; i < n_lc; ++i) {
				mg_lchain_t *p = &lc[i];
				int32_t cnt = p->cnt, off = p->off;
				mm_fix_bad_ends(a, opt->lc_max_occ, opt->lc_max_trim, &off, &cnt);
				mm_fix_bad_ends_alt(a, p->score, opt->bw, 100, &off, &cnt);
				mm_filter_bad_seeds(b->km, off, cnt, a, 10, 40, opt->max_gap>>1, 10);
				mm_filter_bad_seeds_alt(b->km, off, cnt, a, 30, opt->max_gap>>1);
				p->off = off, p->cnt = cnt;
				if (cnt >= opt->min_lc_cnt) {
					int32_t q_span = a[p->off].y>>32 & 0xff;
					p->rs = (int32_t)a[p->off].x + 1 - q_span;
					p->qs = (int32_t)a[p->off].y + 1 - q_span;
					p->re = (int32_t)a[p->off + p->cnt - 1].x + 1;
					p->qe = (int32_t)a[p->off + p->cnt - 1].y + 1;
					lc[n_lc_new++] = *p;
				}
			}
			n_lc = n_lc_new;
		}
		for (i = 0; i < n_lc; ++i)
			mg_update_anchors(lc[i].cnt, &a[lc[i].off], n_mini_pos, mini_pos);
	} else lc = 0;
	kfree(b->km, mini_pos);
	kfree(b->km, u);
	if (mg_dbg_flag & MG_DBG_QNAME) t = print_time(t, 2, qname);

	if (mg_dbg_flag & MG_DBG_LCHAIN)
		mg_print_lchain(stdout, gi, n_lc, lc, a, qname);

	/*KMALLOC(b->km, seq_cat, qlen_sum);
	for (i = l = 0; i < n_segs; ++i) {
		strncpy(&seq_cat[l], seqs[i], qlens[i]);
		l += qlens[i];
	}*/
	n_gc = mg_gchain1_dp(b->km, gi->g, &n_lc, lc, qlen_sum, opt->bw_long, opt->bw_long, opt->bw_long, opt->max_gc_skip, opt->ref_bonus,
						 chn_pen_gap, chn_pen_skip, opt->mask_level, a, &u);
	if (mg_dbg_flag & MG_DBG_QNAME) t = print_time(t, 3, qname);
	gcs[0] = mg_gchain_gen(0, b->km, gi->g, gi->es, n_gc, u, lc, a, hash, opt->min_gc_cnt, opt->min_gc_score, opt->gdp_max_ed, n_segs, seq_cat);
// ---------- 验证最终锚点（种子） ----------
/*if (gcs[0] && gcs[0]->n_gc > 0 && g_qname && gi) {
    const char *method = "ws";
    match_check_ctx_t ctx;
    memset(&ctx, 0, sizeof(ctx));
    ctx.method = method;
    ctx.gi = gi;
    ctx.k = gi->k;

    // 加载 read 采样文件
    char path[2048];
    snprintf(path, sizeof(path),
             "/home/lvjinshuo/BIBM_paper/scripts/new_centromere_fraction/reads_index/%s_m/%s_%s_sample.txt",
             ctx.method, ctx.method, g_qname);
    FILE *fp = fopen(path, "r");
    if (fp) {
        int32_t tmp;
        while (fscanf(fp, "%d", &tmp) == 1) {
            ctx.read_starts = (int32_t*)realloc(ctx.read_starts, (ctx.n_read + 1) * sizeof(int32_t));
            ctx.read_starts[ctx.n_read++] = tmp;
        }
        fclose(fp);
        if (ctx.n_read > 0) {
            qsort(ctx.read_starts, ctx.n_read, sizeof(int32_t), cmp_int32);
            int u = 1;
            for (size_t i = 1; i < ctx.n_read; ++i)
                if (ctx.read_starts[i] != ctx.read_starts[i-1])
                    ctx.read_starts[u++] = ctx.read_starts[i];
            ctx.n_read = u;
        }
    } else {
        fprintf(stderr, "WARNING: cannot open reads_index file %s, skip anchor validation\n", path);
    }

    int n_seg = gi->g->n_seg;
    ctx.seg_pos = (int32_t**)calloc(n_seg, sizeof(int32_t*));
    ctx.seg_n = (int*)calloc(n_seg, sizeof(int));

    int anchor_count = 0;
    int error_count = 0;
    for (int j = 0; j < gcs[0]->n_gc; ++j) {
        const mg_gchain_t *gc = &gcs[0]->gc[j];
        if (gc->cnt == 0) continue;
        for (int k = 0; k < gc->cnt; ++k) {
            const mg_llchain_t *llc = &gcs[0]->lc[gc->off + k];
            if (llc->cnt == 0) continue;
            for (int m = 0; m < llc->cnt; ++m) {
                int idx = llc->off + m;
                if (idx >= gcs[0]->n_a) {
                    fprintf(stderr, "ERROR: anchor index %d out of range (n_a=%d) for read %s\n",
                            idx, gcs[0]->n_a, g_qname);
                    abort();
                }
                mg128_t anchor = gcs[0]->a[idx];

                // 正确的坐标提取
                int32_t q_end = (int32_t)((uint32_t)anchor.y >> 1);      // 查询结束位置 (0‑based)
                int32_t r_end = (int32_t)(uint32_t)anchor.x;              // 参考结束位置 (0‑based)
                uint32_t seg_id = (uint32_t)(anchor.x >> 32);

                if (q_end < 0 || r_end < 0 || seg_id >= (uint32_t)n_seg) {
                    fprintf(stderr, "ERROR: invalid anchor idx=%d for read %s: q_end=%d, r_end=%d, seg_id=%u, x=0x%016llx, y=0x%016llx\n",
                            idx, g_qname, q_end, r_end, seg_id,
                            (unsigned long long)anchor.x, (unsigned long long)anchor.y);
                    error_count++;
                    continue;
                }

                int32_t q_start = q_end - gi->k + 1;
                int32_t r_start = r_end - gi->k + 1;
                if (q_start >= 0 && r_start >= 0 && ctx.read_starts) {
                        ctx.anchor_y = anchor.y;          // <-- 增加这一行

                    validate_kmer_match(&ctx, q_end, seg_id, r_end);
                }
                anchor_count++;
            }
        }
    }

    if (ctx.read_starts) {
        if (error_count > 0) {
            fprintf(stderr, "Read %s: %d anchors validated, %d anchors had errors.\n", g_qname, anchor_count, error_count);
        } else {
            fprintf(stderr, "Read %s: all %d final anchors validated.\n", g_qname, anchor_count);
        }
    }

    free(ctx.read_starts);
    if (ctx.seg_pos) {
        for (int j = 0; j < n_seg; ++j) free(ctx.seg_pos[j]);
        free(ctx.seg_pos);
        free(ctx.seg_n);
    }
}*/
//abort();
	if (gcs[0] == NULL) {
		gcs[0] = (mg_gchains_t*)calloc(1, sizeof(mg_gchains_t));
		gcs[0]->km = NULL;
		gcs[0]->n_a = n_a;
		gcs[0]->rep_len = rep_len;
		gcs[0]->n_mini = n_mini_total;
		gcs[0]->n_filt = n_filt;
	}

	if (all_ref_pos && n_mini_total > 0) {
		int64_t *saved_ref_pos = (int64_t*)malloc(n_mini_total * sizeof(int64_t));
		memcpy(saved_ref_pos, all_ref_pos, n_mini_total * sizeof(int64_t));
		gcs[0]->all_ref_pos = saved_ref_pos;
		kfree(b->km, all_ref_pos);
	} else {
		gcs[0]->all_ref_pos = NULL;
	}

	if (all_mini_pos && n_mini_total > 0) {
		int32_t *saved_mini_pos = (int32_t*)malloc(n_mini_total * sizeof(int32_t));
		memcpy(saved_mini_pos, all_mini_pos, n_mini_total * sizeof(int32_t));
		gcs[0]->all_mini_pos = saved_mini_pos;
		kfree(b->km, all_mini_pos);
	} else {
		gcs[0]->all_mini_pos = NULL;
	}

	if (all_mini_info && n_mini_total > 0) {
		int32_t *saved_mini_info = (int32_t*)malloc(n_mini_total * sizeof(int32_t));
		memcpy(saved_mini_info, all_mini_info, n_mini_total * sizeof(int32_t));
		gcs[0]->all_mini_info = saved_mini_info;
		kfree(b->km, all_mini_info);
	} else {
		gcs[0]->all_mini_info = NULL;
	}
	if (mg_dbg_flag & MG_DBG_QNAME) t = print_time(t, 4, qname);

	if (gcs[0]) {
		if (filt_pos && n_filt > 0) {
			gcs[0]->filt_pos = (int32_t*)malloc(n_filt * sizeof(int32_t));
			if (gcs[0]->filt_pos) {
				memcpy(gcs[0]->filt_pos, filt_pos, n_filt * sizeof(int32_t));
				gcs[0]->n_filt = n_filt;
			} else {
				gcs[0]->n_filt = 0;
			}
			kfree(b->km, filt_pos);
		} else {
			if (filt_pos) kfree(b->km, filt_pos);
			gcs[0]->filt_pos = NULL;
			gcs[0]->n_filt = 0;
		}
	} else {
		if (filt_pos) kfree(b->km, filt_pos);
	}
	filt_pos = NULL;
	if (gcs[0]) gcs[0]->rep_len = rep_len;
	if (gcs[0]) gcs[0]->n_mini = n_mini_total;
	kfree(b->km, a);
	kfree(b->km, lc);
	kfree(b->km, u);

	mg_gchain_set_parent(b->km, opt->mask_level, gcs[0]->n_gc, gcs[0]->gc, opt->sub_diff, 0);
	mg_gchain_flt_sub(opt->pri_ratio, gi->k * 2, opt->best_n, gcs[0]->n_gc, gcs[0]->gc);
	mg_gchain_drop_flt(b->km, gcs[0]);
	mg_gchain_set_mapq(b->km, gcs[0], qlen_sum, mv.n, opt->min_gc_score);
	if ((opt->flag&MG_M_CIGAR) && n_segs == 1) {
		mg_gchain_cigar(b->km, gi->g, gi->es, seq_cat, gcs[0], qname);
		mg_gchain_gen_ds(b->km, gi->g, gi->es, seq_cat, gcs[0]);
	}
	kfree(b->km, seq_cat);
	if (mg_dbg_flag & MG_DBG_QNAME) t = print_time(t, 5, qname);

	if (b->km) {
		km_stat(b->km, &kmst);
		if (mg_dbg_flag & MG_DBG_QNAME)
			fprintf(stderr, "QM\t%s\t%d\tcap=%ld,nCore=%ld,largest=%ld\n", qname, qlen_sum, kmst.capacity, kmst.n_cores, kmst.largest);
		if (kmst.n_blocks != kmst.n_cores) {
			fprintf(stderr, "[E::%s] memory leak at %s\n", __func__, qname);
			abort();
		}
		if (kmst.largest > 1U<<28 || (opt->cap_kalloc > 0 && kmst.capacity > (size_t)(opt->cap_kalloc))) {
			km_destroy(b->km);
			b->km = km_init();
		}
	}
}

mg_gchains_t *mg_map(const mg_idx_t *gi, int qlen, const char *seq, mg_tbuf_t *b, const mg_mapopt_t *opt, const char *qname)
{
	mg_gchains_t *gcs;
	mg_map_frag(gi, 1, &qlen, &seq, &gcs, b, opt, qname, NULL);
	return gcs;
}