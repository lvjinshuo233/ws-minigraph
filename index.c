#include <assert.h>
#include "mgpriv.h"
#include "khashl.h"
#include "kthread.h"
#include "kvec-km.h"
#include "sys.h"
#include "minigraph.h"
#include <fstream>
#include <iostream>
#include <inttypes.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <collect-I.h>
#include <ws-strobe.h>
#define idx_hash(a) ((a)>>1)
#define idx_eq(a, b) ((a)>>1 == (b)>>1)
KHASHL_MAP_INIT(KH_LOCAL, idxhash_t, mg_hidx, uint64_t, uint64_t, idx_hash, idx_eq)
extern __thread const char *g_qname;

typedef struct mg_idx_bucket_s {
	mg128_v a;   // (minimizer, position) array
	int32_t n;   // size of the _p_ array
	uint64_t *p; // position array for minimizers appearing >1 times
	void *h;     // hash table indexing _p_ and minimizers appearing once
} mg_idx_bucket_t;

mg_idx_t *mg_idx_init(int k, int w, int b, int s, int t, int syncmer)
{
	mg_idx_t *gi;
	if (k*2 < b) b = k * 2;
	if (w < 1) w = 1;
	KCALLOC(0, gi, 1);
	gi->w = w, gi->k = k, gi->b = b, gi->s = s, gi->t = t, gi->syncmer = syncmer;
	KCALLOC(0, gi->B, 1<<b);
	return gi;
}

void mg_idx_destroy(mg_idx_t *gi)
{
	uint32_t i;
	if (gi == 0) return;
	if (gi->B) {
		for (i = 0; i < 1U<<gi->b; ++i) {
			free(gi->B[i].p);
			free(gi->B[i].a.a);
			mg_hidx_destroy((idxhash_t*)gi->B[i].h);
		}
		free(gi->B);
	}
	if (gi->strobe_idx) {
    	ws_strobe_destroy((strobe_hash_t*)gi->strobe_idx);
    	gi->strobe_idx = NULL;
	}
	if (gi->all_syncs.a) {
		free(gi->all_syncs.a);
		gi->all_syncs.a = NULL;
		gi->all_syncs.n = gi->all_syncs.m = 0;
	}
	gfa_edseq_destroy(gi->n_seg, gi->es);
	free(gi);
}

/****************
 * Index access *
 ****************/

const uint64_t *mg_idx_hget(const void *h_, const uint64_t *q, int suflen, uint64_t minier, int *n)
{
	khint_t k;
	const idxhash_t *h = (const idxhash_t*)h_;
	*n = 0;
	if (h == 0) return 0;
	k = mg_hidx_get(h, minier>>suflen<<1);
	if (k == kh_end(h)) return 0;
	if (kh_key(h, k)&1) { // special casing when there is only one k-mer
		*n = 1;
		return &kh_val(h, k);
	} else {
		*n = (uint32_t)kh_val(h, k);
		return &q[kh_val(h, k)>>32];
	}
}

const uint64_t *mg_idx_get(const mg_idx_t *gi, uint64_t minier, int *n)
{
	int mask = (1<<gi->b) - 1;
	mg_idx_bucket_t *b = &gi->B[minier&mask];
	return mg_idx_hget(b->h, b->p, gi->b, minier, n);
}

void mg_idx_cal_quantile(const mg_idx_t *gi, int32_t m, float f[], int32_t q[])
{
	int32_t i;
	uint64_t n = 0;
	khint_t *a, k;
	for (i = 0; i < 1<<gi->b; ++i)
		if (gi->B[i].h) n += kh_size((idxhash_t*)gi->B[i].h);
	a = (uint32_t*)malloc(n * 4);
	for (i = 0, n = 0; i < 1<<gi->b; ++i) {
		idxhash_t *h = (idxhash_t*)gi->B[i].h;
		if (h == 0) continue;
		for (k = 0; k < kh_end(h); ++k) {
			if (!kh_exist(h, k)) continue;
			a[n++] = kh_key(h, k)&1? 1 : (uint32_t)kh_val(h, k);
		}
	}
	for (i = 0; i < m; ++i)
		q[i] = ks_ksmall_uint32_t(n, a, (size_t)((1.0 - (double)f[i]) * n));
	free(a);
}

/***************
 * Index build *
 ***************/

static void mg_idx_add(mg_idx_t *gi, int n, const mg128_t *a)
{
	int i, mask = (1<<gi->b) - 1;
	for (i = 0; i < n; ++i) {
		mg128_v *p = &gi->B[a[i].x>>8&mask].a;
		kv_push(mg128_t, 0, *p, a[i]);
	}
}

void mg_idx_hfree(void *h_)
{
	idxhash_t *h = (idxhash_t*)h_;
	if (h == 0) return;
	mg_hidx_destroy(h);
}

void *mg_idx_a2h(void *km, int32_t n_a, mg128_t *a, int suflen, uint64_t **q_, int32_t *n_)
{
	int32_t N, n, n_keys;
	int32_t j, start_a, start_q;
	idxhash_t *h;
	uint64_t *q;

	*q_ = 0, *n_ = 0;
	if (n_a == 0) return 0;

	// sort by minimizer
	radix_sort_128x(a, a + n_a);

	// count and preallocate
	for (j = 1, n = 1, n_keys = 0, N = 0; j <= n_a; ++j) {
		if (j == n_a || a[j].x>>8 != a[j-1].x>>8) {
			++n_keys;
			if (n > 1) N += n;
			n = 1;
		} else ++n;
	}
	h = mg_hidx_init2(km);
	mg_hidx_resize(h, n_keys);
	KCALLOC(km, q, N);
	*q_ = q, *n_ = N;

	// create the hash table
	for (j = 1, n = 1, start_a = start_q = 0; j <= n_a; ++j) {
		if (j == n_a || a[j].x>>8 != a[j-1].x>>8) {
			khint_t itr;
			int absent;
			mg128_t *p = &a[j-1];
			itr = mg_hidx_put(h, p->x>>8>>suflen<<1, &absent);
			assert(absent && j == start_a + n);
			if (n == 1) {
				kh_key(h, itr) |= 1;
				kh_val(h, itr) = p->y;
			} else {
				int k;
				for (k = 0; k < n; ++k)
					q[start_q + k] = a[start_a + k].y;
				radix_sort_gfa64(&q[start_q], &q[start_q + n]); // sort by position; needed as in-place radix_sort_128x() is not stable
				kh_val(h, itr) = (uint64_t)start_q<<32 | n;
				start_q += n;
			}
			start_a = j, n = 1;
		} else ++n;
	}
	assert(N == start_q);
	return h;
}

static void worker_post(void *g, long i, int tid)
{
	mg_idx_t *gi = (mg_idx_t*)g;
	mg_idx_bucket_t *b = &gi->B[i];
	if (b->a.n == 0) return;
	b->h = (idxhash_t*)mg_idx_a2h(0, b->a.n, b->a.a, gi->b, &b->p, &b->n);
	kfree(0, b->a.a);
	b->a.n = b->a.m = 0, b->a.a = 0;
}

int mg_gfa_overlap(const gfa_t *g)
{
	uint64_t i;
	for (i = 0; i < g->n_arc; ++i) // non-zero overlap
		if (g->arc[i].ov != 0 || g->arc[i].ow != 0)
			return 1;
	return 0;
}

uint64_t encodeKmer(const std::string &str)
{
	uint64_t kmer[2] = {0,0};
	int k = str.length();
	uint64_t shift1 = 2 * (k - 1);

	for (int i = 0; i < k; ++i) 
	{
		int c = seq_nt4_table[(uint8_t)str[i]];
		kmer[0] = kmer[0] << 2 | c;  
		kmer[1] = (kmer[1] >> 2) | (3ULL^c) << shift1;
	}

	return kmer[0] < kmer[1]? kmer[0] : kmer[1];
}

/**
 * kmer_freq_filename参数传递
 */
mg_idx_t *mg_index_core(gfa_t *g, int k, int w, int s, int t, int syncmer, int b, 
	int n_threads, const char* kmer_freq_filename, int flag, collect_rep_t *collect_rep, int min_rep_len)
{
	mg_idx_t *gi;
	mg128_v a = {0,0,0};
	mg128_v all_syncs = {0,0,0};
	int i;
	const char *out_dir = "/home/lvjinshuo/BIBM_paper/test/find_problem_centro/sketchs";   // 你的绝对路径
	mkdir(out_dir, 0755);  // 若目录已存在则忽略错误
	if (mg_gfa_overlap(g)) {
		if (mg_verbose >= 1)
			fprintf(stderr, "[E::%s] minigraph doesn't work with graphs containing overlapping segments\n", __func__);
		return 0;
	}
	gi = mg_idx_init(k, w, b, s, t, syncmer);
	gi->g = g;
	//Syncmer parameter error
    assert(syncmer);

	//Read kmer frequency file
	/*------------------------*/
	std::ifstream idt (kmer_freq_filename);

	fprintf(stderr, "[M::%s::%.3f*%.2f] reading downweighted kmers\n", __func__, realtime() - mg_realtime0, cputime() / (realtime() - mg_realtime0));
	std::string kmer;
	uint64_t cnt = 0;
	uint64_t freq;
	while(idt >> kmer >> freq)
		cnt++;

	//kmer length used for kmer counting and mapping must be consistent
	if(cnt > 0)
	{
		if(kmer.length() != (uint64_t)s)
		{
			fprintf(stderr, "ERROR: input list of k-mers and winnowmap parameter k are inconsistent : kmer = %ld, smer = %d\n",kmer.length(), s);
			abort();
		}
	}

	//set up bloom filter
	bloom_parameters parameters;
	parameters.projected_element_count = std::max(cnt, (uint64_t)1000);
	parameters.false_positive_probability = 0.001; 
	parameters.maximum_number_of_hashes = 2;

	if (!parameters)
	{
		std::cout << "Error - Invalid set of bloom filter parameters!" << std::endl;
		exit(1);
	}

	parameters.compute_optimal_parameters();
	gi->downFilter = new bloom_filter(parameters);

	//read the file again
	idt.clear();
	idt.seekg(0);
	while(idt >> kmer >> freq)
	{
		gi->downFilter->insert(encodeKmer(kmer));
	}
	
	assert(cnt == gi->downFilter->element_count());

	fprintf(stderr, "[M::%s::%.3f*%.2f] collected downweighted kmers, no. of kmers read=%" PRIu64"\n", __func__, realtime() - mg_realtime0, cputime() / (realtime() - mg_realtime0), cnt);
	fprintf(stderr, "[M::%s::%.3f*%.2f] saved the kmers in a bloom filter: hash functions=%u and size=%llu \n", __func__, realtime() - mg_realtime0, cputime() / (realtime() - mg_realtime0), parameters.optimal_parameters.number_of_hashes, parameters.optimal_parameters.table_size);
	/*------------------------*/
	for (i = 0; (uint32_t)i < g->n_seg; ++i) {
		gfa_seg_t *seg = &g->seg[i];
		seg_rep_t *seg_rep_info = NULL;
		if(collect_rep){
			seg_rep_info = &collect_rep->segs[i];
		}
		//fprintf(stderr, "[segs] name %s snid %d\n", g->seg[i].name, g->seg[i].snid);

		a.n = 0;
		//mg_sketch(0, s->seq, s->len, w, k, i, &a); // TODO: this can be parallelized
		open_sync_sketch(0, seg->seq, seg->len, k, s, t, gi->flag&MM_I_HPC,i, &a, gi, flag, seg_rep_info, min_rep_len);
		//SW.n_sketch++;
		// ---- 修改点：仅当构建索引时输出采样位置 ----
        if (flag) {
            char path[2048];
            snprintf(path, sizeof(path),
                     "/home/lvjinshuo/BIBM_paper/scripts/new_centromere_fraction/ref_index/ws_m/ws_%s_sample.txt",
                     seg->name);
            FILE *fp = fopen(path, "w");
            if (fp) {
                for (size_t j = 0; j < a.n; ++j) {
                    uint32_t end_pos = (uint32_t)(a.a[j].y) >> 1;   // 0‑based 结束位置
                    int start_pos = end_pos - k + 1;                // 0‑based 起始位置
                    if (start_pos >= 0)
                        fprintf(fp, "%d\n", start_pos);
                }
                fclose(fp);
            } else {
                fprintf(stderr, "[WARNING] cannot write %s\n", path);
            }
        }
        // --------------------------------------------
		// 在遍历 segment 的循环内，open_sync_sketch 之后，添加写入文件的代码：
/*if (1) {  // 始终输出，或可加条件判断
    char path[2048];
    snprintf(path, sizeof(path), "%s/%s_ws_sketch.txt", out_dir, seg->name);
    FILE *fp = fopen(path, "w");
    if (fp) {
        for (size_t j = 0; j < a.n; ++j) {
            uint32_t pos = (uint32_t)(a.a[j].y) >> 1;  // 0‑based 结束位置
            fprintf(fp, "%u\n", pos);
        }
        fclose(fp);
    } else {
        fprintf(stderr, "[WARNING] cannot write %s\n", path);
    }
}*/
		mg_idx_add(gi, a.n, a.a);	
		if (all_syncs.n + a.n > all_syncs.m && (flag==1||flag==3)) {
    		all_syncs.m = all_syncs.m ? all_syncs.m : 1024;
    		while (all_syncs.m < all_syncs.n + a.n) all_syncs.m <<= 1;
    		all_syncs.a = (mg128_t*)realloc(all_syncs.a, all_syncs.m * sizeof(mg128_t));
		}
		//if(flag) n_syncmers += a.n;
		
		if(flag==1||flag==3) memcpy(&all_syncs.a[all_syncs.n], a.a, a.n * sizeof(mg128_t));
		if(flag==1||flag==3) all_syncs.n += a.n;
	}
	/*if(flag){
		bases++;
		SW.seq[bases] = '\0';
        fprintf(f_check, "%ld\t%s\n", window, SW.seq);
	}*/

	free(a.a);
	kt_for(n_threads, worker_post, gi, 1<<gi->b);
	gi->all_syncs = all_syncs;
	return gi;
}

static int cmp_int32_asc(const void *a, const void *b) {
    int32_t x = *(const int32_t*)a;
    int32_t y = *(const int32_t*)b;
    return (x > y) - (x < y);
}

/**
 * 加了kmer_freq_filename参数
 */
mg_idx_t *mg_index(gfa_t *g, const mg_idxopt_t *io, int n_threads, mg_mapopt_t *mo, const char* kmer_freq_filename, int flag, collect_rep_t *collect_rep)
{
	int32_t i, j;
	mg_idx_t *gi;
	for (i = 0; (uint32_t)i < (g->n_seg); ++i) { // uppercase
		gfa_seg_t *s = &g->seg[i];
		for (j = 0; j < s->len; ++j)
			if (s->seq[j] >= 'a' && s->seq[j] <= 'z')
				s->seq[j] -= 32;
	}
	/**
	 * kmer_freq_filename参数传递
	 */
	
	/*char path[1024];
        snprintf(path, sizeof(path),
                 "/home/lvjinshuo/BIBM_paper/scripts/new_centromere_fraction/reads_index/ws_m/ws_%s_sample.txt",
                 g_qname);
        FILE *fp = fopen(path, "w");
		if (fp) {
            fprintf(fp, "");
            fclose(fp);
        }*/
	gi = mg_index_core(g, io->k, io->w, io->s, io->t, io->syncmer, io->bucket_bits, n_threads, kmer_freq_filename, flag, collect_rep, io->min_rep_len);
	if (gi == 0) return 0;
	gi->es = gfa_edseq_init(gi->g);
	gi->n_seg = g->n_seg;
	gi->min_rep_len = io->min_rep_len;
	gi->max_rep_len = io->max_rep_len;

	if (mo) mg_opt_update(gi, mo, 0);
	/* ============ 输出高频区间长度分布 ============ */
/*{
    size_t total_itvs = 0;
    for (int32_t ii = 0; (uint32_t)ii < collect_rep->len; ++ii) {
        const seg_rep_t *seg = &collect_rep->segs[ii];
        for (const rep_interval_t *itv = seg->high_first; itv; itv = itv->next)
            total_itvs++;
    }

    if (total_itvs > 0) {
        int32_t *lengths = (int32_t*)malloc(total_itvs * sizeof(int32_t));
        if (lengths) {
            size_t k = 0;
            for (int32_t ii = 0; (uint32_t)ii < collect_rep->len; ++ii) {
                const seg_rep_t *seg = &collect_rep->segs[ii];
                for (const rep_interval_t *itv = seg->high_first; itv; itv = itv->next)
                    lengths[k++] = itv->end - itv->start;
            }

            qsort(lengths, total_itvs, sizeof(int32_t), cmp_int32_asc);

            FILE *fp = fopen("high_freq_interval.txt", "w");
            if (fp) {
                size_t j = 0;
                while (j < total_itvs) {
                    size_t j2 = j;
                    while (j2 < total_itvs && lengths[j2] == lengths[j]) j2++;
                    fprintf(fp, "%d\t%zu\n", lengths[j], j2 - j);
                    j = j2;
                }
                fclose(fp);

                fprintf(stderr, "[M::%s] wrote high_freq_interval.txt: %zu intervals, %zu distinct lengths\n",
                        __func__, total_itvs, (size_t)0);   /* distinct 统计下面可选 */
            /*} else {
                fprintf(stderr, "[WARNING] cannot open high_freq_interval.txt for writing\n");
            }
            free(lengths);
        } else {
            fprintf(stderr, "[WARNING] malloc failed for interval lengths\n");
        }
    } else {
        fprintf(stderr, "[M::%s] no high-frequency intervals collected, skip writing high_freq_interval.txt\n", __func__);
    }
}
	/**
	 * 修剪高频k-mer区间，过滤太短的区间，最小阈值偏移量是REP_MIN_OFFSET
	 */
	/*if(!fix_rep(collect_rep, io->s, io->min_rep_len)){
		fprintf(stderr, "ERROR: failed to fix repetitive intervals\n");
		assert(0);
	}*/
	/*if(!show_rep(collect_rep)){
		fprintf(stderr, "ERROR: failed to show repetitive intervals\n");
		assert(0);
	}*/
	/**
	 * 构建基于加权syncmer索引gi的strobemer索引：si
	 */

	//gi->strobe_idx = ws_strobe_build(gi, collect_rep, &gi->all_syncs);
	gi->strobe_idx = build_ref_interval_index(gi, collect_rep);
	free_collect_rep(collect_rep);   /* 释放链表和 segs 数组 */
	/*free(gi->all_syncs.a);
	gi->all_syncs.a = NULL;
	gi->all_syncs.n = gi->all_syncs.m = 0;*/
	if (!gi->strobe_idx) {
    	fprintf(stderr, "WARNING: failed to build strobemer index\n");
	}
	//fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec; Peak RSS: %.3f GB\n", __func__, realtime() - mg_realtime0, cputime(), peakrss() / 1024.0 / 1024.0 / 1024.0);
	//assert(0);
	if (mg_verbose >= 3)
		fprintf(stderr, "[M::%s::%.3f*%.2f] indexed the graph\n", __func__,
				realtime() - mg_realtime0, cputime() / (realtime() - mg_realtime0));
	//assert(0);
	return gi;
}
