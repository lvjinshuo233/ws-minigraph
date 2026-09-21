#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#define __STDC_LIMIT_MACROS
#include "kvec-km.h"
#include "mgpriv.h"
#include "minigraph.h"
#include "collect-I.h"
extern __thread const char *g_qname;
typedef struct { // a simplified version of kdq
	int front, count;
	int a[32];
} tiny_queue_t;

unsigned char seq_nt4_table[256] = {
	0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

uint64_t murmur64(uint64_t h)
{
	h+=1;
	h ^= (h >> 33);
	h *= 0xff51afd7ed558ccdL;
	h ^= (h >> 33);
	h *= 0xc4ceb9fe1a85ec53L;
	h ^= (h >> 33);
	return h;
}

static inline uint64_t murmerhash64(uint64_t key, uint64_t mask)
{
	key+=1;
	key ^= key >> 33;
	key *= 0xff51afd7ed558ccd;
	key ^= key >> 33;
	key *= 0xc4ceb9fe1a85ec53;
	key ^= key >> 33;
	return key & mask;
}


static inline uint64_t hash64(uint64_t key, uint64_t mask)
{
	key = (~key + (key << 21)) & mask; // key = (key << 21) - key - 1;
	key = key ^ key >> 24;
	key = ((key + (key << 3)) + (key << 8)) & mask; // key * 265
	key = key ^ key >> 14;
	key = ((key + (key << 2)) + (key << 4)) & mask; // key * 21
	key = key ^ key >> 28;
	key = (key + (key << 31)) & mask;
	return key;
}

/*
 * Inserts a k-mer representing the window around a syncmer 
 * @param t     1-indexed position offset t for syncmer
 * 
 *
 *
 */
static inline void insert_syncmer(int len, const char *str, int t, uint32_t rid, int i, int s, int k, mg128_v *p, void *km, int strand, int flag){
    int start_pos = i - k + 1;
    int end_pos = start_pos + k;
    //char *temp = (char*)malloc(sizeof(char)*(end_pos-start_pos+1));
    if (end_pos > len){
        fprintf(stderr,"end pos, len, start_pos,i = %d,%d,%d,%d",end_pos,len,start_pos,i);
        //assert(end_pos < len);
        return;
    }
    if (start_pos < 0){
        fprintf(stderr,"end pos, len, start_pos, i = %d,%d,%d,%d",end_pos,len,start_pos,i);
        assert(start_pos > -1);
    }
    mg128_t info = { UINT64_MAX, UINT64_MAX };
    int is_ambiguous = 0;
    //int mark = 0;
    uint64_t shift1 = 2 * (k - 1), mask = (1ULL<<2*k) - 1, kmer[2] = {0,0};
    for (int j = start_pos; j < end_pos; ++j){
		int c = seq_nt4_table[(uint8_t)str[j]];
        //temp[mark] = str[j];
        //mark++;
        if (c >= 4){
            is_ambiguous = 1;
            break;
        }
        kmer[0] = (kmer[0] << 2 | c) & mask;           // forward k-mer
        kmer[1] = (kmer[1] >> 2) | (3ULL^c) << shift1; // reverse k-mer
    }
    //temp[mark] = '\0';
    if (kmer[0] == kmer[1]) return; // skip "symmetric k-mers" as we don't know it strand
    if (is_ambiguous) return;
    //int z = kmer[0] < kmer[1]? 0 : 1; // strand
    info.x = hash64(kmer[strand], mask) << 8 | k;
    info.y = (uint64_t)rid<<32 | (uint32_t)(end_pos-1)<<1 | strand;
    kv_push(mg128_t, km, *p, info);
    /*if (flag) {
        SW.n_kmers++;
        //fprintf(f_samplingInfo, "%s\t%ld\t%d\n", temp, sync_distance, SW.n_windows);
        if(flag==1){
            n_syncmers++;
        }
    }*/
    //free(temp);
}

/**
 * @brief		takes hash value of kmer and adjusts it based on kmer's weight
 *					this value will determine its order for minimizer selection
 * @details	this is inspired from Chum et al.'s min-Hash and tf-idf weighting 
 */
static inline double applyWeight(uint64_t kmer, const mg_idx_t *mi, int flag, int* is_high_freq, int* is_weighted_max, int base_len,
     seg_rep_t *seg_rep_info, int s, int min_rep_len)
{

	uint64_t hash = murmerhash64(kmer, UINT64_MAX);//hash64() 
	double x = hash * 1.0 / UINT64_MAX;  //bring it within [0, 1]
	//assert (x >= 0.0 && x <= 1.0);

	if (mi->downFilter->contains(kmer))
	{
        if((*is_high_freq)==0) *is_high_freq = 1;
        *is_weighted_max += 1;
        double p2 = x*x;
		double p4 = p2 * p2;
        /*if(flag==1&&(window==57||window==58||window==59||window==60||window==61||window==62||window==63)){
            fprintf(f_check, "%ld\t%f\t%f\t<---\n", window, x, -1.0 * (p4 * p4));
            n_weights_centro++;
        }else if(flag==1){
            //fprintf(f_check, "%ld\t%f\t%f", window, x, -1.0 * (p4 * p4));
            n_weights++;
        }*/
        //if(flag) n_w_k++;
		/* downweigting by a factor of 8 */
		/* further aggressive downweigting may affect accuracy */
		/* TODO: Consider making it a user parameter */
		return -1.0 * (p4 * p4);
	}
    /**
     * 收集高频k-mer区间信息
     * @param seg_rep_info  segment的重复区间信息
     * @param is_high_freq  上一个s-mer是否高频
     * @param is_weighted_max  上一个高频k-mer区域的长度
     * @param base_len  当前读取碱基的长度
     * @param s  s-mer的长度
     * 区间的起始位置为：base_len - is_weighted_max - s + 1
     * 区间的结束位置为：base_len
     */
    if((flag==1||flag==3)&&*is_high_freq==1&&*is_weighted_max>=min_rep_len){
        //if(flag==1) fprintf(f_check, "上方高频k-mer区域共长：%d\t区间：[%d, %d)\t1-based index\n", *is_weighted_max, base_len - *is_weighted_max - s + 1, base_len);
        if(!add_rep(seg_rep_info, base_len - *is_weighted_max - s + 1, base_len - s + 1)){
            fprintf(stderr, "Failed to add repetitive interval\n");
            assert(0);
        }
        *is_high_freq = 0;
        *is_weighted_max = 0;
    }
    /*if(flag==1&&(window==57||window==58||window==59||window==60||window==61||window==62||window==63)){
        fprintf(f_check, "%ld\t%f\t%f\n", window, x, -1.0 * x);
    }*/
    
	return -1.0 * x;

	//range of returned value is between [-1,0]
	//we avoid adding one for better double precision 
}

static inline int tq_shift(tiny_queue_t *q)
{
	int x;
	if (q->count == 0) return -1;
	x = q->a[q->front++];
	q->front &= 0x1f;
	--q->count;
	return x;
}

static inline void tq_push(tiny_queue_t *q, int x)
{
	q->a[((q->count++) + q->front) & 0x1f] = x;
}


/**
 * Find symmetric (k,s,t)-open syncmers on a DNA sequence. Note: only takes syncmer
 * if the last minimum in the window is at position t. can modify such that
 * takes syncmer if any minimium is at position t. 
 *
 * @param km     thread-local memory pool; using NULL falls back to malloc()
 * @param str    DNA sequence
 * @param len    length of $str
 * @param w      find a minimizer for every $w consecutive k-mers
 * @param k      k-mer size
 * @param s      s-mer size
 * @param t      t offset
 * @param rid    reference ID; will be copied to the output $p array
 * @param is_hpc homopolymer-compressed or not
 * @param p      minimizers
 *               p->a[i].x = kMer<<8 | kmerSpan
 *               p->a[i].y = rid<<32 | lastPos<<1 | strand
 *               where lastPos is the position of the last base of the i-th minimizer,
 *               and strand indicates whether the minimizer comes from the top or the bottom strand.
 *               Callers may want to set "p->n = 0"; otherwise results are appended to p
 */
void open_sync_sketch(void *km, const char *str, int len, int k, int s, int t, int is_hpc, uint32_t rid, mg128_v *p, const mg_idx_t *gi, 
    int flag, seg_rep_t *seg_rep_info, int min_rep_len)
{
    int is_high_freq = 0;//判断上一个s-mer是否高频。
    int is_weighted_max = 0;
    int base_len = 0;

    int used_t = 0;
    if (t == 0) {
        used_t = (k - s)/2 + 1;
    }
    else{
        used_t = t;
    }
    uint64_t shift1 = 2 * (s - 1), mask = (1ULL<<2*s) - 1, smer[2] = {0,0};
    uint64_t kmer[2] = {0,0}, kmer_mask = (1ULL<<2*k) - 1, shift_kmer = 2*(k-1);
	int i, j, l, buf_pos, min_pos, smer_span = 0;
	mg128_t buf[256], min = { UINT64_MAX, UINT64_MAX };
    int w = k - s + 1;
    double buf_order[256], min_order = 2.0; // 2.0 value is indicating uninitialized
    tiny_queue_t tq;
	assert(len > 0 && (w > 0 && w < 256) && (k > 0 && k <= 28)); // 56 bits for k-mer; could use long k-mers, but 28 enough in practice
    assert(k > s);
    memset(buf, 0xff, w * 16);
    for (i = 0; i < w; i++) buf_order[i] = 2.0;
    memset(&tq, 0, sizeof(tiny_queue_t));
	kv_resize(mg128_t, km, *p, p->n + len/w);
    //窗口工具↓
    /*char *temp = (char*)malloc(sizeof(char) * (k + 1));
    int mark = 0;*/

    //Take k-mer from pos i
    //look into how buffer_pos interacts with t. Itl ooks like a circular buffer sort of thing for a linear window
    //try to understand what /smer_span l is doing. i think l keeps track of the k-mer that's building, and when it's
    //>= s then we hash the k-mer. After it's greater than s, I think it just becomes some variable hwich keeps increasing
    //which becomes 0 afer a bad character, signaling the start of the process again. 
	for (i = l = buf_pos = min_pos = 0; i < len; ++i) {
		int c = seq_nt4_table[(uint8_t)str[i]];
        //if(flag) sync_distance++;//每读取一个碱基，则相邻syncmer距离加一
        if(flag==1||flag==3) base_len++;//计算区间用的碱基长度
        /*if(flag==1&&bases<SUPER_WINDOW_SIZE){
            //SW.seq[bases] = str[i];
            bases++;
            if(bases == SUPER_WINDOW_SIZE){
                //SW.seq[bases] = '\0';
                bases = 0;
                //fprintf(f_check, "%ld\t%s\n", window, SW.seq);
                window++;
            }
        }*/
        //截取k-mer
        /*temp[mark] = str[i];
        mark++;
        if(mark == k){
            temp[mark] = '\0';
            mark = 0;
        }*/
        /*if(flag==1&&SW.n_bases<=SUPER_WINDOW_SIZE) SW.n_bases++;
        if(flag==1&&SW.n_bases>SUPER_WINDOW_SIZE){
            //fprintf(f_freq, "%d\t%ld\n", SW.n_windows, SW.n_kmers);// window : 
            SW.n_bases=1;
            SW.n_kmers=0;
            SW.n_windows++;
        }*/
		mg128_t info = { UINT64_MAX, UINT64_MAX };
        double info_order = 2.0; // 2.0 value is indicating uninitialized
        int strand = 0;
        /**
         * 将字符i转换为无符号8位整数，然后再以这数值为索引，查找seq_nt4_table中的值
         * c = 0,1,2,3,4
         */
		if (c < 4) { // not an ambiguous base

             /**
             * z的值是smer的正反向
             * strand的值是kmer的正反向
             */
			int z;
            if (is_hpc) {
                int skip_len = 1;
                if (i + 1 < len && seq_nt4_table[(uint8_t)str[i + 1]] == c) {
                    for (skip_len = 2; i + skip_len < len; ++skip_len)
                        if (seq_nt4_table[(uint8_t)str[i + skip_len]] != c)
                            break;
                    i += skip_len - 1; // put $i at the end of the current homopolymer run
                }
                tq_push(&tq, skip_len);
                smer_span += skip_len;
                if (tq.count > s) smer_span -= tq_shift(&tq);
            }
            /**
             * l 是一个计数器，用于跟踪当前构建的 s-mer 的长度。
             * 确保smer_span不超过s的值
             */
            else smer_span = l + 1 < s? l + 1 : s;
			smer[0] = (smer[0] << 2 | c) & mask;           // forward s-mer
			kmer[0] = (kmer[0] << 2 | c) & kmer_mask;           // forward k-mer
			smer[1] = (smer[1] >> 2) | (3ULL^c) << shift1; // reverse s-mer
			kmer[1] = (kmer[1] >> 2) | (3ULL^c) << shift_kmer; // reverse k-mer
			if (smer[0] == smer[1]) continue; // skip "symmetric s-mers" as we don't know it strand
			if (kmer[0] == kmer[1]) continue; // skip "symmetric k-mers" as we don't know it strand
			z = smer[0] < smer[1]? 0 : 1; // strand
            strand = kmer[0] < kmer[1]? 0 : 1;
			++l;
			if (l >= s && smer_span < 256) {
				info.x = hash64(smer[z], mask) << 8 | smer_span;
				info.y = (uint64_t)rid<<32 | (uint32_t)i<<1 | z;
                info_order = applyWeight(smer[z], gi, flag, &is_high_freq, &is_weighted_max, base_len, seg_rep_info, s, min_rep_len);
			}
		} else{
            l = 0;
            smer_span = 0;
            tq.count = tq.front = 0;
            is_high_freq = 0;
        } 
        //buf_pos is the last s-mer added. buf_pos+1 is the first k-mer in the window. 
		buf[buf_pos] = info; // need to do this here as appropriate buf_pos and buf[buf_pos] are needed below
        //fprintf(stderr,"info %d,%d, s = %d, t = %d, smer = %s \n",info.x,info.y,s,t);
        buf_order[buf_pos] = info_order;
        int offset = (strand == 0)? used_t : w - used_t + 1;
        //int offset = t;
		if (info_order <= min_order) { // a new minimum; then write the old min
            min = info,min_pos = buf_pos;
            min_order = info_order;
            if (offset == w && l >= w + s - 1 && min.x != UINT64_MAX) {
                insert_syncmer(len, str, offset, rid, i, s, k, p, km, strand, flag);
                /*if(flag){
                    fprintf(f_freq, "%s\n", temp);
                }*/
               /*if(flag==1){
                    fprintf(f_seq,"%ld\t%ld\t%d\n", dist_count, sync_distance, SW.n_windows);
                    sync_distance = 0;
                    dist_count++;
               }*/
              if (flag == 2) {
    int start_pos = i - k + 1;
    /*if (start_pos >= 0 && g_qname) {
        char path[1024];
        snprintf(path, sizeof(path),
                 "/home/lvjinshuo/BIBM_paper/scripts/new_centromere_fraction/reads_index/ws_m/ws_%s_sample.txt",
                 g_qname);
        FILE *fp = fopen(path, "a");
        if (fp) {
            fprintf(fp, "%d\n", start_pos);
            fclose(fp);
        }
    }*/
}
            }
        } else if (buf_pos == min_pos) { // old min has moved outside the window
			for (j = buf_pos + 1, min.x = UINT64_MAX, min_order = 2.0; j < w; ++j) {
                if (min_order >= buf_order[j]) {
                    min = buf[j], min_pos = j;
                    min_order = buf_order[j];
                }
            }
            for (j = 0; j <= buf_pos; ++j) {
                if (min_order >= buf_order[j]) {
                    min = buf[j], min_pos = j;
                    min_order = buf_order[j];
                }
            }
            if (l >= w + s - 1 && min.x != UINT64_MAX) {
                if ((buf_pos + offset) % w == min_pos) {
                    insert_syncmer(len, str, offset, rid, i, s, k, p, km, strand, flag);
                    /*if(flag){
                        //fprintf(f_freq, "%s\n", temp);
                    }*/
                    /*if(flag==1){
                        fprintf(f_seq,"%ld\t%ld\t%d\n", dist_count, sync_distance, SW.n_windows);
                        sync_distance = 0;
                        dist_count++;
                    }*/
                   if (flag == 2) {
    int start_pos = i - k + 1;
    /*if (start_pos >= 0 && g_qname) {
        char path[1024];
        snprintf(path, sizeof(path),
                 "/home/lvjinshuo/BIBM_paper/scripts/new_centromere_fraction/reads_index/ws_m/ws_%s_sample.txt",
                 g_qname);
        FILE *fp = fopen(path, "a");
        if (fp) {
            fprintf(fp, "%d\n", start_pos);
            fclose(fp);
        }
    }*/
}
                }
            }
        } else if ((buf_pos + offset) % w == min_pos && l >= w + s - 1) {
            insert_syncmer(len, str, offset, rid, i, s, k, p, km, min.y % 2, flag);
            /*if(flag){
                //fprintf(f_freq, "%s\n", temp);
            }*/
            /*if(flag==1){
                fprintf(f_seq,"%ld\t%ld\t%d\n", dist_count, sync_distance, SW.n_windows);
                sync_distance = 0;
                dist_count++;
            }*/
           if (flag == 2) {
    int start_pos = i - k + 1;
    /*if (start_pos >= 0 && g_qname) {
        char path[1024];
        snprintf(path, sizeof(path),
                 "/home/lvjinshuo/BIBM_paper/scripts/new_centromere_fraction/reads_index/ws_m/ws_%s_sample.txt",
                 g_qname);
        FILE *fp = fopen(path, "a");
        if (fp) {
            fprintf(fp, "%d\n", start_pos);
            fclose(fp);
        }
    }*/
}
        }
		if (++buf_pos == w) buf_pos = 0;
	}
    if((flag==1||flag==3)&&seg_rep_info&&is_weighted_max>0){
        if(is_weighted_max>=min_rep_len){
            if(!add_rep(seg_rep_info, base_len - is_weighted_max - s + 1, base_len - s + 1)){
                fprintf(stderr, "Failed to add repetitive interval\n");
                assert(0);
            }
        }
        //if(flag==1) fprintf(f_check, "上方高频k-mer区域共长：%d\t区间：[%d, %d)\n************************以上为segment：%s\t共长：%d\n", is_weighted_max, base_len - is_weighted_max - s + 1, base_len - s + 1, seg_rep_info->sname, base_len);
    } 
    //if(flag==1&&seg_rep_info&&is_weighted_max==0) fprintf(f_check, "************************以上为segment：%s\t共长：%d\n",seg_rep_info->sname, base_len);
    /*if(flag==1&&(window==59||window==60||window==61)){
        fprintf(stderr, "centromere:seg:%s", seg_rep_info->sname);
    }*/
}
