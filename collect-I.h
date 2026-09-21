#ifndef COLLECTI_H
#define COLLECTI_H

#include <stdint.h>
#define REP_MIN_OFFSET 20


#ifdef __cplusplus
extern "C" {
#endif

typedef struct rep_interval_t{
    int start;  // Start position of the repetitive interval
    int end;    // 实际上是最后一个高频k-mer的起始位置
    struct rep_interval_t *next; // Pointer to the next interval in the linked list

}rep_interval_t;

typedef struct {
    rep_interval_t *high_reps_list;//高频k-mer区间信息头指针。空指针
    rep_interval_t *low_reps_list;//低频的头指针。空指针
    int high_count;
    int low_count;
    rep_interval_t *high_first;//高频第一个区间
    rep_interval_t *high_last;//高频最后一个区间
    rep_interval_t *low_first;//低频第一个区间
    rep_interval_t *low_last;//低频最后一个区间
    char *sname; // Segment name or reads name

} seg_rep_t;//这个也用做reads的统计，一个reads相当于一个segment

typedef struct {
    seg_rep_t *segs;  
    int len;//seg数量

    int rep_num;//区间数量
    double cov_length;//区间平均长度
    int rep_free_num;//过滤掉的短区间
} collect_rep_t;


int init_collect_rep(uint32_t n_seg, collect_rep_t *collect_rep);
int init_seg_rep(seg_rep_t *seg_rep, char *sname);
int add_rep(seg_rep_t *seg_rep, int start, int end);
int show_rep(collect_rep_t *collect_rep);
int show_seg_rep(seg_rep_t *seg_rep);
int fix_rep(collect_rep_t *collect_rep, int s, int min_offset);
int fix_seg_rep(seg_rep_t *seg_rep, int s, int min_offset);
void free_collect_rep(collect_rep_t *collect_rep);
#ifdef __cplusplus
}
#endif

#endif