#include<stdio.h>
#include<stdlib.h>
#include<collect-I.h>
#include<minigraph.h>
int init_collect_rep(uint32_t n_seg, collect_rep_t *collect_rep) {
    collect_rep->segs = (seg_rep_t *)malloc(n_seg * sizeof(seg_rep_t));
    collect_rep->len = n_seg;
    collect_rep->rep_num = 0;
    collect_rep->cov_length = 0.0;
    collect_rep->rep_free_num = 0;
    if (collect_rep->segs == NULL) {
        fprintf(stderr, "Memory allocation failed for segs_rep\n");
        return 0;
    }
    return 1;
}

int init_seg_rep(seg_rep_t *seg_rep, char *sname) {
    seg_rep->high_reps_list = (rep_interval_t *)malloc(sizeof(rep_interval_t));
    seg_rep->low_reps_list = (rep_interval_t *)malloc(sizeof(rep_interval_t));
    seg_rep->high_reps_list->next = NULL;
    seg_rep->low_reps_list->next = NULL;
    seg_rep->high_reps_list->start = 0;
    seg_rep->low_reps_list->start = 0;
    seg_rep->high_reps_list->end = 0;
    seg_rep->low_reps_list->end = 0;
    seg_rep->high_count = 0;
    seg_rep->low_count = 0;
    seg_rep->sname = sname; // Assign the segment name
    
    seg_rep->high_first = seg_rep->high_reps_list->next;
    seg_rep->high_last = seg_rep->high_reps_list->next;
    seg_rep->low_first = seg_rep->low_reps_list->next;
    seg_rep->low_last = seg_rep->low_reps_list->next;

    if(seg_rep->high_reps_list == NULL || seg_rep->low_reps_list == NULL) {
        fprintf(stderr, "Memory allocation failed for seg_rep_t\n");
        return 0;
    }
    return 1;
}

int add_rep(seg_rep_t *seg_rep, int start, int end) {
    if(seg_rep == NULL){
        return 1;
    }
    /*if(seg_rep->high_reps_list->next&&seg_rep->high_last->end >= start){
        seg_rep->high_last->end = end;
        return 1;
    }*/
    rep_interval_t *new_interval = (rep_interval_t *)malloc(sizeof(rep_interval_t));
    if (new_interval == NULL) {
        fprintf(stderr, "Memory allocation failed for new interval\n");
        return 0;
    }
    new_interval->start = start;
    new_interval->end = end;
    new_interval->next = NULL;

    if(seg_rep->high_reps_list->next){
        seg_rep->high_last->next = new_interval;
        seg_rep->high_last = new_interval;
    }else{
        seg_rep->high_reps_list->next = new_interval;
        seg_rep->high_first = new_interval;
        seg_rep->high_last = new_interval;
    }
    seg_rep->high_count++;
    return 1;
}

int show_rep(collect_rep_t *collect_rep){
    fprintf(f_rep, "******************Collective Repetitive Intervals : 1-based index\n");
    //fprintf(f_rep, "Legal intervals: %d\n", collect_rep->rep_num);
    //fprintf(f_rep, "Average length: %.2f\n", collect_rep->cov_length);
    //fprintf(f_rep, "Filtered intervals: %d\n", collect_rep->rep_free_num);
    for(int i = 0; i < collect_rep->len; i++){
        seg_rep_t *seg_rep = &collect_rep->segs[i];
        fprintf(f_rep, "******************Segment: %s\n", seg_rep->sname);
        rep_interval_t *interval = seg_rep->high_first;
        while(interval){
            fprintf(f_rep, "[%d, %d)\n", interval->start, interval->end);
            interval = interval->next;
        }
    }
    return 1;
}
//弃用。区间改为高频k-mer起始位置的区间，而不是所有序列覆盖区间。
int fix_rep(collect_rep_t *collect_rep, int s, int min_offset) {
    int total_len = 0;
    collect_rep->rep_num = 0;   // 重置以便重新计数
    collect_rep->rep_free_num = 0;
    for (int i = 0; i < collect_rep->len; i++) {
        seg_rep_t *seg = &collect_rep->segs[i];
        rep_interval_t *prev = seg->high_reps_list;
        rep_interval_t *cur = prev->next;
        seg->high_count = 0;   // 重新计数，因为下面会重建链表
        while (cur) {
            int t = cur->end - cur->start;
            if (t < s + min_offset) {
                // 删除节点
                rep_interval_t *tmp = cur;
                prev->next = cur->next;
                cur = cur->next;
                free(tmp);
                collect_rep->rep_free_num++;
            } else {
                // 保留节点
                collect_rep->rep_num++;
                total_len += t;
                seg->high_count++;
                prev = cur;
                cur = cur->next;
            }
        }
        seg->high_first = seg->high_reps_list->next;
        seg->high_last = (seg->high_first == NULL) ? NULL : prev; // prev 是最后一个保留节点
    }
    if (collect_rep->rep_num > 0)
        collect_rep->cov_length = (double)total_len / collect_rep->rep_num;
    else
        collect_rep->cov_length = 0.0;  // 避免除零
    return 1;
}

int show_seg_rep(seg_rep_t *seg_rep){
    fprintf(f_rep, "******************Segment: %s\n", seg_rep->sname);
    rep_interval_t *interval = seg_rep->high_first;
    while(interval){
        fprintf(f_rep, "[%d, %d)\n", interval->start, interval->end);
        interval = interval->next;
    }
    return 1;
}
//弃用
int fix_seg_rep(seg_rep_t *seg_rep, int s, int min_offset) {
    int total_len = 0;
    seg_rep->high_count = 0;   // 重置以便重新计数
    rep_interval_t *prev = seg_rep->high_reps_list;
    rep_interval_t *cur = prev->next;
    while (cur) {
        int t = cur->end - cur->start;
        if (t < s + min_offset) {
            // 删除节点
            rep_interval_t *tmp = cur;
            prev->next = cur->next;
            cur = cur->next;
            free(tmp);
        } else {
            // 保留节点
            total_len += t;
            seg_rep->high_count++;
            prev = cur;
            cur = cur->next;
        }
    }
    seg_rep->high_first = seg_rep->high_reps_list->next;
    seg_rep->high_last = (seg_rep->high_first == NULL) ? NULL : prev; // prev 是最后一个保留节点
    return 1;
}

void free_collect_rep(collect_rep_t *collect_rep) {
    if (!collect_rep) return;

    if (collect_rep->segs) {
        for (int i = 0; i < collect_rep->len; ++i) {
            seg_rep_t *seg = &collect_rep->segs[i];

            /* 释放 high_reps_list 链表（含哨兵节点） */
            if (seg->high_reps_list) {
                rep_interval_t *cur = seg->high_reps_list->next;
                while (cur) {
                    rep_interval_t *next = cur->next;
                    free(cur);
                    cur = next;
                }
                free(seg->high_reps_list);
                seg->high_reps_list = NULL;
            }

            /* 释放 low_reps_list 链表（含哨兵节点） */
            if (seg->low_reps_list) {
                rep_interval_t *cur = seg->low_reps_list->next;
                while (cur) {
                    rep_interval_t *next = cur->next;
                    free(cur);
                    cur = next;
                }
                free(seg->low_reps_list);
                seg->low_reps_list = NULL;
            }

            /* high_first / high_last / low_first / low_last
             * 只是指向链表节点的裸指针，不拥有内存，不需要单独释放。
             * sname 指向外部内存（g->seg[i].name 或 read name），
             * 不是本模块 malloc 出来的，也不能在这里释放。 */
            seg->high_first = seg->high_last = NULL;
            seg->low_first  = seg->low_last  = NULL;
        }
        free(collect_rep->segs);
        collect_rep->segs = NULL;
    }

    collect_rep->len = 0;
    collect_rep->rep_num = 0;
    collect_rep->rep_free_num = 0;
    collect_rep->cov_length = 0.0;
}