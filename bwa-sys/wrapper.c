#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "ext/bwa/bwamem.h"

// Function to create a new bseq1_t instance
bseq1_t *new_bseq1_t(int l_seq, int id, const char *name, const char *seq, const char *qual) {
    // Allocate memory for the structure itself
    bseq1_t *new_seq = (bseq1_t *)malloc(sizeof(bseq1_t));
    if (new_seq == NULL) {
        printf("Memory allocation failed!\n");
        return NULL;
    }

    // Initialize the id and length of the sequence
    new_seq->id = id;
    new_seq->l_seq = l_seq;

    // Dynamically allocate and copy strings for name, comment, seq, qual, and sam
    new_seq->name = strdup(name);
    new_seq->comment = 0;
    new_seq->sam = 0;

	new_seq->seq = (char *)malloc(l_seq);
	new_seq->qual = (char *)malloc(l_seq);
	memcpy(new_seq->seq, seq, l_seq);
	memcpy(new_seq->qual, qual, l_seq);

    return new_seq;
}

// Function to free the memory allocated for a bseq1_t instance
void free_bseq1_t(bseq1_t *seq) {
    if (seq) {
        free(seq->name);
        free(seq->comment);
        free(seq->seq);
        free(seq->qual);
        free(seq->sam);
        free(seq);
    }
}

typedef struct {
	bwtintv_v mem, mem1, *tmpv[2];
} smem_aux_t;

smem_aux_t *bwa_smem_aux_init()
{
	smem_aux_t *a;
	a = calloc(1, sizeof(smem_aux_t));
	a->tmpv[0] = calloc(1, sizeof(bwtintv_v));
	a->tmpv[1] = calloc(1, sizeof(bwtintv_v));
	return a;
}

void bwa_smem_aux_destroy(smem_aux_t *a)
{	
	free(a->tmpv[0]->a); free(a->tmpv[0]);
	free(a->tmpv[1]->a); free(a->tmpv[1]);
	free(a->mem.a); free(a->mem1.a);
	free(a);
}

char* mem_align(const mem_opt_t* opt, const bwaidx_t* idx, const char *seq_name, const char *seq_seq, const char *seq_qual, int l_seq)
{
	extern mem_alnreg_v mem_align1_core(const mem_opt_t *opt, const bwt_t *bwt, const bntseq_t *bns, const uint8_t *pac, int l_seq, char *seq, void *buf);
	extern int mem_mark_primary_se(const mem_opt_t *opt, int n, mem_alnreg_t *a, int64_t id);
	extern void mem_reorder_primary5(int T, mem_alnreg_v *a);
	extern void mem_reg2sam(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, bseq1_t *s, mem_alnreg_v *a, int extra_flag, const mem_aln_t *m);

	bseq1_t *seq = new_bseq1_t(l_seq, 0, seq_name, seq_seq, seq_qual);

	smem_aux_t *buf = bwa_smem_aux_init();
	mem_alnreg_v reg = mem_align1_core(opt, idx->bwt, idx->bns, idx->pac, seq->l_seq, seq->seq, buf);
	bwa_smem_aux_destroy(buf);

	mem_mark_primary_se(opt, reg.n, reg.a, 0);
	if (opt->flag & MEM_F_PRIMARY5) {
		mem_reorder_primary5(opt->T, &reg);
	}
	mem_reg2sam(opt, idx->bns, idx->pac, seq, &reg, 0, 0);
	char *sam = strdup(seq->sam);

	free(reg.a);
	free_bseq1_t(seq);
	return sam;
}