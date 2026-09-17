/* =============================================================================
 * Tarea 1 Problema 2 iii. - Prof. Alex Di Genova
 * Juan I. Riquelme & Martín Azúa
 * ============================================================================= */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <pthread.h>
#include <getopt.h>
#include <time.h>
#include <sys/resource.h>

/* ---------- Hash: identico a mapper.c (ver justificacion alli) ---------- */
#define HASH_P  4294967291ULL
#define HASH_B1 1000003ULL
#define HASH_B2 999983ULL

static int encode_kmer(const char *s, size_t k, uint8_t *code_out) {
    for (size_t j = 0; j < k; j++) {
        switch (s[j]) {
            case 'A': case 'a': code_out[j] = 0; break;
            case 'C': case 'c': code_out[j] = 1; break;
            case 'G': case 'g': code_out[j] = 2; break;
            case 'T': case 't': code_out[j] = 3; break;
            default: return 0;
        }
    }
    return 1;
}

static void hash_kmer(const uint8_t *code, size_t k, size_t M, size_t *index, uint16_t *fingerprint) {
    uint64_t h1 = 0, h2 = 0;
    for (size_t j = 0; j < k; j++) {
        h1 = (h1 * HASH_B1 + code[j]) % HASH_P;
        h2 = (h2 * HASH_B2 + code[j]) % HASH_P;
    }
    *index = (size_t)(h1 % (uint64_t)M);
    *fingerprint = (uint16_t)(h2 % 65536ULL);
}

/* ---------- Estructuras compartidas (indice, de solo lectura durante el mapeo) ---------- */
typedef struct {
    uint32_t count;
    uint64_t position;
    uint32_t reference_id;
    uint16_t fingerprint;
    uint8_t  collided;
    uint8_t  used;
} index_cell;

typedef struct {
    char *name;
    char *seq;
    size_t len;
} ref_seq;

static ref_seq *load_fasta(const char *path, size_t *nrefs_out, size_t *G_out) {
    FILE *f = fopen(path, "r");
    if (!f) { perror("fopen referencia"); exit(1); }
    size_t cap = 8, n = 0;
    ref_seq *refs = malloc(cap * sizeof(ref_seq));
    char *line = NULL; size_t linecap = 0; ssize_t len;
    size_t G = 0;
    while ((len = getline(&line, &linecap, f)) != -1) {
        while (len > 0 && (line[len-1] == '\n' || line[len-1] == '\r')) line[--len] = '\0';
        if (len == 0) continue;
        if (line[0] == '>') {
            if (n == cap) { cap *= 2; refs = realloc(refs, cap * sizeof(ref_seq)); }
            refs[n].name = strdup(line + 1);
            refs[n].seq = malloc(1); refs[n].seq[0] = '\0';
            refs[n].len = 0;
            n++;
        } else {
            if (n == 0) { fprintf(stderr, "FASTA invalido\n"); exit(1); }
            ref_seq *r = &refs[n-1];
            size_t oldlen = r->len;
            r->seq = realloc(r->seq, oldlen + (size_t)len + 1);
            memcpy(r->seq + oldlen, line, (size_t)len);
            r->len += (size_t)len;
            r->seq[r->len] = '\0';
            G += (size_t)len;
        }
    }
    free(line); fclose(f);
    *nrefs_out = n; *G_out = G;
    return refs;
}

static void build_index(ref_seq *refs, size_t nrefs, size_t K, size_t L,
                         index_cell *table, size_t M, int use_fp,
                         size_t *occupied_out, size_t *discarded_by_c_out, size_t C) {
    size_t occupied = 0;
    for (size_t r = 0; r < nrefs; r++) {
        size_t n = refs[r].len;
        if (n < K) continue;
        uint8_t *codebuf = malloc(K);
        for (size_t p = 0; p + K <= n; p += L) {
            if (!encode_kmer(refs[r].seq + p, K, codebuf)) continue;
            size_t idx; uint16_t fp;
            hash_kmer(codebuf, K, M, &idx, &fp);
            index_cell *cell = &table[idx];
            if (!cell->used) {
                cell->used = 1; cell->count = 0; cell->position = p;
                cell->reference_id = (uint32_t)r; cell->fingerprint = fp; cell->collided = 0;
                occupied++;
            } else if (use_fp) {
                if (cell->fingerprint != fp) cell->collided = 1;
            }
            cell->count++;
        }
        free(codebuf);
    }
    *occupied_out = occupied;
    size_t discarded = 0;
    for (size_t i = 0; i < M; i++)
        if (table[i].used && (table[i].count < 1 || table[i].count > C)) discarded++;
    *discarded_by_c_out = discarded;
}

typedef struct { uint32_t ref; int64_t start; uint32_t votes; } candidate;

static void process_read(const char *rseq, size_t rlen, size_t K, size_t L,
                          size_t C, int use_fp, index_cell *table, size_t M,
                          candidate *cands, size_t cand_cap, size_t *ncands_out) {
    size_t ncands = 0;
    uint8_t *codebuf = malloc(K);
    for (size_t rp = 0; rp + K <= rlen; rp += L) {
        if (!encode_kmer(rseq + rp, K, codebuf)) continue;
        size_t idx; uint16_t fp;
        hash_kmer(codebuf, K, M, &idx, &fp);
        index_cell *cell = &table[idx];
        if (!cell->used) continue;
        if (use_fp) {
            if (cell->collided) continue;
            if (cell->fingerprint != fp) continue;
        }
        if (cell->count < 1 || cell->count > C) continue;
        int64_t cand_start = (int64_t)cell->position - (int64_t)rp;
        int found = 0;
        for (size_t i = 0; i < ncands; i++) {
            if (cands[i].ref == cell->reference_id && cands[i].start == cand_start) {
                cands[i].votes++; found = 1; break;
            }
        }
        if (!found && ncands < cand_cap) {
            cands[ncands].ref = cell->reference_id;
            cands[ncands].start = cand_start;
            cands[ncands].votes = 1;
            ncands++;
        }
    }
    free(codebuf);
    *ncands_out = ncands;
}

static int pick_best(candidate *cands, size_t ncands, size_t *best_idx) {
    if (ncands == 0) return 0;
    size_t best = 0;
    for (size_t i = 1; i < ncands; i++) {
        if (cands[i].votes > cands[best].votes) best = i;
        else if (cands[i].votes == cands[best].votes) {
            if (cands[i].start < cands[best].start) best = i;
            else if (cands[i].start == cands[best].start && cands[i].ref < cands[best].ref) best = i;
        }
    }
    *best_idx = best;
    return 1;
}

/* ---------- Bloques de lecturas y cola acotada productor/consumidor ---------- */
typedef struct { char *id; char *seq; } fastq_rec;

typedef struct {
    long block_id;
    fastq_rec *reads;
    size_t nreads;
    char **out_lines;      /* llenado por el trabajador que procesa el bloque */
} block_t;

#define QCAP_MAX 4096
typedef struct {
    block_t *buf[QCAP_MAX];
    size_t cap, head, tail, count;
    int done;              /* el productor ya termino de leer el FASTQ */
    pthread_mutex_t mutex;
    pthread_cond_t not_empty, not_full;
} bounded_queue;

static void q_init(bounded_queue *q, size_t cap) {
    q->cap = cap; q->head = q->tail = q->count = 0; q->done = 0;
    pthread_mutex_init(&q->mutex, NULL);
    pthread_cond_init(&q->not_empty, NULL);
    pthread_cond_init(&q->not_full, NULL);
}
static void q_push(bounded_queue *q, block_t *b) {
    pthread_mutex_lock(&q->mutex);
    while (q->count == q->cap) pthread_cond_wait(&q->not_full, &q->mutex);
    q->buf[q->tail] = b; q->tail = (q->tail + 1) % q->cap; q->count++;
    pthread_cond_signal(&q->not_empty);
    pthread_mutex_unlock(&q->mutex);
}
static block_t *q_pop(bounded_queue *q) {
    pthread_mutex_lock(&q->mutex);
    while (q->count == 0 && !q->done) pthread_cond_wait(&q->not_empty, &q->mutex);
    if (q->count == 0 && q->done) { pthread_mutex_unlock(&q->mutex); return NULL; }
    block_t *b = q->buf[q->head]; q->head = (q->head + 1) % q->cap; q->count--;
    pthread_cond_signal(&q->not_full);
    pthread_mutex_unlock(&q->mutex);
    return b;
}
static void q_mark_done(bounded_queue *q) {
    pthread_mutex_lock(&q->mutex);
    q->done = 1;
    pthread_cond_broadcast(&q->not_empty);
    pthread_mutex_unlock(&q->mutex);
}

/* ---------- Contexto compartido ---------- */
typedef struct {
    ref_seq *refs;
    index_cell *table;
    size_t M, K, L, C;
    int use_fp;
    bounded_queue *q;

    /* impresion ordenada de bloques ya resueltos */
    pthread_mutex_t out_mutex;
    block_t **finished;      /* finished[block_id] = bloque listo, o NULL */
    size_t finished_cap;
    long next_to_print;

    /* estadisticas globales (protegidas por out_mutex) */
    size_t total_reads, hit_reads;
    size_t distcap; size_t *dist;
} shared_ctx;

static int read_fastq_record(FILE *f, char **line, size_t *linecap, fastq_rec *rec) {
    ssize_t len;
    if ((len = getline(line, linecap, f)) == -1) return 0;
    while (len > 0 && ((*line)[len-1]=='\n' || (*line)[len-1]=='\r')) (*line)[--len]='\0';
    char *id = (*line)[0]=='@' ? *line + 1 : *line;
    rec->id = strdup(id);
    if ((len = getline(line, linecap, f)) == -1) return 0;
    while (len > 0 && ((*line)[len-1]=='\n' || (*line)[len-1]=='\r')) (*line)[--len]='\0';
    rec->seq = strdup(*line);
    if (getline(line, linecap, f) == -1) return 0;
    if (getline(line, linecap, f) == -1) return 0;
    return 1;
}

typedef struct { FILE *fin; size_t X; bounded_queue *q; } producer_args;

static void *producer_main(void *arg) {
    producer_args *pa = (producer_args *)arg;
    char *line = NULL; size_t linecap = 0;
    long block_id = 0;
    for (;;) {
        block_t *b = malloc(sizeof(block_t));
        b->reads = malloc(pa->X * sizeof(fastq_rec));
        b->nreads = 0;
        while (b->nreads < pa->X) {
            fastq_rec rec;
            if (!read_fastq_record(pa->fin, &line, &linecap, &rec)) break;
            b->reads[b->nreads++] = rec;
        }
        if (b->nreads == 0) { free(b->reads); free(b); break; }
        b->block_id = block_id++;
        b->out_lines = NULL;
        q_push(pa->q, b);
        if (b->nreads < pa->X) break; /* fue el ultimo bloque (parcial) */
    }
    free(line);
    q_mark_done(pa->q);
    return NULL;
}

/* Intenta imprimir todos los bloques consecutivos ya resueltos a partir de
 * next_to_print. Debe llamarse con out_mutex tomado. */
static void drain_ready_blocks(shared_ctx *ctx) {
    while (ctx->next_to_print < (long)ctx->finished_cap &&
           ctx->finished[ctx->next_to_print] != NULL) {
        block_t *b = ctx->finished[ctx->next_to_print];
        for (size_t i = 0; i < b->nreads; i++) fputs(b->out_lines[i], stdout);
        for (size_t i = 0; i < b->nreads; i++) { free(b->out_lines[i]); free(b->reads[i].id); free(b->reads[i].seq); }
        free(b->out_lines); free(b->reads);
        ctx->finished[ctx->next_to_print] = NULL;
        free(b);
        ctx->next_to_print++;
    }
}

static void *worker_main(void *arg) {
    shared_ctx *ctx = (shared_ctx *)arg;
    size_t cand_cap = 4096;
    candidate *cands = malloc(cand_cap * sizeof(candidate));

    for (;;) {
        block_t *b = q_pop(ctx->q);
        if (!b) break;

        b->out_lines = malloc(b->nreads * sizeof(char *));
        size_t local_hits = 0;
        size_t local_distcap = 64;
        size_t *local_dist = calloc(local_distcap, sizeof(size_t));

        for (size_t i = 0; i < b->nreads; i++) {
            size_t rlen = strlen(b->reads[i].seq);
            size_t ncands;
            process_read(b->reads[i].seq, rlen, ctx->K, ctx->L, ctx->C, ctx->use_fp,
                         ctx->table, ctx->M, cands, cand_cap, &ncands);
            size_t best; char buf[512];
            if (pick_best(cands, ncands, &best)) {
                uint32_t score = cands[best].votes;
                int64_t pos1 = cands[best].start + 1; if (pos1 < 1) pos1 = 1;
                snprintf(buf, sizeof(buf), "%s\thit\t%s\t%ld\t%u\n",
                         b->reads[i].id, ctx->refs[cands[best].ref].name, (long)pos1, score);
                local_hits++;
                if (score >= local_distcap) {
                    size_t newcap = score + 1;
                    local_dist = realloc(local_dist, newcap * sizeof(size_t));
                    memset(local_dist + local_distcap, 0, (newcap - local_distcap) * sizeof(size_t));
                    local_distcap = newcap;
                }
                local_dist[score]++;
            } else {
                snprintf(buf, sizeof(buf), "%s\tno-hit\t*\t0\t0\n", b->reads[i].id);
                local_dist[0]++;
            }
            b->out_lines[i] = strdup(buf);
        }

        pthread_mutex_lock(&ctx->out_mutex);
        if ((size_t)b->block_id >= ctx->finished_cap) {
            size_t newcap = (size_t)b->block_id + 16;
            ctx->finished = realloc(ctx->finished, newcap * sizeof(block_t *));
            for (size_t i = ctx->finished_cap; i < newcap; i++) ctx->finished[i] = NULL;
            ctx->finished_cap = newcap;
        }
        ctx->finished[b->block_id] = b;
        ctx->total_reads += b->nreads;
        ctx->hit_reads += local_hits;
        for (size_t s = 0; s < local_distcap; s++) {
            if (!local_dist[s]) continue;
            if (s >= ctx->distcap) {
                size_t newcap = s + 1;
                ctx->dist = realloc(ctx->dist, newcap * sizeof(size_t));
                memset(ctx->dist + ctx->distcap, 0, (newcap - ctx->distcap) * sizeof(size_t));
                ctx->distcap = newcap;
            }
            ctx->dist[s] += local_dist[s];
        }
        free(local_dist);
        drain_ready_blocks(ctx);
        pthread_mutex_unlock(&ctx->out_mutex);
    }
    free(cands);
    return NULL;
}

int main(int argc, char **argv) {
    const char *ref_path = NULL, *fastq_path = NULL;
    size_t K = 31, L = 10, C = 1, M = 0;
    int use_fp = 0;
    int nthreads = 4;
    size_t X = 256, Y = 8;
    int opt;
    while ((opt = getopt(argc, argv, "r:q:k:l:c:m:ft:x:y:")) != -1) {
        switch (opt) {
            case 'r': ref_path = optarg; break;
            case 'q': fastq_path = optarg; break;
            case 'k': K = (size_t)atoi(optarg); break;
            case 'l': L = (size_t)atoi(optarg); break;
            case 'c': C = (size_t)atoi(optarg); break;
            case 'm': M = (size_t)atol(optarg); break;
            case 'f': use_fp = 1; break;
            case 't': nthreads = atoi(optarg); break;
            case 'x': X = (size_t)atoi(optarg); break;
            case 'y': Y = (size_t)atoi(optarg); break;
            default:
                fprintf(stderr, "Uso: %s -r ref.fasta -q lecturas.fastq [-k 31] [-l 10] [-c 1] [-m M] [-f] [-t hilos] [-x X] [-y Y]\n", argv[0]);
                return 1;
        }
    }
    if (!ref_path || !fastq_path) { fprintf(stderr, "Faltan -r y/o -q\n"); return 1; }
    if (Y > QCAP_MAX) Y = QCAP_MAX;
    if (nthreads < 1) nthreads = 1;

    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);

    size_t nrefs, G;
    ref_seq *refs = load_fasta(ref_path, &nrefs, &G);
    if (M == 0) M = G > 0 ? G : 1;

    index_cell *table = calloc(M, sizeof(index_cell));
    if (!table) { fprintf(stderr, "No hay memoria para el indice (M=%zu)\n", M); return 1; }

    size_t occupied = 0, discarded = 0;
    build_index(refs, nrefs, K, L, table, M, use_fp, &occupied, &discarded, C);
    fprintf(stderr, "# Indice: nrefs=%zu G=%zu K=%zu L=%zu M=%zu C=%zu fingerprint=%d hilos=%d X=%zu Y=%zu\n",
            nrefs, G, K, L, M, C, use_fp, nthreads, X, Y);
    fprintf(stderr, "# Celdas ocupadas=%zu ocupacion=%.4f descartadas_por_c=%zu\n",
            occupied, (double)occupied / (double)M, discarded);

    FILE *fin = strcmp(fastq_path, "-") == 0 ? stdin : fopen(fastq_path, "r");
    if (!fin) { perror("fopen fastq"); return 1; }

    bounded_queue q; q_init(&q, Y);
    shared_ctx ctx = {0};
    ctx.refs = refs; ctx.table = table; ctx.M = M; ctx.K = K; ctx.L = L; ctx.C = C;
    ctx.use_fp = use_fp; ctx.q = &q;
    pthread_mutex_init(&ctx.out_mutex, NULL);
    ctx.finished_cap = 16;
    ctx.finished = calloc(ctx.finished_cap, sizeof(block_t *));
    ctx.next_to_print = 0;
    ctx.distcap = 64;
    ctx.dist = calloc(ctx.distcap, sizeof(size_t));

    producer_args pa = { fin, X, &q };
    pthread_t producer;
    pthread_create(&producer, NULL, producer_main, &pa);

    pthread_t *workers = malloc(nthreads * sizeof(pthread_t));
    for (int i = 0; i < nthreads; i++) pthread_create(&workers[i], NULL, worker_main, &ctx);

    pthread_join(producer, NULL);
    for (int i = 0; i < nthreads; i++) pthread_join(workers[i], NULL);

    clock_gettime(CLOCK_MONOTONIC, &t1);
    double secs = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1e9;
    struct rusage ru; getrusage(RUSAGE_SELF, &ru);

    fprintf(stderr, "# --- Resumen ---\n");
    fprintf(stderr, "# total_reads=%zu hits=%zu no_hits=%zu pct_hits=%.2f%%\n",
            ctx.total_reads, ctx.hit_reads, ctx.total_reads - ctx.hit_reads,
            ctx.total_reads ? (100.0*ctx.hit_reads/ctx.total_reads) : 0.0);
    fprintf(stderr, "# tiempo_real_seg=%.3f reads_por_seg=%.1f mem_residente_kb=%ld nucleos_hilos=%d\n",
            secs, secs > 0 ? ctx.total_reads/secs : 0.0, ru.ru_maxrss, nthreads);
    fprintf(stderr, "# distribucion_de_scores:");
    for (size_t s = 0; s < ctx.distcap; s++) if (ctx.dist[s]) fprintf(stderr, " %zu:%zu", s, ctx.dist[s]);
    fprintf(stderr, "\n");

    return 0;
}
