/* =============================================================================
 * Tarea 1 Problema 2 ii. - Prof. Alex Di Genova
 * Juan I. Riquelme & Martín Azúa
 * ============================================================================= */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <ctype.h>
#include <getopt.h>
#include <time.h>
#include <sys/resource.h>


#define HASH_P  4294967291ULL
#define HASH_B1 1000003ULL
#define HASH_B2 999983ULL

/* Codifica k bases desde s (sin '\0' en medio). Retorna 1 si son todas
 * validas (A,C,G,T), 0 si alguna base es distinta. Deja el codigo de bases
 * (0..3 cada una) en code_out[0..k-1] para reusar en las dos pasadas del hash. */
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

/* ---------- Estructuras ---------- */
typedef struct {
    uint32_t count;
    uint64_t position;      /* 0-based, primera posicion observada */
    uint32_t reference_id;
    uint16_t fingerprint;
    uint8_t  collided;
    uint8_t  used;          /* celda ocupada? */
} index_cell;

typedef struct {
    char *name;
    char *seq;
    size_t len;
} ref_seq;

/* ---------- Lectura de FASTA (sencilla, guarda todo en memoria) ---------- */
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
            if (n == 0) { fprintf(stderr, "FASTA invalido: secuencia sin encabezado\n"); exit(1); }
            ref_seq *r = &refs[n-1];
            size_t oldlen = r->len;
            r->seq = realloc(r->seq, oldlen + (size_t)len + 1);
            memcpy(r->seq + oldlen, line, (size_t)len);
            r->len += (size_t)len;
            r->seq[r->len] = '\0';
            G += (size_t)len;
        }
    }
    free(line);
    fclose(f);
    *nrefs_out = n;
    *G_out = G;
    return refs;
}

/* ---------- Construccion del indice ---------- */
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
                cell->used = 1;
                cell->count = 0;
                cell->position = p;         /* 0-based, primera posicion */
                cell->reference_id = (uint32_t)r;
                cell->fingerprint = fp;
                cell->collided = 0;
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

/* ---------- Votos por candidato durante la consulta ---------- */
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

        int64_t cand_start = (int64_t)cell->position - (int64_t)rp; /* 0-based */

        /* buscar si ya existe ese candidato (ref, start) y acumular voto */
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
        if (cands[i].votes > cands[best].votes) {
            best = i;
        } else if (cands[i].votes == cands[best].votes) {
            /* desempate: posicion menor, luego orden de referencia (indice menor) */
            if (cands[i].start < cands[best].start) best = i;
            else if (cands[i].start == cands[best].start && cands[i].ref < cands[best].ref) best = i;
        }
    }
    *best_idx = best;
    return 1;
}

/* ---------- Lectura de FASTQ registro a registro ---------- */
typedef struct { char *id; char *seq; } fastq_rec;

/* Lee el siguiente registro (4 lineas). Retorna 1 si ok, 0 si EOF. */
static int read_fastq_record(FILE *f, char **line, size_t *linecap, fastq_rec *rec) {
    ssize_t len;
    if ((len = getline(line, linecap, f)) == -1) return 0;             /* @id */
    while (len > 0 && ((*line)[len-1]=='\n' || (*line)[len-1]=='\r')) (*line)[--len]='\0';
    char *id = (*line)[0]=='@' ? *line + 1 : *line;
    free(rec->id); rec->id = strdup(id);

    if ((len = getline(line, linecap, f)) == -1) return 0;             /* seq */
    while (len > 0 && ((*line)[len-1]=='\n' || (*line)[len-1]=='\r')) (*line)[--len]='\0';
    free(rec->seq); rec->seq = strdup(*line);

    if (getline(line, linecap, f) == -1) return 0;                     /* + */
    if (getline(line, linecap, f) == -1) return 0;                     /* calidad */
    return 1;
}

int main(int argc, char **argv) {
    const char *ref_path = NULL, *fastq_path = NULL;
    size_t K = 31, L = 10, C = 1, M = 0;
    int use_fp = 0;
    int opt;
    while ((opt = getopt(argc, argv, "r:q:k:l:c:m:f")) != -1) {
        switch (opt) {
            case 'r': ref_path = optarg; break;
            case 'q': fastq_path = optarg; break;
            case 'k': K = (size_t)atoi(optarg); break;
            case 'l': L = (size_t)atoi(optarg); break;
            case 'c': C = (size_t)atoi(optarg); break;
            case 'm': M = (size_t)atol(optarg); break;
            case 'f': use_fp = 1; break;
            default:
                fprintf(stderr, "Uso: %s -r ref.fasta -q lecturas.fastq [-k 31] [-l 10] [-c 1] [-m M] [-f]\n", argv[0]);
                return 1;
        }
    }
    if (!ref_path || !fastq_path) {
        fprintf(stderr, "Faltan -r y/o -q\n");
        return 1;
    }

    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);

    size_t nrefs, G;
    ref_seq *refs = load_fasta(ref_path, &nrefs, &G);
    if (M == 0) M = G > 0 ? G : 1;

    index_cell *table = calloc(M, sizeof(index_cell));
    if (!table) { fprintf(stderr, "No hay memoria para el indice (M=%zu)\n", M); return 1; }

    size_t occupied = 0, discarded = 0;
    build_index(refs, nrefs, K, L, table, M, use_fp, &occupied, &discarded, C);

    fprintf(stderr, "# Indice: nrefs=%zu G=%zu K=%zu L=%zu M=%zu C=%zu fingerprint=%d\n",
            nrefs, G, K, L, M, C, use_fp);
    fprintf(stderr, "# Celdas ocupadas=%zu ocupacion=%.4f descartadas_por_c=%zu\n",
            occupied, (double)occupied / (double)M, discarded);

    FILE *fin = strcmp(fastq_path, "-") == 0 ? stdin : fopen(fastq_path, "r");
    if (!fin) { perror("fopen fastq"); return 1; }

    char *line = NULL; size_t linecap = 0;
    fastq_rec rec = {0};
    size_t cand_cap = 4096;
    candidate *cands = malloc(cand_cap * sizeof(candidate));

    size_t total = 0, hits = 0;
    /* distribucion de score: usamos un arreglo dinamico simple indexado por score */
    size_t distcap = 64; size_t *dist = calloc(distcap, sizeof(size_t));

    while (read_fastq_record(fin, &line, &linecap, &rec)) {
        total++;
        size_t rlen = strlen(rec.seq);
        size_t ncands;
        process_read(rec.seq, rlen, K, L, C, use_fp, table, M, cands, cand_cap, &ncands);
        size_t best;
        if (pick_best(cands, ncands, &best)) {
            uint32_t score = cands[best].votes;
            int64_t pos0 = cands[best].start;
            int64_t pos1 = pos0 + 1; if (pos1 < 1) pos1 = 1; /* 0-based -> 1-based, clamp defensivo */
            printf("%s\thit\t%s\t%ld\t%u\n", rec.id, refs[cands[best].ref].name, (long)pos1, score);
            hits++;
            if (score >= distcap) { size_t newcap = score + 1; dist = realloc(dist, newcap * sizeof(size_t)); memset(dist+distcap, 0, (newcap-distcap)*sizeof(size_t)); distcap = newcap; }
            dist[score]++;
        } else {
            printf("%s\tno-hit\t*\t0\t0\n", rec.id);
            dist[0]++;
        }
    }

    clock_gettime(CLOCK_MONOTONIC, &t1);
    double secs = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1e9;
    struct rusage ru; getrusage(RUSAGE_SELF, &ru);

    fprintf(stderr, "# --- Resumen ---\n");
    fprintf(stderr, "# total_reads=%zu hits=%zu no_hits=%zu pct_hits=%.2f%%\n",
            total, hits, total - hits, total ? (100.0*hits/total) : 0.0);
    fprintf(stderr, "# tiempo_real_seg=%.3f reads_por_seg=%.1f mem_residente_kb=%ld\n",
            secs, secs > 0 ? total/secs : 0.0, ru.ru_maxrss);
    fprintf(stderr, "# distribucion_de_scores:");
    for (size_t s = 0; s < distcap; s++) if (dist[s]) fprintf(stderr, " %zu:%zu", s, dist[s]);
    fprintf(stderr, "\n");

    return 0;
}
