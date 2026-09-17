# =============================================================================
# Tarea 1 Problema 2 i. - Prof. Alex Di Genova
# Juan I. Riquelme & Martín Azúa
# =============================================================================

BEGIN {
    HASH_P  = 4294967291
    HASH_B1 = 1000003
    HASH_B2 = 999983
}

function base_code(ch,    c) {
    c = tolower(ch)
    if (c == "a") return 0
    if (c == "c") return 1
    if (c == "g") return 2
    if (c == "t") return 3
    return -1
}

# Calcula index=hash(kmer) mod M y fingerprint de 16 bits (0..65535).
# Retorna 1 si el kmer es valido (solo A,C,G,T), 0 si no.
function hash_kmer(seq, pos, k, M, out,    i, ch, bc, h1, h2) {
    h1 = 0; h2 = 0
    for (i = 0; i < k; i++) {
        ch = substr(seq, pos + i, 1)
        bc = base_code(ch)
        if (bc < 0) return 0
        h1 = (h1 * HASH_B1 + bc) % HASH_P
        h2 = (h2 * HASH_B2 + bc) % HASH_P
    }
    out["index"] = h1 % M
    out["fp"] = h2 % 65536
    return 1
}

BEGIN {
    if (REF == "") { print "ERROR: debe indicar REF=referencia.fasta" > "/dev/stderr"; exit 1 }
    if (K == "") K = 31
    if (L == "") L = 10
    if (C == "") C = 1
    USE_FP = (USE_FP=="")? 0 : USE_FP   # 0/1: usar firma de 16 bits para detectar colisiones

    # ---- 1) Leer FASTA de referencia y concatenar secuencias ----
    nrefs = 0
    seqname[0] = ""
    seq[0] = ""
    curid = -1
    while ((getline line < REF) > 0) {
        if (substr(line,1,1) == ">") {
            curid++
            nrefs++
            seqname[curid] = substr(line,2)
            seq[curid] = ""
        } else {
            gsub(/[ \t\r]/, "", line)
            seq[curid] = seq[curid] line
        }
    }
    close(REF)

    G = 0
    for (r = 0; r < nrefs; r++) G += length(seq[r])
    if (M == "") M = G
    if (M < 1) M = 1

    printf("# Indice: nrefs=%d G=%d K=%d L=%d M=%d C=%d fingerprint=%d\n", nrefs, G, K, L, M, C, USE_FP) > "/dev/stderr"

    # ---- 2) Construir indice: un k-mer cada L bases por referencia ----
    # NOTA de convencion: p y rp (mas abajo) son posiciones 0-based
    # internamente (0,L,2L,...). hash_kmer recibe una posicion 1-based para
    # poder usar substr() de awk, por eso se le pasa p+1. Al reportar la
    # posicion final se suma 1 para volver a 1-based, tal como exige el
    # formato de salida. Esta misma convencion (0-based interno, +1 al
    # imprimir) se usa identica en las versiones C y C+pthreads.
    for (r = 0; r < nrefs; r++) {
        n = length(seq[r])
        for (p = 0; p + K <= n; p += L) {
            if (hash_kmer(seq[r], p+1, K, M, hk) == 1) {
                idx = hk["index"]
                if (!(idx in cellCount)) {
                    cellCount[idx] = 0
                    cellPos[idx] = p          # se conserva la PRIMERA posicion observada (0-based)
                    cellRef[idx] = r
                    cellFp[idx]  = hk["fp"]
                    cellCol[idx] = 0
                } else if (USE_FP) {
                    if (cellFp[idx] != hk["fp"]) cellCol[idx] = 1
                }
                cellCount[idx]++
            }
        }
    }

    occupied = 0
    for (idx in cellCount) occupied++
    printf("# Celdas ocupadas=%d  ocupacion=%.4f\n", occupied, occupied/M) > "/dev/stderr"

    FS = "\n"
    lineno = 0
    totalReads = 0; hitReads = 0
}

# ---- 3) Leer lecturas FASTQ desde stdin (linea a linea: 4 lineas/registro) ----
{
    lineno++
    m = (lineno - 1) % 4
    if (m == 0) { rid = $0; sub(/^@/, "", rid) }
    else if (m == 1) {
        rseq = $0
        totalReads++
        delete votes

        rn = length(rseq)
        for (rp = 0; rp + K <= rn; rp += L) {
            if (hash_kmer(rseq, rp+1, K, M, hk) == 1) {
                idx = hk["index"]
                if (idx in cellCount) {
                    if (USE_FP && cellCol[idx] == 1) continue          # celda colisionada: descartar
                    if (USE_FP && cellFp[idx] != hk["fp"]) continue    # firma no coincide: descartar
                    cnt = cellCount[idx]
                    if (cnt >= 1 && cnt <= C) {
                        cand_start = cellPos[idx] - rp                 # posicion propuesta (0-based)
                        key = cellRef[idx] SUBSEP cand_start
                        votes[key]++
                    }
                }
            }
        }

        bestVotes = 0; bestKey = ""
        for (key in votes) {
            split(key, parts, SUBSEP)
            v = votes[key]
            if (v > bestVotes) {
                bestVotes = v; bestKey = key
            } else if (v == bestVotes && bestVotes > 0) {
                split(bestKey, bp, SUBSEP)
                if ( (parts[2]+0 < bp[2]+0) || (parts[2]+0 == bp[2]+0 && parts[1]+0 < bp[1]+0) )
                    bestKey = key
            }
        }

        if (bestVotes > 0) {
            split(bestKey, bp, SUBSEP)
            refid = bp[1]+0
            pos1  = bp[2]+0 + 1        # convertir de 0-based a 1-based para la salida
            if (pos1 < 1) pos1 = 1     # clamp defensivo (lectura "cuelga" antes del inicio)
            print rid "\t" "hit" "\t" seqname[refid] "\t" pos1 "\t" bestVotes
            hitReads++
            scoreDist[bestVotes]++
        } else {
            print rid "\t" "no-hit" "\t" "*" "\t" 0 "\t" 0
            scoreDist[0]++
        }
    }
}

END {
    print "# --- Resumen ---" > "/dev/stderr"
    printf("# total_reads=%d hits=%d no_hits=%d pct_hits=%.2f%%\n", totalReads, hitReads, totalReads-hitReads, (totalReads>0? hitReads/totalReads*100:0)) > "/dev/stderr"
    printf("# distribucion_de_scores:") > "/dev/stderr"
    for (s in scoreDist) printf(" %d:%d", s, scoreDist[s]) > "/dev/stderr"
    printf("\n") > "/dev/stderr"
}
