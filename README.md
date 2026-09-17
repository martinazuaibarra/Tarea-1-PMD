# Tarea 1 — Procesamiento Masivo de Datos

Este paquete contiene el desarrollo completo de la tarea: Parte 1 (comandos
Linux/awk) y Parte 2 (alineador probabilístico por k-mers en awk, C y
C+pthreads).

## Estructura

```
tarea1/
├── README.md              
├── informe.pdf            
├── part1/
│   └── P1.sh              todos los comandos awk de la Parte 1 (a-j)
└── part2/
    ├── mapper.awk              versión secuencial en AWK 
    ├── mapper.c                versión secuencial en C
    ├── mapper_pthreads.c       versión paralela en C con pthreads
    └── Makefile

```

## Parte 1 — Comandos Linux/awk

Ejecución
```bash
cd part1
bash P1.sh data/worldcitiespop.csv.gz data/matrix.txt
```

Genera un archivo `out_X_*.{tsv,txt}` por cada apartado (a–j), con el
comando awk correspondiente documentado en línea dentro del script.

## Parte 2 — Alineador por k-mers

### Compilación

```bash
cd part2
make            # genera ./mapper y ./mapper_pthreads
```

La versión awk no requiere compilación (usa gawk/awk directamente).

### Ejecución para E. coli - AWK

```bash
# k = 15
awk -v REF=data/ecoli-k12-ref.fna -v K=15 -v L=10 -v C=1 -f mapper.awk \
    <(gzip -dc data/EC.50X.R1.fastq.gz) > ecoli_awk_k15.tsv 2> ecoli.awk_k15.log

# k = 31
awk -v REF=data/ecoli-k12-ref.fna -v K=31 -v L=10 -v C=1 -f mapper.awk \
    <(gzip -dc data/EC.50X.R1.fastq.gz) > ecoli_awk_k31.tsv 2> ecoli.awk_k31.log    
```

### Ejecución para E. coli - C secuencial

```bash
# k = 15
./mapper -r data/ecoli-k12-ref.fna -q <(gzip -dc data/EC.50X.R1.fastq.gz) -k 15 -l 10 -c 1 \
    > ecoli_c_k15.tsv 2> ecoli_c_k15.log

# k = 31
./mapper -r data/ecoli-k12-ref.fna -q <(gzip -dc data/EC.50X.R1.fastq.gz) -k 31 -l 10 -c 1 \
    > ecoli_c_k31.tsv 2> ecoli_c_k31.log    
```

Flags: `-r` referencia FASTA, `-q` lecturas FASTQ (`-` para stdin), `-k`
largo de k-mer, `-l` salto L, `-c` límite de frecuencia c, `-m` tamaño de
índice M (default: G), `-f` activa la variante con fingerprint.

### Ejecución para E. coli - C con pthreads

```bash
# k = 15
./mapper_pthreads -r data/ecoli-k12-ref.fna -q <(gzip -dc data/EC.50X.R1.fastq.gz) \
    -k 15 -l 10 -c 1 -t 4 > ecoli_pt_k15.tsv 2> ecoli_pt_k15.log

# k = 31
./mapper_pthreads -r data/ecoli-k12-ref.fna -q <(gzip -dc data/EC.50X.R1.fastq.gz) \
    -k 31 -l 10 -c 1 -t 4 > ecoli_pt_k31.tsv 2> ecoli_pt_k31.log    
```

Flags adicionales: `-t` número de hilos trabajadores, `-x` tamaño de bloque
X (lecturas por bloque que arma el productor), `-y` capacidad Y de la cola
acotada (en bloques).

### Verificar que las tres versiones dan el mismo resultado

```bash
awk -v REF=data/ecoli-k12-ref.fna -v K=31 -v L=10 -v C=1 -f mapper.awk \
    <(gzip -dc data/EC.50X.R1.fastq.gz) > ecoli_awk_k31.tsv 2> ecoli.awk_k31.log  
./mapper -r data/ecoli-k12-ref.fna -q <(gzip -dc data/EC.50X.R1.fastq.gz) -k 31 -l 10 -c 1 \
    > ecoli_c_k31.tsv 2> ecoli_c_k31.log    
./mapper_pthreads -r data/ecoli-k12-ref.fna -q <(gzip -dc data/EC.50X.R1.fastq.gz) \
    -k 31 -l 10 -c 1 -t 4 > ecoli_pt_k31.tsv 2> ecoli_pt_k31.log  
diff out_awk.tsv out_c.tsv && diff out_c.tsv out_pt.tsv && echo "Las tres versiones coinciden"
```

### Ejecución para Humano - C secuencial y C con pthreads 

```bash
# Tomar una muestra reproducible de N lecturas antes de correr, por ejemplo:
zcat data/HG002-MGISEQ-L03-1.fq.gz | head -n 400000 | gzip > data/HG002_sample_100k.fastq.gz  

timeout 3600 ./mapper_pthreads -r data/GRCh38-full-analysis-set.fna \
    -q <(gzip -dc data/HG002_sample_100k.fastq.gz) \
    -k 31 -l 10 -c 1 -m 100000000 -t 8 \         # Tamaño reducido por memoria RAM del ordenador 
    > human_pt.tsv 2> human_pt.log
```

Si no termina dentro de 1 hora, se debe interrumpir (`timeout 3600 ...`) y
reportar el avance (líneas de salida ya escritas) como indica el enunciado.

## Formato de salida

Cada ejecución produce, por stdout, una línea por lectura:

```
read_id  status  reference  position  score
```

y por stderr, un resumen con parámetros usados, ocupación del índice,
tiempo, memoria y distribución de scores — pensado para parsear
directamente hacia las tablas del informe (grep de líneas que empiezan con
`#`).
