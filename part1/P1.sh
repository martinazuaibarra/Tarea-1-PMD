# Tarea 1 Problema 1
# Procesamiento Masivo de Datos - Prof. Alex Di Genova
# Juan I. Riquelme & Martín Azúa
# Ejecución: bash P1.sh data/worldcities.csv.gz data/matrix.txt

echo "a) Promedio de habitantes por ciudad y país"

zcat worldcitiespop.csv.gz | awk -F "," '
NR > 1 && $5 != "" {
    lugar = $1 "," $2
    suma[lugar] += $5
    cantidad[lugar]++
}

END {
    print "País, Ciudad, Promedio_habitantes"

    for (lugar in suma) {
        print lugar "," suma[lugar] / cantidad[lugar]
    }    
}
'

echo "b) Top 10 ciudades con mayor población"

zcat worldcitiespop.csv.gz |
awk -F "," 'NR > 1 && &5 != "" {print $5, $3, $1}' |
sort -nr |
head -n 10

echo "c) Porcentaje de ciudades sin población de cada país"

zcat worldcitiespop.csv.gz |
awk -F "," '
NR > 1 {
    total[$1]++
    if ($5 == "") vacias[$1]++
}
END {
    for (pais in total)
        print pais, 100 * vacias[pais] / total[pais] "%"
}'

echo "d) Número de habitantes y ciudad más poblada de América del Sur"

zcat worldcitiespop.csv.gz |
awk -F "," '
NR > 1 &&
tolower($1) - /^(ar|bo|br|cl|co|ec|fk|gf|gy|pe|py|sr|uy|ve)$/ &&
$5 > maxima {
    maxima = $5
    ciudad = $3
    pais = $1
}
END {
    print ciudad, pais, maximo
}'

echo "e) Máximo y mínimo de cada columna de matrix.txt"

awk '
{
    for (i = 1; i <= NF; i++) {
        if (NR == 1 || $i < minimo[i])
            minimo[i] = $i
            
        if (NR == 1 || $i > maximo[i])
            maximo[i] = $i    
    }
}

END {
    for (i = 1; i <= NF; i++)
        print i, minimo[i], maximo[i]
}
' matrix.txt

echo "f) Invertir orden de las columnas de matrix.txt"

awk '
{
    fot (i = NF; i > 1; i--)
        print $i
        
    print $1
}
' matrix.txt

echo "g) Filas con promedio mayor al promedio general"

awk '
FNR = NR {
    for (i = 1; i <= NF; i++) {
        total += $i
        cantidad++
    } next
}
{suma = 0
    for (i = 1; i <= NF; i++)
        suma += $i
        
    if (suma / NF > total / cantidad)
        print
}
' matrix.txt

echo "h) Cantidad de número impares en matrix.txt"

awk '
{for (i = 1; i <= NF; i++) {
    if ($i % 2 != 0)
        impares++
    }}
END {
    print "Cantidad de números impares:", impares
}' matrix.txt

echo "i) Cantidad de números divisibles por 4 y 9"

awk '
{for (i = 1; i <= NF; i++) {
        if ($i % 4 == 0)
            divisibles_4++
        if ($i % 9 == 0)
            divisibles_9++
        if ($i % 4 == 0 && $i % 9 == 0)
            divisibles_ambos++}}
END {
    print "Divisibles por 4:", divisibles_4
    print "Divisibles por 9:", divisibles_9
    print "Divisibles por 4 y 9:", divisibles_ambos
}
' matrix.txt

echo "j) Fila con mayor suma acumulada"

awk '
{   suma_fila = 0
    for (i = 1; i <= NF; i++)
        suma_fila += $i
    if (NR == 1 || suma_fila > suma_maxima) {
        suma_maxima = suma_fila
        numero_fila = NR
        fila_maxima = $0
    }}
END {
    print "Numero de fila:", numero_fila
    print "Suma maxima:", suma_maxima
    print "Contenidas de la fila:", fila_maxima
}
' matrix.txt


