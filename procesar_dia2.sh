#!/bin/bash
FECHA=$1

if [ -z "$FECHA" ]; then 
    echo "Uso: ./procesar_dia.sh 20251122"
    exit 1
fi

DIR_SALIDA="analisis_$FECHA"
mkdir -p "$DIR_SALIDA"

echo "[1/3] Compilando analizador..."
# Compilamos con un nombre claro
g++ analyze_waveforms2.cpp -o analyze_exe `root-config --cflags --libs`

echo "[2/3] Procesando archivos de la carpeta $FECHA..."
for raw in "$FECHA"/acq_*.root; do
    if [ -f "$raw" ]; then
        # Saltamos archivos vacíos (error de cuota previo)
        if [ ! -s "$raw" ]; then
            echo "Saltando $raw (archivo vacío)"
            continue
        fi

        base=$(basename "$raw" .root)
        salida="$DIR_SALIDA/${base/acq/analisis}.root"
        
        echo "Analizando: $raw"
        # Usamos el nombre correcto del ejecutable
        ./analyze_exe "$raw" "$salida"
    fi
done

echo "[3/3] Uniendo archivos con hadd..."
ARCHIVO_FINAL="total_$FECHA.root"

# Unimos todo en el directorio actual o en el de salida
hadd -f "$ARCHIVO_FINAL" "$DIR_SALIDA"/analisis_*.root

echo "------------------------------------------------"
echo "¡Listo! Archivo final: $ARCHIVO_FINAL"
echo "------------------------------------------------"
