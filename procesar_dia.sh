#!/bin/bash
FECHA=$1

# 1. Validación de entrada
if [ -z "$FECHA" ]; then 
    echo "Uso: ./procesar_dia.sh 20251122"
    exit 1
fi

# 2. Preparación
DIR_SALIDA="analisis_$FECHA"
mkdir -p "$DIR_SALIDA"

echo "[1/3] Compilando analizador..."
g++ analyze_waveforms.cpp -o analyze `root-config --cflags --libs`

# 3. Procesamiento individual
echo "[2/3] Procesando archivos individuales de la carpeta $FECHA..."
for raw in "$FECHA"/acq_*.root; do
    if [ -f "$raw" ]; then
        base=$(basename "$raw" .root)
        # Cambiamos acq_ por analisis_ para el nombre de salida
        salida="$DIR_SALIDA/${base/acq/analisis}.root"
        echo "Analizando: $raw -> $salida"
        ./analyze "$raw" "$salida"
    fi
done

# 4. Unión de archivos (hadd)
echo "[3/3] Uniendo todos los archivos procesados en uno solo..."
ARCHIVO_FINAL="total_$FECHA.root"

# -f obliga a sobrescribir si el archivo ya existe
hadd -f "$DIR_SALIDA/$ARCHIVO_FINAL" "$DIR_SALIDA"/analisis_*.root

echo "------------------------------------------------"
echo "✅ Proceso completado con éxito."
echo "📂 Resultados en: $DIR_SALIDA/"
echo "📊 Archivo unificado: $DIR_SALIDA/$ARCHIVO_FINAL"
echo "------------------------------------------------"
