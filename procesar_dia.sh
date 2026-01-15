#!/bin/bash

# Comprobar si se pasó el argumento de la fecha
if [ -z "$1" ]; then
    echo "Uso: ./procesar_dia.sh AAAAMMDD"
    echo "Ejemplo: ./procesar_dia.sh 20260113"
    exit 1
fi

DIA=$1
FOLDER_DATA="testLargeHodoscope/$DIA"
FOLDER_OUT="analisis_$DIA"

# 1. Verificar si la carpeta de datos originales existe
if [ ! -d "$FOLDER_DATA" ]; then
    echo "Error: La carpeta de datos $FOLDER_DATA no existe."
    exit 1
fi

# 2. Crear la carpeta de salida si no existe
if [ ! -d "$FOLDER_OUT" ]; then
    echo "Creando carpeta de destino: $FOLDER_OUT"
    mkdir -p "$FOLDER_OUT"
fi

echo "===================================================="
echo "Iniciando procesamiento del día: $DIA"
echo "Origen: $FOLDER_DATA"
echo "Destino: $FOLDER_OUT"
echo "===================================================="

# 3. Iterar sobre archivos de adquisición
# Usamos un contador para informar el progreso
count=0
for file in "$FOLDER_DATA"/acq_*.root; do
    # Verificar si existen archivos para evitar errores de loop vacío
    [ -e "$file" ] || continue
    
    count=$((count + 1))
    echo ""
    echo "[$count] Procesando: $(basename "$file")"
    
    # Ejecutar el binario automático
    ./bin/analyze_auto "$file"
    
    # El binario genera un archivo llamado 'analisis_pulsos_acq_XXXX.root' en la raíz.
    # Lo movemos inmediatamente a la carpeta de destino.
    mv analisis_pulsos_*.root "$FOLDER_OUT/" 2>/dev/null
done

echo ""
echo "===================================================="
echo "PROCESO FINALIZADO"
echo "Total de archivos procesados: $count"
echo "Los resultados están en: $FOLDER_OUT/"
echo "===================================================="