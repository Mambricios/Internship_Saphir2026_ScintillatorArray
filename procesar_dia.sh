#!/bin/bash

# Comprobar si se pasó el argumento de la fecha
if [ -z "$1" ]; then
    echo "Uso: ./procesar_dia.sh AAAAMMDD"
    exit 1
fi

DIA=$1
# RUTA ABSOLUTA DE LOS DATOS (Solo lectura)
FOLDER_DATA="/data/user/mayala/testLargeHodoscope/$DIA"
# RUTA EN TU HOME DE USUARIO
FOLDER_OUT="/user/santibanez/analisis_$DIA"

# 1. Verificar si la carpeta de datos originales existe
if [ ! -d "$FOLDER_DATA" ]; then
    echo "Error: La carpeta de datos $FOLDER_DATA no existe."
    exit 1
fi

# 2. Crear la carpeta de salida en tu usuario
mkdir -p "$FOLDER_OUT"

echo "===================================================="
echo "Iniciando procesamiento en Clúster: $DIA"
echo "Leyendo de: $FOLDER_DATA"
echo "Guardando en: $FOLDER_OUT"
echo "===================================================="

count=0
for file in "$FOLDER_DATA"/acq_*.root; do
    [ -e "$file" ] || continue
    
    count=$((count + 1))
    FILENAME=$(basename "$file")
    echo "[$count] Procesando: $FILENAME"
    
    # BUENA PRÁCTICA: Copiar a TMPDIR para procesamiento rápido [cite: 229]
    cp "$file" "$TMPDIR/$FILENAME"
    
    # Ejecutar el binario (asegúrate de que esté en tu carpeta bin)
    ./bin/analyze_auto "$TMPDIR/$FILENAME"
    
    # Mover el resultado a tu carpeta de usuario
    # El binario genera 'analisis_pulsos_acq_XXXX.root' en el directorio actual de ejecución
    mv analisis_pulsos_*.root "$FOLDER_OUT/" 2>/dev/null
    
    # Limpiar el archivo temporal para no llenar el nodo [cite: 229]
    rm "$TMPDIR/$FILENAME"
done

echo "===================================================="
echo "PROCESO FINALIZADO. Archivos en: $FOLDER_OUT/"
echo "===================================================="
