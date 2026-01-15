#include <iostream>
#include <vector>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <cfloat>
#include <iomanip>
#include <ctime>
#include <sstream>

// ROOT Includes
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TSystem.h>

// Constantes
const int NUM_CHANNELS = 32;
const int SAMPLES = 40;
const float BIN_WIDTH = 8.0;

// Estructuras de Datos
class Pulse {
public:
    float f_max;    
    float f_t0;     
    float f_int;    
    double f_acqTime; 

    void Reset() {
        f_max = 0; f_t0 = 0; f_int = 0; f_acqTime = 0;
    }
};

// Convierte microsegundos a formato legible
std::string FormatearTiempoGlobal(double total_us) {
    time_t segundos = static_cast<time_t>(total_us / 1e6);
    long microsegundos = static_cast<long>(std::fmod(total_us, 1e6));
    struct tm *timeinfo = std::localtime(&segundos);
    
    char buffer[80];
    std::strftime(buffer, sizeof(buffer), "%Y-%m-%d %H:%M:%S", timeinfo);
    
    std::stringstream ss;
    ss << buffer << "." << std::setfill('0') << std::setw(6) << microsegundos;
    return ss.str();
}

// --- Bloques de Procesamiento ---

void ProcesarPSA(std::string inputFile, std::string outputFile) {
    // TFile::Open es más robusto que new TFile
    TFile *fIn = TFile::Open(inputFile.c_str(), "READ");
    if (!fIn || fIn->IsZombie()) {
        std::cerr << "Error al abrir el archivo de entrada o archivo corrupto: " << inputFile << std::endl;
        return;
    }
    
    TTree *events = (TTree*)fIn->Get("events");
    if (!events) {
        std::cerr << "Error: No se encontró el árbol 'events' en " << inputFile << std::endl;
        fIn->Close();
        return;
    }

    // IMPORTANTE: ULong64_t evita el Segmentation Violation en SetBranchAddress
    ULong64_t event_time_us;
    unsigned short waveforms[NUM_CHANNELS][SAMPLES];
    
    events->SetBranchAddress("event_time_us", &event_time_us);
    events->SetBranchAddress("waveforms", waveforms);

    TFile *fOut = new TFile(outputFile.c_str(), "RECREATE");
    TTree* pTree[NUM_CHANNELS];
    TH2F* hPersistencia[NUM_CHANNELS];
    Pulse pData[NUM_CHANNELS];

    for (int ch = 0; ch < NUM_CHANNELS; ++ch) {
        pTree[ch] = new TTree(Form("p%d", ch), Form("Análisis Canal %d", ch));
        pTree[ch]->Branch("f_max", &pData[ch].f_max, "f_max/F");
        pTree[ch]->Branch("f_t0", &pData[ch].f_t0, "f_t0/F");
        pTree[ch]->Branch("f_int", &pData[ch].f_int, "f_int/F");
        pTree[ch]->Branch("f_acqTime", &pData[ch].f_acqTime, "f_acqTime/D");

        hPersistencia[ch] = new TH2F(Form("hPersist_Ch%d", ch), 
                                     Form("Persistencia Canal %d;Time (ns);Amplitude (ADC)", ch),
                                     SAMPLES, 0, SAMPLES * BIN_WIDTH, 
                                     550, -50, 1500);
    }

    Long64_t nEntries = events->GetEntries();
    if (nEntries == 0) {
        std::cerr << "Archivo vacío." << std::endl;
        fOut->Close(); fIn->Close(); return;
    }

    events->GetEntry(0);
    double startTime = static_cast<double>(event_time_us);
    events->GetEntry(nEntries - 1);
    double endTime = static_cast<double>(event_time_us);

    std::cout << "Procesando " << nEntries << " eventos de " << inputFile << "..." << std::endl;

    for (Long64_t i = 0; i < nEntries; ++i) {
        events->GetEntry(i);
        for (int ch = 0; ch < NUM_CHANNELS; ++ch) {
            pData[ch].Reset();
            float baseline = 0;
            for (int s = 0; s < 5; ++s) baseline += waveforms[ch][s];
            baseline /= 5.0;

            float current_max = -FLT_MAX;
            float current_int = 0;
            int peak_sample = 0;

            for (int s = 0; s < SAMPLES; ++s) {
                float val = -(static_cast<float>(waveforms[ch][s]) - baseline);
                hPersistencia[ch]->Fill(s * BIN_WIDTH, val);
                if (val > current_max) {
                    current_max = val;
                    peak_sample = s;
                }
                current_int += val;
            }
            pData[ch].f_max = current_max;
            pData[ch].f_int = current_int * BIN_WIDTH;
            pData[ch].f_t0  = peak_sample * BIN_WIDTH;
            pData[ch].f_acqTime = static_cast<double>(event_time_us);
            pTree[ch]->Fill();
        }
    }

    std::cout << "\n===========================================" << std::endl;
    std::cout << "INFORMACIÓN DE LA ADQUISICIÓN (TIEMPO REAL)" << std::endl;
    std::cout << "===========================================" << std::endl;
    std::cout << "Archivo de salida: " << outputFile << std::endl;
    std::cout << "Fecha/Hora Inicio: " << FormatearTiempoGlobal(startTime) << std::endl;
    std::cout << "Fecha/Hora Fin:    " << FormatearTiempoGlobal(endTime) << std::endl;
    std::cout << "-------------------------------------------" << std::endl;
    std::cout << "Duración total:    " << (endTime - startTime) / 1e6 << " s" << std::endl;
    std::cout << "===========================================\n" << std::endl;

    fOut->cd();
    for (int ch = 0; ch < NUM_CHANNELS; ++ch) {
        pTree[ch]->Write();
        hPersistencia[ch]->Write();
    }
    fOut->Close();
    fIn->Close();
}

void CalcularRate(std::string fileName, int canal_objetivo = 0, float threshold = 150.0) {
    TFile *f = TFile::Open(fileName.c_str(), "UPDATE");
    if (!f || f->IsZombie()) return;

    TTree *t = (TTree*)f->Get(Form("p%d", canal_objetivo));
    if(!t) { f->Close(); return; }

    Pulse p;
    t->SetBranchAddress("f_max", &p.f_max);
    t->SetBranchAddress("f_acqTime", &p.f_acqTime);

    Long64_t nEntries = t->GetEntries();
    if (nEntries < 2) { f->Close(); return; }

    t->GetEntry(0);
    double t_start = p.f_acqTime;
    t->GetEntry(nEntries - 1);
    double t_end = p.f_acqTime;

    double duracion_total_s = (t_end - t_start) / 1e6;
    if (duracion_total_s <= 0) duracion_total_s = 1.0; 
    
    int nBins10 = std::max(1, (int)std::ceil(duracion_total_s / 600.0));
    TH1F *hRate10 = new TH1F("hRate10", Form("Rate Ch%d (10 min bins);Tiempo (min);Rate (Hz)", canal_objetivo), nBins10, 0, nBins10 * 10);

    for (Long64_t i = 0; i < nEntries; ++i) {
        t->GetEntry(i);
        if (p.f_max > threshold) {
            double tiempo_relativo_min = ((p.f_acqTime - t_start) / 1e6) / 60.0;
            hRate10->Fill(tiempo_relativo_min);
        }
    }

    for (int b = 1; b <= nBins10; ++b) {
        double conteo = hRate10->GetBinContent(b);
        hRate10->SetBinContent(b, conteo / 600.0); 
        hRate10->SetBinError(b, std::sqrt(conteo) / 600.0);
    }

    hRate10->SetMarkerStyle(20); 
    hRate10->SetMarkerColor(kBlue+2);
    hRate10->Write("", TObject::kOverwrite);
    f->Close();
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cout << "Uso: ./analyze_waveforms <ruta/al/archivo.root>" << std::endl;
        return 1;
    }

    std::string fullPath = argv[1];
    
    // Extraer solo el nombre del archivo para evitar errores de carpetas inexistentes
    // Ejemplo: testLargeHodoscope/20251121/acq_XYZ.root -> acq_XYZ.root
    size_t lastSlash = fullPath.find_last_of("/\\");
    std::string fileName = (lastSlash == std::string::npos) ? fullPath : fullPath.substr(lastSlash + 1);

    // Evitar procesar archivos que ya son de salida
    if (fileName.find("analisis_pulsos_") != std::string::npos) {
        std::cout << "Saltando archivo de salida: " << fileName << std::endl;
        return 0;
    }

    std::string outputFile = "analisis_pulsos_" + fileName;

    ProcesarPSA(fullPath, outputFile);
    CalcularRate(outputFile, 0, 150.0); 

    return 0;
}