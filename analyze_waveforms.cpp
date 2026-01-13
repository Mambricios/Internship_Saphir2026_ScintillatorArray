#include <iostream>
#include <vector>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <cfloat>
#include <iomanip>
#include <ctime>

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>

// --- Función mejorada para formato jerárquico Año-Mes-Día Hora:Min:Seg.us ---
std::string FormatearTiempoGlobal(double total_us) {
    // 1. Extraer segundos y microsegundos
    time_t segundos = static_cast<time_t>(total_us / 1e6);
    long microsegundos = static_cast<long>(std::fmod(total_us, 1e6));

    // 2. Convertir segundos a estructura de tiempo local (año, mes, día, etc.)
    struct tm *timeinfo = std::localtime(&segundos);

    // 3. Crear el buffer para la fecha/hora hasta los segundos
    char buffer[80];
    // Formato: %Y (Año), %m (Mes), %d (Día), %H:%M:%S (Hora:Min:Seg)
    std::strftime(buffer, sizeof(buffer), "%Y-%m-%d %H:%M:%S", timeinfo);

    // 4. Unir todo en un stream para añadir los microsegundos con precisión
    std::stringstream ss;
    ss << buffer << "." << std::setfill('0') << std::setw(6) << microsegundos;
    
    return ss.str();
}

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

const int NUM_CHANNELS = 32;
const int SAMPLES = 40;
const float BIN_WIDTH = 8.0;

void ProcesarPSA() {
    TFile *fIn = new TFile("acq_20251231_192759.root", "READ");
    if (!fIn || fIn->IsZombie()) return;
    
    TTree *events = (TTree*)fIn->Get("events");
    
    unsigned long long event_time_us;
    unsigned short waveforms[NUM_CHANNELS][SAMPLES];
    
    events->SetBranchAddress("event_time_us", &event_time_us);
    events->SetBranchAddress("waveforms", waveforms);

    TFile *fOut = new TFile("analisis_pulsos.root", "RECREATE");
    
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
                                     550, -50, 500);
    }

    Long64_t nEntries = events->GetEntries();
    
    events->GetEntry(0);
    double startTime = static_cast<double>(event_time_us);
    
    events->GetEntry(nEntries - 1);
    double endTime = static_cast<double>(event_time_us);

    std::cout << "Procesando " << nEntries << " eventos..." << std::endl;

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

    // --- SALIDA ORDENADA JERÁRQUICAMENTE ---
    std::cout << "\n===========================================" << std::endl;
    std::cout << "INFORMACIÓN DE LA ADQUISICIÓN (TIEMPO REAL)" << std::endl;
    std::cout << "===========================================" << std::endl;
    std::cout << "Fecha/Hora Inicio: " << FormatearTiempoGlobal(startTime) << std::endl;
    std::cout << "Fecha/Hora Fin:    " << FormatearTiempoGlobal(endTime) << std::endl;
    std::cout << "-------------------------------------------" << std::endl;
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Duración total:    " << (endTime - startTime) / 1e6 << " s" << std::endl;
    std::cout << "===========================================\n" << std::endl;

    fOut->cd();
    for (int ch = 0; ch < NUM_CHANNELS; ++ch) {
        pTree[ch]->Write();
        hPersistencia[ch]->SetMarkerStyle(1);
        hPersistencia[ch]->SetMarkerColor(kBlack);
        hPersistencia[ch]->SetOption("SCAT");
        hPersistencia[ch]->Write();
    }
    
    fOut->Close();
    fIn->Close();
}

int main() {
    ProcesarPSA();
    return 0;
}