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
#include <TCanvas.h>

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

// Convierte microsegundos a formato: Año-Mes-Día Hora:Min:Seg.us
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
                                     550, -50, 1500);
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

    std::cout << "\n===========================================" << std::endl;
    std::cout << "INFORMACIÓN DE LA ADQUISICIÓN (TIEMPO REAL)" << std::endl;
    std::cout << "===========================================" << std::endl;
    std::cout << "Fecha/Hora Inicio: " << FormatearTiempoGlobal(startTime) << std::endl;
    std::cout << "Fecha/Hora Fin:    " << FormatearTiempoGlobal(endTime) << std::endl;
    std::cout << "-------------------------------------------" << std::endl;
    std::cout << "Duración total:    " << (endTime - startTime) / 1e6 << " s" << std::endl;
    std::cout << "===========================================\n" << std::endl;

    // Histograma de Persistencia
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

// Rate of Events
void CalcularRate(int canal_objetivo = 0, float threshold = 150.0) {
    TFile *f = new TFile("analisis_pulsos.root", "UPDATE");
    if (!f || f->IsZombie()) return;

    TTree *t = (TTree*)f->Get(Form("p%d", canal_objetivo));
    Pulse p;
    t->SetBranchAddress("f_max", &p.f_max);
    t->SetBranchAddress("f_acqTime", &p.f_acqTime);

    Long64_t nEntries = t->GetEntries();
    t->GetEntry(0);
    double t_start = p.f_acqTime;
    t->GetEntry(nEntries - 1);
    double t_end = p.f_acqTime;

    double duracion_total_s = (t_end - t_start) / 1e6;
    
    // Histograma 10 min
    int nBins10 = std::ceil(duracion_total_s / 600.0);
    TH1F *hRate10 = new TH1F("hRate10", Form("Rate Ch%d (10 min bins);Tiempo (min);Rate (Hz)", canal_objetivo), nBins10, 0, nBins10 * 10);

    // Histograma 5 min
    int nBins5 = std::ceil(duracion_total_s / 300.0);
    TH1F *hRate5 = new TH1F("hRate5", Form("Rate Ch%d (5 min bins);Tiempo (min);Rate (Hz)", canal_objetivo), nBins5, 0, nBins5 * 5);

    for (Long64_t i = 0; i < nEntries; ++i) {
        t->GetEntry(i);
        if (p.f_max > threshold) {
            double tiempo_relativo_min = ((p.f_acqTime - t_start) / 1e6) / 60.0;
            hRate10->Fill(tiempo_relativo_min);
            hRate5->Fill(tiempo_relativo_min);
        }
    }

    // Normalización y cálculo de error (10 min)
    for (int b = 1; b <= nBins10; ++b) {
        double conteo = hRate10->GetBinContent(b);
        hRate10->SetBinContent(b, conteo / 600.0); 
        hRate10->SetBinError(b, std::sqrt(conteo) / 600.0);
    }

    // Normalización y cálculo de error (5 min)
    for (int b = 1; b <= nBins5; ++b) {
        double conteo = hRate5->GetBinContent(b);
        hRate5->SetBinContent(b, conteo / 300.0); 
        hRate5->SetBinError(b, std::sqrt(conteo) / 300.0);
    }

    hRate10->SetMarkerStyle(20); hRate10->SetMarkerColor(kBlue+2);
    hRate5->SetMarkerStyle(21); hRate5->SetMarkerColor(kRed+2);
    
    hRate10->Write("", TObject::kOverwrite);
    hRate5->Write("", TObject::kOverwrite);
    
    f->Close();
}

// Función para generar la correlación 
void GenerarCorrelacionCanales(int canal_objetivo1 = 0, int canal_objetivo2 = 1) {
    TFile *f = new TFile("analisis_pulsos.root", "UPDATE");
    if (!f || f->IsZombie()) return;

    TTree *t1 = (TTree*)f->Get(Form("p%d", canal_objetivo1));
    TTree *t2 = (TTree*)f->Get(Form("p%d", canal_objetivo2));

    if (!t1 || !t2) {
        std::cout << "Error: No se encontraron los árboles de los canales seleccionados." << std::endl;
        f->Close();
        return;
    }

    // Variables para f_int
    float f_int1, f_int2;
    t1->SetBranchAddress("f_int", &f_int1);
    t2->SetBranchAddress("f_int", &f_int2);

    // Variables para f_max
    float f_max1, f_max2;
    t1->SetBranchAddress("f_max", &f_max1);
    t2->SetBranchAddress("f_max", &f_max2);

    // Histograma de correlación para f_int
    TH2F *hCorrInt = new TH2F(Form("hCorr_Int_Ch%d_Ch%d", canal_objetivo1, canal_objetivo2), 
                              Form("Correlacion Integral: Ch%d vs Ch%d;Integral Ch%d;Integral Ch%d", canal_objetivo1, canal_objetivo2, canal_objetivo1, canal_objetivo2), 
                              200, 0, 50000, 
                              200, 0, 50000);

    // Histograma de correlación para f_max
    TH2F *hCorrMax = new TH2F(Form("hCorr_Max_Ch%d_Ch%d", canal_objetivo1, canal_objetivo2), 
                              Form("Correlacion Amplitud Max: Ch%d vs Ch%d;f_max Ch%d (ADC);f_max Ch%d (ADC)", canal_objetivo1, canal_objetivo2, canal_objetivo1, canal_objetivo2), 
                              200, 0, 1500, 
                              200, 0, 1500);

    Long64_t nEntries = t1->GetEntries();
    
    for (Long64_t i = 0; i < nEntries; ++i) {
        t1->GetEntry(i);
        t2->GetEntry(i);
        
        hCorrInt->Fill(f_int1, f_int2);
        hCorrMax->Fill(f_max1, f_max2);
    }

    // Guardar ambos con escala de colores
    hCorrInt->SetOption("P");
    hCorrMax->SetOption("P");
    
    hCorrInt->Write("", TObject::kOverwrite);
    hCorrMax->Write("", TObject::kOverwrite);
    
    f->Close();
}

// Ejecución actualizada
int main() {
    ProcesarPSA();
    CalcularRate(0, 150.0);
    //GenerarCorrelacionCanales(0, 25); 
}

