#include <iostream>
#include <vector>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <cfloat>
#include <iomanip>
#include <cstdlib>
#include <string>

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TSpectrum.h>
#include <TF1.h>

// --- Constantes Globales ---
const int NUM_CHANNELS = 32;
const int SAMPLES = 40; 
const float BIN_WIDTH = 8.0;

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cout << "Uso: ./analisis <archivo_raw.root> <canal_excluido>" << std::endl;
        return 1;
    }

    std::string inputPath = argv[1];
    int ch_excluido = std::atoi(argv[2]);
    float trigger_threshold = 150.0; 

    std::string fileNameOnly = inputPath;
    size_t lastSlash = fileNameOnly.find_last_of("/\\");
    if (lastSlash != std::string::npos) fileNameOnly = fileNameOnly.substr(lastSlash + 1);
    
    std::string datePart = "output"; 
    size_t firstUnderscore = fileNameOnly.find("_");
    size_t secondUnderscore = fileNameOnly.find("_", firstUnderscore + 1);
    
    if (firstUnderscore != std::string::npos && secondUnderscore != std::string::npos) {
        datePart = fileNameOnly.substr(firstUnderscore + 1, secondUnderscore - firstUnderscore - 1);
    }
    
    std::string outputFileName = datePart + "_ch" + std::to_string(ch_excluido) + ".root";

    std::vector<int> canales_trigger = {0, 26, 1, 27, 4, 30, 5, 31, 6, 28, 7, 29};

    TFile *fIn = TFile::Open(inputPath.c_str(), "READ");
    if (!fIn || fIn->IsZombie()) return 1;
    TTree *events = (TTree*)fIn->Get("events");
    
    unsigned short waveforms[NUM_CHANNELS][SAMPLES];
    events->SetBranchAddress("waveforms", waveforms);

    TFile *fOut = new TFile(outputFileName.c_str(), "RECREATE");

    TH1F* hChargeDist = new TH1F("hChargeDist_Excluido", 
                                 Form("Distribucion de Carga Canal %d;ADC*ns;Entries", ch_excluido), 
                                 100, 0, 3000);

    TH1F* hMaxAmplitude = new TH1F("hMaxAmp_Excluido", 
                                   Form("Amplitud Maxima SPE Canal %d;ADC;Entries", ch_excluido), 
                                   100, 0, 800);

    TH2F* hPersistencia = new TH2F("hPersistencia_Excluido", 
                                   Form("Persistencia Canal %d;Time (ns);ADC", ch_excluido),
                                   SAMPLES, 0, SAMPLES * BIN_WIDTH, 400, -50, 800);

    TH1F* hTotalEventCharge = new TH1F("hTotalEventCharge_Excluido",
                                       Form("Carga Total Acumulada por Evento Canal %d;ADC*ns;Entries", ch_excluido),
                                       100, 0, 5000);

    std::vector<float> event_sums_vector;

    Long64_t nEntries = events->GetEntries();

    for (Long64_t i = 0; i < nEntries; ++i) {
        events->GetEntry(i);

        bool condicion_disparo = false;
        for (int ch_t : canales_trigger) {
            float bl = 0; for (int s = 0; s < 5; ++s) bl += waveforms[ch_t][s]; bl /= 5.0;
            float max_val = -FLT_MAX;
            for (int s = 0; s < SAMPLES; ++s) {
                float val = -(static_cast<float>(waveforms[ch_t][s]) - bl);
                if (val > max_val) max_val = val;
            }
            if (max_val > trigger_threshold) { condicion_disparo = true; break; }
        }

        if (condicion_disparo) {
            float baseline = 0; for (int s = 0; s < 5; ++s) baseline += waveforms[ch_excluido][s]; baseline /= 5.0;
            float rms = 0;
            for (int s = 0; s < 5; ++s) {
                float diff = waveforms[ch_excluido][s] - baseline;
                rms += diff * diff;
            }
            float sigma = std::sqrt(rms / 5.0);
            if (sigma < 1.0) sigma = 1.0; 

            std::vector<float> data_adj(SAMPLES);
            for (int s = 0; s < SAMPLES; ++s) {
                data_adj[s] = -(static_cast<float>(waveforms[ch_excluido][s]) - baseline);
                hPersistencia->Fill(s * BIN_WIDTH, data_adj[s]);
            }

            float total_event_sum = 0.0; 

            for (int s = 2; s < SAMPLES - 5; ++s) {
                float amp = data_adj[s];
                if (amp >= 5 * sigma && amp <= 20 * sigma) {
                    float sum_charge = 0.0;
                    float max_amp_peak = 0.0;
                    for (int j = s - 2; j <= s + 4; ++j) {
                        sum_charge += data_adj[j];
                        if (data_adj[j] > max_amp_peak) max_amp_peak = data_adj[j];
                    }
                    float pulse_charge = sum_charge * BIN_WIDTH;
                    hChargeDist->Fill(pulse_charge);
                    hMaxAmplitude->Fill(max_amp_peak);
                    total_event_sum += pulse_charge; 
                    s += 4; 
                }
            }
            if (total_event_sum > 0) {
                hTotalEventCharge->Fill(total_event_sum);
                event_sums_vector.push_back(total_event_sum); 
            }
        }
    }

    // --- Cálculo Automático de MU (SPE) Forzado entre 500 y 700 ---
    float mu = 0;
    // Forzamos el rango del histograma antes de buscar el pico
    hTotalEventCharge->GetXaxis()->SetRangeUser(500, 700);
    
    TSpectrum *s_spec = new TSpectrum(1);
    // Buscamos picos solo en el rango visible
    Int_t nfound = s_spec->Search(hTotalEventCharge, 1, "nobackground", 0.1); 
    
    if (nfound > 0) {
        Double_t *xpeaks = s_spec->GetPositionX();
        Double_t peak_x = xpeaks[0];
        
        // Ajustamos la Gaussiana específicamente en la zona del pico detectado
        TF1 *fGaus = new TF1("fGaus", "gaus", peak_x - 200, peak_x + 200);
        hTotalEventCharge->Fit(fGaus, "RQ"); // Q para modo silencioso, R para usar el rango del TF1
        mu = fGaus->GetParameter(1);

        std::cout << "\n==========================================" << std::endl;
        std::cout << " CALIBRACIÓN SPE FORZADA (500-700 ADC*ns)" << std::endl;
        std::cout << " Canal:            " << ch_excluido << std::endl;
        std::cout << " mu (Carga 1 spe): " << std::fixed << std::setprecision(2) << mu << " ADC*ns" << std::endl;
        std::cout << "==========================================\n" << std::endl;

        // Resetear el rango para guardar el histograma completo pero procesar la normalización
        hTotalEventCharge->GetXaxis()->SetRange(0, 0); 

        TH1F* hNormalizedEventCharge = new TH1F("hNormalizedEventCharge_Excluido",
                                                Form("Carga Evento Normalizada Canal %d;No. Photoelectrons (Carga/mu);Entries", ch_excluido),
                                                100, 0, 10);
        for (float val : event_sums_vector) {
            hNormalizedEventCharge->Fill(val / mu);
        }
        hNormalizedEventCharge->Write();
    } else {
        std::cout << "\n[!] No se encontró un pico claro entre 500 y 700 ADC*ns." << std::endl;
        hTotalEventCharge->GetXaxis()->SetRange(0, 0); // Resetear rango si falla
    }

    hChargeDist->Write();
    hMaxAmplitude->Write();
    hPersistencia->Write();
    hTotalEventCharge->Write();

    fOut->Close();
    fIn->Close();

    std::cout << "Analisis SPE completado. Resultados en: " << outputFileName << std::endl;
    return 0;
}