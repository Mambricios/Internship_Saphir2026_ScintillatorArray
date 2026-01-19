#include <iostream>
#include <vector>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <cfloat>
#include <iomanip>
#include <ctime>
#include <sstream>

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TAxis.h>

// --- Constantes Globales ---
const int NUM_CHANNELS = 32;
const int SAMPLES = 40;
const float BIN_WIDTH = 8.0;

// --- Estructuras de Datos ---
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

// --- Tiempo Universal a Hora local ---
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

// --- Procesamiento Base ---
void ProcesarPSA(const char* inputFile, const char* outputFile) {
    TFile *fIn = TFile::Open(inputFile, "READ");
    if (!fIn || fIn->IsZombie()) {
        std::cerr << "[!] Error: No se pudo abrir " << inputFile << std::endl;
        return;
    }
    
    TTree *events = (TTree*)fIn->Get("events");
    if (!events) {
        std::cerr << "[!] Error: No se encontró el árbol 'events' en " << inputFile << std::endl;
        fIn->Close();
        return;
    }

    unsigned long long event_time_us;
    unsigned short waveforms[NUM_CHANNELS][SAMPLES];
    events->SetBranchAddress("event_time_us", &event_time_us);
    events->SetBranchAddress("waveforms", waveforms);

    TFile *fOut = new TFile(outputFile, "RECREATE");
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
                                     550, -50, 2000);
    }

    Long64_t nEntries = events->GetEntries();
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
    std::cout << "Procesamiento listo, Guardado en: " << outputFile << std::endl;
}

// --- Cálculo de Rate ---
void CalcularRate(const char* targetFile, int canal_objetivo, float threshold) {
    TFile *f = TFile::Open(targetFile, "UPDATE");
    if (!f || f->IsZombie()) return;

    TTree *t = (TTree*)f->Get(Form("p%d", canal_objetivo));
    if(!t) { f->Close(); return; }

    Pulse p;
    t->SetBranchAddress("f_max", &p.f_max);
    t->SetBranchAddress("f_acqTime", &p.f_acqTime);

    Long64_t nEntries = t->GetEntries();
    if (nEntries < 2) { f->Close(); return; }

    t->GetEntry(0);
    double t_start_sec = p.f_acqTime / 1e6; 
    t->GetEntry(nEntries - 1);
    double t_end_sec = p.f_acqTime / 1e6;
    
    double duracion_total_s = t_end_sec - t_start_sec;
    int nBins = std::ceil(duracion_total_s / 600.0);

    TH1F *hRateGlobal = new TH1F(Form("hRateGlobal_Ch%d", canal_objetivo), 
                                 Form("Rate Canal %d;Hora del dia;Rate (Hz)", canal_objetivo), 
                                 nBins, t_start_sec, t_end_sec);

    for (Long64_t i = 0; i < nEntries; ++i) {
        t->GetEntry(i);
        if (p.f_max > threshold) {
            hRateGlobal->Fill(p.f_acqTime / 1e6);
        }
    }

    hRateGlobal->GetXaxis()->SetTimeDisplay(1);
    hRateGlobal->GetXaxis()->SetTimeFormat("%H:%M");
    hRateGlobal->GetXaxis()->SetTimeOffset(0, "gmt");

    hRateGlobal->Scale(1.0 / 600.0);
    hRateGlobal->Write("", TObject::kOverwrite);
    f->Close();
    std::cout << "[OK] Histograma de Rate generado con eje de tiempo real." << std::endl;
}

// --- Correlación entre 4 canales ---
void GenerarCorrelacionCanales(const char* targetFile, std::vector<int> chs) {
    TFile *f = TFile::Open(targetFile, "UPDATE");
    if (!f || f->IsZombie()) return;

    TCanvas *cCorr = new TCanvas("cCorrelaciones", "Correlaciones Múltiples", 1200, 800);
    cCorr->Divide(3, 2);
    int padCount = 1;
    double threshold = 150.0;

    for (size_t i = 0; i < chs.size(); ++i) {
        for (size_t j = i + 1; j < chs.size(); ++j) {
            int c1 = chs[i];
            int c2 = chs[j];

            TTree *t1 = (TTree*)f->Get(Form("p%d", c1));
            TTree *t2 = (TTree*)f->Get(Form("p%d", c2));

            if (!t1 || !t2) continue;

            float max1, max2;
            t1->SetBranchAddress("f_max", &max1);
            t2->SetBranchAddress("f_max", &max2);

            cCorr->cd(padCount++);
            TH2F *hCorr = new TH2F(Form("hCorr_Ch%d_Ch%d", c1, c2), 
                                   Form("Ch%d vs Ch%d;Ch%d Max;Ch%d Max", c1, c2, c1, c2),
                                   150, 0, 1000, 150, 0, 1000);

            for (Long64_t n = 0; n < t1->GetEntries(); ++n) {
                t1->GetEntry(n);
                t2->GetEntry(n);

                if (max1 > threshold && max2 > threshold) {
                    hCorr->Fill(max1, max2);
                }

            }
            hCorr->Draw("COLZ");
            hCorr->Write("", TObject::kOverwrite);
        }
    }

    cCorr->Write();
    f->Close();
    std::cout << "[OK] Canvas de correlaciones guardado en el archivo." << std::endl;
}

// --- Rate por Coincidencia (AND) ---
void CalcularRateCoincidencia(const char* targetFile, int c1, int c2, float threshold) {
    TFile *f = TFile::Open(targetFile, "UPDATE");
    if (!f || f->IsZombie()) return;

    TTree *t1 = (TTree*)f->Get(Form("p%d", c1));
    TTree *t2 = (TTree*)f->Get(Form("p%d", c2));
    if(!t1 || !t2) { f->Close(); return; }

    float max1, max2;
    double time1; 
    t1->SetBranchAddress("f_max", &max1);
    t1->SetBranchAddress("f_acqTime", &time1);
    t2->SetBranchAddress("f_max", &max2);

    Long64_t nEntries = t1->GetEntries();
    if (nEntries < 2) { f->Close(); return; }

    t1->GetEntry(0);
    double t_start_sec = time1 / 1e6;
    t1->GetEntry(nEntries - 1);
    double t_end_sec = time1 / 1e6;

    double duracion_total_s = t_end_sec - t_start_sec;
    int nBins = std::ceil(duracion_total_s / 600.0);

    TH1F *hRateCoin = new TH1F(Form("hRateCoin_Ch%d_Ch%d", c1, c2), 
                               Form("Rate Coincidencia Ch%d AND Ch%d;Hora del dia;Rate (Hz)", c1, c2), 
                               nBins, t_start_sec, t_end_sec);

    for (Long64_t i = 0; i < nEntries; ++i) {
        t1->GetEntry(i);
        t2->GetEntry(i);
        if (max1 > threshold && max2 > threshold) {
            hRateCoin->Fill(time1 / 1e6);
        }
    }

    hRateCoin->GetXaxis()->SetTimeDisplay(1);
    hRateCoin->GetXaxis()->SetTimeFormat("%H:%M");
    hRateCoin->GetXaxis()->SetTimeOffset(0, "gmt");
    hRateCoin->Scale(1.0 / 600.0);
    hRateCoin->Write("", TObject::kOverwrite);
    f->Close();
    std::cout << "[OK] Rate por coincidencia generado exitosamente." << std::endl;
}

// --- Rate por Coincidencia Cuádruple (1&2&3&4) ---
void CalcularRateCoincidenciaCuadruple(const char* targetFile, std::vector<int> chs, float threshold) {
    TFile *f = TFile::Open(targetFile, "UPDATE");
    if (!f || f->IsZombie() || chs.size() < 4) return;

    TTree *trees[4];
    float maxs[4];
    double time_ref;

    for(int i = 0; i < 4; ++i) {
        trees[i] = (TTree*)f->Get(Form("p%d", chs[i]));
        if(!trees[i]) { f->Close(); return; }
        trees[i]->SetBranchAddress("f_max", &maxs[i]);
        if(i == 0) trees[i]->SetBranchAddress("f_acqTime", &time_ref);
    }

    Long64_t nEntries = trees[0]->GetEntries();
    trees[0]->GetEntry(0);
    double t_start_sec = time_ref / 1e6;
    trees[0]->GetEntry(nEntries - 1);
    double t_end_sec = time_ref / 1e6;

    double duracion_total_s = t_end_sec - t_start_sec;
    int nBins = std::ceil(duracion_total_s / 600.0);

    TH1F *hRateQuad = new TH1F(Form("hRateQuad_%d_%d_%d_%d", chs[0], chs[1], chs[2], chs[3]), 
                               Form("Rate Coincidencia 4 Chs;Hora del dia;Rate (Hz)"), 
                               nBins, t_start_sec, t_end_sec);

    for (Long64_t i = 0; i < nEntries; ++i) {
        bool all_over = true;
        for(int j = 0; j < 4; ++j) {
            trees[j]->GetEntry(i);
            if(maxs[j] <= threshold) {
                all_over = false;
                break;
            }
        }
        if(all_over) hRateQuad->Fill(time_ref / 1e6);
    }

    hRateQuad->GetXaxis()->SetTimeDisplay(1);
    hRateQuad->GetXaxis()->SetTimeFormat("%H:%M");
    hRateQuad->GetXaxis()->SetTimeOffset(0, "gmt");
    hRateQuad->Scale(1.0 / 600.0);
    hRateQuad->Write("", TObject::kOverwrite);
    f->Close();
    std::cout << "[OK] Rate de coincidencia cuádruple generado." << std::endl;
}

// --- Rate OR de Barras (S1 || S2 || S3 || S4) ---
void CalcularRateORBarras(const char* targetFile, float threshold) {
    TFile *f = TFile::Open(targetFile, "UPDATE");
    if (!f || f->IsZombie()) return;

    // Configuración física de las barras del plano superior
    std::vector<std::vector<int>> barras = {
        {0, 26, 1, 27}, // S1
        {3, 24, 2, 25}, // S2
        {4, 30, 5, 31}, // S3
        {6, 28, 7, 29}  // S4
    };

    const int nB = 4; // Número de barras
    TTree* trees[nB][4];
    float maxs[nB][4];
    double time_ref;

    // Vincular todos los árboles y ramas
    for(int i = 0; i < nB; ++i) {
        for(int j = 0; j < 4; ++j) {
            trees[i][j] = (TTree*)f->Get(Form("p%d", barras[i][j]));
            if(!trees[i][j]) {
                std::cerr << "[!] Error: No se halló el canal p" << barras[i][j] << std::endl;
                f->Close(); return;
            }
            trees[i][j]->SetBranchAddress("f_max", &maxs[i][j]);
            // Usamos el tiempo del primer canal de la primera barra como referencia global
            if(i == 0 && j == 0) trees[i][j]->SetBranchAddress("f_acqTime", &time_ref);
        }
    }

    Long64_t nEntries = trees[0][0]->GetEntries();
    trees[0][0]->GetEntry(0);
    double t_start = time_ref / 1e6;
    trees[0][0]->GetEntry(nEntries - 1);
    double t_end = time_ref / 1e6;

    int nBins = std::ceil((t_end - t_start) / 600.0);
    TH1F *hRateOR = new TH1F("hRate_PlanoSuperior_OR", 
                             "Rate Plano Superior (S1||S2||S3||S4);Hora;Rate (Hz)", 
                             nBins, t_start, t_end);

    std::cout << "Analizando " << nEntries << " eventos para OR de Barras..." << std::endl;

    for (Long64_t i = 0; i < nEntries; ++i) {
        bool hitPlano = false;
        
        // Iteramos sobre cada barra
        for(int b = 0; b < nB; ++b) {
            bool coincidenciaBarra = true;
            // AND de los 4 SiPMs de la barra 'b'
            for(int c = 0; c < 4; ++c) {
                trees[b][c]->GetEntry(i);
                if(maxs[b][c] <= threshold) {
                    coincidenciaBarra = false;
                    break;
                }
            }
            // Si la barra b tuvo coincidencia, el plano se activó
            if(coincidenciaBarra) {
                hitPlano = true;
                break; // OR: con una barra basta
            }
        }
        
        if(hitPlano) hRateOR->Fill(time_ref / 1e6);
    }

    hRateOR->GetXaxis()->SetTimeDisplay(1);
    hRateOR->GetXaxis()->SetTimeFormat("%H:%M");
    hRateOR->GetXaxis()->SetTimeOffset(0, "gmt");
    hRateOR->Scale(1.0 / 600.0);
    hRateOR->Write("", TObject::kOverwrite);
    
    f->Close();
    std::cout << "[OK] Histograma de Rate OR guardado exitosamente." << std::endl;
}

// --- MAIN ---
int main(int argc, char** argv) {
    if (argc >= 3) {
        ProcesarPSA(argv[1], argv[2]);
        return 0;
    }

    int opcion;
    std::cout << "\n--- SOFTWARE DE ANÁLISIS DE PULSOS ---" << std::endl;
    std::cout << "1. Procesamiento PSA (Raw -> Root)" << std::endl;
    std::cout << "2. Análisis Individual (Rate - Umbral 150)" << std::endl;
    std::cout << "3. Correlación entre 4 canales" << std::endl;
    std::cout << "4. Rate por Coincidencia AND (2 canales)" << std::endl;
    std::cout << "5. Rate por Coincidencia AND (4 canales)" << std::endl;
    std::cout << "6. Rate OR de Barras (S1 || S2 || S3 || S4)" << std::endl;
    std::cout << "0. Salir" << std::endl;
    std::cout << "Seleccione: ";
    if (!(std::cin >> opcion)) return 0;

    switch(opcion) {
        case 1: {
            std::string in, out;
            std::cout << "Archivo entrada: "; std::cin >> in;
            std::cout << "Archivo salida: "; std::cin >> out;
            ProcesarPSA(in.c_str(), out.c_str());
            break;
        }
        case 2: {
            int ch; std::string file;
            std::cout << "Archivo .root procesado: "; std::cin >> file;
            std::cout << "Canal: "; std::cin >> ch;
            CalcularRate(file.c_str(), ch, 150.0);
            break;
        }
        case 3: {
            std::vector<int> chs(4);
            std::string file;
            std::cout << "Archivo .root procesado: "; std::cin >> file;
            for(int i=0; i<4; i++) {
                std::cout << "Ingrese Canal " << i+1 << ": ";
                std::cin >> chs[i];
            }
            GenerarCorrelacionCanales(file.c_str(), chs);
            break;
        }
        case 4: {
            int c1, c2; std::string file;
            std::cout << "Archivo .root procesado: "; std::cin >> file;
            std::cout << "Canal 1: "; std::cin >> c1;
            std::cout << "Canal 2: "; std::cin >> c2;
            CalcularRateCoincidencia(file.c_str(), c1, c2, 150.0);
            break;
        }
        case 5: {
            std::vector<int> chs(4);
            std::string file;
            std::cout << "Archivo .root procesado: "; std::cin >> file;
            for(int i=0; i<4; i++) {
                std::cout << "Ingrese Canal " << i+1 << ": ";
                std::cin >> chs[i];
            }
            CalcularRateCoincidenciaCuadruple(file.c_str(), chs, 150.0);
            break;
        }
        case 6: {
            std::string file;
            std::cout << "Archivo .root procesado: "; std::cin >> file;
            CalcularRateORBarras(file.c_str(), 150.0);
            break;
        }
    }
    return 0;
}
