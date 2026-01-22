#include <iostream>
#include <vector>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <cfloat>
#include <iomanip>
#include <ctime>
#include <sstream>
#include <cstdlib>

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

// --- Estructura de Datos Empaqueta ---
#pragma pack(push, 1)
struct Pulse {
    float fMax;    
    float ft0;     
    float fInt;    
    double fAcqTime; 
    void Reset() {
        fMax = 0; ft0 = 0; fInt = 0; fAcqTime = 0;
    }
};
#pragma pack(pop)

// --- Tiempo Universal a Hora local ---
std::string FormatearTiempoGlobal(double total_s_local) {  // Ahora total_s_local es segundos en hora local
    time_t segundos = static_cast<time_t>(total_s_local);
    long microsegundos = static_cast<long>((total_s_local - segundos) * 1e6);
    struct tm *timeinfo = localtime(&segundos);  // Ya usa zona de Chile forzada
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
    // UN SOLO TTREE PARA TODOS LOS CANALES
    TTree* T = new TTree("T", "Datos Procesados Hodoscopio");
    
    Pulse pData[NUM_CHANNELS];
    TH2F* hPersistencia[NUM_CHANNELS];

    for (int ch = 0; ch < NUM_CHANNELS; ++ch) {
        // CADA CANAL ES UNA BRANCH CON SUS LEAFS
        T->Branch(Form("p%d", ch), &pData[ch], "fMax/F:ft0/F:fInt/F:fAcqTime/D");
        
        hPersistencia[ch] = new TH2F(Form("hPersist_Ch%d", ch), 
                                     Form("Persistencia Ch %d;Time (ns);ADC", ch),
                                     SAMPLES, 0, SAMPLES * BIN_WIDTH, 550, -50, 2000);
    }

    const double OFFSET_CHILE = -10800.0; 


    Long64_t nEntries = events->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        events->GetEntry(i);
        double utc_seconds = static_cast<double>(event_time_us) / 1e6;  // Asumiendo microsegundos UTC
        double local_seconds = utc_seconds + OFFSET_CHILE;  // Convertir a segundos en hora local de Chile
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
            pData[ch].fMax = current_max;
            pData[ch].fInt = current_int * BIN_WIDTH;
            pData[ch].ft0  = peak_sample * BIN_WIDTH;
            pData[ch].fAcqTime = local_seconds;
        }
        T->Fill(); // Se llena una vez por evento (contiene los 32 canales)
    }

    T->Write(); 
    for (int ch = 0; ch < NUM_CHANNELS; ++ch) {
        hPersistencia[ch]->SetMarkerStyle(1);
        hPersistencia[ch]->SetOption("SCAT");
        hPersistencia[ch]->Write(); 
    }
    fOut->Close(); 
    fIn->Close();
    std::cout << "Procesamiento listo, Guardado en: " << outputFile << std::endl;
}

// --- Funciones auxiliares de tiempo ---
void ConfigurarEjeTiempo(TAxis *axis) {
    axis->SetTimeDisplay(1);
    axis->SetTimeFormat("%H:%M");  // Solo hora:minuto
    axis->SetTimeOffset(0);  // Sin offset, ya que el tiempo está en local
}

// --- Análisis Individual (Rate) ---
void CalcularRate(const char* targetFile, int canal_objetivo, float threshold) {
    TFile *f = TFile::Open(targetFile, "UPDATE");
    TTree *T = (TTree*)f->Get("T");
    if(!T) { f->Close(); return; }

    Pulse p;
    T->SetBranchAddress(Form("p%d", canal_objetivo), &p);

    Long64_t nEntries = T->GetEntries();
    T->GetEntry(0); double t_start = p.fAcqTime;
    T->GetEntry(nEntries - 1); double t_end = p.fAcqTime;
    if(t_end <= t_start) t_end = t_start + 1.0;

    int nBins = std::max(1, (int)std::ceil((t_end - t_start) / 600.0));
    TH1F *hRate = new TH1F(Form("hRateGlobal_Ch%d", canal_objetivo), 
                           Form("Rate Canal %d;Hora;Hz", canal_objetivo), nBins, t_start, t_end);

    for (Long64_t i = 0; i < nEntries; ++i) {
        T->GetEntry(i);
        if (p.fMax > threshold) hRate->Fill(p.fAcqTime);
    }
    ConfigurarEjeTiempo(hRate->GetXaxis());
    hRate->Scale(1.0 / 600.0);
    hRate->Write("", TObject::kOverwrite);
    f->Close();
}

// --- Correlación ---
void GenerarCorrelacionCanales(const char* targetFile, std::vector<int> chs) {
    TFile *f = TFile::Open(targetFile, "UPDATE");
    TTree *T = (TTree*)f->Get("T");
    if(!T) return;

    Pulse p1, p2;
    TCanvas *cCorr = new TCanvas("cCorrelaciones", "Correlaciones", 1200, 800);
    cCorr->Divide(3, 2);
    int pad = 1;

    for (size_t i = 0; i < chs.size(); ++i) {
        for (size_t j = i + 1; j < chs.size(); ++j) {
            if(pad > 6) break;
            cCorr->cd(pad++);
            
            // Acceso a las branches específicas del mismo TTree
            T->SetBranchAddress(Form("p%d", chs[i]), &p1);
            T->SetBranchAddress(Form("p%d", chs[j]), &p2);
            
            TH2F *hCorr = new TH2F(Form("hCorr_Ch%d_Ch%d", chs[i], chs[j]), 
                                   Form("Ch%d vs Ch%d;Max;Max", chs[i], chs[j]), 150, 150, 1000, 150, 150, 1000);
            
            for (Long64_t n = 0; n < T->GetEntries(); ++n) {
                T->GetEntry(n);
                if (p1.fMax > 150 || p2.fMax > 150) hCorr->Fill(p1.fMax, p2.fMax);
            }
            hCorr->Draw("COLZ");
            hCorr->Write("", TObject::kOverwrite);
        }
    }
    cCorr->Write(); f->Close();
}

// --- Rate Coincidencia AND (2 canales) ---
void CalcularRateCoincidencia(const char* targetFile, int c1, int c2, float threshold) {
    TFile *f = TFile::Open(targetFile, "UPDATE");
    TTree *T = (TTree*)f->Get("T");
    Pulse p1, p2;
    T->SetBranchAddress(Form("p%d", c1), &p1);
    T->SetBranchAddress(Form("p%d", c2), &p2);

    T->GetEntry(0); double ts = p1.fAcqTime;
    T->GetEntry(T->GetEntries()-1); double te = p1.fAcqTime;

    TH1F *h = new TH1F(Form("hRateCoin_Ch%d_Ch%d", c1, c2), "AND;Hora;Hz", std::max(1, (int)std::ceil((te-ts)/600.0)), ts, te);
    for(Long64_t i=0; i<T->GetEntries(); i++) {
        T->GetEntry(i);
        if(p1.fMax > threshold && p2.fMax > threshold) h->Fill(p1.fAcqTime);
    }
    ConfigurarEjeTiempo(h->GetXaxis());
    h->Scale(1.0/600.0); h->Write("", TObject::kOverwrite); f->Close();
}

// --- Rate Coincidencia AND (4 canales) ---
void CalcularRateCoincidenciaCuadruple(const char* targetFile, std::vector<int> chs, float threshold) {
    TFile *f = TFile::Open(targetFile, "UPDATE");
    TTree *T = (TTree*)f->Get("T");
    Pulse p[4];
    for(int i=0; i<4; i++) T->SetBranchAddress(Form("p%d", chs[i]), &p[i]);

    T->GetEntry(0); double ts = p[0].fAcqTime;
    T->GetEntry(T->GetEntries()-1); double te = p[0].fAcqTime;

    TH1F *h = new TH1F(Form("hRateQuad_%d_%d_%d_%d", chs[0], chs[1], chs[2], chs[3]), "AND 4;Hora;Hz", std::max(1,(int)std::ceil((te-ts)/600.0)), ts, te);
    for(Long64_t i=0; i<T->GetEntries(); i++){
        T->GetEntry(i);
        if(p[0].fMax>threshold && p[1].fMax>threshold && p[2].fMax>threshold && p[3].fMax>threshold) h->Fill(p[0].fAcqTime);
    }
    ConfigurarEjeTiempo(h->GetXaxis());
    h->Scale(1.0/600.0); h->Write("", TObject::kOverwrite); f->Close();
}

// --- Rate OR Planos ---
void CalcularRateORBarrasSuperior(const char* targetFile, float threshold) {
    TFile *f = TFile::Open(targetFile, "UPDATE");
    TTree *T = (TTree*)f->Get("T");
    std::vector<std::vector<int>> barras = {{0,26,1,27}, {3,24,2,25}, {4,30,5,31}, {6,28,7,29}};
    Pulse p[4][4];
    for(int i=0; i<4; i++) for(int j=0; j<4; j++) T->SetBranchAddress(Form("p%d", barras[i][j]), &p[i][j]);

    T->GetEntry(0); double ts = p[0][0].fAcqTime;
    T->GetEntry(T->GetEntries()-1); double te = p[0][0].fAcqTime;
    TH1F *h = new TH1F("hRate_PlanoSuperior_OR", "OR Superior;Hora;Hz", std::max(1,(int)std::ceil((te-ts)/600.0)), ts, te);
    for(Long64_t i=0; i<T->GetEntries(); i++) {
        T->GetEntry(i);
        bool hit = false;
        for(int b=0; b<4; b++) {
            if(p[b][0].fMax>threshold && p[b][1].fMax>threshold && p[b][2].fMax>threshold && p[b][3].fMax>threshold) { hit=true; break; }
        }
        if(hit) h->Fill(p[0][0].fAcqTime);
    }
    ConfigurarEjeTiempo(h->GetXaxis());
    h->Scale(1.0/600.0); h->Write("", TObject::kOverwrite); f->Close();
}

void CalcularRateORBarrasInferior(const char* targetFile, float threshold) {
    TFile *f = TFile::Open(targetFile, "UPDATE");
    TTree *T = (TTree*)f->Get("T");
    std::vector<std::vector<int>> barras = {{14,22,15,23}, {10,20,11,21}, {12,18,13,19}, {8,17,9,16}};
    Pulse p[4][4];
    for(int i=0; i<4; i++) for(int j=0; j<4; j++) T->SetBranchAddress(Form("p%d", barras[i][j]), &p[i][j]);

    T->GetEntry(0); double ts = p[0][0].fAcqTime;
    T->GetEntry(T->GetEntries()-1); double te = p[0][0].fAcqTime;
    TH1F *h = new TH1F("hRate_PlanoInferior_OR", "OR Inferior;Hora;Hz", std::max(1,(int)std::ceil((te-ts)/600.0)), ts, te);
    for(Long64_t i=0; i<T->GetEntries(); i++) {
        T->GetEntry(i);
        bool hit = false;
        for(int b=0; b<4; b++) {
            if(p[b][0].fMax>threshold && p[b][1].fMax>threshold && p[b][2].fMax>threshold && p[b][3].fMax>threshold) { hit=true; break; }
        }
        if(hit) h->Fill(p[0][0].fAcqTime);
    }
    ConfigurarEjeTiempo(h->GetXaxis());
    h->Scale(1.0/600.0); h->Write("", TObject::kOverwrite); f->Close();
}

void CalcularRateCoincidenciaPlanos(const char* targetFile, float threshold) {
    TFile *f = TFile::Open(targetFile, "UPDATE");
    if (!f || f->IsZombie()) return;
    TTree *T = (TTree*)f->Get("T");
    if (!T) { f->Close(); return; }

    std::vector<std::vector<int>> bS = {{0,26,1,27}, {3,24,2,25}, {4,30,5,31}, {6,28,7,29}};
    std::vector<std::vector<int>> bI = {{14,22,15,23}, {10,20,11,21}, {12,18,13,19}, {8,17,9,16}};
    
    Pulse pS[4][4], pI[4][4];
    for(int i=0; i<4; i++) {
        for(int j=0; j<4; j++) {
            T->SetBranchAddress(Form("p%d", bS[i][j]), &pS[i][j]);
            T->SetBranchAddress(Form("p%d", bI[i][j]), &pI[i][j]);
        }
    }

    T->GetEntry(0); double ts = pS[0][0].fAcqTime;
    T->GetEntry(T->GetEntries()-1); double te = pS[0][0].fAcqTime;

    // Histograma 1: Rate (Hz)
    TH1F *hRate = new TH1F("hRate_Coincidencia_Planos_AND", "AND Planos;Hora;Hz", 
                           std::max(1,(int)std::ceil((te-ts)/600.0)), ts, te);

    // Histograma 2: Distribución de Amplitudes (fMax)
    TH1F *hAmp = new TH1F("hAmp_Distribucion_Coincidencia", 
                          "Distribucion Amplitud en Coincidencia;fMax (ADC);Entries", 
                          500, 0, 2000);

    for(Long64_t i=0; i<T->GetEntries(); i++) {
        T->GetEntry(i);
        bool hS = false, hI = false;
        
        // Verificamos hit en plano superior
        for(int b=0; b<4; b++) {
            if(pS[b][0].fMax>threshold && pS[b][1].fMax>threshold && 
               pS[b][2].fMax>threshold && pS[b][3].fMax>threshold) { hS=true; break; }
        }
        // Verificamos hit en plano inferior
        for(int b=0; b<4; b++) {
            if(pI[b][0].fMax>threshold && pI[b][1].fMax>threshold && 
               pI[b][2].fMax>threshold && pI[b][3].fMax>threshold) { hI=true; break; }
        }

        // Si hay coincidencia en ambos planos
        if(hS && hI) {
            hRate->Fill(pS[0][0].fAcqTime);
            
            // Llenamos la distribución de amplitud con los fMax de los canales disparados
            for(int b=0; b<4; b++) {
                for(int j=0; j<4; j++) {
                    if(pS[b][j].fMax > threshold) hAmp->Fill(pS[b][j].fMax);
                    if(pI[b][j].fMax > threshold) hAmp->Fill(pI[b][j].fMax);
                }
            }
        }
    }

    ConfigurarEjeTiempo(hRate->GetXaxis());
    hRate->Scale(1.0/600.0);
    
    hRate->Write("", TObject::kOverwrite);
    hAmp->Write("", TObject::kOverwrite);
    
    std::cout << "Se generaron 'hRate_Coincidencia_Planos_AND' y 'hAmp_Distribucion_Coincidencia'" << std::endl;
    f->Close();
}


// --- MAIN ---
int main(int argc, char** argv) {
    setenv("TZ", "America/Santiago", 1);
    tzset();
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
    std::cout << "7. Rate OR de Barras (I1 || I2 || I3 || I4)" << std::endl;
    std::cout << "8. Rate Coincidencia AND (Superior && Inferior)" << std::endl;
    std::cout << "0. Salir" << std::endl;
    std::cout << "Seleccione: ";
    
    if (!(std::cin >> opcion) || opcion == 0) return 0;

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
            CalcularRateORBarrasSuperior(file.c_str(), 150.0);
            break;
        }
        case 7: {
            std::string file;
            std::cout << "Archivo .root procesado: "; std::cin >> file;
            CalcularRateORBarrasInferior(file.c_str(), 150.0);
            break;
        }
        case 8: {
            std::string file;
            std::cout << "Archivo .root procesado: "; std::cin >> file;
            CalcularRateCoincidenciaPlanos(file.c_str(), 150.0);
            break;
        }
        default:
            std::cout << "Opción no válida." << std::endl;
            break;
    }
    return 0;
}