#include <iostream>
#include <vector>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <cfloat>
#include <iomanip>
#include <ctime>
#include <sstream>
#include <TLegend.h>
#include <cstdlib>
#include <TGraph2D.h> 
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TProfile.h>
#include <TAxis.h>
#include <TGraph.h>
#include <TLine.h>
#include <TLatex.h>
#include <TDirectory.h>

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
std::string FormatearTiempoGlobal(double total_s_local) {
    time_t segundos = static_cast<time_t>(total_s_local);
    long microsegundos = static_cast<long>((total_s_local - segundos) * 1e6);
    struct tm *timeinfo = localtime(&segundos);
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
    
    TDirectory *dirPersist = fOut->mkdir("Persistencia");
    TDirectory *dirGlobal = fOut->mkdir("MaxAmp_Global");

    TTree* T = new TTree("T", "Datos Procesados Hodoscopio");
    
    Pulse pData[NUM_CHANNELS];
    TH2F* hPersistencia[NUM_CHANNELS];

    for (int ch = 0; ch < NUM_CHANNELS; ++ch) {
        T->Branch(Form("p%d", ch), &pData[ch], "fMax/F:ft0/F:fInt/F:fAcqTime/D");
        hPersistencia[ch] = new TH2F(Form("hPersist_Ch%d", ch), 
                                     Form("Persistencia Ch %d;Time (ns);ADC", ch),
                                     SAMPLES, 0, SAMPLES * BIN_WIDTH, 550, -50, 2000);
    }

    const double OFFSET_CHILE = -10800.0; 

    Long64_t nEntries = events->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        events->GetEntry(i);
        double utc_seconds = static_cast<double>(event_time_us) / 1e6;
        double local_seconds = utc_seconds + OFFSET_CHILE;
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
        T->Fill();
    }

    fOut->cd();
    T->Write(); 

    dirPersist->cd();
    for (int ch = 0; ch < NUM_CHANNELS; ++ch) {
        hPersistencia[ch]->SetMarkerStyle(1);
        hPersistencia[ch]->SetOption("SCAT");
        hPersistencia[ch]->Write(); 
    }

    dirGlobal->cd();
    for (int ch = 0; ch < NUM_CHANNELS; ++ch) {
        TH1F *hMaxGlobal = new TH1F(Form("hMax_Ch%d_Global", ch), 
                                    Form("Max Amplitude Ch %d Global;ADC;Entries", ch), 500, 0, 2000);
        T->Project(hMaxGlobal->GetName(), Form("p%d.fMax", ch));
        hMaxGlobal->Write();
        delete hMaxGlobal;
    }

    fOut->Close(); 
    fIn->Close();
    std::cout << "Procesamiento listo, Guardado en: " << outputFile << std::endl;
}

void ConfigurarEjeTiempo(TAxis *axis) {
    axis->SetTimeDisplay(1);
    axis->SetTimeFormat("%H:%M");
    axis->SetTimeOffset(0);
}

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

void CalcularRateCoincidenciaTriple(const char* targetFile, std::vector<int> chs, float threshold) {
    TFile *f = TFile::Open(targetFile, "UPDATE");
    if (!f || f->IsZombie()) return;
    TTree *T = (TTree*)f->Get("T");
    if (!T) { f->Close(); return; }

    Pulse p[3];
    for(int i=0; i<3; i++) T->SetBranchAddress(Form("p%d", chs[i]), &p[i]);

    T->GetEntry(0); double ts = p[0].fAcqTime;
    T->GetEntry(T->GetEntries()-1); double te = p[0].fAcqTime;

    TH1F *h = new TH1F(Form("hRateTriple_%d_%d_%d", chs[0], chs[1], chs[2]), 
                       "AND 3;Hora;Hz", 
                       std::max(1,(int)std::ceil((te-ts)/600.0)), ts, te);

    for(Long64_t i=0; i<T->GetEntries(); i++){
        T->GetEntry(i);
        if(p[0].fMax > threshold && p[1].fMax > threshold && p[2].fMax > threshold) {
            h->Fill(p[0].fAcqTime);
        }
    }

    ConfigurarEjeTiempo(h->GetXaxis());
    h->Scale(1.0/600.0); 
    h->Write("", TObject::kOverwrite); 
    f->Close();
}

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

    TH1F *hRate = new TH1F("hRate_Coincidencia_Planos_AND", "AND Planos;Hora;Hz", std::max(1,(int)std::ceil((te-ts)/600.0)), ts, te);
    TH1F *hAmp = new TH1F("hAmp_Distribucion_Coincidencia", "Distribucion Amplitud en Coincidencia;fMax (ADC);Entries", 500, 0, 2000);

    for(Long64_t i=0; i<T->GetEntries(); i++) {
        T->GetEntry(i);
        bool hS = false, hI = false;
        for(int b=0; b<4; b++) {
            if(pS[b][0].fMax>threshold && pS[b][1].fMax>threshold && pS[b][2].fMax>threshold && pS[b][3].fMax>threshold) { hS=true; break; }
        }
        for(int b=0; b<4; b++) {
            if(pI[b][0].fMax>threshold && pI[b][1].fMax>threshold && pI[b][2].fMax>threshold && pI[b][3].fMax>threshold) { hI=true; break; }
        }

        if(hS && hI) {
            hRate->Fill(pS[0][0].fAcqTime);
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
    f->Close();
}
struct ParUmbral {
    int chA;
    int chB;
    float thrA; // Umbral que debe superar chB para registrar chA
    float thrB; // Umbral que debe superar chA para registrar chB
};
void GenerarHistogramasMaxAmpCondicionados(const char* targetFile) {
    TFile *f = TFile::Open(targetFile, "UPDATE");
    if (!f || f->IsZombie()) return;

    TDirectory *dirCond = f->GetDirectory("MaxAmp_Condicionado");
    if (!dirCond) dirCond = f->mkdir("MaxAmp_Condicionado");

    TTree *T = (TTree*)f->Get("T");
    if (!T) { f->Close(); return; }

    std::vector<ParUmbral> pares = {
        {0, 26, 20.0, 20.0}, {1, 27, 22.0, 18.0}, {3, 24, 20.0, 17.0}, {2, 25, 20.0, 19.0},
        {6, 28, 21.0, 18.0}, {7, 29, 19.0, 19.0}, {4, 30, 19.0, 20.0}, {5, 31, 20.0, 20.0},
        {8, 17, 23.0, 20.0}, {9, 16, 20.0, 20.0}, {12, 18, 18.0, 21.0}, {13, 19, 21.0, 22.0},
        {10, 20, 20.0, 16.0}, {11, 21, 19.0, 21.0}, {14, 22, 18.0, 20.0}, {15, 23, 17.0, 22.0}
    };

    Pulse pA, pB;
    for (auto const& par : pares) {
        T->ResetBranchAddresses();
        dirCond->cd(); 
        // Importante: Usamos nombres únicos para evitar conflictos al recuperar
        TH1F *hMaxA_cond = new TH1F(Form("hMax_Ch%d_c%d", par.chA, par.chB), Form("Ch%d cond Ch%d", par.chA, par.chB), 200, 0, 200);
        TH1F *hMaxB_cond = new TH1F(Form("hMax_Ch%d_c%d", par.chB, par.chA), Form("Ch%d cond Ch%d", par.chB, par.chA), 200, 0, 200);

        T->SetBranchAddress(Form("p%d", par.chA), &pA);
        T->SetBranchAddress(Form("p%d", par.chB), &pB);

        for (Long64_t i = 0; i < T->GetEntries(); ++i) {
            T->GetEntry(i);
            if (pB.fMax > par.thrB) hMaxA_cond->Fill(pA.fMax);
            if (pA.fMax > par.thrA) hMaxB_cond->Fill(pB.fMax);
        }
        hMaxA_cond->Write("", TObject::kOverwrite);
        hMaxB_cond->Write("", TObject::kOverwrite);
        delete hMaxA_cond; delete hMaxB_cond;
    }

auto DibujarCascada = [&](const std::string& nombrePlano, int inicioPar, int finPar) {
        TCanvas *c = new TCanvas(Form("c_%s", nombrePlano.c_str()), nombrePlano.c_str(), 1100, 900);
        c->SetLeftMargin(0.12);
        gStyle->SetOptStat(0);
        
        TLegend *leg = new TLegend(0.13, 0.45, 0.30, 0.89); 
        leg->SetTextSize(0.018);
        leg->SetHeader(Form("Canales %s", nombrePlano.c_str()));

        // --- NUEVA LÓGICA DE ORDENAMIENTO ---
        // Recolectamos todos los canales involucrados en este rango de pares
        struct InfoCanal { int id; int idCond; };
        std::vector<InfoCanal> listaOrdenada;
        for (int p = inicioPar; p <= finPar; ++p) {
            listaOrdenada.push_back({pares[p].chA, pares[p].chB});
            listaOrdenada.push_back({pares[p].chB, pares[p].chA});
        }
        // Ordenamos de menor a mayor ID de canal
        std::sort(listaOrdenada.begin(), listaOrdenada.end(), [](const InfoCanal& a, const InfoCanal& b) {
            return a.id < b.id;
        });

        double dx = 8.0;  
        double dy = 12.0; 
        int totalHists = listaOrdenada.size();
        
        TLine *lDepth0 = new TLine(0, 0, (totalHists-1)*dx, (totalHists-1)*dy);
        lDepth0->SetLineColor(kGray+2);
        lDepth0->SetLineStyle(2);

        int colores[] = {kBlue+1, kRed+1, kGreen+2, kOrange+2, kMagenta+1, kCyan+2, kAzure+7, kViolet-3,
                         kPink+9, kSpring-6, kTeal+3, kYellow+2, kGray+2, kOrange-3, kRed-7, kBlue-9};

        for (int j = 0; j < totalHists; ++j) {
            int canalActual = listaOrdenada[j].id;
            int canalCondicion = listaOrdenada[j].idCond;

            TH1F *h = (TH1F*)dirCond->Get(Form("hMax_Ch%d_c%d", canalActual, canalCondicion));
            if (!h) continue;

            TGraph *g = new TGraph();
            double xShift = j * dx;
            double yShift = j * dy;

            for (int b = 1; b <= h->GetNbinsX(); b++) {
                g->SetPoint(b-1, h->GetBinCenter(b) + xShift, h->GetBinContent(b) + yShift);
            }

            g->SetLineColor(colores[j % 16]);
            g->SetLineWidth(2);
            leg->AddEntry(g, Form("Ch%d", canalActual), "l");

            if (j == 0) {
                g->Draw("AL"); 
                g->GetXaxis()->SetTitle("ADC");
                g->GetYaxis()->SetTitle("Entries");
                g->GetXaxis()->SetRangeUser(-10, 200 + (totalHists * dx) + 20);
                g->GetYaxis()->SetRangeUser(-10, h->GetMaximum() + (totalHists * dy) + 50);
                lDepth0->Draw("SAME");
            } else {
                g->Draw("L SAME");
            }

            TLine *base = new TLine(xShift, yShift, 200 + xShift, yShift);
            base->SetLineColor(kGray+1);
            base->SetLineStyle(3);
            base->Draw();
        }
        
        TLatex *title = new TLatex();
        title->SetNDC();
        title->SetTextSize(0.035);
        title->DrawLatex(0.4, 0.93, Form("Amplitudes: Plano %s", nombrePlano.c_str()));

        leg->Draw();
        c->Write("", TObject::kOverwrite);
    };

    dirCond->cd();
    DibujarCascada("Superior", 0, 7);  // Pares 0-7 -> 16 canales
    DibujarCascada("Inferior", 8, 15); // Pares 8-15 -> 16 canales

    T->ResetBranchAddresses();
    f->Close();
}
void AnalizarEstabilidadTemperatura(const char* targetFile, int canal) {
    TFile *f = TFile::Open(targetFile, "UPDATE");
    if (!f || f->IsZombie()) return;
    TTree *T = (TTree*)f->Get("T");
    if (!T) { f->Close(); return; }

    Pulse p;
    T->SetBranchAddress(Form("p%d", canal), &p);

    T->GetEntry(0); double t_start = p.fAcqTime;
    T->GetEntry(T->GetEntries() - 1); double t_end = p.fAcqTime;

    int nBins = std::max(1, (int)std::ceil((t_end - t_start) / 1800.0));
    
    TProfile *hProf = new TProfile(Form("hTempProf_Ch%d", canal), 
                                   Form("Estabilidad Ch %d;Hora;Media fMax (ADC)", canal), 
                                   nBins, t_start, t_end);

    for (Long64_t i = 0; i < T->GetEntries(); ++i) {
        T->GetEntry(i);
        if (p.fMax > 20.0) hProf->Fill(p.fAcqTime, p.fMax);
    }

    // Configurar el eje del perfil original ANTES de clonar
    ConfigurarEjeTiempo(hProf->GetXaxis());

    TCanvas *c = new TCanvas(Form("cTemp_Ch%d", canal), "Analisis Temperatura", 1000, 600);
    
    // Clonamos el perfil ya configurado
    TProfile *hRed = (TProfile*)hProf->Clone(Form("%s_red", hProf->GetName()));
    TProfile *hBlue = (TProfile*)hProf->Clone(Form("%s_blue", hProf->GetName()));

    for (int b = 1; b <= nBins; b++) {
        double binCenter = hProf->GetBinCenter(b);
        time_t rawtime = (time_t)binCenter;
        struct tm *timeinfo = localtime(&rawtime);
        int hora = timeinfo->tm_hour;

        if (hora >= 9 && hora < 19) {
            hRed->SetBinContent(b, 0); 
            hRed->SetBinEntries(b, 0);
        } else {
            hBlue->SetBinContent(b, 0); 
            hBlue->SetBinEntries(b, 0);
        }
    }

    // Configuración estética
    hRed->SetMarkerStyle(20);
    hRed->SetMarkerSize(0.8);
    hRed->SetMarkerColor(kRed+1);
    hRed->SetLineColor(kRed+1);

    hBlue->SetMarkerStyle(20);
    hBlue->SetMarkerSize(0.8);
    hBlue->SetMarkerColor(kBlue+1);
    hBlue->SetLineColor(kBlue+1);
    
    // Re-asegurar formato tiempo en los clones para el Canvas
    ConfigurarEjeTiempo(hRed->GetXaxis());
    ConfigurarEjeTiempo(hBlue->GetXaxis());
    
    hRed->Draw("E1"); 
    hBlue->Draw("E1 SAME");

    // Guardar el perfil original (que ya tiene ConfigurarEjeTiempo arriba)
    hProf->Write("", TObject::kOverwrite);
    c->Write("", TObject::kOverwrite);
    
    std::cout << "Análisis de estabilidad generado para Canal " << canal << std::endl;
    f->Close();
}
void AnalizarDeltaTime(const char* targetFile, int canal, float threshold) {
    TFile *f = TFile::Open(targetFile, "UPDATE");
    if (!f || f->IsZombie()) return;
    TTree *T = (TTree*)f->Get("T");
    if (!T) { f->Close(); return; }

    Pulse p;
    T->SetBranchAddress(Form("p%d", canal), &p);

    Long64_t nEntries = T->GetEntries();
    if (nEntries < 2) {
        std::cout << "[!] No hay suficientes entradas para calcular Delta T." << std::endl;
        f->Close();
        return;
    }

    T->GetEntry(0); double t_start = p.fAcqTime;
    T->GetEntry(nEntries - 1); double t_end = p.fAcqTime;
    if (t_end <= t_start) t_end = t_start + 3600.0; 

    int nBins = std::max(1, (int)std::ceil((t_end - t_start) / 1800));

    TProfile *hProfDT = new TProfile(Form("hDeltaT_Evol_Ch%d", canal), 
                                     Form("Tiempo entre eventos Ch %d (Thr > %.1f);Hora;#Delta t promedio [ms]", canal, threshold), 
                                     nBins, t_start, t_end, "s");
    
    double lastTime = -1.0;
    for (Long64_t i = 0; i < nEntries; ++i) {
        T->GetEntry(i);
        if (p.fMax > threshold) {
            if (lastTime > 0) {
                double dt_ms = (p.fAcqTime - lastTime) * 1000.0;
                if (dt_ms > 0 && dt_ms < 3600000) hProfDT->Fill(p.fAcqTime, dt_ms);
            }
            lastTime = p.fAcqTime;
        }
    }

    TCanvas *c = new TCanvas(Form("cDeltaT_Evol_Ch%d", canal), "Evolucion Delta T", 1200, 700);
    hProfDT->SetMarkerStyle(20);
    hProfDT->SetMarkerSize(0.8);
    hProfDT->SetMarkerColor(kBlue+1);
    hProfDT->SetLineColor(kBlue+1);
    
    ConfigurarEjeTiempo(hProfDT->GetXaxis());
    hProfDT->Draw("E1");

    hProfDT->Write("", TObject::kOverwrite);
    c->Write("", TObject::kOverwrite);
    f->Close();
    std::cout << "[+] Análisis Delta T (ms) completado para el canal " << canal << std::endl;
}

void AnalizarDeltaTimeCoincidencia(const char* targetFile, int c1, int c2, float thr1, float thr2) {
    TFile *f = TFile::Open(targetFile, "UPDATE");
    TTree *T = (TTree*)f->Get("T");
    Pulse p1, p2;
    T->SetBranchAddress(Form("p%d", c1), &p1);
    T->SetBranchAddress(Form("p%d", c2), &p2);

    T->GetEntry(0); double ts = p1.fAcqTime;
    T->GetEntry(T->GetEntries()-1); double te = p1.fAcqTime;
    int nBins = std::max(1, (int)std::ceil((te-ts)/1800.0));

    TProfile *h = new TProfile(Form("hDeltaT_Coin_%d_%d", c1, c2), "Delta T Coincidencia;Hora;#Delta t promedio [ms]", nBins, ts, te, "s");

    double lastTime = -1.0;
    for(Long64_t i=0; i<T->GetEntries(); i++) {
        T->GetEntry(i);
        if(p1.fMax > thr1 && p2.fMax > thr2) {
            if(lastTime > 0) h->Fill(p1.fAcqTime, (p1.fAcqTime - lastTime)*1000.0);
            lastTime = p1.fAcqTime;
        }
    }
    ConfigurarEjeTiempo(h->GetXaxis());
    h->SetMarkerStyle(20); h->SetMarkerColor(kBlue+1); h->SetLineColor(kBlue+1);
    h->Write("", TObject::kOverwrite); f->Close();
}

void AnalizarDeltaTimeCoincidenciaTriple(const char* targetFile, std::vector<int> chs, std::vector<float> thrs) {
    TFile *f = TFile::Open(targetFile, "UPDATE");
    TTree *T = (TTree*)f->Get("T");
    Pulse p[3];
    for(int i=0; i<3; i++) T->SetBranchAddress(Form("p%d", chs[i]), &p[i]);

    T->GetEntry(0); double ts = p[0].fAcqTime;
    T->GetEntry(T->GetEntries()-1); double te = p[0].fAcqTime;
    int nBins = std::max(1, (int)std::ceil((te-ts)/1800.0));

    TProfile *h = new TProfile(Form("hDeltaT_Triple_%d_%d_%d", chs[0], chs[1], chs[2]), "Delta T Triple;Hora;#Delta t promedio [ms]", nBins, ts, te, "s");

    double lastTime = -1.0;
    for(Long64_t i=0; i<T->GetEntries(); i++) {
        T->GetEntry(i);
        if(p[0].fMax > thrs[0] && p[1].fMax > thrs[1] && p[2].fMax > thrs[2]) {
            if(lastTime > 0) h->Fill(p[0].fAcqTime, (p[0].fAcqTime - lastTime)*1000.0);
            lastTime = p[0].fAcqTime;
        }
    }
    ConfigurarEjeTiempo(h->GetXaxis());
    h->SetMarkerStyle(20); h->SetMarkerColor(kBlue+1); h->SetLineColor(kBlue+1);
    h->Write("", TObject::kOverwrite); f->Close();
}

void AnalizarDeltaTimeCoincidenciaCuadruple(const char* targetFile, std::vector<int> chs, std::vector<float> thrs) {
    TFile *f = TFile::Open(targetFile, "UPDATE");
    TTree *T = (TTree*)f->Get("T");
    Pulse p[4];
    for(int i=0; i<4; i++) T->SetBranchAddress(Form("p%d", chs[i]), &p[i]);

    T->GetEntry(0); double ts = p[0].fAcqTime;
    T->GetEntry(T->GetEntries()-1); double te = p[0].fAcqTime;
    int nBins = std::max(1, (int)std::ceil((te-ts)/1800.0));

    TProfile *h = new TProfile(Form("hDeltaT_Quad_%d_%d_%d_%d", chs[0], chs[1], chs[2], chs[3]), "Delta T Cuadruple;Hora;#Delta t promedio [ms]", nBins, ts, te, "s");

    double lastTime = -1.0;
    for(Long64_t i=0; i<T->GetEntries(); i++) {
        T->GetEntry(i);
        bool allPass = true;
        for(int j=0; j<4; j++) if(p[j].fMax <= thrs[j]) { allPass = false; break; }
        if(allPass) {
            if(lastTime > 0) h->Fill(p[0].fAcqTime, (p[0].fAcqTime - lastTime)*1000.0);
            lastTime = p[0].fAcqTime;
        }
    }
    ConfigurarEjeTiempo(h->GetXaxis());
    h->SetMarkerStyle(20); h->SetMarkerColor(kBlue+1); h->SetLineColor(kBlue+1);
    h->Write("", TObject::kOverwrite); f->Close();
}

void AnalizarVoltajeTemporal(const char* targetFile, int canal, float threshold) {
    TFile *f = TFile::Open(targetFile, "UPDATE");
    TTree *T = (TTree*)f->Get("T");
    Pulse p;
    T->SetBranchAddress(Form("p%d", canal), &p);
    T->GetEntry(0); double ts = p.fAcqTime;
    T->GetEntry(T->GetEntries() - 1); double te = p.fAcqTime;

    TProfile *hV = new TProfile(Form("hVolt_mV_Ch%d", canal), "Amplitud Promedio;Hora;Voltaje [mV]", std::max(1, (int)std::ceil((te - ts) / 1200.0)), ts, te);

    for (Long64_t i = 0; i < T->GetEntries(); ++i) {
        T->GetEntry(i);
        if (p.fMax > threshold) hV->Fill(p.fAcqTime, p.fMax * 0.091);
    }

    ConfigurarEjeTiempo(hV->GetXaxis());
    hV->SetMarkerStyle(20); hV->SetMarkerColor(kBlue+1); hV->SetLineColor(kBlue+1);
    hV->Write("", TObject::kOverwrite); f->Close();
}
// --- Nuevas funciones de Amplitud en Coincidencia ---

void AnalizarVoltajeCoincidenciaDoble(const char* targetFile, int chA, float thrA, int chB, float thrB) {
    TFile *f = TFile::Open(targetFile, "UPDATE");
    TTree *T = (TTree*)f->Get("T");
    Pulse pA, pB;
    T->SetBranchAddress(Form("p%d", chA), &pA);
    T->SetBranchAddress(Form("p%d", chB), &pB);
    T->GetEntry(0); double ts = pA.fAcqTime;
    T->GetEntry(T->GetEntries()-1); double te = pA.fAcqTime;

    TProfile *hV = new TProfile(Form("hVolt_Ch%d_en_Coin_Ch%d", chA, chB), 
                                Form("Amplitud Ch%d (Coin con Ch%d);Hora;Voltaje [mV]", chA, chB), 
                                std::max(1, (int)std::ceil((te-ts)/1200.0)), ts, te);

    for(Long64_t i=0; i<T->GetEntries(); i++) {
        T->GetEntry(i);
        // Cada canal compara contra su propio threshold
        if(pA.fMax > thrA && pB.fMax > thrB) hV->Fill(pA.fAcqTime, pA.fMax * 0.091);
    }
    ConfigurarEjeTiempo(hV->GetXaxis());
    hV->SetMarkerStyle(20); hV->SetMarkerColor(kBlue+1); hV->SetLineColor(kBlue+1);
    hV->Write("", TObject::kOverwrite); f->Close();
}

void AnalizarVoltajeCoincidenciaTriple(const char* targetFile, int targetCh, float targetThr, std::vector<int> others, std::vector<float> otherThrs) {
    TFile *f = TFile::Open(targetFile, "UPDATE");
    TTree *T = (TTree*)f->Get("T");
    Pulse pTarget, pOthers[2];
    T->SetBranchAddress(Form("p%d", targetCh), &pTarget);
    for(int i=0; i<2; i++) T->SetBranchAddress(Form("p%d", others[i]), &pOthers[i]);
    T->GetEntry(0); double ts = pTarget.fAcqTime;
    T->GetEntry(T->GetEntries()-1); double te = pTarget.fAcqTime;

    TProfile *hV = new TProfile(Form("hVolt_Ch%d_en_Coin_Triple", targetCh), 
                                Form("Amplitud Ch%d (Coin Triple);Hora;Voltaje [mV]", targetCh), 
                                std::max(1, (int)std::ceil((te-ts)/1200.0)), ts, te);

    for(Long64_t i=0; i<T->GetEntries(); i++) {
        T->GetEntry(i);
        if(pTarget.fMax > targetThr && pOthers[0].fMax > otherThrs[0] && pOthers[1].fMax > otherThrs[1])
            hV->Fill(pTarget.fAcqTime, pTarget.fMax * 0.091);
    }
    ConfigurarEjeTiempo(hV->GetXaxis());
    hV->SetMarkerStyle(20); hV->SetMarkerColor(kBlue+1); hV->SetLineColor(kBlue+1);
    hV->Write("", TObject::kOverwrite); f->Close();
}

void AnalizarVoltajeCoincidenciaCuadruple(const char* targetFile, int targetCh, float targetThr, std::vector<int> others, std::vector<float> otherThrs) {
    TFile *f = TFile::Open(targetFile, "UPDATE");
    TTree *T = (TTree*)f->Get("T");
    Pulse pTarget, pOthers[3];
    T->SetBranchAddress(Form("p%d", targetCh), &pTarget);
    for(int i=0; i<3; i++) T->SetBranchAddress(Form("p%d", others[i]), &pOthers[i]);
    T->GetEntry(0); double ts = pTarget.fAcqTime;
    T->GetEntry(T->GetEntries()-1); double te = pTarget.fAcqTime;

    TProfile *hV = new TProfile(Form("hVolt_Ch%d_en_Coin_Quad", targetCh), 
                                Form("Amplitud Ch%d (Coin Cuadruple);Hora;Voltaje [mV]", targetCh), 
                                std::max(1, (int)std::ceil((te-ts)/1200.0)), ts, te);

    for(Long64_t i=0; i<T->GetEntries(); i++) {
        T->GetEntry(i);
        bool coin = (pTarget.fMax > targetThr);
        for(int j=0; j<3; j++) if(pOthers[j].fMax <= otherThrs[j]) coin = false;
        if(coin) hV->Fill(pTarget.fAcqTime, pTarget.fMax * 0.091);
    }
    ConfigurarEjeTiempo(hV->GetXaxis());
    hV->SetMarkerStyle(20); hV->SetMarkerColor(kBlue+1); hV->SetLineColor(kBlue+1);
    hV->Write("", TObject::kOverwrite); f->Close();
}

void AnalizarCargaTemporal(const char* targetFile, int canal, float threshold) {
    TFile *f = TFile::Open(targetFile, "UPDATE");
    TTree *T = (TTree*)f->Get("T");
    Pulse p;
    T->SetBranchAddress(Form("p%d", canal), &p);
    T->GetEntry(0); double ts = p.fAcqTime;
    T->GetEntry(T->GetEntries() - 1); double te = p.fAcqTime;

    TProfile *hQ = new TProfile(Form("hCharge_Evol_Ch%d", canal), "Evolucion Carga;Hora;Carga [mV ns]", std::max(1, (int)std::ceil((te - ts) / 1200.0)), ts, te);

    for (Long64_t i = 0; i < T->GetEntries(); ++i) {
        T->GetEntry(i);
        if (p.fMax > threshold) hQ->Fill(p.fAcqTime, p.fInt * 0.091); 
    }

    ConfigurarEjeTiempo(hQ->GetXaxis());
    hQ->SetMarkerStyle(20); hQ->SetMarkerColor(kBlue+1); hQ->SetLineColor(kBlue+1);
    hQ->Write("", TObject::kOverwrite); f->Close();
}
// --- FUNCIONES DE CARGA EN COINCIDENCIA ---

void AnalizarCargaCoincidenciaDoble(const char* targetFile, int chA, float thrA, int chB, float thrB) {
    TFile *f = TFile::Open(targetFile, "UPDATE");
    if (!f || f->IsZombie()) { std::cout << "Error al abrir archivo." << std::endl; return; }
    TTree *T = (TTree*)f->Get("T");
    if (!T) { f->Close(); return; }

    Pulse pA, pB;
    T->SetBranchAddress(Form("p%d", chA), &pA);
    T->SetBranchAddress(Form("p%d", chB), &pB);
    T->GetEntry(0); double ts = pA.fAcqTime;
    T->GetEntry(T->GetEntries()-1); double te = pA.fAcqTime;

    TProfile *hQ = new TProfile(Form("hCharge_Ch%d_en_Coin_Ch%d", chA, chB), 
                                Form("Carga Ch%d (Coin con Ch%d);Hora;Carga [mV ns]", chA, chB), 
                                std::max(1, (int)std::ceil((te-ts)/1200.0)), ts, te);

    for(Long64_t i=0; i<T->GetEntries(); i++) {
        T->GetEntry(i);
        if(pA.fMax > thrA && pB.fMax > thrB) hQ->Fill(pA.fAcqTime, pA.fInt * 0.091);
    }
    ConfigurarEjeTiempo(hQ->GetXaxis());
    hQ->SetMarkerStyle(20); hQ->SetMarkerColor(kBlue+1); hQ->SetLineColor(kBlue+1);
    hQ->Write("", TObject::kOverwrite); f->Close();
}

void AnalizarCargaCoincidenciaTriple(const char* targetFile, int targetCh, float targetThr, std::vector<int> others, std::vector<float> otherThrs) {
    TFile *f = TFile::Open(targetFile, "UPDATE");
    if (!f || f->IsZombie()) { std::cout << "Error al abrir archivo." << std::endl; return; }
    TTree *T = (TTree*)f->Get("T");
    if (!T) { f->Close(); return; }

    Pulse pTarget, pOthers[2];
    T->SetBranchAddress(Form("p%d", targetCh), &pTarget);
    for(int i=0; i<2; i++) T->SetBranchAddress(Form("p%d", others[i]), &pOthers[i]);
    T->GetEntry(0); double ts = pTarget.fAcqTime;
    T->GetEntry(T->GetEntries()-1); double te = pTarget.fAcqTime;

    TProfile *hQ = new TProfile(Form("hCharge_Ch%d_en_Coin_Triple", targetCh), 
                                Form("Carga Ch%d (Coin Triple);Hora;Carga [mV ns]", targetCh), 
                                std::max(1, (int)std::ceil((te-ts)/1200.0)), ts, te);

    for(Long64_t i=0; i<T->GetEntries(); i++) {
        T->GetEntry(i);
        if(pTarget.fMax > targetThr && pOthers[0].fMax > otherThrs[0] && pOthers[1].fMax > otherThrs[1])
            hQ->Fill(pTarget.fAcqTime, pTarget.fInt * 0.091);
    }
    ConfigurarEjeTiempo(hQ->GetXaxis());
    hQ->SetMarkerStyle(20); hQ->SetMarkerColor(kBlue+1); hQ->SetLineColor(kBlue+1);
    hQ->Write("", TObject::kOverwrite); f->Close();
}

void AnalizarCargaCoincidenciaCuadruple(const char* targetFile, int targetCh, float targetThr, std::vector<int> others, std::vector<float> otherThrs) {
    TFile *f = TFile::Open(targetFile, "UPDATE");
    if (!f || f->IsZombie()) { std::cout << "Error al abrir archivo." << std::endl; return; }
    TTree *T = (TTree*)f->Get("T");
    if (!T) { f->Close(); return; }

    Pulse pTarget, pOthers[3];
    T->SetBranchAddress(Form("p%d", targetCh), &pTarget);
    for(int i=0; i<3; i++) T->SetBranchAddress(Form("p%d", others[i]), &pOthers[i]);
    T->GetEntry(0); double ts = pTarget.fAcqTime;
    T->GetEntry(T->GetEntries()-1); double te = pTarget.fAcqTime;

    TProfile *hQ = new TProfile(Form("hCharge_Ch%d_en_Coin_Quad", targetCh), 
                                Form("Carga Ch%d (Coin Cuadruple);Hora;Carga [mV ns]", targetCh), 
                                std::max(1, (int)std::ceil((te-ts)/1200.0)), ts, te);

    for(Long64_t i=0; i<T->GetEntries(); i++) {
        T->GetEntry(i);
        bool coin = (pTarget.fMax > targetThr);
        for(int j=0; j<3; j++) if(pOthers[j].fMax <= otherThrs[j]) coin = false;
        if(coin) hQ->Fill(pTarget.fAcqTime, pTarget.fInt * 0.091);
    }
    ConfigurarEjeTiempo(hQ->GetXaxis());
    hQ->SetMarkerStyle(20); hQ->SetMarkerColor(kBlue+1); hQ->SetLineColor(kBlue+1);
    hQ->Write("", TObject::kOverwrite); f->Close();
}
// --- Constantes de Flujo (Áreas en cm2) ---
const float AREA_BARRA = 500.0;      // 125cm * 4cm
const float AREA_PLANO_OR = 2000.0;  // 125cm * 16cm
const float AREA_PIXEL_AND = 16.0;   // 4cm * 4cm (Cruce de 2 barras)
const float AREA_GRILLA_AND = 256.0; // 16cm * 16cm (Cruce de 2 planos)
const float BIN_TIME = 600.0;        // Tiempo por bin en segundos (10 min)

// Estructura para organizar los 4 canales de una barra y sus umbrales
struct BarraCalibrada {
    int ch[4];
    float thr[4];
};

// --- FUNCIONES DE FLUJO CORREGIDAS ---
// 1. FLUJO BARRA INDIVIDUAL (4 canales en coincidencia)
void CalcularFlujoBarra(const char* targetFile, int barraID, bool esSuperior) {
    TFile *f = TFile::Open(targetFile, "UPDATE");
    TTree *T = (TTree*)f->Get("T");
    
    // Definimos las estructuras de canales y umbrales igual que en tus otras funciones
    std::vector<BarraCalibrada> mS = {{{0,26,1,27},{20,20,22,18}}, {{3,24,2,25},{20,17,20,19}}, {{6,28,7,29},{21,18,19,19}}, {{4,30,5,31},{19,20,20,20}}};
    std::vector<BarraCalibrada> mI = {{{8,17,9,16},{23,20,20,20}}, {{12,18,13,19},{18,21,21,22}}, {{10,20,11,21},{20,16,19,21}}, {{14,22,15,23},{18,20,17,22}}};
    
    // Seleccionamos la barra específica
    BarraCalibrada b = esSuperior ? mS[barraID] : mI[barraID];
    
    // Mapeamos los 4 canales de la barra
    Pulse p[4];
    for(int i=0; i<4; i++) T->SetBranchAddress(Form("p%d", b.ch[i]), &p[i]);

    T->GetEntry(0); double ts = p[0].fAcqTime;
    T->GetEntry(T->GetEntries()-1); double te = p[0].fAcqTime;

    // Nombre dinámico para el histograma
    TString name = Form("hFlux_Barra_%s_%d", esSuperior ? "Sup" : "Inf", barraID);
    TH1F *h = new TH1F(name, name + ";Hora;Flujo [cm^{-2}s^{-1}]", std::max(1, (int)std::ceil((te-ts)/BIN_TIME)), ts, te);

    for(Long64_t i=0; i<T->GetEntries(); i++) {
        T->GetEntry(i);
        // Lógica de coincidencia de los 4 canales de LA barra
        if(p[0].fMax > b.thr[0] && p[1].fMax > b.thr[1] && 
           p[2].fMax > b.thr[2] && p[3].fMax > b.thr[3]) {
            h->Fill(p[0].fAcqTime);
        }
    }

    ConfigurarEjeTiempo(h->GetXaxis());
    h->Scale(1.0 / (BIN_TIME * AREA_BARRA)); // AREA_BARRA = 500 cm2
    h->Write("", TObject::kOverwrite); 
    f->Close();
    
    std::cout << "[+] Flujo de Barra " << (esSuperior ? "Superior " : "Inferior ") << barraID << " calculado." << std::endl;
}

// 2. FLUJO PLANO OR (Máxima Estadística)
void CalcularFlujoPlanoOR(const char* targetFile, bool esSuperior) {
    TFile *f = TFile::Open(targetFile, "UPDATE");
    TTree *T = (TTree*)f->Get("T");
    Pulse p[32]; 
    for(int i=0; i<32; i++) T->SetBranchAddress(Form("p%d", i), &p[i]);

    // Definición de barras (4 canales cada una) con umbrales de tu lista de pares
    std::vector<BarraCalibrada> mS = {{{0,26,1,27},{20,20,22,18}}, {{3,24,2,25},{20,17,20,19}}, {{6,28,7,29},{21,18,19,19}}, {{4,30,5,31},{19,20,20,20}}};
    std::vector<BarraCalibrada> mI = {{{8,17,9,16},{23,20,20,20}}, {{12,18,13,19},{18,21,21,22}}, {{10,20,11,21},{20,16,19,21}}, {{14,22,15,23},{18,20,17,22}}};
    std::vector<BarraCalibrada> &m = esSuperior ? mS : mI;

    T->GetEntry(0); double ts = p[0].fAcqTime;
    T->GetEntry(T->GetEntries()-1); double te = p[0].fAcqTime;
    TH1F *h = new TH1F(esSuperior ? "hFlux_Sup_OR" : "hFlux_Inf_OR", "Flujo Plano OR;Hora;Flujo", std::max(1,(int)std::ceil((te-ts)/BIN_TIME)), ts, te);

    for(Long64_t i=0; i<T->GetEntries(); i++) {
        T->GetEntry(i);
        bool hitPlano = false;
        for(int b=0; b<4; b++) {
            if(p[m[b].ch[0]].fMax > m[b].thr[0] && p[m[b].ch[1]].fMax > m[b].thr[1] &&
               p[m[b].ch[2]].fMax > m[b].thr[2] && p[m[b].ch[3]].fMax > m[b].thr[3]) {
                hitPlano = true; break;
            }
        }
        if(hitPlano) h->Fill(p[0].fAcqTime);
    }
    ConfigurarEjeTiempo(h->GetXaxis());
    h->Scale(1.0 / (BIN_TIME * AREA_PLANO_OR));
    h->Write("", TObject::kOverwrite); f->Close();
}

// 3. FLUJO PÍXEL (AND 2 barras cruzadas - 4 canales por barra)
void CalcularFlujoPixel(const char* targetFile, int idBarraSup, int idBarraInf) {
    TFile *f = TFile::Open(targetFile, "UPDATE");
    TTree *T = (TTree*)f->Get("T");
    Pulse p[32];
    for(int i=0; i<32; i++) T->SetBranchAddress(Form("p%d", i), &p[i]);
    
    std::vector<BarraCalibrada> mS = {{{0,26,1,27},{20,20,22,18}}, {{3,24,2,25},{20,17,20,19}}, {{6,28,7,29},{21,18,19,19}}, {{4,30,5,31},{19,20,20,20}}};
    std::vector<BarraCalibrada> mI = {{{8,17,9,16},{23,20,20,20}}, {{12,18,13,19},{18,21,21,22}}, {{10,20,11,21},{20,16,19,21}}, {{14,22,15,23},{18,20,17,22}}};

    T->GetEntry(0); double ts = p[0].fAcqTime;
    T->GetEntry(T->GetEntries()-1); double te = p[0].fAcqTime;
    TH1F *h = new TH1F(Form("hFluxPixel_S%d_I%d", idBarraSup, idBarraInf), "Flujo Pixel;Hora;Flujo", std::max(1,(int)std::ceil((te-ts)/BIN_TIME)), ts, te);

    for(Long64_t i=0; i<T->GetEntries(); i++) {
        T->GetEntry(i);
        bool hitS = (p[mS[idBarraSup].ch[0]].fMax > mS[idBarraSup].thr[0] && p[mS[idBarraSup].ch[1]].fMax > mS[idBarraSup].thr[1] &&
                     p[mS[idBarraSup].ch[2]].fMax > mS[idBarraSup].thr[2] && p[mS[idBarraSup].ch[3]].fMax > mS[idBarraSup].thr[3]);
        bool hitI = (p[mI[idBarraInf].ch[0]].fMax > mI[idBarraInf].thr[0] && p[mI[idBarraInf].ch[1]].fMax > mI[idBarraInf].thr[1] &&
                     p[mI[idBarraInf].ch[2]].fMax > mI[idBarraInf].thr[2] && p[mI[idBarraInf].ch[3]].fMax > mI[idBarraInf].thr[3]);
        if(hitS && hitI) h->Fill(p[0].fAcqTime);
    }
    ConfigurarEjeTiempo(h->GetXaxis());
    h->Scale(1.0 / (BIN_TIME * AREA_PIXEL_AND));
    h->Write("", TObject::kOverwrite); f->Close();
}

// 4. FLUJO VERTICAL AND (Grilla Completa - AND de Planos)
void CalcularFlujoVerticalAND(const char* targetFile) {
    TFile *f = TFile::Open(targetFile, "UPDATE");
    TTree *T = (TTree*)f->Get("T");
    Pulse p[32];
    for(int i=0; i<32; i++) T->SetBranchAddress(Form("p%d", i), &p[i]);

    std::vector<BarraCalibrada> mS = {{{0,26,1,27},{20,20,22,18}}, {{3,24,2,25},{20,17,20,19}}, {{6,28,7,29},{21,18,19,19}}, {{4,30,5,31},{19,20,20,20}}};
    std::vector<BarraCalibrada> mI = {{{8,17,9,16},{23,20,20,20}}, {{12,18,13,19},{18,21,21,22}}, {{10,20,11,21},{20,16,19,21}}, {{14,22,15,23},{18,20,17,22}}};

    T->GetEntry(0); double ts = p[0].fAcqTime;
    T->GetEntry(T->GetEntries()-1); double te = p[0].fAcqTime;
    TH1F *h = new TH1F("hFlux_Grilla_AND", "Flujo Vertical Grilla;Hora;Flujo", std::max(1,(int)std::ceil((te-ts)/BIN_TIME)), ts, te);

    for(Long64_t i=0; i<T->GetEntries(); i++) {
        T->GetEntry(i);
        bool hitS = false, hitI = false;
        for(int b=0; b<4; b++) {
            if(p[mS[b].ch[0]].fMax > mS[b].thr[0] && p[mS[b].ch[1]].fMax > mS[b].thr[1] && p[mS[b].ch[2]].fMax > mS[b].thr[2] && p[mS[b].ch[3]].fMax > mS[b].thr[3]) { hitS = true; break; }
        }
        for(int b=0; b<4; b++) {
            if(p[mI[b].ch[0]].fMax > mI[b].thr[0] && p[mI[b].ch[1]].fMax > mI[b].thr[1] && p[mI[b].ch[2]].fMax > mI[b].thr[2] && p[mI[b].ch[3]].fMax > mI[b].thr[3]) { hitI = true; break; }
        }
        if(hitS && hitI) h->Fill(p[0].fAcqTime);
    }
    ConfigurarEjeTiempo(h->GetXaxis());
    h->Scale(1.0 / (BIN_TIME * AREA_GRILLA_AND));
    h->Write("", TObject::kOverwrite); f->Close();
}
// --- MAIN ACTUALIZADO ---

int main(int argc, char** argv) {
    setenv("TZ", "America/Santiago", 1);
    tzset();
    if (argc >= 3) {
        ProcesarPSA(argv[1], argv[2]);
        return 0;
    }
    int opcion;
    std::cout << "\n--- SOFTWARE DE ANÁLISIS DE PULSOS ---" << std::endl;
    std::cout << "1.  Procesamiento PSA (Raw -> Root)" << std::endl;
    std::cout << "2.  Análisis Individual (Rate)" << std::endl;
    std::cout << "3.  Correlación entre 4 canales" << std::endl;
    std::cout << "4.  Rate Coincidencia AND (2 canales)" << std::endl;
    std::cout << "5.  Rate Coincidencia AND (3 canales)" << std::endl;
    std::cout << "6.  Rate Coincidencia AND (4 canales)" << std::endl;
    std::cout << "7.  Rate OR Barras Superior (S1||S2||S3||S4)" << std::endl;
    std::cout << "8.  Rate OR Barras Inferior (I1||I2||I3||I4)" << std::endl;
    std::cout << "9.  Rate Coincidencia Planos (Sup && Inf)" << std::endl;
    std::cout << "10. Histogramas de Amplitud Máxima Condicionados" << std::endl;
    std::cout << "11. Estabilidad vs Temperatura (fMax vs Tiempo)" << std::endl;
    std::cout << "12. Delta T: Un solo canal" << std::endl;
    std::cout << "13. Delta T: Coincidencia AND (2 canales)" << std::endl;
    std::cout << "14. Delta T: Coincidencia AND (3 canales)" << std::endl;
    std::cout << "15. Delta T: Coincidencia AND (4 canales)" << std::endl;
    std::cout << "16. Evolución de Voltaje (mV vs Tiempo)" << std::endl;
    std::cout << "17. Amplitud Ch Objetivo en Coincidencia DOBLE" << std::endl;
    std::cout << "18. Amplitud Ch Objetivo en Coincidencia TRIPLE" << std::endl;
    std::cout << "19. Amplitud Ch Objetivo en Coincidencia CUADRUPLE" << std::endl;
    std::cout << "20. Evolución de Carga (mV ns vs Tiempo)" << std::endl;
    std::cout << "21. Carga Ch Objetivo en Coincidencia DOBLE" << std::endl;
    std::cout << "22. Carga Ch Objetivo en Coincidencia TRIPLE" << std::endl;
    std::cout << "23. Carga Ch Objetivo en Coincidencia CUADRUPLE" << std::endl;
    std::cout << "--- OPCIONES DE FLUJO (NORMALIZADO POR ÁREA) ---" << std::endl;
    std::cout << "24. Flujo Individual [cm-2 s-1]" << std::endl;
    std::cout << "25. Flujo Plano OR (Superior e Inferior)" << std::endl;
    std::cout << "26. Flujo Píxel (AND 2 barras cruzadas 4ch+4ch)" << std::endl;
    std::cout << "27. Flujo Vertical AND (Grilla 16x16)" << std::endl;
    std::cout << "0.  Salir" << std::endl;
    std::cout << "--------------------------------------" << std::endl;
    std::cout << "Seleccione: ";
    
    if (!(std::cin >> opcion) || opcion == 0) return 0;

    std::string file;
    if (opcion >= 2 && opcion <= 27) {
        std::cout << "Archivo .root a analizar: ";
        std::cin >> file;
    }

    switch(opcion) {
        case 1: {
            std::string in, out;
            std::cout << "Archivo binario entrada: "; std::cin >> in;
            std::cout << "Archivo .root salida: "; std::cin >> out;
            ProcesarPSA(in.c_str(), out.c_str());
            break;
        }
        case 2: {
            int ch; float thr;
            std::cout << "Canal: "; std::cin >> ch;
            std::cout << "Umbral para Ch" << ch << ": "; std::cin >> thr;
            CalcularRate(file.c_str(), ch, thr);
            break;
        }
        case 3: {
            std::vector<int> chs(4);
            for(int i=0; i<4; i++) {
                std::cout << "Canal " << i+1 << ": "; std::cin >> chs[i];
            }
            GenerarCorrelacionCanales(file.c_str(), chs);
            break;
        }
        case 4: {
            int c1, c2; float thr;
            std::cout << "Canal 1: "; std::cin >> c1;
            std::cout << "Canal 2: "; std::cin >> c2;
            std::cout << "Umbral común (ADC): "; std::cin >> thr;
            CalcularRateCoincidencia(file.c_str(), c1, c2, thr);
            break;
        }
        case 5: { 
            std::vector<int> chs(3); float thr;
            for(int i=0; i<3; i++) { std::cout << "Canal " << i+1 << ": "; std::cin >> chs[i]; }
            std::cout << "Umbral común (ADC): "; std::cin >> thr;
            CalcularRateCoincidenciaTriple(file.c_str(), chs, thr);
            break;
        }
        case 6: {
            std::vector<int> chs(4); float thr;
            for(int i=0; i<4; i++) { std::cout << "Canal " << i+1 << ": "; std::cin >> chs[i]; }
            std::cout << "Umbral común (ADC): "; std::cin >> thr;
            CalcularRateCoincidenciaCuadruple(file.c_str(), chs, thr);
            break;
        }
        case 7: {
            float thr; std::cout << "Umbral común para OR Superior: "; std::cin >> thr;
            CalcularRateORBarrasSuperior(file.c_str(), thr);
            break;
        }
        case 8: {
            float thr; std::cout << "Umbral común para OR Inferior: "; std::cin >> thr;
            CalcularRateORBarrasInferior(file.c_str(), thr);
            break;
        }
        case 9: {
            float thr; std::cout << "Umbral común para Coincidencia Planos: "; std::cin >> thr;
            CalcularRateCoincidenciaPlanos(file.c_str(), thr);
            break;
        }
        case 10: { 
            GenerarHistogramasMaxAmpCondicionados(file.c_str());
            break;
        }
        case 11: {
            int ch; std::cout << "Canal a monitorear: "; std::cin >> ch;
            AnalizarEstabilidadTemperatura(file.c_str(), ch);
            break;
        }
        case 12: {
            int ch; float thr;
            std::cout << "Canal: "; std::cin >> ch;
            std::cout << "Umbral para Ch" << ch << ": "; std::cin >> thr;
            AnalizarDeltaTime(file.c_str(), ch, thr);
            break;
        }
        case 13: { 
            int c1, c2; float t1, t2;
            std::cout << "Canal 1: "; std::cin >> c1;
            std::cout << "Umbral para Ch" << c1 << ": "; std::cin >> t1;
            std::cout << "Canal 2: "; std::cin >> c2;
            std::cout << "Umbral para Ch" << c2 << ": "; std::cin >> t2;
            AnalizarDeltaTimeCoincidencia(file.c_str(), c1, c2, t1, t2);
            break;
        }
        case 14: { 
            std::vector<int> chs(3); std::vector<float> thrs(3);
            for(int i=0; i<3; i++) {
                std::cout << "Canal " << i+1 << ": "; std::cin >> chs[i];
                std::cout << "Umbral para Ch" << chs[i] << ": "; std::cin >> thrs[i];
            }
            AnalizarDeltaTimeCoincidenciaTriple(file.c_str(), chs, thrs);
            break;
        }
        case 15: { 
            std::vector<int> chs(4); std::vector<float> thrs(4);
            for(int i=0; i<4; i++) {
                std::cout << "Canal " << i+1 << ": "; std::cin >> chs[i];
                std::cout << "Umbral para Ch" << chs[i] << ": "; std::cin >> thrs[i];
            }
            AnalizarDeltaTimeCoincidenciaCuadruple(file.c_str(), chs, thrs);
            break;
        }
        case 16: {
            int ch; float thr;
            std::cout << "Canal: "; std::cin >> ch;
            std::cout << "Umbral para Ch" << ch << ": "; std::cin >> thr;
            AnalizarVoltajeTemporal(file.c_str(), ch, thr);
            break;
        }
        case 17: {
            int target, coin; float tTarget, tCoin;
            std::cout << "Canal Objetivo: "; std::cin >> target;
            std::cout << "Umbral para Ch" << target << ": "; std::cin >> tTarget;
            std::cout << "Canal Coincidencia: "; std::cin >> coin;
            std::cout << "Umbral para Ch" << coin << ": "; std::cin >> tCoin;
            AnalizarVoltajeCoincidenciaDoble(file.c_str(), target, tTarget, coin, tCoin);
            break;
        }
        case 18: {
            int target; std::vector<int> others(2); 
            float tTarget; std::vector<float> tOthers(2);
            std::cout << "Canal Objetivo: "; std::cin >> target;
            std::cout << "Umbral para Ch" << target << ": "; std::cin >> tTarget;
            for(int i=0; i<2; i++) {
                std::cout << "Canal Coincidencia " << i+1 << ": "; std::cin >> others[i];
                std::cout << "Umbral para Ch" << others[i] << ": "; std::cin >> tOthers[i];
            }
            AnalizarVoltajeCoincidenciaTriple(file.c_str(), target, tTarget, others, tOthers);
            break;
        }
        case 19: {
            int target; std::vector<int> others(3); 
            float tTarget; std::vector<float> tOthers(3);
            std::cout << "Canal Objetivo: "; std::cin >> target;
            std::cout << "Umbral para Ch" << target << ": "; std::cin >> tTarget;
            for(int i=0; i<3; i++) {
                std::cout << "Canal Coincidencia " << i+1 << ": "; std::cin >> others[i];
                std::cout << "Umbral para Ch" << others[i] << ": "; std::cin >> tOthers[i];
            }
            AnalizarVoltajeCoincidenciaCuadruple(file.c_str(), target, tTarget, others, tOthers);
            break;
        }
        case 20: {
            int ch; float thr;
            std::cout << "Canal: "; std::cin >> ch;
            std::cout << "Umbral para Ch" << ch << ": "; std::cin >> thr;
            AnalizarCargaTemporal(file.c_str(), ch, thr);
            break;
        }
        case 21: {
            int target, coin; float tTarget, tCoin;
            std::cout << "Canal Objetivo: "; std::cin >> target;
            std::cout << "Umbral para Ch" << target << ": "; std::cin >> tTarget;
            std::cout << "Canal Coincidencia: "; std::cin >> coin;
            std::cout << "Umbral para Ch" << coin << ": "; std::cin >> tCoin;
            AnalizarCargaCoincidenciaDoble(file.c_str(), target, tTarget, coin, tCoin);
            break;
        }
        case 22: {
            int target; std::vector<int> others(2); 
            float tTarget; std::vector<float> tOthers(2);
            std::cout << "Canal Objetivo: "; std::cin >> target;
            std::cout << "Umbral para Ch" << target << ": "; std::cin >> tTarget;
            for(int i=0; i<2; i++) {
                std::cout << "Canal Coincidencia " << i+1 << ": "; std::cin >> others[i];
                std::cout << "Umbral para Ch" << others[i] << ": "; std::cin >> tOthers[i];
            }
            AnalizarCargaCoincidenciaTriple(file.c_str(), target, tTarget, others, tOthers);
            break;
        }
        case 23: {
            int target; std::vector<int> others(3); 
            float tTarget; std::vector<float> tOthers(3);
            std::cout << "Canal Objetivo: "; std::cin >> target;
            std::cout << "Umbral para Ch" << target << ": "; std::cin >> tTarget;
            for(int i=0; i<3; i++) {
                std::cout << "Canal Coincidencia " << i+1 << ": "; std::cin >> others[i];
                std::cout << "Umbral para Ch" << others[i] << ": "; std::cin >> tOthers[i];
            }
            AnalizarCargaCoincidenciaCuadruple(file.c_str(), target, tTarget, others, tOthers);
            break;
        }
        case 24: { // Flujo de una Barra Individual (Coincidencia de sus 4 canales)
            int barraID, plano;
            std::cout << "Plano (1: Superior, 0: Inferior): "; std::cin >> plano;
            std::cout << "ID de Barra (0-3): "; std::cin >> barraID;
            if(barraID < 0 || barraID > 3) {
                std::cout << "ID no válido." << std::endl;
            } else {
                CalcularFlujoBarra(file.c_str(), barraID, (plano == 1));
            }
            break;
        }
        case 25: { // Flujo Plano OR (Máxima Estadística)
            std::cout << "Calculando Flujo Plano OR para ambos planos..." << std::endl;
            // Ejecuta ambos automáticamente o puedes preguntar cuál
            CalcularFlujoPlanoOR(file.c_str(), true);  // Superior
            CalcularFlujoPlanoOR(file.c_str(), false); // Inferior
            break;
        }
        case 26: { // Flujo Píxel (AND entre dos barras cruzadas)
            int bS, bI;
            std::cout << "ID Barra Superior (0-3): "; std::cin >> bS;
            std::cout << "ID Barra Inferior (0-3): "; std::cin >> bI;
            if(bS >= 0 && bS <= 3 && bI >= 0 && bI <= 3) {
                CalcularFlujoPixel(file.c_str(), bS, bI);
            } else {
                std::cout << "IDs fuera de rango." << std::endl;
            }
            break;
        }
        case 27: { // Flujo Vertical AND (Grilla completa)
            std::cout << "Calculando Flujo Vertical Grilla (AND de Planos)..." << std::endl;
            CalcularFlujoVerticalAND(file.c_str()); 
            break;
        }
        default:
            std::cout << "Opción no válida o saliendo." << std::endl;
            break;
    }
    return 0;
}
