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

    // --- (Tu lista de pares y thresholds se mantiene igual) ---
    std::vector<ParUmbral> pares = {
        {0, 26, 20.0, 20.0}, {1, 27, 22.0, 18.0}, {3, 24, 20.0, 17.0}, {2, 25, 20.0, 19.0},
        {6, 28, 21.0, 18.0}, {7, 29, 19.0, 19.0}, {4, 30, 19.0, 20.0}, {5, 31, 20.0, 20.0},
        {8, 17, 23.0, 20.0}, {9, 16, 20.0, 20.0}, {12, 18, 18.0, 21.0}, {13, 19, 21.0, 22.0},
        {10, 20, 20.0, 16.0}, {11, 21, 19.0, 21.0}, {14, 22, 18.0, 20.0}, {15, 23, 17.0, 22.0}
    };

    Pulse pA, pB;
    // --- Lógica de llenado original (Sin cambios) ---
    for (auto const& par : pares) {
        T->ResetBranchAddresses();
        dirCond->cd(); 
        TH1F *hMaxA_cond = new TH1F(Form("hMax_Ch%d_cond_Ch%d", par.chA, par.chB), 
                                    Form("Max %d cond %d;ADC;Entries", par.chA, par.chB), 200, 0, 200);
        TH1F *hMaxB_cond = new TH1F(Form("hMax_Ch%d_cond_Ch%d", par.chB, par.chA), 
                                    Form("Max %d cond %d;ADC;Entries", par.chB, par.chA), 200, 0, 200);

        T->SetBranchAddress(Form("p%d", par.chA), &pA);
        T->SetBranchAddress(Form("p%d", par.chB), &pB);

        Long64_t nEntries = T->GetEntries();
        for (Long64_t i = 0; i < nEntries; ++i) {
            T->GetEntry(i);
            if (pB.fMax > par.thrA) hMaxA_cond->Fill(pA.fMax);
            if (pA.fMax > par.thrB) hMaxB_cond->Fill(pB.fMax);
        }
        hMaxA_cond->Write("", TObject::kOverwrite);
        hMaxB_cond->Write("", TObject::kOverwrite);
        delete hMaxA_cond; delete hMaxB_cond;
    }

    // =========================================================================
    // NUEVA FUNCIONALIDAD: Cascada Diagonal (Efecto Isométrico)
    // =========================================================================
    
auto DibujarCascada = [&](const std::string& nombrePlano, int inicio, int fin) {
    TCanvas *c = new TCanvas(Form("c_%s", nombrePlano.c_str()), nombrePlano.c_str(), 1000, 900);
    c->SetLeftMargin(0.12);
    c->SetBottomMargin(0.12);
    gStyle->SetOptStat(0);
    
    TLegend *leg = new TLegend(0.15, 0.60, 0.45, 0.88); 
    leg->SetBorderSize(1);
    leg->SetFillColor(0);
    leg->SetTextSize(0.022);
    leg->SetHeader("Condiciones de Coincidencia", "C");

    double dx = 12.0; 
    double dy = 15.0; 
    int colores[] = {kBlue+1, kRed+1, kGreen+2, kOrange+2, kMagenta+1, kCyan+2, kAzure+7, kViolet-3};

    TLine *lDepth0 = new TLine(0, 0, (fin-inicio)*dx, (fin-inicio)*dy);
    lDepth0->SetLineColor(kGray+2);
    lDepth0->SetLineStyle(2);

    for (int i = inicio; i <= fin; ++i) {
        TH1F *h = (TH1F*)dirCond->Get(Form("hMax_Ch%d_cond_Ch%d", pares[i].chA, pares[i].chB));
        if (!h) continue;

        TGraph *g = new TGraph();
        double xShift = (i - inicio) * dx;
        double yShift = (i - inicio) * dy;

        int pt = 0;
        for (int b = 1; b <= h->GetNbinsX(); b++) {
            double x = h->GetBinCenter(b);
            double y = h->GetBinContent(b);
            if (b == 1) g->SetPoint(pt++, x + xShift, yShift); 
            g->SetPoint(pt++, x + xShift, y + yShift);
            if (b == h->GetNbinsX()) g->SetPoint(pt++, x + xShift, yShift);
        }

        int colorActual = colores[(i - inicio) % 8];
        g->SetLineColor(colorActual);
        g->SetLineWidth(2);
        leg->AddEntry(g, Form("Ch%d (si Ch%d > %.1f)", pares[i].chA, pares[i].chB, pares[i].thrA), "l");

        if (i == inicio) {
            g->Draw("AL"); 
            g->GetXaxis()->SetTitle("ADC (con offset)");
            g->GetYaxis()->SetTitle("Entries (con offset)");
            
            double xMaxVisual = 200 + (fin - inicio) * dx;
            g->GetXaxis()->SetRangeUser(-10, xMaxVisual + 20);
            
            double yMaxVisual = h->GetMaximum() + (fin - inicio) * dy;
            g->GetYaxis()->SetRangeUser(-5, yMaxVisual + 50); 
            
            lDepth0->Draw();
        } else {
            g->Draw("L SAME");
        }

        TLine *base = new TLine(xShift, yShift, 200 + xShift, yShift);
        base->SetLineColor(kGray + 1);
        base->SetLineStyle(2);
        base->Draw();
    }
    
    leg->Draw(); 
    // --- EL CAMBIO ESTÁ AQUÍ ---
    // kOverwrite evita que se creen múltiples versiones (;1, ;2, etc.)
    c->Write("", TObject::kOverwrite); 
    
    // Opcional: También podrías querer borrar el canvas de memoria 
    // después de escribirlo para que no se acumule si llamas a la función muchas veces
    // delete c; 
};

    dirCond->cd();
    DibujarCascada("Cascada_Superior", 0, 7);
    DibujarCascada("Cascada_Inferior", 8, 15);

    T->ResetBranchAddresses();
    f->Close();
    std::cout << "Gráficas en cascada diagonal generadas con éxito." << std::endl;
}

void AnalizarEstabilidadTemperatura(const char* targetFile, int canal) {
    TFile *f = TFile::Open(targetFile, "UPDATE");
    if (!f || f->IsZombie()) return;
    TTree *T = (TTree*)f->Get("T");
    if (!T) { f->Close(); return; }

    Pulse p;
    T->SetBranchAddress(Form("p%d", canal), &p);

    // Obtener rango de tiempo
    T->GetEntry(0); double t_start = p.fAcqTime;
    T->GetEntry(T->GetEntries() - 1); double t_end = p.fAcqTime;

    // Bins de 30 minutos (1800 segundos) para tener buena estadística por punto
    int nBins = std::max(1, (int)std::ceil((t_end - t_start) / 1800.0));
    
    // Perfil: Calcula automáticamente la media de Y para cada bin de X
    TProfile *hProf = new TProfile(Form("hTempProf_Ch%d", canal), 
                                   Form("Estabilidad Ch %d;Hora;Media fMax (ADC)", canal), 
                                   nBins, t_start, t_end);

    for (Long64_t i = 0; i < T->GetEntries(); ++i) {
        T->GetEntry(i);
        // Filtramos ruidos bajos para que la media no se sesgue
        if (p.fMax > 20.0) hProf->Fill(p.fAcqTime, p.fMax);
    }

    TCanvas *c = new TCanvas(Form("cTemp_Ch%d", canal), "Analisis Temperatura", 1000, 600);
    hProf->SetMarkerStyle(20);
    hProf->SetMarkerSize(0.8);
    hProf->SetMarkerColor(kBlue+1);
    hProf->SetLineColor(kBlue+1);
    
    ConfigurarEjeTiempo(hProf->GetXaxis());
    hProf->Draw("E1"); // Dibuja puntos con barras de error (error de la media)

    // Dibujar líneas verticales que marcan el horario del Aire Acondicionado (aprox)
    // Nota: Esto asume que tus datos cubren esas horas
    TLine *lOn = new TLine(); lOn->SetLineColor(kGreen+2); lOn->SetLineStyle(2);
    TLine *lOff = new TLine(); lOff->SetLineColor(kRed+2); lOff->SetLineStyle(2);

    // Aquí podrías automatizar la detección del día, pero por ahora marcamos
    // visualmente las 10:00 y las 19:00 si están en el rango.
    
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

    // Obtener rango de tiempo real para el eje X
    T->GetEntry(0); double t_start = p.fAcqTime;
    T->GetEntry(nEntries - 1); double t_end = p.fAcqTime;
    if (t_end <= t_start) t_end = t_start + 3600.0; 

    // Bins de 600 segundos (10 minutos) para seguir la estructura del Rate
    int nBins = std::max(1, (int)std::ceil((t_end - t_start) / 600.0));

    // Usamos TProfile con opción "s" para que el error sea la Desviación Estándar
    TProfile *hProfDT = new TProfile(Form("hDeltaT_Evol_Ch%d", canal), 
                                     Form("Tiempo entre eventos Ch %d (Thr > %.1f);Hora;#Delta t promedio [ms]", canal, threshold), 
                                     nBins, t_start, t_end, "s");
    
    if (!hProfDT) { f->Close(); return; }

    double lastTime = -1.0;
    for (Long64_t i = 0; i < nEntries; ++i) {
        T->GetEntry(i);
        // Solo consideramos pulsos que superen el umbral definido
        if (p.fMax > threshold) {
            if (lastTime > 0) {
                // Calculamos la diferencia de tiempo y convertimos a ms
                double dt_ms = (p.fAcqTime - lastTime) * 1000.0;
                
                // Filtro de seguridad para evitar saltos temporales incoherentes
                if (dt_ms > 0 && dt_ms < 3600000) { 
                    hProfDT->Fill(p.fAcqTime, dt_ms);
                }
            }
            lastTime = p.fAcqTime;
        }
    }

    // Configuración estética similar a los otros histogramas
    TCanvas *c = new TCanvas(Form("cDeltaT_Evol_Ch%d", canal), "Evolucion Delta T", 1200, 700);
    hProfDT->SetMarkerStyle(20);
    hProfDT->SetMarkerSize(0.8);
    hProfDT->SetMarkerColor(kSpring+9);
    hProfDT->SetLineColor(kSpring+9);
    
    ConfigurarEjeTiempo(hProfDT->GetXaxis());
    
    // Dibujamos con barras de error
    hProfDT->Draw("E1");

    // Guardado y cierre
    hProfDT->Write("", TObject::kOverwrite);
    c->Write("", TObject::kOverwrite);
    
    std::cout << "[+] Análisis Delta T (ms) completado para el canal " << canal << std::endl;
    f->Close();
}

void AnalizarDeltaTimeCoincidencia(const char* targetFile, int c1, int c2, float thr1, float thr2) {
    TFile *f = TFile::Open(targetFile, "UPDATE");
    TTree *T = (TTree*)f->Get("T");
    if(!T) { f->Close(); return; }

    Pulse p1, p2;
    T->SetBranchAddress(Form("p%d", c1), &p1);
    T->SetBranchAddress(Form("p%d", c2), &p2);

    T->GetEntry(0); double ts = p1.fAcqTime;
    T->GetEntry(T->GetEntries()-1); double te = p1.fAcqTime;
    int nBins = std::max(1, (int)std::ceil((te-ts)/1800.0));

    TProfile *h = new TProfile(Form("hDeltaT_Coin_%d_%d", c1, c2), 
                               Form("Delta T Coincidencia %d (Thr>%.1f) & %d (Thr>%.1f);Hora;#Delta t promedio [ms]", c1, thr1, c2, thr2), nBins, ts, te, "s");

    double lastTime = -1.0;
    for(Long64_t i=0; i<T->GetEntries(); i++) {
        T->GetEntry(i);
        if(p1.fMax > thr1 && p2.fMax > thr2) {
            if(lastTime > 0) h->Fill(p1.fAcqTime, (p1.fAcqTime - lastTime)*1000.0);
            lastTime = p1.fAcqTime;
        }
    }
    ConfigurarEjeTiempo(h->GetXaxis());
    h->SetLineColor(kAzure+1); h->SetMarkerStyle(20);
    h->Write("", TObject::kOverwrite); f->Close();
    std::cout << "[+] Delta T Coincidencia Doble guardado." << std::endl;
}

void AnalizarDeltaTimeCoincidenciaTriple(const char* targetFile, std::vector<int> chs, std::vector<float> thrs) {
    TFile *f = TFile::Open(targetFile, "UPDATE");
    TTree *T = (TTree*)f->Get("T");
    if(!T) { f->Close(); return; }

    Pulse p[3];
    for(int i=0; i<3; i++) T->SetBranchAddress(Form("p%d", chs[i]), &p[i]);

    T->GetEntry(0); double ts = p[0].fAcqTime;
    T->GetEntry(T->GetEntries()-1); double te = p[0].fAcqTime;
    int nBins = std::max(1, (int)std::ceil((te-ts)/1800.0));

    TProfile *h = new TProfile(Form("hDeltaT_Triple_%d_%d_%d", chs[0], chs[1], chs[2]), 
                               "Delta T Triple (Umbrales Indep.);Hora;#Delta t promedio [ms]", nBins, ts, te, "s");

    double lastTime = -1.0;
    for(Long64_t i=0; i<T->GetEntries(); i++) {
        T->GetEntry(i);
        if(p[0].fMax > thrs[0] && p[1].fMax > thrs[1] && p[2].fMax > thrs[2]) {
            if(lastTime > 0) h->Fill(p[0].fAcqTime, (p[0].fAcqTime - lastTime)*1000.0);
            lastTime = p[0].fAcqTime;
        }
    }
    ConfigurarEjeTiempo(h->GetXaxis());
    h->SetLineColor(kViolet+1); h->SetMarkerStyle(20);
    h->Write("", TObject::kOverwrite); f->Close();
}

void AnalizarDeltaTimeCoincidenciaCuadruple(const char* targetFile, std::vector<int> chs, std::vector<float> thrs) {
    TFile *f = TFile::Open(targetFile, "UPDATE");
    TTree *T = (TTree*)f->Get("T");
    if(!T) { f->Close(); return; }

    Pulse p[4];
    for(int i=0; i<4; i++) T->SetBranchAddress(Form("p%d", chs[i]), &p[i]);

    T->GetEntry(0); double ts = p[0].fAcqTime;
    T->GetEntry(T->GetEntries()-1); double te = p[0].fAcqTime;
    int nBins = std::max(1, (int)std::ceil((te-ts)/1800.0));

    TProfile *h = new TProfile(Form("hDeltaT_Quad_%d_%d_%d_%d", chs[0], chs[1], chs[2], chs[3]), 
                               "Delta T Cuadruple (Umbrales Indep.);Hora;#Delta t promedio [ms]", nBins, ts, te, "s");

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
    h->SetLineColor(kOrange+7); h->SetMarkerStyle(20);
    h->Write("", TObject::kOverwrite); f->Close();
}
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
    std::cout << "2.  Análisis Individual (Rate - Umbral 150)" << std::endl;
    std::cout << "3.  Correlación entre 4 canales" << std::endl;
    std::cout << "4.  Rate por Coincidencia AND (2 canales)" << std::endl;
    std::cout << "5.  Rate por Coincidencia AND (3 canales)" << std::endl;
    std::cout << "6.  Rate por Coincidencia AND (4 canales)" << std::endl;
    std::cout << "7.  Rate OR de Barras (S1 || S2 || S3 || S4)" << std::endl;
    std::cout << "8.  Rate OR de Barras (I1 || I2 || I3 || I4)" << std::endl;
    std::cout << "9.  Rate Coincidencia AND (Superior && Inferior)" << std::endl;
    std::cout << "10. Análisis Amplitud Condicionada (Pares de Canales)" << std::endl;
    std::cout << "11. Analizar estabilidad vs Temperatura (Media fMax vs Tiempo)" << std::endl;
    std::cout << "12. Delta T: Un solo canal" << std::endl;
    std::cout << "13. Delta T: Coincidencia AND (2 canales + Umbrales Indep.)" << std::endl;
    std::cout << "14. Delta T: Coincidencia AND (3 canales + Umbrales Indep.)" << std::endl;
    std::cout << "15. Delta T: Coincidencia AND (4 canales + Umbrales Indep.)" << std::endl;
    std::cout << "0.  Salir" << std::endl;
    std::cout << "--------------------------------------" << std::endl;
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
            std::vector<int> chs(3);
            std::string file;
            std::cout << "Archivo .root procesado: "; std::cin >> file;
            for(int i=0; i<3; i++) {
                std::cout << "Ingrese Canal " << i+1 << ": ";
                std::cin >> chs[i];
            }
            CalcularRateCoincidenciaTriple(file.c_str(), chs, 150.0);
            break;
        }
        case 6: {
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
        case 7: {
            std::string file;
            std::cout << "Archivo .root procesado: "; std::cin >> file;
            CalcularRateORBarrasSuperior(file.c_str(), 150.0);
            break;
        }
        case 8: {
            std::string file;
            std::cout << "Archivo .root procesado: "; std::cin >> file;
            CalcularRateORBarrasInferior(file.c_str(), 150.0);
            break;
        }
        case 9: {
            std::string file;
            std::cout << "Archivo .root procesado: "; std::cin >> file;
            CalcularRateCoincidenciaPlanos(file.c_str(), 150.0);
            break;
        }
        case 10: { 
            std::string file;
            std::cout << "Archivo .root procesado: "; std::cin >> file;
            GenerarHistogramasMaxAmpCondicionados(file.c_str());
            break;
        }
        // ... dentro del switch(opcion) ...
        case 11: {
            int ch; std::string file;
            std::cout << "Archivo .root: "; std::cin >> file;
            std::cout << "Canal a monitorear: "; std::cin >> ch;
            AnalizarEstabilidadTemperatura(file.c_str(), ch);
            break;
        }
        case 12: {
            int ch; std::string file; float thr;
            std::cout << "Archivo .root: "; std::cin >> file;
            std::cout << "Canal: "; std::cin >> ch;
            std::cout << "Threshold para el tiempo (ej 150): "; std::cin >> thr;
            AnalizarDeltaTime(file.c_str(), ch, thr);
            break;
        }
        case 13: { // Delta T 2 canales con thresholds individuales
            int c1, c2; float t1, t2; std::string file;
            std::cout << "Archivo .root: "; std::cin >> file;
            std::cout << "Canal A: "; std::cin >> c1;
            std::cout << "Threshold para Ch" << c1 << ": "; std::cin >> t1;
            std::cout << "Canal B: "; std::cin >> c2;
            std::cout << "Threshold para Ch" << c2 << ": "; std::cin >> t2;
            AnalizarDeltaTimeCoincidencia(file.c_str(), c1, c2, t1, t2);
            break;
        }
        case 14: { // Delta T 3 canales con thresholds individuales
            std::vector<int> chs(3); std::vector<float> thrs(3); std::string file;
            std::cout << "Archivo .root: "; std::cin >> file;
            for(int i=0; i<3; i++) {
                std::cout << "Canal " << i+1 << ": "; std::cin >> chs[i];
                std::cout << "Threshold para Ch" << chs[i] << ": "; std::cin >> thrs[i];
            }
            AnalizarDeltaTimeCoincidenciaTriple(file.c_str(), chs, thrs);
            break;
        }
        case 15: { // Delta T 4 canales con thresholds individuales
            std::vector<int> chs(4); std::vector<float> thrs(4); std::string file;
            std::cout << "Archivo .root: "; std::cin >> file;
            for(int i=0; i<4; i++) {
                std::cout << "Canal " << i+1 << ": "; std::cin >> chs[i];
                std::cout << "Threshold para Ch" << chs[i] << ": "; std::cin >> thrs[i];
            }
            AnalizarDeltaTimeCoincidenciaCuadruple(file.c_str(), chs, thrs);
            break;
        }
// ...
        default:
            std::cout << "Opción no válida." << std::endl;
            break;
    }
    return 0;
}