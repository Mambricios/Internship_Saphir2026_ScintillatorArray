#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <iostream>

void plot_evolucion(string fecha = "20251121") {
    // Configuración del Canvas
    TCanvas *c1 = new TCanvas("c1", "Evolucion del Rate", 1200, 600);
    c1->SetGrid();

    // El TGraphErrors guardará todos los puntos de todos los archivos del día
    TGraphErrors *gRate = new TGraphErrors();
    gRate->SetTitle(Form("Evolucion del Rate (Ch0) - Dia %s;Hora del Dia;Rate (Hz)", fecha.c_str()));

    // Definir la carpeta según la estructura del nuevo procesar_dia.sh
    string folder = "analisis_" + fecha;
    
    // Comando para listar archivos, el 2>/dev/null evita mensajes feos si no hay archivos
    TString comando = Form("ls %s/analisis_pulsos_acq_*.root 2>/dev/null", folder.c_str());
    TString lista_archivos = gSystem->GetFromPipe(comando);
    TObjArray *files = lista_archivos.Tokenize("\n");

    if (files->GetEntries() == 0) {
        std::cerr << "!!! ERROR: No se encontraron archivos en la carpeta: " << folder << std::endl;
        std::cerr << "Asegurate de haber corrido ./procesar_dia.sh " << fecha << std::endl;
        return;
    }

    int pointIdx = 0;
    for (int i = 0; i < files->GetEntries(); ++i) {
        string fname = ((TObjString*)files->At(i))->GetString().Data();
        TFile *f = TFile::Open(fname.c_str());
        
        if (!f || f->IsZombie()) continue;

        // Extraemos el histograma de 10 min y el árbol del canal 0
        TH1F *h = (TH1F*)f->Get("hRate10");
        TTree *t = (TTree*)f->Get("p0");
        
        if (!h || !t) {
            std::cout << "Saltando archivo (faltan objetos): " << fname << std::endl;
            f->Close(); 
            continue; 
        }

        // Obtener el tiempo de inicio real del archivo (Unix Time en microsegundos)
        double f_acqTime;
        t->SetBranchAddress("f_acqTime", &f_acqTime);
        t->GetEntry(0);
        double startTimeSec = f_acqTime / 1e6;

        // Pasar los bins del histograma al TGraph
        for (int b = 1; b <= h->GetNbinsX(); ++b) {
            double rate = h->GetBinContent(b);
            double errY = h->GetBinError(b);
            
            if (rate > 0) {
                // Eje X: Tiempo absoluto (Inicio archivo + posición del bin en segundos)
                double binCenterSec = (h->GetBinCenter(b) * 60.0) + startTimeSec;
                
                // Error X: 10 minutos son 600 segundos, la barra mide 300 a cada lado
                double errX = (h->GetBinWidth(b) * 60.0) / 2.0;

                gRate->SetPoint(pointIdx, binCenterSec, rate);
                gRate->SetPointError(pointIdx, errX, errY);
                pointIdx++;
            }
        }
        f->Close();
    }

    // --- Estética del gráfico ---
    gRate->SetMarkerStyle(20);
    gRate->SetMarkerSize(0.8);
    gRate->SetMarkerColor(kAzure+2);
    gRate->SetLineColor(kAzure+2);
    gRate->SetLineWidth(1);

    // Ajustar el eje Y para que siempre empiece en 0
    gRate->SetMinimum(0.0);

    // --- Configuración del Eje X de Tiempo ---
    gRate->GetXaxis()->SetTimeDisplay(1);
    gRate->GetXaxis()->SetTimeFormat("%H:%M");
    gRate->GetXaxis()->SetTimeOffset(0, "local"); // Ajusta a tu zona horaria local
    gRate->GetXaxis()->SetLabelSize(0.035);
    gRate->GetXaxis()->SetTitleOffset(1.2);

    // Dibujar
    gRate->Draw("AP"); 

    // Guardar resultado
    c1->SaveAs(Form("Plot_Rate_%s_Final.png", fecha.c_str()));
    
    std::cout << "Done! Se procesaron " << pointIdx << " puntos de rate." << std::endl;
}