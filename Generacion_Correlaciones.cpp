#include <TFile.h>
#include <TTree.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <iostream>

void generar_correlaciones(int canal_ref) {
    const int NUM_CHANNELS = 32;
    
    if (canal_ref < 0 || canal_ref >= NUM_CHANNELS) {
        std::cout << "Error: Canal de referencia inválido." << std::endl;
        return;
    }

    // 1. Abrimos en modo READ porque ya no vamos a escribir nada
    TFile *f = new TFile("analisis_pulsos.root", "READ");
    if (!f || f->IsZombie()) {
        std::cout << "Error: No se pudo abrir analisis_pulsos.root." << std::endl;
        return;
    }

    // Preparamos los Canvas antes del loop
    TCanvas *cInt = new TCanvas(Form("cInt_Ch%d", canal_ref), "Correlaciones Integrales", 1800, 1200);
    cInt->Divide(6, 6);

    TCanvas *cMax = new TCanvas(Form("cMax_Ch%d", canal_ref), "Correlaciones Amplitud Max", 1800, 1200);
    cMax->Divide(6, 6);

    int pad = 1;

    for (int i = 0; i < NUM_CHANNELS; ++i) {
        if (i == canal_ref) continue;

        TTree *t1 = (TTree*)f->Get(Form("p%d", canal_ref));
        TTree *t2 = (TTree*)f->Get(Form("p%d", i));

        if (!t1 || !t2) continue;

        float f_int1, f_int2, f_max1, f_max2;
        t1->SetBranchAddress("f_int", &f_int1);
        t2->SetBranchAddress("f_int", &f_int2);
        t1->SetBranchAddress("f_max", &f_max1);
        t2->SetBranchAddress("f_max", &f_max2);

        // 2. Creamos los histogramas
        TH2F *hCorrInt = new TH2F(Form("hCorr_Int_Ch%d_Ch%d", canal_ref, i), 
                                  Form("Ch%d vs Ch%d (Int);Ch%d;Ch%d", canal_ref, i, canal_ref, i), 
                                  200, 0, 50000, 200, 0, 50000);

        TH2F *hCorrMax = new TH2F(Form("hCorr_Max_Ch%d_Ch%d", canal_ref, i), 
                                  Form("Ch%d vs Ch%d (Max);Ch%d;Ch%d", canal_ref, i, canal_ref, i), 
                                  200, 0, 1500, 200, 0, 1500);

        // 3. ¡IMPORTANTE! Desvincular del archivo .root para que no se guarden
        hCorrInt->SetDirectory(0);
        hCorrMax->SetDirectory(0);

        Long64_t nEntries = t1->GetEntries();
        for (Long64_t j = 0; j < nEntries; ++j) {
            t1->GetEntry(j);
            t2->GetEntry(j);
            hCorrInt->Fill(f_int1, f_int2);
            hCorrMax->Fill(f_max1, f_max2);
        }

        // 4. Dibujar en los pads correspondientes
        cInt->cd(pad);
        hCorrInt->Draw("COLZ");

        cMax->cd(pad);
        hCorrMax->Draw("COLZ");

        pad++;
    }

    // 5. Actualizar los canvas
    cInt->Update();
    cMax->Update();

    // Nota: No cerramos el archivo 'f' inmediatamente si los árboles son necesarios 
    // para el dibujado (aunque en TH2F al estar llenos no debería haber problema).
    // Si los canvas se quedan en blanco, mueve f->Close() al final de tu programa.
    std::cout << "Correlaciones generadas en memoria y mostradas en pantalla." << std::endl;
}