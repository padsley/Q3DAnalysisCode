{

    gROOT->ProcessLine(".L Kinematics.cpp+");

    TFile *f232 = TFile::Open("/home/padsley/data/MunichQ3D_30Si_3He_d/run232_an.root");
    TTree *t232 = (TTree*)f232->Get("readout");
    
    gROOT->ProcessLine("/home/padsley/data/MunichQ3D_30Si_3He_d/CUT_anodeanode1_run232.C");
    gROOT->ProcessLine("/home/padsley/data/MunichQ3D_30Si_3He_d/CUT_anodepm_run232.C");
    
    t232->Draw("BrhoToEx(PosToBrho(pos))>>hEx(2000,6,10)","CUT_anodeanode1_run232 && CUT_anodepm_run232 && pos>0");
    
}