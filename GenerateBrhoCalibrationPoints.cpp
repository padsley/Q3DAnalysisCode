{
    //mass 3He  = 2809.413242
//mass d    = 1876.123740
//mass 30Si = 27920.38734
//mass 31P  = 28851.87339
    
    gROOT->ProcessLine(".L Kinematics.cpp+");
    
    vector<double> Ex;
    Ex.push_back(6.922);
    Ex.push_back(6.943);
    Ex.push_back(7.141);
    Ex.push_back(7.214);
    Ex.push_back(7.310);
    Ex.push_back(7.340);
    Ex.push_back(7.430);
    Ex.push_back(7.682);
    Ex.push_back(7.718);
    Ex.push_back(7.738);
    Ex.push_back(7.780);
    
    //ExToBrho(double TBeam, double Q3DAngle, double TargetAngle, double TargetThickness, double CarbonThickness, double MassBeam, double MassTarget, double MassEjectile, double MassRecoil, double Ex)
    
    for(unsigned int i=0;i<Ex.size();i++)
    {
        cout << "Ex = " << Ex.at(i) << " MeV \t Brho = " << ExToBrho(25., 20., 0., 40./1000/2.65, 20./1000/2.2530, 2809.413242, 27920.38734, 1876.123740, 28851.87339, Ex.at(i)) << " Tm" << endl;
        
        TestKinematics(25., 20., 0., 40./1000/2.65, 20./1000/2.2530, 2809.413242, 27920.38734, 1876.123740, 28851.87339, Ex.at(i));
    }
}