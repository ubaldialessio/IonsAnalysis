#include "definition.h"
#include "binning.h"
#include "utils.h"
#include "StreamUtility.h"

int Z;
double A,mass;
const double prod_min = 0.9957, prod_max = 2001;
double RigToEkn(double Rig, double Z, double A, double M) {
  return (sqrt((Rig * Z) * (Rig * Z) + M * M) - M) / A;
}
double EknToRig(double Ekn, double Z, double A, double M) {
  // return sqrt((Ekn * A + M) * (Ekn * A + M) - (M * M)) / Z;
  return sqrt(A * Ekn * (A * Ekn + 2. * M)) / Z;
}
//original
double force_field(double* x, double* p) {
  double Rig      = x[0];
  double Norm     = p[0];
  double phi      = p[1];
  double exponent = p[2];

  double Etot = RigToEkn(Rig, Z, A, mass) + mass;

  double SM_factor = (Etot * Etot - mass * mass) / ((Etot + phi) * (Etot + phi) - mass * mass);

  double flux = pow(EknToRig(Etot - mass + phi, Z, A, mass), exponent);

  return Norm * flux * SM_factor;
}

int pidToNucleusCharge(int pID);
TH1D *BuildPublishedFlux(unsigned int charge);

int main(int argc, char *argv[]) {
    TH1::SetDefaultSumw2();
    // this mask will check for charge > 2  according to both tracker or tof
	NAIA::Category cat = NAIA::Category::ChargeGT2_Trk | NAIA::Category::ChargeGT2_Tof;

    const int nRbins_Track = 54;
	const double Rbins_Track[nRbins_Track] = {0.1,0.3,0.7,1.,1.5,2.15,2.40,2.67,2.97,3.29,3.64,
	4.02,4.43,4.88,5.37,5.90, 6.47,7.09,7.76,8.48,9.26, 10.1,11.0,12.0,13.0,
	14.1,15.3,16.6,18.0,19.5,21.1,22.8,24.7,26.7,28.8,33.5,38.9,45.1,
	52.2,60.3,69.7,80.5,93.0,108.,125.,147.,175.,211,259.,330.,441.,
	660.,1200.,3000.};
	auto hist_rig_Track= new TH1D("","", nRbins_Track-1,Rbins_Track);
	auto firstTrack_Zinn	= new TH1D("",";Z_{inn};Counts",550,0,28);
	auto secondTrack_Zinn   = new TH1D("",";Z_{inn};Counts",550,0,28);
	auto firstVsSecond_Zinn = new TH2D("",";Z1_{inn};Z2_{inn}",550,0,28,550,0,28);
    auto firstTrackVsSecond                = new TH2D("", ";R1_{Inn} (GV); R2_{INN} (GV)", nRbins_Track-1, Rbins_Track, nRbins_Track-1, Rbins_Track);
    auto firstTrackVsSecond_2ndGood        = new TH2D("", ";R1_{Inn} (GV); R2_{INN} (GV)", nRbins_Track-1, Rbins_Track, nRbins_Track-1, Rbins_Track);
    auto firstTrackVsSecond_2ndGood_ZGT1_8 = new TH2D("", ";R1_{Inn} (GV); R2_{INN} (GV)", nRbins_Track-1, Rbins_Track, nRbins_Track-1, Rbins_Track);

    //All
    std::vector<TH2D*> clusterAllX_vs_R(9);
    std::vector<TH2D*> clusterAllY_vs_R(9);
    std::vector<TH2D*> clusterAllXY_vs_R(9);
    std::vector<TH2D*> clusterAllX_vs_R_secTrack(9);
    std::vector<TH2D*> clusterAllY_vs_R_secTrack(9);
    std::vector<TH2D*> clusterAllXY_vs_R_secTrack(9);
    std::vector<TH2D*> clusterAllX_vs_2R(9);
    std::vector<TH2D*> clusterAllY_vs_2R(9);
    std::vector<TH2D*> clusterAllXY_vs_2R(9);
    std::vector<TH2D*> clusterAllX_vs_2R_secTrack(9);
    std::vector<TH2D*> clusterAllY_vs_2R_secTrack(9);
    std::vector<TH2D*> clusterAllXY_vs_2R_secTrack(9);

    //Ten cm
    std::vector<TH2D*> clusterTenX_vs_R(9);
    std::vector<TH2D*> clusterTenY_vs_R(9);
    std::vector<TH2D*> clusterTenXY_vs_R(9);
    std::vector<TH2D*> clusterTenX_vs_R_secTrack(9);
    std::vector<TH2D*> clusterTenY_vs_R_secTrack(9);
    std::vector<TH2D*> clusterTenXY_vs_R_secTrack(9);
    std::vector<TH2D*> clusterTenX_vs_2R(9);
    std::vector<TH2D*> clusterTenY_vs_2R(9);
    std::vector<TH2D*> clusterTenXY_vs_2R(9);
    std::vector<TH2D*> clusterTenX_vs_2R_secTrack(9);
    std::vector<TH2D*> clusterTenY_vs_2R_secTrack(9);
    std::vector<TH2D*> clusterTenXY_vs_2R_secTrack(9);

    //Five cm
    std::vector<TH2D*> clusterFiveX_vs_R(9);
    std::vector<TH2D*> clusterFiveY_vs_R(9);
    std::vector<TH2D*> clusterFiveXY_vs_R(9);
    std::vector<TH2D*> clusterFiveX_vs_R_secTrack(9);
    std::vector<TH2D*> clusterFiveY_vs_R_secTrack(9);
    std::vector<TH2D*> clusterFiveXY_vs_R_secTrack(9);
    std::vector<TH2D*> clusterFiveX_vs_2R(9);
    std::vector<TH2D*> clusterFiveY_vs_2R(9);
    std::vector<TH2D*> clusterFiveXY_vs_2R(9);
    std::vector<TH2D*> clusterFiveX_vs_2R_secTrack(9);
    std::vector<TH2D*> clusterFiveY_vs_2R_secTrack(9);
    std::vector<TH2D*> clusterFiveXY_vs_2R_secTrack(9);

    //Two cm
    std::vector<TH2D*> clusterTwoX_vs_R(9);
    std::vector<TH2D*> clusterTwoY_vs_R(9);
    std::vector<TH2D*> clusterTwoXY_vs_R(9);
    std::vector<TH2D*> clusterTwoX_vs_R_secTrack(9);
    std::vector<TH2D*> clusterTwoY_vs_R_secTrack(9);
    std::vector<TH2D*> clusterTwoXY_vs_R_secTrack(9);
    std::vector<TH2D*> clusterTwoX_vs_2R(9);
    std::vector<TH2D*> clusterTwoY_vs_2R(9);
    std::vector<TH2D*> clusterTwoXY_vs_2R(9);
    std::vector<TH2D*> clusterTwoX_vs_2R_secTrack(9);
    std::vector<TH2D*> clusterTwoY_vs_2R_secTrack(9);
    std::vector<TH2D*> clusterTwoXY_vs_2R_secTrack(9);

    //One cm
    std::vector<TH2D*> clusterOneX_vs_R(9);
    std::vector<TH2D*> clusterOneY_vs_R(9);
    std::vector<TH2D*> clusterOneXY_vs_R(9);
    std::vector<TH2D*> clusterOneX_vs_R_secTrack(9);
    std::vector<TH2D*> clusterOneY_vs_R_secTrack(9);
    std::vector<TH2D*> clusterOneXY_vs_R_secTrack(9);
    std::vector<TH2D*> clusterOneX_vs_2R(9);
    std::vector<TH2D*> clusterOneY_vs_2R(9);
    std::vector<TH2D*> clusterOneXY_vs_2R(9);
    std::vector<TH2D*> clusterOneX_vs_2R_secTrack(9);
    std::vector<TH2D*> clusterOneY_vs_2R_secTrack(9);
    std::vector<TH2D*> clusterOneXY_vs_2R_secTrack(9);

    //One mm
    std::vector<TH2D*> clusterOneMMX_vs_R(9);
    std::vector<TH2D*> clusterOneMMY_vs_R(9);
    std::vector<TH2D*> clusterOneMMXY_vs_R(9);
    std::vector<TH2D*> clusterOneMMX_vs_R_secTrack(9);
    std::vector<TH2D*> clusterOneMMY_vs_R_secTrack(9);
    std::vector<TH2D*> clusterOneMMXY_vs_R_secTrack(9);
    std::vector<TH2D*> clusterOneMMX_vs_2R(9);
    std::vector<TH2D*> clusterOneMMY_vs_2R(9);
    std::vector<TH2D*> clusterOneMMXY_vs_2R(9);
    std::vector<TH2D*> clusterOneMMX_vs_2R_secTrack(9);
    std::vector<TH2D*> clusterOneMMY_vs_2R_secTrack(9);
    std::vector<TH2D*> clusterOneMMXY_vs_2R_secTrack(9);

    std::vector<TH2D*> edep_vs_R(20);
    std::vector<TH2D*> edep_vs_2R(20);
    std::vector<TH2D*> edep_vs_R_secTrack(20);
    std::vector<TH2D*> edep_vs_2R_secTrack(20);

    for (int i = 0; i < 20; i++) {
        edep_vs_R[i] = new TH2D(Form("l%dedep_vs_R", i+1),
                            Form(";R1 (GV); L%d edep", i+1),
                            nRbins_HighZ-1, Rbins_HighZ, 750, 0, 10000);
        edep_vs_2R[i] = new TH2D(Form("l%dedep_vs_2R", i+1),
                            Form(";R2 (GV); L%d edep", i+1),
                            nRbins_HighZ-1, Rbins_HighZ, 750, 0, 10000);    
        edep_vs_R_secTrack[i] = new TH2D(Form("l%dedep_vs_R_secTrack", i+1),
                            Form(";R1 (GV); L%d edep", i+1),
                            nRbins_HighZ-1, Rbins_HighZ, 750, 0, 10000);
        edep_vs_2R_secTrack[i] = new TH2D(Form("l%dedep_vs_2R_secTrack", i+1),
                            Form(";R2 (GV); L%d edep", i+1),
                            nRbins_HighZ-1, Rbins_HighZ, 750, 0, 10000);                          
    }

    for (int i = 0; i < 9; ++i) {

        //All
        clusterAllX_vs_R[i] = new TH2D(Form("l%dClusterAllX_vs_R", i+1),
                                    Form(";R1 (GV); L%d clusters", i+1),
                                    nRbins_HighZ-1, Rbins_HighZ, 150, 0, 250);

        clusterAllY_vs_R[i] = new TH2D(Form("l%dClusterAllY_vs_R", i+1),
                                    Form(";R1 (GV); L%d clusters", i+1),
                                    nRbins_HighZ-1, Rbins_HighZ, 150, 0, 250);

        clusterAllXY_vs_R[i] = new TH2D(Form("l%dClusterAllXY_vs_R", i+1),
                                        Form(";R1 (GV); L%d clusters", i+1),
                                        nRbins_HighZ-1, Rbins_HighZ, 150, 0, 250);
        clusterAllX_vs_R_secTrack[i] = new TH2D(Form("l%dClusterAllX_vs_R_secTrack", i+1),
                                    Form(";R1 (GV); L%d clusters", i+1),
                                    nRbins_HighZ-1, Rbins_HighZ, 150, 0, 250);

        clusterAllY_vs_R_secTrack[i] = new TH2D(Form("l%dClusterAllY_vs_R_secTrack", i+1),
                                    Form(";R1 (GV); L%d clusters", i+1),
                                    nRbins_HighZ-1, Rbins_HighZ, 150, 0, 250);

        clusterAllXY_vs_R_secTrack[i] = new TH2D(Form("l%dClusterAllXY_vs_R_secTrack", i+1),
                                        Form(";R1 (GV); L%d clusters", i+1),
                                        nRbins_HighZ-1, Rbins_HighZ, 150, 0, 250);

        clusterAllX_vs_2R[i] = new TH2D(Form("l%dClusterAllX_vs_2R", i+1),
                                    Form(";R2 (GV); L%d clusters", i+1),
                                    nRbins_HighZ-1, Rbins_HighZ, 150, 0, 250);

        clusterAllY_vs_2R[i] = new TH2D(Form("l%dClusterAllY_vs_2R", i+1),
                                    Form(";R2 (GV); L%d clusters", i+1),
                                    nRbins_HighZ-1, Rbins_HighZ, 150, 0, 250);

        clusterAllXY_vs_2R[i] = new TH2D(Form("l%dClusterAllXY_vs_2R", i+1),
                                        Form(";R2 (GV); L%d clusters", i+1),
                                        nRbins_HighZ-1, Rbins_HighZ, 150, 0, 250);
        clusterAllX_vs_2R_secTrack[i] = new TH2D(Form("l%dClusterAllX_vs_2R_secTrack", i+1),
                                    Form(";R2 (GV); L%d clusters", i+1),
                                    nRbins_HighZ-1, Rbins_HighZ, 150, 0, 250);

        clusterAllY_vs_2R_secTrack[i] = new TH2D(Form("l%dClusterAllY_vs_2R_secTrack", i+1),
                                    Form(";R2 (GV); L%d clusters", i+1),
                                    nRbins_HighZ-1, Rbins_HighZ, 150, 0, 250);

        clusterAllXY_vs_2R_secTrack[i] = new TH2D(Form("l%dClusterAllXY_vs_2R_secTrack", i+1),
                                        Form(";R2 (GV); L%d clusters", i+1),
                                        nRbins_HighZ-1, Rbins_HighZ, 150, 0, 250);

        //10 cm
        clusterTenX_vs_R[i] = new TH2D(Form("l%dClusterTenX_vs_R", i+1),
                                    Form(";R1 (GV); L%d clusters", i+1),
                                    nRbins_HighZ-1, Rbins_HighZ, 150, 0, 250);

        clusterTenY_vs_R[i] = new TH2D(Form("l%dClusterTenY_vs_R", i+1),
                                    Form(";R1 (GV); L%d clusters", i+1),
                                    nRbins_HighZ-1, Rbins_HighZ, 150, 0, 250);

        clusterTenXY_vs_R[i] = new TH2D(Form("l%dClusterTenXY_vs_R", i+1),
                                        Form(";R1 (GV); L%d clusters", i+1),
                                        nRbins_HighZ-1, Rbins_HighZ, 150, 0, 250);
        clusterTenX_vs_R_secTrack[i] = new TH2D(Form("l%dClusterTenX_vs_R_secTrack", i+1),
                                    Form(";R1 (GV); L%d clusters", i+1),
                                    nRbins_HighZ-1, Rbins_HighZ, 150, 0, 250);

        clusterTenY_vs_R_secTrack[i] = new TH2D(Form("l%dClusterTenY_vs_R_secTrack", i+1),
                                    Form(";R1 (GV); L%d clusters", i+1),
                                    nRbins_HighZ-1, Rbins_HighZ, 150, 0, 250);

        clusterTenXY_vs_R_secTrack[i] = new TH2D(Form("l%dClusterTenXY_vs_R_secTrack", i+1),
                                        Form(";R1 (GV); L%d clusters", i+1),
                                        nRbins_HighZ-1, Rbins_HighZ, 150, 0, 250);

        clusterTenX_vs_2R[i] = new TH2D(Form("l%dClusterTenX_vs_2R", i+1),
                                    Form(";R2 (GV); L%d clusters", i+1),
                                    nRbins_HighZ-1, Rbins_HighZ, 150, 0, 250);

        clusterTenY_vs_2R[i] = new TH2D(Form("l%dClusterTenY_vs_2R", i+1),
                                    Form(";R2 (GV); L%d clusters", i+1),
                                    nRbins_HighZ-1, Rbins_HighZ, 150, 0, 250);

        clusterTenXY_vs_2R[i] = new TH2D(Form("l%dClusterTenXY_vs_2R", i+1),
                                        Form(";R2 (GV); L%d clusters", i+1),
                                        nRbins_HighZ-1, Rbins_HighZ, 150, 0, 250);
        clusterTenX_vs_2R_secTrack[i] = new TH2D(Form("l%dClusterTenX_vs_2R_secTrack", i+1),
                                    Form(";R2 (GV); L%d clusters", i+1),
                                    nRbins_HighZ-1, Rbins_HighZ, 150, 0, 250);

        clusterTenY_vs_2R_secTrack[i] = new TH2D(Form("l%dClusterTenY_vs_2R_secTrack", i+1),
                                    Form(";R2 (GV); L%d clusters", i+1),
                                    nRbins_HighZ-1, Rbins_HighZ, 150, 0, 250);

        clusterTenXY_vs_2R_secTrack[i] = new TH2D(Form("l%dClusterTenXY_vs_2R_secTrack", i+1),
                                        Form(";R2 (GV); L%d clusters", i+1),
                                        nRbins_HighZ-1, Rbins_HighZ, 150, 0, 250);    
    
        //5 cm
        clusterFiveX_vs_R[i] = new TH2D(Form("l%dClusterFiveX_vs_R", i+1),
                                    Form(";R1 (GV); L%d clusters", i+1),
                                    nRbins_HighZ-1, Rbins_HighZ, 150, 0, 250);

        clusterFiveY_vs_R[i] = new TH2D(Form("l%dClusterFiveY_vs_R", i+1),
                                    Form(";R1 (GV); L%d clusters", i+1),
                                    nRbins_HighZ-1, Rbins_HighZ, 150, 0, 250);

        clusterFiveXY_vs_R[i] = new TH2D(Form("l%dClusterFiveXY_vs_R", i+1),
                                        Form(";R1 (GV); L%d clusters", i+1),
                                        nRbins_HighZ-1, Rbins_HighZ, 150, 0, 250);
        clusterFiveX_vs_R_secTrack[i] = new TH2D(Form("l%dClusterFiveX_vs_R_secTrack", i+1),
                                    Form(";R1 (GV); L%d clusters", i+1),
                                    nRbins_HighZ-1, Rbins_HighZ, 150, 0, 250);

        clusterFiveY_vs_R_secTrack[i] = new TH2D(Form("l%dClusterFiveY_vs_R_secTrack", i+1),
                                    Form(";R1 (GV); L%d clusters", i+1),
                                    nRbins_HighZ-1, Rbins_HighZ, 150, 0, 250);

        clusterFiveXY_vs_R_secTrack[i] = new TH2D(Form("l%dClusterFiveXY_vs_R_secTrack", i+1),
                                        Form(";R1 (GV); L%d clusters", i+1),
                                        nRbins_HighZ-1, Rbins_HighZ, 150, 0, 250);

        clusterFiveX_vs_2R[i] = new TH2D(Form("l%dClusterFiveX_vs_2R", i+1),
                                    Form(";R2 (GV); L%d clusters", i+1),
                                    nRbins_HighZ-1, Rbins_HighZ, 150, 0, 250);

        clusterFiveY_vs_2R[i] = new TH2D(Form("l%dClusterFiveY_vs_2R", i+1),
                                    Form(";R2 (GV); L%d clusters", i+1),
                                    nRbins_HighZ-1, Rbins_HighZ, 150, 0, 250);

        clusterFiveXY_vs_2R[i] = new TH2D(Form("l%dClusterFiveXY_vs_2R", i+1),
                                        Form(";R2 (GV); L%d clusters", i+1),
                                        nRbins_HighZ-1, Rbins_HighZ, 150, 0, 250);
        clusterFiveX_vs_2R_secTrack[i] = new TH2D(Form("l%dClusterFiveX_vs_2R_secTrack", i+1),
                                    Form(";R2 (GV); L%d clusters", i+1),
                                    nRbins_HighZ-1, Rbins_HighZ, 150, 0, 250);

        clusterFiveY_vs_2R_secTrack[i] = new TH2D(Form("l%dClusterFiveY_vs_2R_secTrack", i+1),
                                    Form(";R2 (GV); L%d clusters", i+1),
                                    nRbins_HighZ-1, Rbins_HighZ, 150, 0, 250);

        clusterFiveXY_vs_2R_secTrack[i] = new TH2D(Form("l%dClusterFiveXY_vs_2R_secTrack", i+1),
                                        Form(";R2 (GV); L%d clusters", i+1),
                                        nRbins_HighZ-1, Rbins_HighZ, 150, 0, 250);    
        //2 cm
        clusterTwoX_vs_R[i] = new TH2D(Form("l%dClusterTwoX_vs_R", i+1),
                                    Form(";R1 (GV); L%d clusters", i+1),
                                    nRbins_HighZ-1, Rbins_HighZ, 150, 0, 250);

        clusterTwoY_vs_R[i] = new TH2D(Form("l%dClusterTwoY_vs_R", i+1),
                                    Form(";R1 (GV); L%d clusters", i+1),
                                    nRbins_HighZ-1, Rbins_HighZ, 150, 0, 250);

        clusterTwoXY_vs_R[i] = new TH2D(Form("l%dClusterTwoXY_vs_R", i+1),
                                        Form(";R1 (GV); L%d clusters", i+1),
                                        nRbins_HighZ-1, Rbins_HighZ, 150, 0, 250);
        clusterTwoX_vs_R_secTrack[i] = new TH2D(Form("l%dClusterTwoX_vs_R_secTrack", i+1),
                                    Form(";R1 (GV); L%d clusters", i+1),
                                    nRbins_HighZ-1, Rbins_HighZ, 150, 0, 250);

        clusterTwoY_vs_R_secTrack[i] = new TH2D(Form("l%dClusterTwoY_vs_R_secTrack", i+1),
                                    Form(";R1 (GV); L%d clusters", i+1),
                                    nRbins_HighZ-1, Rbins_HighZ, 150, 0, 250);

        clusterTwoXY_vs_R_secTrack[i] = new TH2D(Form("l%dClusterTwoXY_vs_R_secTrack", i+1),
                                        Form(";R1 (GV); L%d clusters", i+1),
                                        nRbins_HighZ-1, Rbins_HighZ, 150, 0, 250);

        clusterTwoX_vs_2R[i] = new TH2D(Form("l%dClusterTwoX_vs_2R", i+1),
                                    Form(";R2 (GV); L%d clusters", i+1),
                                    nRbins_HighZ-1, Rbins_HighZ, 150, 0, 250);

        clusterTwoY_vs_2R[i] = new TH2D(Form("l%dClusterTwoY_vs_2R", i+1),
                                    Form(";R2 (GV); L%d clusters", i+1),
                                    nRbins_HighZ-1, Rbins_HighZ, 150, 0, 250);

        clusterTwoXY_vs_2R[i] = new TH2D(Form("l%dClusterTwoXY_vs_2R", i+1),
                                        Form(";R2 (GV); L%d clusters", i+1),
                                        nRbins_HighZ-1, Rbins_HighZ, 150, 0, 250);
        clusterTwoX_vs_2R_secTrack[i] = new TH2D(Form("l%dClusterTwoX_vs_2R_secTrack", i+1),
                                    Form(";R2 (GV); L%d clusters", i+1),
                                    nRbins_HighZ-1, Rbins_HighZ, 150, 0, 250);

        clusterTwoY_vs_2R_secTrack[i] = new TH2D(Form("l%dClusterTwoY_vs_2R_secTrack", i+1),
                                    Form(";R2 (GV); L%d clusters", i+1),
                                    nRbins_HighZ-1, Rbins_HighZ, 150, 0, 250);

        clusterTwoXY_vs_2R_secTrack[i] = new TH2D(Form("l%dClusterTwoXY_vs_2R_secTrack", i+1),
                                        Form(";R2 (GV); L%d clusters", i+1),
                                        nRbins_HighZ-1, Rbins_HighZ, 150, 0, 250);   
        //1 cm
        clusterOneX_vs_R[i] = new TH2D(Form("l%dClusterOneX_vs_R", i+1),
                                    Form(";R1 (GV); L%d clusters", i+1),
                                    nRbins_HighZ-1, Rbins_HighZ, 150, 0, 250);

        clusterOneY_vs_R[i] = new TH2D(Form("l%dClusterOneY_vs_R", i+1),
                                    Form(";R1 (GV); L%d clusters", i+1),
                                    nRbins_HighZ-1, Rbins_HighZ, 150, 0, 250);

        clusterOneXY_vs_R[i] = new TH2D(Form("l%dClusterOneXY_vs_R", i+1),
                                        Form(";R1 (GV); L%d clusters", i+1),
                                        nRbins_HighZ-1, Rbins_HighZ, 150, 0, 250);
        clusterOneX_vs_R_secTrack[i] = new TH2D(Form("l%dClusterOneX_vs_R_secTrack", i+1),
                                    Form(";R1 (GV); L%d clusters", i+1),
                                    nRbins_HighZ-1, Rbins_HighZ, 150, 0, 250);

        clusterOneY_vs_R_secTrack[i] = new TH2D(Form("l%dClusterOneY_vs_R_secTrack", i+1),
                                    Form(";R1 (GV); L%d clusters", i+1),
                                    nRbins_HighZ-1, Rbins_HighZ, 150, 0, 250);

        clusterOneXY_vs_R_secTrack[i] = new TH2D(Form("l%dClusterOneXY_vs_R_secTrack", i+1),
                                        Form(";R1 (GV); L%d clusters", i+1),
                                        nRbins_HighZ-1, Rbins_HighZ, 150, 0, 250);

        clusterOneX_vs_2R[i] = new TH2D(Form("l%dClusterOneX_vs_2R", i+1),
                                    Form(";R2 (GV); L%d clusters", i+1),
                                    nRbins_HighZ-1, Rbins_HighZ, 150, 0, 250);

        clusterOneY_vs_2R[i] = new TH2D(Form("l%dClusterOneY_vs_2R", i+1),
                                    Form(";R2 (GV); L%d clusters", i+1),
                                    nRbins_HighZ-1, Rbins_HighZ, 150, 0, 250);

        clusterOneXY_vs_2R[i] = new TH2D(Form("l%dClusterOneXY_vs_2R", i+1),
                                        Form(";R2 (GV); L%d clusters", i+1),
                                        nRbins_HighZ-1, Rbins_HighZ, 150, 0, 250);
        clusterOneX_vs_2R_secTrack[i] = new TH2D(Form("l%dClusterOneX_vs_2R_secTrack", i+1),
                                    Form(";R2 (GV); L%d clusters", i+1),
                                    nRbins_HighZ-1, Rbins_HighZ, 150, 0, 250);

        clusterOneY_vs_2R_secTrack[i] = new TH2D(Form("l%dClusterOneY_vs_2R_secTrack", i+1),
                                    Form(";R2 (GV); L%d clusters", i+1),
                                    nRbins_HighZ-1, Rbins_HighZ, 150, 0, 250);

        clusterOneXY_vs_2R_secTrack[i] = new TH2D(Form("l%dClusterOneXY_vs_2R_secTrack", i+1),
                                        Form(";R2 (GV); L%d clusters", i+1),
                                        nRbins_HighZ-1, Rbins_HighZ, 150, 0, 250);
        //1 mm
        clusterOneMMX_vs_R[i] = new TH2D(Form("l%dClusterOneMMX_vs_R", i+1),
                                    Form(";R1 (GV); L%d clusters", i+1),
                                    nRbins_HighZ-1, Rbins_HighZ, 150, 0, 250);

        clusterOneMMY_vs_R[i] = new TH2D(Form("l%dClusterOneMMY_vs_R", i+1),
                                    Form(";R1 (GV); L%d clusters", i+1),
                                    nRbins_HighZ-1, Rbins_HighZ, 150, 0, 250);

        clusterOneMMXY_vs_R[i] = new TH2D(Form("l%dClusterOneMMXY_vs_R", i+1),
                                        Form(";R1 (GV); L%d clusters", i+1),
                                        nRbins_HighZ-1, Rbins_HighZ, 150, 0, 250);
        clusterOneMMX_vs_R_secTrack[i] = new TH2D(Form("l%dClusterOneMMX_vs_R_secTrack", i+1),
                                    Form(";R1 (GV); L%d clusters", i+1),
                                    nRbins_HighZ-1, Rbins_HighZ, 150, 0, 250);

        clusterOneMMY_vs_R_secTrack[i] = new TH2D(Form("l%dClusterOneMMY_vs_R_secTrack", i+1),
                                    Form(";R1 (GV); L%d clusters", i+1),
                                    nRbins_HighZ-1, Rbins_HighZ, 150, 0, 250);

        clusterOneMMXY_vs_R_secTrack[i] = new TH2D(Form("l%dClusterOneMMXY_vs_R_secTrack", i+1),
                                        Form(";R1 (GV); L%d clusters", i+1),
                                        nRbins_HighZ-1, Rbins_HighZ, 150, 0, 250);

        clusterOneMMX_vs_2R[i] = new TH2D(Form("l%dClusterOneMMX_vs_2R", i+1),
                                    Form(";R2 (GV); L%d clusters", i+1),
                                    nRbins_HighZ-1, Rbins_HighZ, 150, 0, 250);

        clusterOneMMY_vs_2R[i] = new TH2D(Form("l%dClusterOneMMY_vs_2R", i+1),
                                    Form(";R2 (GV); L%d clusters", i+1),
                                    nRbins_HighZ-1, Rbins_HighZ, 150, 0, 250);

        clusterOneMMXY_vs_2R[i] = new TH2D(Form("l%dClusterOneMMXY_vs_2R", i+1),
                                        Form(";R2 (GV); L%d clusters", i+1),
                                        nRbins_HighZ-1, Rbins_HighZ, 150, 0, 250);
        clusterOneMMX_vs_2R_secTrack[i] = new TH2D(Form("l%dClusterOneMMX_vs_2R_secTrack", i+1),
                                    Form(";R2 (GV); L%d clusters", i+1),
                                    nRbins_HighZ-1, Rbins_HighZ, 150, 0, 250);

        clusterOneMMY_vs_2R_secTrack[i] = new TH2D(Form("l%dClusterOneMMY_vs_2R_secTrack", i+1),
                                    Form(";R2 (GV); L%d clusters", i+1),
                                    nRbins_HighZ-1, Rbins_HighZ, 150, 0, 250);

        clusterOneMMXY_vs_2R_secTrack[i] = new TH2D(Form("l%dClusterOneMMXY_vs_2R_secTrack", i+1),
                                        Form(";R2 (GV); L%d clusters", i+1),
                                        nRbins_HighZ-1, Rbins_HighZ, 150, 0, 250);
    }

    auto utofCluster_vs_R         = new TH2D("", ";R1 (GV); utof clusters", nRbins_HighZ-1,Rbins_HighZ, 100, 0 , 100);
    auto ltofCluster_vs_R         = new TH2D("", ";R1 (GV); ltof clusters", nRbins_HighZ-1,Rbins_HighZ, 100, 0 , 100);
    auto tofClusterNoTrack_vs_R   = new TH2D("", ";R1 (GV); utof clusters", nRbins_HighZ-1,Rbins_HighZ, 100, 0 , 100);

    auto utofCluster_vs_R_secTrack         = new TH2D("", ";R1 (GV); utof clusters", nRbins_HighZ-1,Rbins_HighZ, 100, 0 , 100);
    auto ltofCluster_vs_R_secTrack         = new TH2D("", ";R1 (GV); ltof clusters", nRbins_HighZ-1,Rbins_HighZ, 100, 0 , 100);
    auto tofClusterNoTrack_vs_R_secTrack   = new TH2D("", ";R1 (GV); utof clusters", nRbins_HighZ-1,Rbins_HighZ, 100, 0 , 100);

    auto utofCluster_vs_2R         = new TH2D("", ";R2 (GV); utof clusters", nRbins_HighZ-1,Rbins_HighZ, 100, 0 , 100);
    auto ltofCluster_vs_2R         = new TH2D("", ";R2 (GV); ltof clusters", nRbins_HighZ-1,Rbins_HighZ, 100, 0 , 100);
    auto tofClusterNoTrack_vs_2R   = new TH2D("", ";R2 (GV); utof clusters", nRbins_HighZ-1,Rbins_HighZ, 100, 0 , 100);

    auto utofCluster_vs_2R_secTrack         = new TH2D("", ";R2 (GV); utof clusters", nRbins_HighZ-1,Rbins_HighZ, 100, 0 , 100);
    auto ltofCluster_vs_2R_secTrack         = new TH2D("", ";R2 (GV); ltof clusters", nRbins_HighZ-1,Rbins_HighZ, 100, 0 , 100);
    auto tofClusterNoTrack_vs_2R_secTrack   = new TH2D("", ";R2 (GV); utof clusters", nRbins_HighZ-1,Rbins_HighZ, 100, 0 , 100);

    auto TRDHitTot_vs_R = new TH2D("", ";R1 (GV); TRD hit tot", nRbins_HighZ-1,Rbins_HighZ, 250, 0 , 250);
    auto TRDHitOn_vs_R = new TH2D("", ";R1 (GV); TRD hit on",   nRbins_HighZ-1,Rbins_HighZ, 250, 0 , 250);
    auto TRDHitOff_vs_R = new TH2D("", ";R1 (GV); TRD hit off",  nRbins_HighZ-1,Rbins_HighZ, 250, 0 , 250);
    auto TRDAmpTot_vs_R= new TH2D("", ";R1 (GV); TRD amp tot",  nRbins_HighZ-1,Rbins_HighZ, 7500, 0 , 100000);
    auto TRDAmpOn_vs_R= new TH2D("", ";R1 (GV); TRD amp on",    nRbins_HighZ-1,Rbins_HighZ, 7500, 0 , 100000);
    auto TRDAmpOff_vs_R= new TH2D("", ";R1 (GV); TRD amp off",   nRbins_HighZ-1,Rbins_HighZ, 7500, 0 , 100000);
    auto TRDChargeTot_vs_R = new TH2D("", ";R1 (GV); TRD charge tot",   nRbins_HighZ-1,Rbins_HighZ, 250, 0 , 50);
    auto TRDChargeUp_vs_R = new TH2D("", ";R1 (GV); TRD charge up",   nRbins_HighZ-1,Rbins_HighZ, 250, 0 , 50);
    auto TRDChargeLow_vs_R = new TH2D("", ";R1 (GV); TRD charge low",   nRbins_HighZ-1,Rbins_HighZ, 250, 0 , 50);

    auto TRDHitTot_vs_2R = new TH2D("", ";R2 (GV); TRD hit tot", nRbins_HighZ-1,Rbins_HighZ, 250, 0 , 250);
    auto TRDHitOn_vs_2R = new TH2D("", ";R2 (GV); TRD hit on",   nRbins_HighZ-1,Rbins_HighZ, 250, 0 , 250);
    auto TRDHitOff_vs_2R = new TH2D("", ";R2 (GV); TRD hit off",  nRbins_HighZ-1,Rbins_HighZ, 250, 0 , 250);
    auto TRDAmpTot_vs_2R= new TH2D("", ";R2 (GV); TRD amp tot",  nRbins_HighZ-1,Rbins_HighZ, 7500, 0 , 100000);
    auto TRDAmpOn_vs_2R= new TH2D("", ";R2 (GV); TRD amp on",    nRbins_HighZ-1,Rbins_HighZ, 7500, 0 , 100000);
    auto TRDAmpOff_vs_2R= new TH2D("", ";R2 (GV); TRD amp off",   nRbins_HighZ-1,Rbins_HighZ, 7500, 0 , 100000);
    auto TRDChargeTot_vs_2R = new TH2D("", ";R2 (GV); TRD charge tot",   nRbins_HighZ-1,Rbins_HighZ, 250, 0 , 50);
    auto TRDChargeUp_vs_2R = new TH2D("", ";R2 (GV); TRD charge up",   nRbins_HighZ-1,Rbins_HighZ, 250, 0 , 50);
    auto TRDChargeLow_vs_2R = new TH2D("", ";R2 (GV); TRD charge low",   nRbins_HighZ-1,Rbins_HighZ, 250, 0 , 50);

    auto TRDHitTot_vs_R_secTrack = new TH2D("", ";R1 (GV); TRD hit tot", nRbins_HighZ-1,Rbins_HighZ, 250, 0 , 250);
    auto TRDHitOn_vs_R_secTrack = new TH2D("", ";R1 (GV); TRD hit on",   nRbins_HighZ-1,Rbins_HighZ, 250, 0 , 250);
    auto TRDHitOff_vs_R_secTrack = new TH2D("", ";R1 (GV); TRD hit off",  nRbins_HighZ-1,Rbins_HighZ, 250, 0 , 250);
    auto TRDAmpTot_vs_R_secTrack= new TH2D("", ";R1 (GV); TRD amp tot",  nRbins_HighZ-1,Rbins_HighZ, 7500, 0 , 100000);
    auto TRDAmpOn_vs_R_secTrack= new TH2D("", ";R1 (GV); TRD amp on",    nRbins_HighZ-1,Rbins_HighZ, 7500, 0 , 100000);
    auto TRDAmpOff_vs_R_secTrack= new TH2D("", ";R1 (GV); TRD amp off",   nRbins_HighZ-1,Rbins_HighZ, 7500, 0 , 100000);
    auto TRDChargeTot_vs_R_secTrack = new TH2D("", ";R1 (GV); TRD charge tot",   nRbins_HighZ-1,Rbins_HighZ, 250, 0 , 50);
    auto TRDChargeUp_vs_R_secTrack = new TH2D("", ";R1 (GV); TRD charge up",   nRbins_HighZ-1,Rbins_HighZ, 250, 0 , 50);
    auto TRDChargeLow_vs_R_secTrack = new TH2D("", ";R1 (GV); TRD charge low",   nRbins_HighZ-1,Rbins_HighZ, 250, 0 , 50);

    auto TRDHitTot_vs_2R_secTrack = new TH2D("", ";R2 (GV); TRD hit tot", nRbins_HighZ-1,Rbins_HighZ, 250, 0 , 250);
    auto TRDHitOn_vs_2R_secTrack = new TH2D("", ";R2 (GV); TRD hit on",   nRbins_HighZ-1,Rbins_HighZ, 250, 0 , 250);
    auto TRDHitOff_vs_2R_secTrack = new TH2D("", ";R2 (GV); TRD hit off",  nRbins_HighZ-1,Rbins_HighZ, 250, 0 , 250);
    auto TRDAmpTot_vs_2R_secTrack= new TH2D("", ";R2 (GV); TRD amp tot",  nRbins_HighZ-1,Rbins_HighZ, 7500, 0 , 100000);
    auto TRDAmpOn_vs_2R_secTrack= new TH2D("", ";R2 (GV); TRD amp on",    nRbins_HighZ-1,Rbins_HighZ, 7500, 0 , 100000);
    auto TRDAmpOff_vs_2R_secTrack= new TH2D("", ";R2 (GV); TRD amp off",   nRbins_HighZ-1,Rbins_HighZ, 7500, 0 , 100000);
    auto TRDChargeTot_vs_2R_secTrack = new TH2D("", ";R2 (GV); TRD charge tot",   nRbins_HighZ-1,Rbins_HighZ, 250, 0 , 50);
    auto TRDChargeUp_vs_2R_secTrack = new TH2D("", ";R2 (GV); TRD charge up",   nRbins_HighZ-1,Rbins_HighZ, 250, 0 , 50);
    auto TRDChargeLow_vs_2R_secTrack = new TH2D("", ";R2 (GV); TRD charge low",   nRbins_HighZ-1,Rbins_HighZ, 250, 0 , 50);

    auto TRDChargeTot_vs_R_offCut = new TH2D("", ";R1 (GV); TRD charge tot_offCut",   nRbins_HighZ-1,Rbins_HighZ, 250, 0 , 50);
    auto TRDChargeUp_vs_R_offCut = new TH2D("", ";R1 (GV); TRD charge up_offCut",   nRbins_HighZ-1,Rbins_HighZ, 250, 0 , 50);
    auto TRDChargeLow_vs_R_offCut = new TH2D("", ";R1 (GV); TRD charge low_offCut",   nRbins_HighZ-1,Rbins_HighZ, 250, 0 , 50);

    if (argc < 3) {
		printf("Usage: \n");
	    printf("%s <charge> <path/to/input.root> <output.root> \n", argv[0]);
	    return 1;
	}
	unsigned int charge=atoi(argv[1]);
	TString exename = argv[0],
		    infilename=argv[2],
		    outname=argv[3],
		    out, ionPath=getIonPath(charge);

    NAIA::NAIAChain chain;
    if(infilename.Contains(".root") /* && filesystem::exists(infilename.Data())*/ ){
		chain.Add(infilename.Data());
	}else if (infilename.Contains(".txt") /* && filesystem::exists(infilename.Data())*/ ){
		ifstream infilelist(infilename.Data());
		TString bufname;
		while(infilelist >> bufname) 
		    if(bufname.Contains(".root" /* && filesystem::exists(bufname.Data()) */))
		        chain.Add(bufname.Data());
		}
	chain.SetupBranches();
	bool isMC = chain.IsMC();
	out = "/storage/gpfs_ams/ams/users/aubaldi/IonsAnalysis/ControlVariables/ClustersNumber/"+getIonPath(charge)+"/"+outname;
	unsigned int utime;
    unsigned int NutofCluster, NltofCluster;
    short NtofClusterMatchTrack;
	double cutoff,reco_il1,sec_track_reco_il1 ;
    std::vector<int> NClusterX(9, 0), NClusterTenX(9,0), NClusterFiveX(9,0), NClusterTwoX(9,0), NClusterOneX(9,0), NClusterOneMMX(9,0);
    std::vector<int> NClusterY(9, 0), NClusterTenY(9,0), NClusterFiveY(9,0), NClusterTwoY(9,0), NClusterOneY(9,0), NClusterOneMMY(9,0);
    std::vector<float> Edep(20,0);
    int trdHit,trdHitOn,trdHitOff;
    float ampTot, ampOn, ampOff, z1_inn, z2_inn;
    float chargeTot, chargeUp, chargeLow;
    double gen;

    auto myselStoP = MyselStoP(charge);
    auto hasGoodSecondTrack = GoodSecTrack();
    auto secTrackOnDiagonal = MySel::IsSecondTrackOnDiagonal(FIT,IL1,0.3);

    for(Event& event : chain) {
		utime = event.header->UTCTime;
		if (!isMC) {
			if (utime>=1620025528 && utime<1635856717) // skipping photon-polarization period
				continue;
            if (isRunBad(utime)==true) //skipping trigger study and wrong configuration
				continue;
			if (!ns::DefaultRTISelection(chain.GetEventRTIInfo() ) || chain.GetEventRTIInfo().IsInSAA() )  //RTI check
				continue;
			cutoff 		= chain.GetEventRTIInfo().MaxIGRFCutoff[2][1];
            ////MY SELECTION////
            if (myselStoP(event)) {
                reco_il1  = event.trTrackBase->Rigidity[FIT][IL1];
                int rbin = hist_rig_highZ->FindBin(reco_il1);
                double rlowedge = hist_rig_highZ->GetBinLowEdge(rbin);
                if(rlowedge > 1.2*cutoff) {
                    sec_track_reco_il1 = event.secondTrTrackBase->Rigidity[FIT][IL1];

                    firstTrackVsSecond->Fill(reco_il1,sec_track_reco_il1);

                    trdHitOn = event.trdKBase->NHits[TRD_ONTRACK];
                    trdHitOff= event.trdKBase->NHits[TRD_OFFTRACK];;
                    trdHit = trdHitOn+trdHitOff;
                    ampOn  = event.trdKBase->Amps[TRD_ONTRACK];
                    ampOff = event.trdKBase->Amps[TRD_OFFTRACK];
                    ampTot = ampOn+ampOff;
                    chargeUp = event.trdKBase->Charge[TRD_UP];
                    chargeLow= event.trdKBase->Charge[TRD_LOW];
                    chargeTot = event.trdKBase->Charge[TRD_TOT];
                    for (int i = 0; i < 20; i++) {
                        Edep[i] = event.trdKBase->Edep[i];
                        edep_vs_R[i]->Fill(reco_il1,Edep[i]);
                        edep_vs_2R[i]->Fill(sec_track_reco_il1,Edep[i]);
                    }
                    for (int i = 0; i < 9; ++i) {
                        NClusterX[i] = event.trTrackPlus->NClusters[i][ALL_LAYER][XSD];
                        NClusterY[i] = event.trTrackPlus->NClusters[i][ALL_LAYER][YSD];
                        NClusterTenX[i] = event.trTrackPlus->NClusters[i][TEN_CM][XSD];
                        NClusterTenY[i] = event.trTrackPlus->NClusters[i][TEN_CM][YSD];
                        NClusterFiveX[i] = event.trTrackPlus->NClusters[i][FIVE_CM][XSD];
                        NClusterFiveY[i] = event.trTrackPlus->NClusters[i][FIVE_CM][YSD];
                        NClusterTwoX[i] = event.trTrackPlus->NClusters[i][TWO_CM][XSD];
                        NClusterTwoY[i] = event.trTrackPlus->NClusters[i][TWO_CM][YSD];
                        NClusterOneX[i] = event.trTrackPlus->NClusters[i][ONE_CM][XSD];
                        NClusterOneY[i] = event.trTrackPlus->NClusters[i][ONE_CM][YSD];
                        NClusterOneMMX[i] = event.trTrackPlus->NClusters[i][ONE_MM][XSD];
                        NClusterOneMMY[i] = event.trTrackPlus->NClusters[i][ONE_MM][YSD];
                        //All vs R1
                        clusterAllX_vs_R[i]->Fill(reco_il1, NClusterX[i]);
                        clusterAllY_vs_R[i]->Fill(reco_il1, NClusterY[i]);
                        clusterAllXY_vs_R[i]->Fill(reco_il1, NClusterX[i] + NClusterY[i]);
                        clusterAllX_vs_2R[i]->Fill(sec_track_reco_il1, NClusterX[i]);
                        clusterAllY_vs_2R[i]->Fill(sec_track_reco_il1, NClusterY[i]);
                        clusterAllXY_vs_2R[i]->Fill(sec_track_reco_il1, NClusterX[i] + NClusterY[i]);

                        //10 cm vs R1
                        clusterTenX_vs_R[i]->Fill(reco_il1, NClusterTenX[i]);
                        clusterTenY_vs_R[i]->Fill(reco_il1, NClusterTenY[i]);
                        clusterTenXY_vs_R[i]->Fill(reco_il1, NClusterTenX[i] + NClusterTenY[i]);
                        clusterTenX_vs_2R[i]->Fill(sec_track_reco_il1, NClusterTenX[i]);
                        clusterTenY_vs_2R[i]->Fill(sec_track_reco_il1, NClusterTenY[i]);
                        clusterTenXY_vs_2R[i]->Fill(sec_track_reco_il1, NClusterTenX[i] + NClusterTenY[i]);

                        //5 cm vs R1
                        clusterFiveX_vs_R[i]->Fill(reco_il1, NClusterFiveX[i]);
                        clusterFiveY_vs_R[i]->Fill(reco_il1, NClusterFiveY[i]);
                        clusterFiveXY_vs_R[i]->Fill(reco_il1, NClusterFiveX[i] + NClusterFiveY[i]);
                        clusterFiveX_vs_2R[i]->Fill(sec_track_reco_il1, NClusterFiveX[i]);
                        clusterFiveY_vs_2R[i]->Fill(sec_track_reco_il1, NClusterFiveY[i]);
                        clusterFiveXY_vs_2R[i]->Fill(sec_track_reco_il1, NClusterFiveX[i] + NClusterFiveY[i]);

                        //2 cm vs R1
                        clusterTwoX_vs_R[i]->Fill(reco_il1, NClusterTwoX[i]);
                        clusterTwoY_vs_R[i]->Fill(reco_il1, NClusterTwoY[i]);
                        clusterTwoXY_vs_R[i]->Fill(reco_il1, NClusterTwoX[i] + NClusterTwoY[i]);
                        clusterTwoX_vs_2R[i]->Fill(sec_track_reco_il1, NClusterTwoX[i]);
                        clusterTwoY_vs_2R[i]->Fill(sec_track_reco_il1, NClusterTwoY[i]);
                        clusterTwoXY_vs_2R[i]->Fill(sec_track_reco_il1, NClusterTwoX[i] + NClusterTwoY[i]);

                        //1 cm vs R1
                        clusterOneX_vs_R[i]->Fill(reco_il1, NClusterOneX[i]);
                        clusterOneY_vs_R[i]->Fill(reco_il1, NClusterOneY[i]);
                        clusterOneXY_vs_R[i]->Fill(reco_il1, NClusterOneX[i] + NClusterOneY[i]);
                        clusterOneX_vs_2R[i]->Fill(sec_track_reco_il1, NClusterOneX[i]);
                        clusterOneY_vs_2R[i]->Fill(sec_track_reco_il1, NClusterOneY[i]);
                        clusterOneXY_vs_2R[i]->Fill(sec_track_reco_il1, NClusterOneX[i] + NClusterOneY[i]);

                        //1 mm vs R1
                        clusterOneMMX_vs_R[i]->Fill(reco_il1, NClusterOneMMX[i]);
                        clusterOneMMY_vs_R[i]->Fill(reco_il1, NClusterOneMMY[i]);
                        clusterOneMMXY_vs_R[i]->Fill(reco_il1, NClusterOneMMX[i] + NClusterOneMMY[i]);
                        clusterOneMMX_vs_2R[i]->Fill(sec_track_reco_il1, NClusterOneMMX[i]);
                        clusterOneMMY_vs_2R[i]->Fill(sec_track_reco_il1, NClusterOneMMY[i]);
                        clusterOneMMXY_vs_2R[i]->Fill(sec_track_reco_il1, NClusterOneMMX[i] + NClusterOneMMY[i]);
                    }
                    NutofCluster = event.tofPlus->NClusters[0][ON] + event.tofPlus->NClusters[0][OFF] +
                                event.tofPlus->NClusters[1][ON] + event.tofPlus->NClusters[1][OFF];
                    NltofCluster = event.tofPlus->NClusters[2][ON] + event.tofPlus->NClusters[2][OFF] +
                                event.tofPlus->NClusters[3][ON] + event.tofPlus->NClusters[3][OFF];
                    NtofClusterMatchTrack = event.tofPlus->NTrkClusters;
                    utofCluster_vs_R->Fill(reco_il1,NutofCluster); 
                    ltofCluster_vs_R->Fill(reco_il1,NltofCluster);
                    tofClusterNoTrack_vs_R->Fill(reco_il1,NutofCluster+NltofCluster-NtofClusterMatchTrack);            
                    utofCluster_vs_2R->Fill(sec_track_reco_il1,NutofCluster); 
                    ltofCluster_vs_2R->Fill(sec_track_reco_il1,NltofCluster);
                    tofClusterNoTrack_vs_2R->Fill(sec_track_reco_il1,NutofCluster+NltofCluster-NtofClusterMatchTrack);
                    TRDHitTot_vs_R->Fill(reco_il1,trdHit);
                    TRDHitOn_vs_R->Fill(reco_il1,trdHitOn);
                    TRDHitOff_vs_R->Fill(reco_il1,trdHitOff);
                    TRDAmpTot_vs_R->Fill(reco_il1,ampTot);
                    TRDAmpOn_vs_R->Fill(reco_il1,ampOn);
                    TRDAmpOff_vs_R->Fill(reco_il1,ampOff);
                    TRDChargeTot_vs_R->Fill(reco_il1,chargeTot);
                    TRDChargeUp_vs_R->Fill(reco_il1,chargeUp);
                    TRDChargeLow_vs_R->Fill(reco_il1,chargeLow);
                    TRDHitTot_vs_2R->Fill(sec_track_reco_il1,trdHit);
                    TRDHitOn_vs_2R->Fill(sec_track_reco_il1,trdHitOn);
                    TRDHitOff_vs_2R->Fill(sec_track_reco_il1,trdHitOff); 
                    TRDAmpTot_vs_2R->Fill(sec_track_reco_il1,ampTot);
                    TRDAmpOn_vs_2R->Fill(sec_track_reco_il1,ampOn);
                    TRDAmpOff_vs_2R->Fill(sec_track_reco_il1,ampOff);
                    TRDChargeTot_vs_2R->Fill(sec_track_reco_il1,chargeTot);
                    TRDChargeUp_vs_2R->Fill(sec_track_reco_il1,chargeUp);
                    TRDChargeLow_vs_2R->Fill(sec_track_reco_il1,chargeLow);
                    if (trdHitOff < 50 || trdHitOff > 180) {
                        TRDChargeTot_vs_R_offCut->Fill(reco_il1,chargeTot);
                        TRDChargeUp_vs_R_offCut->Fill(reco_il1,chargeUp);
                        TRDChargeLow_vs_R_offCut->Fill(reco_il1,chargeLow);
                    }
                    if (hasGoodSecondTrack(event)) {

                        firstTrackVsSecond_2ndGood->Fill(reco_il1,sec_track_reco_il1);
						z1_inn = event.trTrackBase->InnerCharge[CRT];
						z2_inn = event.secondTrTrackBase->InnerCharge[CRT];
						firstTrack_Zinn->Fill(z1_inn);
						secondTrack_Zinn->Fill(z2_inn);
						firstVsSecond_Zinn->Fill(z1_inn,z2_inn);
                        if (event.secondTrTrackBase->InnerCharge[CRT]>1.8) {
							firstTrackVsSecond_2ndGood_ZGT1_8->Fill(reco_il1,sec_track_reco_il1);
						}

                        if (secTrackOnDiagonal(event)) {
                            for (int i = 0; i < 20; i++) {
                                edep_vs_R_secTrack[i]->Fill(reco_il1,Edep[i]);
                                edep_vs_2R_secTrack[i]->Fill(sec_track_reco_il1,Edep[i]);
                            }
                            for (int i = 0; i < 9; ++i) {
                                //All
                                clusterAllX_vs_R_secTrack[i]->Fill(reco_il1, NClusterX[i]);
                                clusterAllY_vs_R_secTrack[i]->Fill(reco_il1, NClusterY[i]);
                                clusterAllXY_vs_R_secTrack[i]->Fill(reco_il1, NClusterX[i] + NClusterY[i]);
                                clusterAllX_vs_2R_secTrack[i]->Fill(sec_track_reco_il1, NClusterX[i]);
                                clusterAllY_vs_2R_secTrack[i]->Fill(sec_track_reco_il1, NClusterY[i]);
                                clusterAllXY_vs_2R_secTrack[i]->Fill(sec_track_reco_il1, NClusterX[i] + NClusterY[i]);

                                //10 cm vs R1
                                clusterTenX_vs_R_secTrack[i]->Fill(reco_il1, NClusterTenX[i]);
                                clusterTenY_vs_R_secTrack[i]->Fill(reco_il1, NClusterTenY[i]);
                                clusterTenXY_vs_R_secTrack[i]->Fill(reco_il1, NClusterTenX[i] + NClusterTenY[i]);
                                clusterTenX_vs_2R_secTrack[i]->Fill(sec_track_reco_il1, NClusterTenX[i]);
                                clusterTenY_vs_2R_secTrack[i]->Fill(sec_track_reco_il1, NClusterTenY[i]);
                                clusterTenXY_vs_2R_secTrack[i]->Fill(sec_track_reco_il1, NClusterTenX[i] + NClusterTenY[i]);

                                //5 cm vs R1
                                clusterFiveX_vs_R_secTrack[i]->Fill(reco_il1, NClusterFiveX[i]);
                                clusterFiveY_vs_R_secTrack[i]->Fill(reco_il1, NClusterFiveY[i]);
                                clusterFiveXY_vs_R_secTrack[i]->Fill(reco_il1, NClusterFiveX[i] + NClusterFiveY[i]);
                                clusterFiveX_vs_2R_secTrack[i]->Fill(sec_track_reco_il1, NClusterFiveX[i]);
                                clusterFiveY_vs_2R_secTrack[i]->Fill(sec_track_reco_il1, NClusterFiveY[i]);
                                clusterFiveXY_vs_2R_secTrack[i]->Fill(sec_track_reco_il1, NClusterFiveX[i] + NClusterFiveY[i]);

                                //2 cm vs R1
                                clusterTwoX_vs_R_secTrack[i]->Fill(reco_il1, NClusterTwoX[i]);
                                clusterTwoY_vs_R_secTrack[i]->Fill(reco_il1, NClusterTwoY[i]);
                                clusterTwoXY_vs_R_secTrack[i]->Fill(reco_il1, NClusterTwoX[i] + NClusterTwoY[i]);
                                clusterTwoX_vs_2R_secTrack[i]->Fill(sec_track_reco_il1, NClusterTwoX[i]);
                                clusterTwoY_vs_2R_secTrack[i]->Fill(sec_track_reco_il1, NClusterTwoY[i]);
                                clusterTwoXY_vs_2R_secTrack[i]->Fill(sec_track_reco_il1, NClusterTwoX[i] + NClusterTwoY[i]);

                                //1 cm vs R1
                                clusterOneX_vs_R_secTrack[i]->Fill(reco_il1, NClusterOneX[i]);
                                clusterOneY_vs_R_secTrack[i]->Fill(reco_il1, NClusterOneY[i]);
                                clusterOneXY_vs_R_secTrack[i]->Fill(reco_il1, NClusterOneX[i] + NClusterOneY[i]);
                                clusterOneX_vs_2R_secTrack[i]->Fill(sec_track_reco_il1, NClusterOneX[i]);
                                clusterOneY_vs_2R_secTrack[i]->Fill(sec_track_reco_il1, NClusterOneY[i]);
                                clusterOneXY_vs_2R_secTrack[i]->Fill(sec_track_reco_il1, NClusterOneX[i] + NClusterOneY[i]);

                                //1 mm vs R1
                                clusterOneMMX_vs_R_secTrack[i]->Fill(reco_il1, NClusterOneMMX[i]);
                                clusterOneMMY_vs_R_secTrack[i]->Fill(reco_il1, NClusterOneMMY[i]);
                                clusterOneMMXY_vs_R_secTrack[i]->Fill(reco_il1, NClusterOneMMX[i] + NClusterOneMMY[i]);
                                clusterOneMMX_vs_2R_secTrack[i]->Fill(sec_track_reco_il1, NClusterOneMMX[i]);
                                clusterOneMMY_vs_2R_secTrack[i]->Fill(sec_track_reco_il1, NClusterOneMMY[i]);
                                clusterOneMMXY_vs_2R_secTrack[i]->Fill(sec_track_reco_il1, NClusterOneMMX[i] + NClusterOneMMY[i]);
                                
                            }
                            utofCluster_vs_R_secTrack->Fill(reco_il1,NutofCluster); 
                            ltofCluster_vs_R_secTrack->Fill(reco_il1,NltofCluster);
                            tofClusterNoTrack_vs_R_secTrack->Fill(reco_il1,NutofCluster+NltofCluster-NtofClusterMatchTrack);
                            utofCluster_vs_2R_secTrack->Fill(sec_track_reco_il1,NutofCluster); 
                            ltofCluster_vs_2R_secTrack->Fill(sec_track_reco_il1,NltofCluster);
                            tofClusterNoTrack_vs_2R_secTrack->Fill(sec_track_reco_il1,NutofCluster+NltofCluster-NtofClusterMatchTrack);

                            TRDHitTot_vs_R_secTrack->Fill(reco_il1,trdHit);
                            TRDHitOn_vs_R_secTrack->Fill(reco_il1,trdHitOn);
                            TRDHitOff_vs_R_secTrack->Fill(reco_il1,trdHitOff);
                            TRDAmpTot_vs_R_secTrack->Fill(reco_il1,ampTot);
                            TRDAmpOn_vs_R_secTrack->Fill(reco_il1,ampOn);
                            TRDAmpOff_vs_R_secTrack->Fill(reco_il1,ampOff);
                            TRDChargeTot_vs_R_secTrack->Fill(reco_il1,chargeTot);
                            TRDChargeUp_vs_R_secTrack->Fill(reco_il1,chargeUp);
                            TRDChargeLow_vs_R_secTrack->Fill(reco_il1,chargeLow);
                            TRDHitTot_vs_2R_secTrack->Fill(sec_track_reco_il1,trdHit);
                            TRDHitOn_vs_2R_secTrack->Fill(sec_track_reco_il1,trdHitOn);
                            TRDHitOff_vs_2R_secTrack->Fill(sec_track_reco_il1,trdHitOff); 
                            TRDAmpTot_vs_2R_secTrack->Fill(sec_track_reco_il1,ampTot);
                            TRDAmpOn_vs_2R_secTrack->Fill(sec_track_reco_il1,ampOn);
                            TRDAmpOff_vs_2R_secTrack->Fill(sec_track_reco_il1,ampOff);
                            TRDChargeTot_vs_2R_secTrack->Fill(sec_track_reco_il1,chargeTot);
                            TRDChargeUp_vs_2R_secTrack->Fill(sec_track_reco_il1,chargeUp);
                            TRDChargeLow_vs_2R_secTrack->Fill(sec_track_reco_il1,chargeLow);
                        }
                    }
                }
            }
		} else {
            gen = event.mcTruthBase->Primary.GetGenMomentum()/event.mcTruthBase->Primary.Z;
            //Check that the generated is S, is still S at L1 and at the L2 I have P
            if (NAIA::ContainsKeys(event.mcTruthPlus->TrackMCHits, 0)) { //If there is a cluster at L1
                if (NAIA::ContainsKeys(event.mcTruthPlus->TrackMCHits, 1)) { //If there is a cluster at L2
                    NAIA::MCParticle primary = event.mcTruthBase->Primary;
                    auto generated_charge    = primary.Z; //primary true charge
                    NAIA::MCTruthPlusData::TrMCHit mchit_l1 = event.mcTruthPlus->TrackMCHits[0]; //Take the cluster at L1
                    NAIA::MCTruthPlusData::TrMCHit mchit_l2 = event.mcTruthPlus->TrackMCHits[1]; //Take the cluster at L2
                    auto charge_at_layer1 = pidToNucleusCharge(mchit_l1.pID); 
                    auto charge_at_layer2 = pidToNucleusCharge(mchit_l2.pID); 
                    if (charge_at_layer2 != 15) continue;
                    bool is_charge_at_l1_as_generated = (charge_at_layer1 == generated_charge);
                    bool good_interaction = (charge_at_layer1 != charge_at_layer2);
                    if (good_interaction) {
                        if (is_charge_at_l1_as_generated) {
                            auto pub_flux = BuildPublishedFlux(generated_charge);
                            //-------Creating model-----------------
                            Z    = charge;
                            A    = charge*2;
                            mass = 0.938272075*A;
                            auto flux_model = new TF1("flux_exp", force_field, prod_min, prod_max, 3);
                            double exponent = -2.7;
                            double phi      = 0.6;
                            flux_model->SetParameters(1., phi, exponent);
                            flux_model->SetParNames("norm","phi","exp");
                            flux_model->SetRange(prod_min,1000);
                            flux_model->SetNpx(100000);
                            pub_flux->Fit(flux_model,"QN");
                            double norm = flux_model->Integral(prod_min,prod_max);
                            double weight = flux_model->Eval(gen)*gen*TMath::Log(prod_max/prod_min)/norm;
                            delete pub_flux;
                            delete flux_model;

                            reco_il1  = event.trTrackBase->Rigidity[FIT][IL1];
                            int rbin = hist_rig_highZ->FindBin(reco_il1);
                            double rlowedge = hist_rig_highZ->GetBinLowEdge(rbin);
                            sec_track_reco_il1 = event.secondTrTrackBase->Rigidity[FIT][IL1];
                            trdHitOn = event.trdKBase->NHits[TRD_ONTRACK];
                            trdHitOff= event.trdKBase->NHits[TRD_OFFTRACK];;
                            trdHit = trdHitOn+trdHitOff;
                            ampOn  = event.trdKBase->Amps[TRD_ONTRACK];
                            ampOff = event.trdKBase->Amps[TRD_OFFTRACK];
                            ampTot = ampOn+ampOff;
                            chargeUp = event.trdKBase->Charge[TRD_UP];
                            chargeLow= event.trdKBase->Charge[TRD_LOW];
                            chargeTot = event.trdKBase->Charge[TRD_TOT];
                            for (int i = 0; i < 20; i++) {
                                Edep[i] = event.trdKBase->Edep[i];
                                edep_vs_R[i]->Fill(reco_il1,Edep[i],weight);
                                edep_vs_2R[i]->Fill(sec_track_reco_il1,Edep[i]);
                            }
                            for (int i = 0; i < 9; ++i) {
                                NClusterX[i] = event.trTrackPlus->NClusters[i][ALL_LAYER][XSD];
                                NClusterY[i] = event.trTrackPlus->NClusters[i][ALL_LAYER][YSD];
                                clusterAllX_vs_R[i]->Fill(reco_il1, NClusterX[i],weight);
                                clusterAllY_vs_R[i]->Fill(reco_il1, NClusterY[i],weight);
                                clusterAllXY_vs_R[i]->Fill(reco_il1, NClusterX[i] + NClusterY[i],weight);
                                clusterAllX_vs_2R[i]->Fill(sec_track_reco_il1, NClusterX[i]);
                                clusterAllY_vs_2R[i]->Fill(sec_track_reco_il1, NClusterY[i]);
                                clusterAllXY_vs_2R[i]->Fill(sec_track_reco_il1, NClusterX[i] + NClusterY[i]);
                            }
                            NutofCluster = event.tofPlus->NClusters[0][ON] + event.tofPlus->NClusters[0][OFF] +
                                        event.tofPlus->NClusters[1][ON] + event.tofPlus->NClusters[1][OFF];
                            NltofCluster = event.tofPlus->NClusters[2][ON] + event.tofPlus->NClusters[2][OFF] +
                                        event.tofPlus->NClusters[3][ON] + event.tofPlus->NClusters[3][OFF];
                            NtofClusterMatchTrack = event.tofPlus->NTrkClusters;
                            utofCluster_vs_R->Fill(reco_il1,NutofCluster,weight); 
                            ltofCluster_vs_R->Fill(reco_il1,NltofCluster,weight);
                            tofClusterNoTrack_vs_R->Fill(reco_il1,NutofCluster+NltofCluster-NtofClusterMatchTrack,weight);            
                            utofCluster_vs_2R->Fill(sec_track_reco_il1,NutofCluster); 
                            ltofCluster_vs_2R->Fill(sec_track_reco_il1,NltofCluster);
                            tofClusterNoTrack_vs_2R->Fill(sec_track_reco_il1,NutofCluster+NltofCluster-NtofClusterMatchTrack);
                            TRDHitTot_vs_R->Fill(reco_il1,trdHit,weight);
                            TRDHitOn_vs_R->Fill(reco_il1,trdHitOn,weight);
                            TRDHitOff_vs_R->Fill(reco_il1,trdHitOff,weight);
                            TRDAmpTot_vs_R->Fill(reco_il1,ampTot,weight);
                            TRDAmpOn_vs_R->Fill(reco_il1,ampOn,weight);
                            TRDAmpOff_vs_R->Fill(reco_il1,ampOff,weight);
                            TRDChargeTot_vs_R->Fill(reco_il1,chargeTot,weight);
                            TRDChargeUp_vs_R->Fill(reco_il1,chargeUp,weight);
                            TRDChargeLow_vs_R->Fill(reco_il1,chargeLow,weight);
                            TRDHitTot_vs_2R->Fill(sec_track_reco_il1,trdHit);
                            TRDHitOn_vs_2R->Fill(sec_track_reco_il1,trdHitOn);
                            TRDHitOff_vs_2R->Fill(sec_track_reco_il1,trdHitOff); 
                            TRDAmpTot_vs_2R->Fill(sec_track_reco_il1,ampTot);
                            TRDAmpOn_vs_2R->Fill(sec_track_reco_il1,ampOn);
                            TRDAmpOff_vs_2R->Fill(sec_track_reco_il1,ampOff);
                            TRDChargeTot_vs_2R->Fill(sec_track_reco_il1,chargeTot);
                            TRDChargeUp_vs_2R->Fill(sec_track_reco_il1,chargeUp);
                            TRDChargeLow_vs_2R->Fill(sec_track_reco_il1,chargeLow);
                            if (trdHitOff < 50 || trdHitOff > 180) {
                                TRDChargeTot_vs_R_offCut->Fill(reco_il1,chargeTot,weight);
                                TRDChargeUp_vs_R_offCut->Fill(reco_il1,chargeUp,weight);
                                TRDChargeLow_vs_R_offCut->Fill(reco_il1,chargeLow,weight);
                            }
                            if (hasGoodSecondTrack(event)) {
                                if (secTrackOnDiagonal(event)) {
                                    for (int i = 0; i < 20; i++) {
                                        edep_vs_R_secTrack[i]->Fill(reco_il1,Edep[i],weight);
                                        edep_vs_2R_secTrack[i]->Fill(sec_track_reco_il1,Edep[i]);
                                    }
                                    for (int i = 0; i < 9; ++i) {
                                        clusterAllX_vs_R_secTrack[i]->Fill(reco_il1, NClusterX[i],weight);
                                        clusterAllY_vs_R_secTrack[i]->Fill(reco_il1, NClusterY[i],weight);
                                        clusterAllXY_vs_R_secTrack[i]->Fill(reco_il1, NClusterX[i] + NClusterY[i],weight);

                                        clusterAllX_vs_2R_secTrack[i]->Fill(sec_track_reco_il1, NClusterX[i]);
                                        clusterAllY_vs_2R_secTrack[i]->Fill(sec_track_reco_il1, NClusterY[i]);
                                        clusterAllXY_vs_2R_secTrack[i]->Fill(sec_track_reco_il1, NClusterX[i] + NClusterY[i]);
                                    }
                                    utofCluster_vs_R_secTrack->Fill(reco_il1,NutofCluster,weight); 
                                    ltofCluster_vs_R_secTrack->Fill(reco_il1,NltofCluster,weight);
                                    tofClusterNoTrack_vs_R_secTrack->Fill(reco_il1,NutofCluster+NltofCluster-NtofClusterMatchTrack,weight);
                                    utofCluster_vs_2R_secTrack->Fill(sec_track_reco_il1,NutofCluster); 
                                    ltofCluster_vs_2R_secTrack->Fill(sec_track_reco_il1,NltofCluster);
                                    tofClusterNoTrack_vs_2R_secTrack->Fill(sec_track_reco_il1,NutofCluster+NltofCluster-NtofClusterMatchTrack);

                                    TRDHitTot_vs_R_secTrack->Fill(reco_il1,trdHit,weight);
                                    TRDHitOn_vs_R_secTrack->Fill(reco_il1,trdHitOn,weight);
                                    TRDHitOff_vs_R_secTrack->Fill(reco_il1,trdHitOff,weight);
                                    TRDAmpTot_vs_R_secTrack->Fill(reco_il1,ampTot,weight);
                                    TRDAmpOn_vs_R_secTrack->Fill(reco_il1,ampOn,weight);
                                    TRDAmpOff_vs_R_secTrack->Fill(reco_il1,ampOff,weight);
                                    TRDChargeTot_vs_R_secTrack->Fill(reco_il1,chargeTot,weight);
                                    TRDChargeUp_vs_R_secTrack->Fill(reco_il1,chargeUp,weight);
                                    TRDChargeLow_vs_R_secTrack->Fill(reco_il1,chargeLow,weight);
                                    TRDHitTot_vs_2R_secTrack->Fill(sec_track_reco_il1,trdHit);
                                    TRDHitOn_vs_2R_secTrack->Fill(sec_track_reco_il1,trdHitOn);
                                    TRDHitOff_vs_2R_secTrack->Fill(sec_track_reco_il1,trdHitOff); 
                                    TRDAmpTot_vs_2R_secTrack->Fill(sec_track_reco_il1,ampTot);
                                    TRDAmpOn_vs_2R_secTrack->Fill(sec_track_reco_il1,ampOn);
                                    TRDAmpOff_vs_2R_secTrack->Fill(sec_track_reco_il1,ampOff);
                                    TRDChargeTot_vs_2R_secTrack->Fill(sec_track_reco_il1,chargeTot);
                                    TRDChargeUp_vs_2R_secTrack->Fill(sec_track_reco_il1,chargeUp);
                                    TRDChargeLow_vs_2R_secTrack->Fill(sec_track_reco_il1,chargeLow);
                                }
                            }
                    }
                    } else {
                        continue;
                    }
                }
            }
        }
    }
	//WRITE//
	auto outfile = new TFile(out, "recreate");
    //All
	for (int i = 0; i < 9; ++i)
        outfile->WriteTObject(clusterAllX_vs_R[i]);
    for (int i = 0; i < 9; ++i)
        outfile->WriteTObject(clusterAllY_vs_R[i]);
    for (int i = 0; i < 9; ++i)
        outfile->WriteTObject(clusterAllXY_vs_R[i]);
    for (int i = 0; i < 9; ++i)
        outfile->WriteTObject(clusterAllX_vs_R_secTrack[i]);
    for (int i = 0; i < 9; ++i)
        outfile->WriteTObject(clusterAllY_vs_R_secTrack[i]);
    for (int i = 0; i < 9; ++i)
        outfile->WriteTObject(clusterAllXY_vs_R_secTrack[i]);
	for (int i = 0; i < 9; ++i)
        outfile->WriteTObject(clusterAllX_vs_2R[i]);
    for (int i = 0; i < 9; ++i)
        outfile->WriteTObject(clusterAllY_vs_2R[i]);
    for (int i = 0; i < 9; ++i)
        outfile->WriteTObject(clusterAllXY_vs_2R[i]);
    for (int i = 0; i < 9; ++i)
        outfile->WriteTObject(clusterAllX_vs_2R_secTrack[i]);
    for (int i = 0; i < 9; ++i)
        outfile->WriteTObject(clusterAllY_vs_2R_secTrack[i]);
    for (int i = 0; i < 9; ++i)
        outfile->WriteTObject(clusterAllXY_vs_2R_secTrack[i]);    
    //10 cm
	for (int i = 0; i < 9; ++i)
        outfile->WriteTObject(clusterTenX_vs_R[i]);
    for (int i = 0; i < 9; ++i)
        outfile->WriteTObject(clusterTenY_vs_R[i]);
    for (int i = 0; i < 9; ++i)
        outfile->WriteTObject(clusterTenXY_vs_R[i]);
    for (int i = 0; i < 9; ++i)
        outfile->WriteTObject(clusterTenX_vs_R_secTrack[i]);
    for (int i = 0; i < 9; ++i)
        outfile->WriteTObject(clusterTenY_vs_R_secTrack[i]);
    for (int i = 0; i < 9; ++i)
        outfile->WriteTObject(clusterTenXY_vs_R_secTrack[i]);
	for (int i = 0; i < 9; ++i)
        outfile->WriteTObject(clusterTenX_vs_2R[i]);
    for (int i = 0; i < 9; ++i)
        outfile->WriteTObject(clusterTenY_vs_2R[i]);
    for (int i = 0; i < 9; ++i)
        outfile->WriteTObject(clusterTenXY_vs_2R[i]);
    for (int i = 0; i < 9; ++i)
        outfile->WriteTObject(clusterTenX_vs_2R_secTrack[i]);
    for (int i = 0; i < 9; ++i)
        outfile->WriteTObject(clusterTenY_vs_2R_secTrack[i]);
    for (int i = 0; i < 9; ++i)
        outfile->WriteTObject(clusterTenXY_vs_2R_secTrack[i]);  
    //5 cm
	for (int i = 0; i < 9; ++i)
        outfile->WriteTObject(clusterFiveX_vs_R[i]);
    for (int i = 0; i < 9; ++i)
        outfile->WriteTObject(clusterFiveY_vs_R[i]);
    for (int i = 0; i < 9; ++i)
        outfile->WriteTObject(clusterFiveXY_vs_R[i]);
    for (int i = 0; i < 9; ++i)
        outfile->WriteTObject(clusterFiveX_vs_R_secTrack[i]);
    for (int i = 0; i < 9; ++i)
        outfile->WriteTObject(clusterFiveY_vs_R_secTrack[i]);
    for (int i = 0; i < 9; ++i)
        outfile->WriteTObject(clusterFiveXY_vs_R_secTrack[i]);
	for (int i = 0; i < 9; ++i)
        outfile->WriteTObject(clusterFiveX_vs_2R[i]);
    for (int i = 0; i < 9; ++i)
        outfile->WriteTObject(clusterFiveY_vs_2R[i]);
    for (int i = 0; i < 9; ++i)
        outfile->WriteTObject(clusterFiveXY_vs_2R[i]);
    for (int i = 0; i < 9; ++i)
        outfile->WriteTObject(clusterFiveX_vs_2R_secTrack[i]);
    for (int i = 0; i < 9; ++i)
        outfile->WriteTObject(clusterFiveY_vs_2R_secTrack[i]);
    for (int i = 0; i < 9; ++i)
        outfile->WriteTObject(clusterFiveXY_vs_2R_secTrack[i]);     
    //2 cm
	for (int i = 0; i < 9; ++i)
        outfile->WriteTObject(clusterTwoX_vs_R[i]);
    for (int i = 0; i < 9; ++i)
        outfile->WriteTObject(clusterTwoY_vs_R[i]);
    for (int i = 0; i < 9; ++i)
        outfile->WriteTObject(clusterTwoXY_vs_R[i]);
    for (int i = 0; i < 9; ++i)
        outfile->WriteTObject(clusterTwoX_vs_R_secTrack[i]);
    for (int i = 0; i < 9; ++i)
        outfile->WriteTObject(clusterTwoY_vs_R_secTrack[i]);
    for (int i = 0; i < 9; ++i)
        outfile->WriteTObject(clusterTwoXY_vs_R_secTrack[i]);
	for (int i = 0; i < 9; ++i)
        outfile->WriteTObject(clusterTwoX_vs_2R[i]);
    for (int i = 0; i < 9; ++i)
        outfile->WriteTObject(clusterTwoY_vs_2R[i]);
    for (int i = 0; i < 9; ++i)
        outfile->WriteTObject(clusterTwoXY_vs_2R[i]);
    for (int i = 0; i < 9; ++i)
        outfile->WriteTObject(clusterTwoX_vs_2R_secTrack[i]);
    for (int i = 0; i < 9; ++i)
        outfile->WriteTObject(clusterTwoY_vs_2R_secTrack[i]);
    for (int i = 0; i < 9; ++i)
        outfile->WriteTObject(clusterTwoXY_vs_2R_secTrack[i]);   
    //1 cm
	for (int i = 0; i < 9; ++i)
        outfile->WriteTObject(clusterOneX_vs_R[i]);
    for (int i = 0; i < 9; ++i)
        outfile->WriteTObject(clusterOneY_vs_R[i]);
    for (int i = 0; i < 9; ++i)
        outfile->WriteTObject(clusterOneXY_vs_R[i]);
    for (int i = 0; i < 9; ++i)
        outfile->WriteTObject(clusterOneX_vs_R_secTrack[i]);
    for (int i = 0; i < 9; ++i)
        outfile->WriteTObject(clusterOneY_vs_R_secTrack[i]);
    for (int i = 0; i < 9; ++i)
        outfile->WriteTObject(clusterOneXY_vs_R_secTrack[i]);
	for (int i = 0; i < 9; ++i)
        outfile->WriteTObject(clusterOneX_vs_2R[i]);
    for (int i = 0; i < 9; ++i)
        outfile->WriteTObject(clusterOneY_vs_2R[i]);
    for (int i = 0; i < 9; ++i)
        outfile->WriteTObject(clusterOneXY_vs_2R[i]);
    for (int i = 0; i < 9; ++i)
        outfile->WriteTObject(clusterOneX_vs_2R_secTrack[i]);
    for (int i = 0; i < 9; ++i)
        outfile->WriteTObject(clusterOneY_vs_2R_secTrack[i]);
    for (int i = 0; i < 9; ++i)
        outfile->WriteTObject(clusterOneXY_vs_2R_secTrack[i]);   
     //1 mm
	for (int i = 0; i < 9; ++i)
        outfile->WriteTObject(clusterOneMMX_vs_R[i]);
    for (int i = 0; i < 9; ++i)
        outfile->WriteTObject(clusterOneMMY_vs_R[i]);
    for (int i = 0; i < 9; ++i)
        outfile->WriteTObject(clusterOneMMXY_vs_R[i]);
    for (int i = 0; i < 9; ++i)
        outfile->WriteTObject(clusterOneMMX_vs_R_secTrack[i]);
    for (int i = 0; i < 9; ++i)
        outfile->WriteTObject(clusterOneMMY_vs_R_secTrack[i]);
    for (int i = 0; i < 9; ++i)
        outfile->WriteTObject(clusterOneMMXY_vs_R_secTrack[i]);
	for (int i = 0; i < 9; ++i)
        outfile->WriteTObject(clusterOneMMX_vs_2R[i]);
    for (int i = 0; i < 9; ++i)
        outfile->WriteTObject(clusterOneMMY_vs_2R[i]);
    for (int i = 0; i < 9; ++i)
        outfile->WriteTObject(clusterOneMMXY_vs_2R[i]);
    for (int i = 0; i < 9; ++i)
        outfile->WriteTObject(clusterOneMMX_vs_2R_secTrack[i]);
    for (int i = 0; i < 9; ++i)
        outfile->WriteTObject(clusterOneMMY_vs_2R_secTrack[i]);
    for (int i = 0; i < 9; ++i)
        outfile->WriteTObject(clusterOneMMXY_vs_2R_secTrack[i]);  

    for (int i = 0; i <20; i++)
        outfile->WriteTObject(edep_vs_R[i]);
    for (int i = 0; i <20; i++)
        outfile->WriteTObject(edep_vs_2R[i]);
    for (int i = 0; i <20; i++)
        outfile->WriteTObject(edep_vs_R_secTrack[i]);
    for (int i = 0; i <20; i++)
        outfile->WriteTObject(edep_vs_2R_secTrack[i]);

    outfile->WriteTObject(utofCluster_vs_R, "utofCluster_vs_R");
	outfile->WriteTObject(ltofCluster_vs_R, "ltofCluster_vs_R");
	outfile->WriteTObject(tofClusterNoTrack_vs_R, "tofClusterNoTrack_vs_R");
    outfile->WriteTObject(utofCluster_vs_R_secTrack, "utofCluster_vs_R_secTrack");
	outfile->WriteTObject(ltofCluster_vs_R_secTrack, "ltofCluster_vs_R_secTrack");
	outfile->WriteTObject(tofClusterNoTrack_vs_R_secTrack, "tofClusterNoTrack_vs_R_secTrack");

    outfile->WriteTObject(utofCluster_vs_2R, "utofCluster_vs_2R");
	outfile->WriteTObject(ltofCluster_vs_2R, "ltofCluster_vs_2R");
	outfile->WriteTObject(tofClusterNoTrack_vs_2R, "tofClusterNoTrack_vs_2R");
    outfile->WriteTObject(utofCluster_vs_2R_secTrack, "utofCluster_vs_2R_secTrack");
	outfile->WriteTObject(ltofCluster_vs_2R_secTrack, "ltofCluster_vs_2R_secTrack");
	outfile->WriteTObject(tofClusterNoTrack_vs_2R_secTrack, "tofClusterNoTrack_vs_2R_secTrack");

    outfile->WriteTObject(TRDHitTot_vs_R,"TRDHitTot_vs_R");
    outfile->WriteTObject(TRDHitOn_vs_R,"TRDHitOn_vs_R");
    outfile->WriteTObject(TRDHitOff_vs_R,"TRDHitOff_vs_R");
    outfile->WriteTObject(TRDAmpTot_vs_R,"TRDAmpTot_vs_R");
    outfile->WriteTObject(TRDAmpOn_vs_R,"TRDAmpOn_vs_R");
    outfile->WriteTObject(TRDAmpOff_vs_R,"TRDAmpOff_vs_R");
    outfile->WriteTObject(TRDChargeTot_vs_R,"TRDChargeTot_vs_R");
    outfile->WriteTObject(TRDChargeUp_vs_R,"TRDChargeUp_vs_R");
    outfile->WriteTObject(TRDChargeLow_vs_R,"TRDChargeLow_vs_R");
    outfile->WriteTObject(TRDHitTot_vs_2R,"TRDHitTot_vs_2R");
    outfile->WriteTObject(TRDHitOn_vs_2R,"TRDHitOn_vs_2R");
    outfile->WriteTObject(TRDHitOff_vs_2R,"TRDHitOff_vs_2R");
    outfile->WriteTObject(TRDAmpTot_vs_2R,"TRDAmpTot_vs_2R");
    outfile->WriteTObject(TRDAmpOn_vs_2R,"TRDAmpOn_vs_2R");
    outfile->WriteTObject(TRDAmpOff_vs_2R,"TRDAmpOff_vs_2R");
    outfile->WriteTObject(TRDChargeTot_vs_2R,"TRDChargeTot_vs_2R");
    outfile->WriteTObject(TRDChargeUp_vs_2R,"TRDChargeUp_vs_2R");
    outfile->WriteTObject(TRDChargeLow_vs_2R,"TRDChargeLow_vs_2R");

    outfile->WriteTObject(TRDHitTot_vs_R_secTrack,"TRDHitTot_vs_R_secTrack");
    outfile->WriteTObject(TRDHitOn_vs_R_secTrack,"TRDHitOn_vs_R_secTrack");
    outfile->WriteTObject(TRDHitOff_vs_R_secTrack,"TRDHitOff_vs_R_secTrack");
    outfile->WriteTObject(TRDAmpTot_vs_R_secTrack,"TRDAmpTot_vs_R_secTrack");
    outfile->WriteTObject(TRDAmpOn_vs_R_secTrack,"TRDAmpOn_vs_R_secTrack");
    outfile->WriteTObject(TRDAmpOff_vs_R_secTrack,"TRDAmpOff_vs_R_secTrack");
    outfile->WriteTObject(TRDChargeTot_vs_R_secTrack,"TRDChargeTot_vs_R_secTrack");
    outfile->WriteTObject(TRDChargeUp_vs_R_secTrack,"TRDChargeUp_vs_R_secTrack");
    outfile->WriteTObject(TRDChargeLow_vs_R_secTrack,"TRDChargeLow_vs_R_secTrack");
    outfile->WriteTObject(TRDHitTot_vs_2R,"TRDHitTot_vs_2R");
    outfile->WriteTObject(TRDHitOn_vs_2R,"TRDHitOn_vs_2R");
    outfile->WriteTObject(TRDHitOff_vs_2R,"TRDHitOff_vs_2R");
    outfile->WriteTObject(TRDAmpTot_vs_2R,"TRDAmpTot_vs_2R");
    outfile->WriteTObject(TRDAmpOn_vs_2R,"TRDAmpOn_vs_2R");
    outfile->WriteTObject(TRDAmpOff_vs_2R,"TRDAmpOff_vs_2R");
    outfile->WriteTObject(TRDChargeTot_vs_2R,"TRDChargeTot_vs_2R");
    outfile->WriteTObject(TRDChargeUp_vs_2R,"TRDChargeUp_vs_2R");
    outfile->WriteTObject(TRDChargeLow_vs_2R,"TRDChargeLow_vs_2R");

    outfile->WriteTObject(TRDChargeTot_vs_R_offCut,"TRDChargeTot_vs_R_offCut");
    outfile->WriteTObject(TRDChargeUp_vs_R_offCut,"TRDChargeUp_vs_R_offCut");
    outfile->WriteTObject(TRDChargeLow_vs_R_offCut,"TRDChargeLow_vs_R_offCut");

    outfile->WriteTObject(firstTrack_Zinn,"firstTrack_Zinn");
	outfile->WriteTObject(secondTrack_Zinn,"secondTrack_Zinn");
	outfile->WriteTObject(firstVsSecond_Zinn,"firstVsSecond_Zinn");
	outfile->WriteTObject(firstTrackVsSecond, "firstTrackVsSecond");
	outfile->WriteTObject(firstTrackVsSecond_2ndGood, "firstTrackVsSecond_2ndGood");
	outfile->WriteTObject(firstTrackVsSecond_2ndGood_ZGT1_8, "firstTrackVsSecond_2ndGood_ZGT1_8");

	outfile->Close();
}
int pidToNucleusCharge(int pID) {
    switch(pID) {
        case 11:   // e-
            return -1;
        case -2112: // anti-neutron
            return 0;
        case 2212: // p
            return 1;
        case +1000020030: // HE3C         Z:   2 A:   3 mass: 2.809230
            return 2;
        case +1000020040: // ALPHA
            return 2;
        case +1000030070: // LI7
            return 3;
        case +1000040070: // BE7
            return 4;
        case +1000040090: // BE9
            return 4;
        case +1000040100: // BE10
            return 4;
        case +1000050100: // B10
            return 5;
        case +1000050110: // B11
            return 5;
        case +1000060120: // C12
            return 6; 
        case +1000060130: // C13
            return 6;
        case +1000070140: // N14
            return 7;
        case +1000070150: // N15
            return 7;
        case +1000080160: // O16
            return 8;     
        case +1000080170: // O17
            return 8;
        case +1000080180: // O18
            return 8;   
        case +1000150310: // P31
            return 15;
        case +1000160320: // S32
            return 16;
        case +1000170350: // CL35
            return 17;
        case +1000180360: // AR36
            return 18;
        case +1000190390: // K39
            return 19;
        case +1000200400: // CA40
            return 20;
        case +1000250550: // MN55
            return 25;
        case +1000260560: // FE56
            return 26;
    }
    return 0;
}
TH1D *BuildPublishedFlux(unsigned int charge) {
    std::vector<double> lower_bounds;
    std::vector<double> upper_bounds;
    std::vector<double> values;
    std::vector<double> syst;
    std::vector<double> errors;
    readFluxTable(lower_bounds,upper_bounds,values,errors,charge,syst);
    auto published = buildPublished(lower_bounds,upper_bounds,values,errors);
    return published;
}