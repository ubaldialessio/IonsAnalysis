#include "ChargeAnalysisManager.h"
#include "definition.h"
/*#include <TFile.h>
#include <TDirectory.h>
#include <stdexcept>
#include <iostream>
#include "binning.h"
#include "StreamUtility.h"*/
const int nLogbinss = 1000;
double Logbinss[nLogbinss] = {
	0.2, 0.202075, 0.204171, 0.206288, 0.208428, 0.21059, 0.212775, 0.214982, 0.217212, 0.219465, 0.221741, 
	0.224041, 0.226365, 0.228713, 0.231085, 0.233482, 0.235904, 0.238351, 0.240823, 0.243321, 0.245845, 
	0.248395, 0.250972, 0.253575, 0.256206, 0.258863, 0.261548, 0.264261, 0.267002, 0.269772, 0.27257, 
	0.275397, 0.278254, 0.28114, 0.284056, 0.287003, 0.28998, 0.292988, 0.296027, 0.299097, 0.3022, 
	0.305334, 0.308502, 0.311702, 0.314935, 0.318201, 0.321502, 0.324837, 0.328206, 0.331611, 0.33505, 
	0.338526, 0.342037, 0.345585, 0.34917, 0.352792, 0.356451, 0.360148, 0.363884, 0.367658, 0.371472, 
	0.375325, 0.379218, 0.383152, 0.387126, 0.391142, 0.395199, 0.399298, 0.40344, 0.407625, 0.411853, 
	0.416125, 0.420441, 0.424802, 0.429209, 0.433661, 0.438159, 0.442704, 0.447296, 0.451936, 0.456623, 
	0.46136, 0.466145, 0.470981, 0.475866, 0.480802, 0.485789, 0.490828, 0.495919, 0.501063, 0.506261, 
	0.511512, 0.516818, 0.522179, 0.527595, 0.533068, 0.538597, 0.544184, 0.549828, 0.555531, 0.561294, 
	0.567116, 0.572999, 0.578942, 0.584947, 0.591015, 0.597145, 0.603339, 0.609597, 0.615921, 0.622309, 
	0.628764, 0.635286, 0.641876, 0.648534, 0.655261, 0.662058, 0.668925, 0.675864, 0.682874, 0.689957, 
	0.697114, 0.704345, 0.711651, 0.719033, 0.726491, 0.734027, 0.741641, 0.749334, 0.757106, 0.764959, 
	0.772894, 0.780911, 0.789011, 0.797195, 0.805464, 0.813819, 0.822261, 0.83079, 0.839407, 0.848114, 
	0.856912, 0.8658, 0.874781, 0.883855, 0.893023, 0.902286, 0.911645, 0.921101, 0.930655, 0.940309, 
	0.950062, 0.959917, 0.969874, 0.979934, 0.990099, 1.00037, 1.01075, 1.02123, 1.03182, 1.04252, 
	1.05334, 1.06426, 1.0753, 1.08646, 1.09773, 1.10911, 1.12062, 1.13224, 1.14399, 1.15585, 
	1.16784, 1.17996, 1.19219, 1.20456, 1.21706, 1.22968, 1.24244, 1.25532, 1.26834, 1.2815, 
	1.29479, 1.30822, 1.32179, 1.3355, 1.34936, 1.36335, 1.37749, 1.39178, 1.40622, 1.42081, 
	1.43554, 1.45043, 1.46548, 1.48068, 1.49604, 1.51156, 1.52723, 1.54308, 1.55908, 1.57525, 
	1.59159, 1.6081, 1.62478, 1.64164, 1.65866, 1.67587, 1.69325, 1.71082, 1.72856, 1.74649, 
	1.76461, 1.78291, 1.80141, 1.82009, 1.83897, 1.85804, 1.87732, 1.89679, 1.91647, 1.93634, 
	1.95643, 1.97672, 1.99723, 2.01794, 2.03888, 2.06002, 2.08139, 2.10298, 2.1248, 2.14683, 
	2.1691, 2.1916, 2.21434, 2.2373, 2.26051, 2.28396, 2.30765, 2.33159, 2.35577, 2.38021, 
	2.4049, 2.42984, 2.45505, 2.48051, 2.50624, 2.53224, 2.5585, 2.58504, 2.61186, 2.63895, 
	2.66632, 2.69398, 2.72192, 2.75015, 2.77868, 2.8075, 2.83663, 2.86605, 2.89578, 2.92581, 
	2.95616, 2.98683, 3.01781, 3.04911, 3.08074, 3.11269, 3.14498, 3.1776, 3.21056, 3.24386, 
	3.27751, 3.31151, 3.34586, 3.38056, 3.41563, 3.45106, 3.48686, 3.52302, 3.55957, 3.59649, 
	3.63379, 3.67149, 3.70957, 3.74805, 3.78693, 3.82621, 3.86589, 3.90599, 3.94651, 3.98745, 
	4.02881, 4.0706, 4.11282, 4.15548, 4.19858, 4.24213, 4.28614, 4.3306, 4.37552, 4.4209, 
	4.46676, 4.51309, 4.5599, 4.6072, 4.65499, 4.70328, 4.75206, 4.80135, 4.85116, 4.90148, 
	4.95232, 5.00369, 5.05559, 5.10803, 5.16101, 5.21455, 5.26863, 5.32328, 5.3785, 5.43429, 
	5.49066, 5.54761, 5.60516, 5.6633, 5.72204, 5.78139, 5.84136, 5.90195, 5.96317, 6.02503, 
	6.08752, 6.15067, 6.21447, 6.27893, 6.34406, 6.40986, 6.47635, 6.54352, 6.6114, 6.67998, 
	6.74927, 6.81927, 6.89001, 6.96148, 7.03369, 7.10664, 7.18036, 7.25484, 7.33009, 7.40612, 
	7.48295, 7.56056, 7.63899, 7.71822, 7.79828, 7.87917, 7.9609, 8.04348, 8.12691, 8.21121, 
	8.29638, 8.38244, 8.46938, 8.55723, 8.646, 8.73568, 8.82629, 8.91784, 9.01035, 9.10381, 
	9.19824, 9.29365, 9.39005, 9.48745, 9.58586, 9.68529, 9.78575, 9.88726, 9.98982, 10.0934, 
	10.1981, 10.3039, 10.4108, 10.5188, 10.6279, 10.7381, 10.8495, 10.9621, 11.0758, 11.1906, 
	11.3067, 11.424, 11.5425, 11.6622, 11.7832, 11.9054, 12.0289, 12.1537, 12.2798, 12.4071, 
	12.5358, 12.6659, 12.7972, 12.93, 13.0641, 13.1996, 13.3365, 13.4749, 13.6146, 13.7558, 
	13.8985, 14.0427, 14.1884, 14.3355, 14.4842, 14.6345, 14.7863, 14.9396, 15.0946, 15.2512, 
	15.4094, 15.5692, 15.7307, 15.8939, 16.0587, 16.2253, 16.3936, 16.5636, 16.7355, 16.9091, 
	17.0844, 17.2617, 17.4407, 17.6216, 17.8044, 17.9891, 18.1757, 18.3642, 18.5547, 18.7471, 
	18.9416, 19.1381, 19.3366, 19.5372, 19.7398, 19.9446, 20.1515, 20.3605, 20.5717, 20.7851, 
	21.0007, 21.2185, 21.4386, 21.661, 21.8856, 22.1127, 22.342, 22.5738, 22.8079, 23.0445, 
	23.2835, 23.525, 23.7691, 24.0156, 24.2647, 24.5164, 24.7707, 25.0277, 25.2873, 25.5496, 
	25.8146, 26.0823, 26.3529, 26.6262, 26.9024, 27.1815, 27.4634, 27.7483, 28.0361, 28.3269, 
	28.6207, 28.9176, 29.2176, 29.5206, 29.8268, 30.1362, 30.4488, 30.7647, 31.0838, 31.4062, 
	31.732, 32.0611, 32.3937, 32.7297, 33.0692, 33.4122, 33.7588, 34.1089, 34.4627, 34.8202, 
	35.1814, 35.5463, 35.915, 36.2876, 36.664, 37.0443, 37.4285, 37.8168, 38.209, 38.6053, 
	39.0058, 39.4104, 39.8192, 40.2322, 40.6495, 41.0712, 41.4972, 41.9276, 42.3625, 42.8019, 
	43.2459, 43.6945, 44.1477, 44.6057, 45.0683, 45.5358, 46.0081, 46.4854, 46.9676, 47.4547, 
	47.947, 48.4443, 48.9468, 49.4545, 49.9675, 50.4858, 51.0095, 51.5386, 52.0732, 52.6133, 
	53.159, 53.7104, 54.2676, 54.8305, 55.3992, 55.9738, 56.5544, 57.1411, 57.7338, 58.3326, 
	58.9377, 59.549, 60.1667, 60.7908, 61.4214, 62.0585, 62.7022, 63.3526, 64.0097, 64.6737, 
	65.3445, 66.0223, 66.7072, 67.3991, 68.0982, 68.8046, 69.5182, 70.2393, 70.9679, 71.704, 
	72.4478, 73.1993, 73.9586, 74.7257, 75.5008, 76.284, 77.0752, 77.8747, 78.6825, 79.4986, 
	80.3232, 81.1564, 81.9982, 82.8488, 83.7081, 84.5764, 85.4537, 86.3401, 87.2357, 88.1405, 
	89.0548, 89.9785, 90.9118, 91.8549, 92.8076, 93.7703, 94.7429, 95.7257, 96.7186, 97.7219, 
	98.7355, 99.7596, 100.794, 101.84, 102.896, 103.964, 105.042, 106.132, 107.232, 108.345, 
	109.469, 110.604, 111.751, 112.91, 114.082, 115.265, 116.461, 117.669, 118.889, 120.122, 
	121.368, 122.627, 123.899, 125.184, 126.483, 127.795, 129.12, 130.46, 131.813, 133.18, 
	134.562, 135.957, 137.368, 138.793, 140.232, 141.687, 143.156, 144.641, 146.142, 147.658, 
	149.189, 150.737, 152.3, 153.88, 155.476, 157.089, 158.718, 160.365, 162.028, 163.709, 
	165.407, 167.123, 168.856, 170.608, 172.377, 174.165, 175.972, 177.797, 179.641, 181.505, 
	183.387, 185.29, 187.212, 189.153, 191.115, 193.098, 195.101, 197.125, 199.169, 201.235, 
	203.323, 205.432, 207.562, 209.715, 211.891, 214.089, 216.309, 218.553, 220.82, 223.11, 
	225.425, 227.763, 230.126, 232.513, 234.924, 237.361, 239.823, 242.311, 244.824, 247.364, 
	249.93, 252.522, 255.141, 257.788, 260.462, 263.163, 265.893, 268.651, 271.438, 274.253, 
	277.098, 279.972, 282.876, 285.811, 288.775, 291.771, 294.797, 297.855, 300.945, 304.066, 
	307.22, 310.407, 313.627, 316.88, 320.167, 323.488, 326.843, 330.233, 333.659, 337.12, 
	340.616, 344.15, 347.719, 351.326, 354.97, 358.652, 362.373, 366.131, 369.929, 373.766, 
	377.643, 381.56, 385.518, 389.517, 393.557, 397.64, 401.764, 405.932, 410.142, 414.397, 
	418.695, 423.038, 427.426, 431.86, 436.339, 440.865, 445.438, 450.058, 454.727, 459.444, 
	464.209, 469.024, 473.889, 478.805, 483.771, 488.789, 493.859, 498.982, 504.158, 509.387, 
	514.671, 520.01, 525.404, 530.853, 536.36, 541.923, 547.544, 553.224, 558.962, 564.76, 
	570.618, 576.537, 582.518, 588.56, 594.665, 600.833, 607.065, 613.362, 619.724, 626.153, 
	632.648, 639.21, 645.84, 652.539, 659.308, 666.147, 673.056, 680.038, 687.092, 694.219, 
	701.42, 708.695, 716.046, 723.474, 730.978, 738.56, 746.221, 753.961, 761.782, 769.684, 
	777.667, 785.734, 793.884, 802.119, 810.439, 818.845, 827.339, 835.921, 844.592, 853.352, 
	862.204, 871.147, 880.183, 889.313, 898.538, 907.858, 917.275, 926.79, 936.403, 946.116, 
	955.93, 965.845, 975.864, 985.986, 996.213, 1006.55, 1016.99, 1027.54, 1038.19, 1048.96, 
	1059.84, 1070.84, 1081.94, 1093.17, 1104.51, 1115.96, 1127.54, 1139.23, 1151.05, 1162.99, 
	1175.05, 1187.24, 1199.56, 1212, 1224.57, 1237.27, 1250.11, 1263.08, 1276.18, 1289.41, 
	1302.79, 1316.3, 1329.96, 1343.75, 1357.69, 1371.77, 1386, 1400.38, 1414.9, 1429.58, 
	1444.41, 1459.39, 1474.53, 1489.82, 1505.28, 1520.89, 1536.67, 1552.61, 1568.71, 1584.98, 
	1601.42, 1618.03, 1634.82, 1651.78, 1668.91, 1686.22, 1703.71, 1721.38, 1739.24, 1757.28, 
	1775.51, 1793.92, 1812.53, 1831.33, 1850.33, 1869.52, 1888.91, 1908.51, 1928.3, 1948.3, 
	1968.51, 1988.93, 2009.56, 2030.41, 2051.47, 2072.75, 2094.25, 2115.97, 2137.92, 2160.09, 
	2182.5, 2205.14, 2228.01, 2251.12, 2274.47, 2298.06, 2321.9, 2345.99, 2370.32, 2394.91, 
	2419.75, 2444.85, 2470.21, 2495.83, 2521.72, 2547.88, 2574.3, 2601.01, 2627.99, 2655.25, 
	2682.79, 2710.61, 2738.73, 2767.14, 2795.84, 2824.84, 2854.14, 2883.75, 2913.66, 2943.88, 
	2974.42, 3005.27, 3036.45, 3067.94, 3099.76, 3131.92, 3164.4, 3197.23, 3230.39, 3263.9, 
	3297.75, 3331.96, 3366.52, 3401.44, 3436.72, 3472.37, 3508.39, 3544.78, 3581.55, 3618.7, 
	3656.24, 3694.16, 3732.48, 3771.2, 3810.31, 3849.84, 3889.77, 3930.12, 3970.88, 4012.07, 
	4053.69, 4095.74, 4138.22, 4181.14, 4224.51, 4268.33, 4312.61, 4357.34, 4402.54, 4448.2, 
	4494.34, 4540.96, 4588.07, 4635.66, 4683.74, 4732.32, 4781.41, 4831.01, 4881.12, 4931.75, 
	4982.9, 5034.59, 5086.81, 5139.57, 5192.89, 5246.75, 5301.17, 5356.16, 5411.72, 5467.85, 
	5524.57, 5581.87, 5639.77, 5698.27, 5757.38, 5817.1, 5877.44, 5938.4, 6000};
double GetDataMasses(unsigned int charge) {
  double A = 0;
  switch (charge) {
  case 1:
    A = 1;
    break;
  case 2:
    A = 4;
    break;
  case 3:
    A = 6.5;
    break;
  case 4:
    A = 7.6;
    break;
  case 5:
    A = 10.65;
    break;
  case 6:
    A = 12;
    break;
  case 7:
    A = 14;
    break;
  case 8:
    A = 16;
    break;
  case 9:
    A = 18;
    break;
  case 10:
    A = 20;
    break;
  case 11:
    A = 22;
    break;
  case 12:
    A = 24;
    break;
  case 13:
    A = 26;
    break;
  case 14:
    A = 28;
    break;
  case 15:
    A = 31;
	  break;
  case 16:
    A = 32;
    break;
  case 17:
    A = 34;
    break;
  case 18:
    A = 36;
	  break;
  case 19:
    A = 38;
    break;
  case 20:
    A = 40;
    break;
  case 21:
    A = 42;
    break;
  case 22:
    A = 44;
    break;
  case 23:
    A = 46;
    break;
  case 24:
    A = 48;
    break;
  case 25:
    A = 50;
    break;
  case 26:
    A = 52;
    break;
  }
  return A;
}
auto hist_logg = new TH1D("", "", nLogbinss - 1, Logbinss);
ChargeAnalysisManager::ChargeAnalysisManager(const std::vector<unsigned int>& charges,
                                             int nTbins, const double* Tbins,
                                             int nRbins, const double* Rbins,
                                             bool isMC, TString exename)

    : fCharges(charges),
      fNTbins(nTbins),
      fNRbins(nRbins),
      fTbins(Tbins, Tbins + nTbins + 1),
      fRbins(Rbins, Rbins + nRbins + 1),
      fIsMC(isMC),
      fExename(exename) {
    for (auto c : charges) {

        if (!fIsMC) {
          if (fExename=="DistributionsSelector") BookDistributionHistograms(c);
          if (fExename=="TemplatesSelector") BookTemplateHistograms(c);
          else BookHistograms(c);
        } else {
          BookMcTrees(c);
          BookMcHist(c);
        }

        ChargeContext ctx;
        ctx.Amass = GetDataMasses(c);
        ctx.daq.fill(0);
        ctx.daq_event_size.fill(0.0);

        ctx.mysel = [sel = Mysel(c)]( NAIA::Event& e) mutable { return sel(e); };
        ctx.myselNoL1 = [sel = MyselNoL1(c)]( NAIA::Event& e) mutable { return sel(e); };

        ctx.num_tof = [sel = Num_tof(c)]( NAIA::Event& e) mutable { return sel(e); };
        ctx.den_tof = [sel = Den_tof(c)]( NAIA::Event& e) mutable { return sel(e); };
        ctx.den_tof_ltof = [sel = Den_tof_ltof(c)]( NAIA::Event& e) mutable { return sel(e); };
        ctx.den_tof_ltof_l9_05 = [sel = Den_tof_ltof_l9_05(c)]( NAIA::Event& e) mutable { return sel(e); };
        ctx.den_tof_ltof_l9_1 = [sel = Den_tof_ltof_l9_1(c)]( NAIA::Event& e) mutable { return sel(e); };

        ctx.num_l1 = [sel = Num_l1(c)]( NAIA::Event& e) mutable { return sel(e); };
        ctx.den_l1 = [sel = Den_l1(c)]( NAIA::Event& e) mutable { return sel(e); };
        ctx.den_l1_ltof = [sel = Den_l1_ltof(c)]( NAIA::Event& e) mutable { return sel(e); };
        ctx.den_l1_ltof_l9_05 = [sel = Den_l1_ltof_l9_05(c)]( NAIA::Event& e) mutable { return sel(e); };
        ctx.den_l1_ltof_l9_1 = [sel = Den_l1_ltof_l9_1(c)]( NAIA::Event& e) mutable { return sel(e); };

        ctx.num_l1Unb = [sel = Num_l1Unb(c)]( NAIA::Event& e) mutable { return sel(e); };
        ctx.den_l1Unb = [sel = Den_l1Unb(c)]( NAIA::Event& e) mutable { return sel(e); };
        ctx.den_l1Unb_ltof = [sel = Den_l1Unb_ltof(c)]( NAIA::Event& e) mutable { return sel(e); };
        ctx.den_l1Unb_ltof_l9_05 = [sel = Den_l1Unb_ltof_l9_05(c)]( NAIA::Event& e) mutable { return sel(e); };
        ctx.den_l1Unb_ltof_l9_1 = [sel = Den_l1Unb_ltof_l9_1(c)]( NAIA::Event& e) mutable { return sel(e); };

        ctx.num_track = [sel = Num_track(c)]( NAIA::Event& e) mutable { return sel(e); };
        ctx.den_track = [sel = Den_track(c)]( NAIA::Event& e) mutable { return sel(e); };

        ctx.num_trackCh = [sel = Num_trackCh(c)]( NAIA::Event& e) mutable { return sel(e); };
        ctx.den_trackCh = [sel = Den_trackCh(c)]( NAIA::Event& e) mutable { return sel(e); };
        ctx.den_trackCh_ltof = [sel = Den_trackCh_ltof(c)]( NAIA::Event& e) mutable { return sel(e); };
        ctx.den_trackCh_ltof_l9_05 = [sel = Den_trackCh_ltof_l9_05(c)]( NAIA::Event& e) mutable { return sel(e); };
        ctx.den_trackCh_ltof_l9_1 = [sel = Den_trackCh_ltof_l9_1(c)]( NAIA::Event& e) mutable { return sel(e); };

        ctx.den_trig = [sel = Den_trig(c)]( NAIA::Event& e) mutable { return sel(e); };

        ctx.l1_temp = [sel = L1_temp(c)]( NAIA::Event& e) mutable { return sel(e); };
        ctx.l2_temp = [sel = L2_temp(c)]( NAIA::Event& e) mutable { return sel(e); };


        contexts[c] = std::move(ctx);
    }
}

void ChargeAnalysisManager::BookHistograms(unsigned int charge) {
    std::vector<std::string> geoms = {"IL1", "FS"};
    std::vector<std::string> names = {
        "acc5_counters_pre", "acc5_counters", "acc8_counters",
        "acc5_clusters_pre", "acc5_clusters", "acc8_clusters",
        "nPhys_1", "nPhys_2",
        "jinj1_aboveBuffer", "jinj1_belowBuffer",
        "jinj2_aboveBuffer", "jinj2_belowBuffer",
        "lvt_25", "rigidity",
      
        "sample_tof_den","sample_tof_den_ltof","sample_tof_den_ltof_l9_05","sample_tof_den_ltof_l9_1",
        "pass_tof_den","pass_tof_den_ltof","pass_tof_den_ltof_l9_05","pass_tof_den_ltof_l9_1",

        "sample_l1_den","sample_l1_den_ltof","sample_l1_den_ltof_l9_05","sample_l1_den_ltof_l9_1",
        "pass_l1_den","pass_l1_den_ltof","pass_l1_den_ltof_l9_05","pass_l1_den_ltof_l9_1",

        "sample_tr", "pass_tr",

        "sample_l1Unb_den","sample_l1Unb_den_ltof","sample_l1Unb_den_ltof_l9_05","sample_l1Unb_den_ltof_l9_1",
        "pass_l1Unb_den","pass_l1Unb_den_ltof","pass_l1Unb_den_ltof_l9_05","pass_l1Unb_den_ltof_l9_1",

        "sample_trch_den","sample_trch_den_ltof","sample_trch_den_ltof_l9_05","sample_trch_den_ltof_l9_1",
        "pass_trch_den","pass_trch_den_ltof","pass_trch_den_ltof_l9_05","pass_trch_den_ltof_l9_1",

        "unbiased_trig",
        "phys_trig",
        "sample_trig",
        "phys_trig_acc4"
        "phys_trig_acc7"
    };

    for (auto& geom : geoms) {
        for (auto& name : names) {
            std::string hname = name;
            histos[charge][geom][name] =
                new TH2D(hname.c_str(), ";Utime;R(GV)",
                         fNTbins-1, fTbins.data(),
                         fNRbins-1, fRbins.data());
        }
    }
}
void ChargeAnalysisManager::BookDistributionHistograms(unsigned int charge) {
    std::vector<std::string> geoms = {"IL1", "FS"};
    std::vector<std::string> names = {
        "qL1_den", "qL1_den_ltof", "qL1_den_ltof_l9_05","qL1_den_ltof_l9_1",
        "qL1Unb_den","qL1Unb_den_ltof","qL1Unb_den_ltof_l9_05","qL1Unb_den_ltof_l9_1",
        "qtof_den","qtof_den_ltof","qtof_den_ltof_l9_05","qtof_den_ltof_l9_1",
        "qInn_R_den","qInn_R_den_ltof","qInn_R_den_ltof_l9_05","qInn_R_den_ltof_l9_1",
        "qUnbL9_R_den","qUnbL9_R_den_ltof","qUnbL9_R_den_ltof_l9_05","qUnbL9_R_den_ltof_l9_1"
    };

    for (auto& geom : geoms) {
        for (auto& name : names) {
            std::string hname = name;
            distributionHistos[charge][geom][name] =
                new TH1D(hname.c_str(), ";Q (c.u.)",550,0,28);
        }
    }
}
void ChargeAnalysisManager::BookTemplateHistograms(unsigned int charge) {
    std::vector<std::string> geoms = {"IL1", "FS"};

    // number of rigidity bins: use fRbins size-1
    int nRtemplateBins = static_cast<int>(fRbins.size()) - 1;
    if (nRtemplateBins <= 0) return;

    for (auto& geom : geoms) {
        // ensure an entry exists
        TemplateSet &set = templateSets[charge][geom];
        // clear any preexisting
        set.Clear();

        // L1 templates: create one histogram + kernel tree per rigidity template bin
        set.L1Templates.reserve(nRtemplateBins);
        set.L1KernelTrees.reserve(nRtemplateBins);
        set.L1charge.reserve(nRtemplateBins);
        set.L1weight.reserve(nRtemplateBins);

        for (int ibin = 0; ibin < nRtemplateBins; ++ibin) {
            double rmin = fRbins[ibin];
            double rmax = fRbins[ibin+1];
            TString hname = Form("L1Template_%03d", ibin);
            TString htitle = Form("%5.3f < R (GV) < %5.3f;Q_{L1X};Counts", rmin, rmax);
            TH1D* h = new TH1D(hname.Data(), htitle.Data(), 600, 0.0, 32.0);
            set.L1Templates.push_back(h);

            // kernel tree
            TTree* kt = new TTree(Form("%s_L1KTree_%03d", geom.c_str(), ibin),
                                  Form("%5.3f < R (GV) < %5.3f", rmin, rmax));
            kt->SetDirectory(0);
            // push default values to vectors so addresses are stable
            set.L1charge.push_back(0.0f);
            set.L1weight.push_back(0.0);
            int idx = static_cast<int>(set.L1KernelTrees.size());
            // create branches using pointers to the storage we just pushed
            kt->Branch("L1charge", &(set.L1charge[idx]), "L1charge/F");
            kt->Branch("L1weight", &(set.L1weight[idx]), "L1weight/D");
            set.L1KernelTrees.push_back(kt);
        }

        // L2 templates: same as L1 but different naming
        set.L2Templates.reserve(nRtemplateBins);
        set.L2KernelTrees.reserve(nRtemplateBins);
        set.L2charge.reserve(nRtemplateBins);
        set.L2weight.reserve(nRtemplateBins);

        for (int ibin = 0; ibin < nRtemplateBins; ++ibin) {
            double rmin = fRbins[ibin];
            double rmax = fRbins[ibin+1];
            TString hname = Form("L2Template_%03d", ibin);
            TString htitle = Form("%5.3f < R (GV) < %5.3f;Q_{L2X};Counts", rmin, rmax);
            TH1D* h = new TH1D(hname.Data(), htitle.Data(), 600, 0.0, 32.0);
            set.L2Templates.push_back(h);

            TTree* kt = new TTree(Form("%s_L2KTree_%03d", geom.c_str(), ibin),
                                  Form("%5.3f < R (GV) < %5.3f", rmin, rmax));
            kt->SetDirectory(0);
            set.L2charge.push_back(0.0f);
            set.L2weight.push_back(0.0);
            int idx = static_cast<int>(set.L2KernelTrees.size());
            kt->Branch("L2charge", &(set.L2charge[idx]), "L2charge/F");
            kt->Branch("L2weight", &(set.L2weight[idx]), "L2weight/D");
            set.L2KernelTrees.push_back(kt);
        }

        // L1 charge distribution histos + trees (same binning)
        set.L1Distributions.reserve(nRtemplateBins);
        set.L1DistKernelTrees.reserve(nRtemplateBins);
        set.L1Dcharge.reserve(nRtemplateBins);
        set.L1Dweight.reserve(nRtemplateBins);

        for (int ibin = 0; ibin < nRtemplateBins; ++ibin) {
            double rmin = fRbins[ibin];
            double rmax = fRbins[ibin+1];
            TString hname = Form("L1Distribution_%03d", ibin);
            TString htitle = Form("%5.3f < R (GV) < %5.3f;Q_{L1};Counts", rmin, rmax);
            TH1D* h = new TH1D(hname.Data(), htitle.Data(), 600, 0.0, 32.0);
            set.L1Distributions.push_back(h);

            TTree* kt = new TTree(Form("%s_L1DKTree_%03d", geom.c_str(), ibin),
                                  Form("%5.3f < R (GV) < %5.3f", rmin, rmax));
            kt->SetDirectory(0);
            set.L1Dcharge.push_back(0.0f);
            set.L1Dweight.push_back(0.0);
            int idx = static_cast<int>(set.L1DistKernelTrees.size());
            kt->Branch("L1Dcharge", &(set.L1Dcharge[idx]), "L1Dcharge/F");
            kt->Branch("L1Dweight", &(set.L1Dweight[idx]), "L1Dweight/D");
            set.L1DistKernelTrees.push_back(kt);
        }
    } // geom loop
}
void ChargeAnalysisManager::BookMcTrees(unsigned int charge) {
    std::vector<std::string> geoms = {"IL1", "FS"};

    for (auto& geom : geoms) {
        // --- Tree 1 ---
        {
            auto& s = mcSamples[charge][geom]["tree1"]; // holds variables
            TTree* t = new TTree(Form("tree1_Z%d_%s", charge, geom.c_str()), "MC tree 1");
            t->Branch("gen", &s.gen, "gen/D");
            t->Branch("den_species", &s.species, "species/I");
            t->Branch("den_p", &s.passed, "passed/S");
            t->Branch("den_inn", &s.reco_inn, "reco_inn/F");
            t->Branch("den_il1", &s.reco_il1, "reco_il1/F");
            t->Branch("den_b", &s.reco_beta, "reco_beta/F");
            mcTrees[charge][geom]["tree1"] = t;
        }

        // --- Tree 2 ---
        {
            auto& s = mcSamples[charge][geom]["tree2"];
            TTree* t = new TTree(Form("tree2_Z%d_%s", charge, geom.c_str()), "MC tree 2");
            t->Branch("den_ltof_species", &s.species, "species/I");
            t->Branch("den_ltof_p", &s.passed, "passed/S");
            t->Branch("den_ltof_inn", &s.reco_inn, "reco_inn/F");
            t->Branch("den_ltof_il1", &s.reco_il1, "reco_il1/F");
            t->Branch("den_ltof_b", &s.reco_beta, "reco_beta/F");
            mcTrees[charge][geom]["tree2"] = t;
        }

        // --- Tree 3 ---
        {
            auto& s = mcSamples[charge][geom]["tree3"];
            TTree* t = new TTree(Form("tree3_Z%d_%s", charge, geom.c_str()), "MC tree 3");
            t->Branch("den_ltof_l9Unb_05_species", &s.species, "species/I");
            t->Branch("den_ltof_l9Unb_05_p", &s.passed, "passed/S");
            t->Branch("den_ltof_l9Unb_05_inn", &s.reco_inn, "reco_inn/F");
            t->Branch("den_ltof_l9Unb_05_il1", &s.reco_il1, "reco_il1/F");
            t->Branch("den_ltof_l9Unb_05_b", &s.reco_beta, "reco_beta/F");
            mcTrees[charge][geom]["tree3"] = t;
        }

        // --- Tree 4 ---
        {
            auto& s = mcSamples[charge][geom]["tree4"];
            TTree* t = new TTree(Form("tree4_Z%d_%s", charge, geom.c_str()), "MC tree 4");
            t->Branch("den_ltof_l9Unb_1_species", &s.species, "species/I");
            t->Branch("den_ltof_l9Unb_1_p", &s.passed, "passed/S");
            t->Branch("den_ltof_l9Unb_1_inn", &s.reco_inn, "reco_inn/F");
            t->Branch("den_ltof_l9Unb_1_il1", &s.reco_il1, "reco_il1/F");
            t->Branch("den_ltof_l9Unb_1_b", &s.reco_beta, "reco_beta/F");
            mcTrees[charge][geom]["tree4"] = t;
        }
    }
}
void ChargeAnalysisManager::BookMcHist(unsigned int charge) {
    std::vector<std::string> geoms = {"IL1", "FS"};
    std::vector<std::string> names = {
        "mc_pass","mc_pass_gen","mc_samp"
    };
    for (auto& geom : geoms) {
        for (auto& name : names) {
            std::string hname = name;
            mcHist[charge][geom][name] =
                (TH1D*)hist_logg->Clone(name.c_str());
        }
    }  
}

void ChargeAnalysisManager::Fill(unsigned int charge,
                                 const std::string& geom,
                                 const std::string& name,
                                 double x, double y, double weight) {
    if (histos.count(charge) &&
        histos.at(charge).count(geom) &&
        histos.at(charge).at(geom).count(name)) 
    {
        histos[charge][geom][name]->Fill(x, y, weight);
    }
}
void ChargeAnalysisManager::Fill(unsigned int charge,
                                 const std::string& geom,
                                 const std::string& name,
                                 double x, double y) {
    ChargeAnalysisManager::Fill(charge,geom,name,x,y,1.);                                
}
void ChargeAnalysisManager::FillDistributions(unsigned int charge,
                                 const std::string& geom,
                                 const std::string& name,
                                 double x) {
    if (distributionHistos.count(charge) &&
        distributionHistos.at(charge).count(geom) &&
        distributionHistos.at(charge).at(geom).count(name)) 
    {
        distributionHistos[charge][geom][name]->Fill(x);
    }
}
void ChargeAnalysisManager::FillMcTrees(NAIA::NAIAChain &chain, NAIA::Event &event, const std::string &geom) {
    unsigned int nucleus = event.mcTruthBase->Primary.Z;
    if (nucleus == 0) return;

    auto ctxIt = contexts.find(nucleus);
    if (ctxIt == contexts.end()) return;
    ChargeContext &ctx = ctxIt->second;

    if (mcSamples.count(nucleus) == 0 ||
        mcSamples[nucleus].count(geom) == 0 ||
        mcTrees.count(nucleus) == 0 ||
        mcTrees[nucleus].count(geom) == 0 ||
        mcHist.count(nucleus) == 0 ||
        mcHist[nucleus].count(geom) == 0) {
        return;
    }

    auto &ms_map   = mcSamples[nucleus][geom];
    auto &tree_map = mcTrees[nucleus][geom];
    auto &hmap     = mcHist[nucleus][geom];

    McSample &den              = ms_map["tree1"];
    McSample &den_ltof         = ms_map["tree2"];
    McSample &den_ltof_l9Unb_05= ms_map["tree3"];
    McSample &den_ltof_l9Unb_1 = ms_map["tree4"];

    // Reset per-event (assumes Reset sets fields to zero but preserves branch bindings)
    for (auto &p : ms_map) p.second.Reset();

    // histo used by FillHistos (or dummy)
    TH2D* sample_tr_ptr = GetHist(nucleus, geom, "sample_tr");
    static TH2D dummy_sample_tr("mc_dummy_sample_tr","dummy;Utime;R", 2, 0, 1, 2, 0, 1);
    TH2D &h_for_fill = (sample_tr_ptr ? *sample_tr_ptr : dummy_sample_tr);

    TH1D* h_mc_pass     = (hmap.count("mc_pass")     ? hmap.at("mc_pass")     : nullptr);
    TH1D* h_mc_pass_gen = (hmap.count("mc_pass_gen") ? hmap.at("mc_pass_gen") : nullptr);

    // precompute common variables
    const float R_inn = event.trTrackBase->Rigidity[NAIA::TrTrack::Fit::GBL][NAIA::TrTrack::Span::InnerOnly];
    const float R_il1 = event.trTrackBase->Rigidity[NAIA::TrTrack::Fit::GBL][NAIA::TrTrack::Span::InnerL1];
    const double gen = event.mcTruthBase->Primary.GetGenMomentum() / double(nucleus);

    // --- Precompute all selection booleans ONCE ---
    const bool pass_myselNoL1 = ctx.myselNoL1(event);
    const bool den_l1            = ctx.den_l1(event);
    const bool den_tof           = ctx.den_tof(event);
    const bool den_trig          = ctx.den_trig(event) && (event.evSummary->PhysBPatt & 65);
    const bool den_track         = ctx.den_track(event);
    const bool den_l1Unb         = ctx.den_l1Unb(event);
    const bool den_trackCh       = ctx.den_trackCh(event);

    const bool den_l1_ltof       = ctx.den_l1_ltof(event);
    const bool den_tof_ltof      = ctx.den_tof_ltof(event);
    const bool den_l1Unb_ltof    = ctx.den_l1Unb_ltof(event);
    const bool den_trackCh_ltof  = ctx.den_trackCh_ltof(event);

    const bool den_l1_ltof_l9_05 = ctx.den_l1_ltof_l9_05(event);
    const bool den_tof_ltof_l9_05= ctx.den_tof_ltof_l9_05(event);
    const bool den_l1Unb_ltof_l9_05 = ctx.den_l1Unb_ltof_l9_05(event);
    const bool den_trackCh_ltof_l9_05 = ctx.den_trackCh_ltof_l9_05(event);

    const bool den_l1_ltof_l9_1  = ctx.den_l1_ltof_l9_1(event);
    const bool den_tof_ltof_l9_1 = ctx.den_tof_ltof_l9_1(event);
    const bool den_l1Unb_ltof_l9_1 = ctx.den_l1Unb_ltof_l9_1(event);
    const bool den_trackCh_ltof_l9_1 = ctx.den_trackCh_ltof_l9_1(event);

    // numerators (subconditions)
    const bool num_l1       = ctx.num_l1(event);
    const bool num_tof      = ctx.num_tof(event);
    const bool num_l1Unb    = ctx.num_l1Unb(event);
    const bool num_track    = ctx.num_track(event);
    const bool num_trackCh  = ctx.num_trackCh(event);

    const bool isPhysicsTrigger = event.evSummary->IsPhysicsTrigger();

    // Inline helper that mirrors your SetSample behavior but ensures species is written whenever passed is set
    auto set_sample = [&](McSample &sample, bool condition, int bit,
                        float &reco_var, float rigidity,
                        bool subcondition = false, int subbit = -1) {
        if (!condition) return;
        reco_var = rigidity;
        sample.passed |= (1 << bit);
        sample.species = nucleus;   // ALWAYS update species when we set a bit
        if (subcondition && subbit >= 0) {
            sample.passed |= (1 << subbit);
            sample.species = nucleus;
        }
    };

    // If base preselection passed, set bit 0 and species/gen for all samples
    if (pass_myselNoL1) {
        for (McSample* s : { &den, &den_ltof, &den_ltof_l9Unb_05, &den_ltof_l9Unb_1 }) {
            s->passed |= 1;         // base bit
            s->species = nucleus;   // important: always set species
            s->gen = gen;
        }
        if (h_mc_pass) h_mc_pass->Fill(R_il1);
        if (h_mc_pass_gen) h_mc_pass_gen->Fill(gen);
    }

    // Track-related heavy call: call FillHistos() at most once per event
    // If any sample requires den_track, we compute it once and reuse
    double track_beta = std::numeric_limits<double>::quiet_NaN();
    bool track_computed = false;
    auto ensure_track = [&](McSample &sample){
        if (!den_track) return;
        if (!track_computed) {
            track_beta = FillHistos(h_for_fill, nucleus, ctx.Amass, /*weight*/1.0, chain, event);
            track_computed = true;
        }
        sample.reco_beta = static_cast<float>(track_beta);
        sample.passed |= (1 << 7);
        sample.species = nucleus;
        if (num_track) sample.passed |= (1 << 8);
    };

    // --- Fill tree 1 (den) ---
    set_sample(den, den_l1, 1, den.reco_inn, R_inn, num_l1, 2);
    set_sample(den, den_tof, 3, den.reco_il1, R_il1, num_tof, 4);
    set_sample(den, den_trig, 5, den.reco_il1, R_il1, isPhysicsTrigger, 6);
    ensure_track(den);
    set_sample(den, den_l1Unb, 9, den.reco_inn, R_inn, num_l1Unb, 10);
    set_sample(den, den_trackCh, 11, den.reco_il1, R_il1, num_trackCh, 12);

    // --- Fill tree 2 (den_ltof) ---
    set_sample(den_ltof, den_l1_ltof, 1, den_ltof.reco_inn, R_inn, num_l1, 2);
    set_sample(den_ltof, den_tof_ltof, 3, den_ltof.reco_il1, R_il1, num_tof, 4);
    set_sample(den_ltof, den_trig, 5, den_ltof.reco_il1, R_il1, isPhysicsTrigger, 6);
    ensure_track(den_ltof);
    set_sample(den_ltof, den_l1Unb_ltof, 9, den_ltof.reco_inn, R_inn, num_l1Unb, 10);
    set_sample(den_ltof, den_trackCh_ltof, 11, den_ltof.reco_il1, R_il1, num_trackCh, 12);

    // --- Fill tree 3 (den_ltof_l9Unb_05) ---
    set_sample(den_ltof_l9Unb_05, den_l1_ltof_l9_05, 1, den_ltof_l9Unb_05.reco_inn, R_inn, num_l1, 2);
    set_sample(den_ltof_l9Unb_05, den_tof_ltof_l9_05, 3, den_ltof_l9Unb_05.reco_il1, R_il1, num_tof, 4);
    set_sample(den_ltof_l9Unb_05, den_trig, 5, den_ltof_l9Unb_05.reco_il1, R_il1, isPhysicsTrigger, 6);
    ensure_track(den_ltof_l9Unb_05);
    set_sample(den_ltof_l9Unb_05, den_l1Unb_ltof_l9_05, 9, den_ltof_l9Unb_05.reco_inn, R_inn, num_l1Unb, 10);
    set_sample(den_ltof_l9Unb_05, den_trackCh_ltof_l9_05, 11, den_ltof_l9Unb_05.reco_il1, R_il1, num_trackCh, 12);

    // --- Fill tree 4 (den_ltof_l9Unb_1) ---
    set_sample(den_ltof_l9Unb_1, den_l1_ltof_l9_1, 1, den_ltof_l9Unb_1.reco_inn, R_inn, num_l1, 2);
    set_sample(den_ltof_l9Unb_1, den_tof_ltof_l9_1, 3, den_ltof_l9Unb_1.reco_il1, R_il1, num_tof, 4);
    set_sample(den_ltof_l9Unb_1, den_trig, 5, den_ltof_l9Unb_1.reco_il1, R_il1, isPhysicsTrigger, 6);
    ensure_track(den_ltof_l9Unb_1);
    set_sample(den_ltof_l9Unb_1, den_l1Unb_ltof_l9_1, 9, den_ltof_l9Unb_1.reco_inn, R_inn, num_l1Unb, 10);
    set_sample(den_ltof_l9Unb_1, den_trackCh_ltof_l9_1, 11, den_ltof_l9Unb_1.reco_il1, R_il1, num_trackCh, 12);

    // --- Fill TTrees if anything passed ---
    if (den.passed)              tree_map["tree1"]->Fill();
    if (den_ltof.passed)         tree_map["tree2"]->Fill();
    if (den_ltof_l9Unb_05.passed) tree_map["tree3"]->Fill();
    if (den_ltof_l9Unb_1.passed) tree_map["tree4"]->Fill();
}

void ChargeAnalysisManager::FillMcSampleFromFileInfo(NAIA::NAIAChain& chain, const std::string &geom) {
    auto mcchain = chain.GetFileInfoTree();
    if (!mcchain) {
        std::cerr << "[ChargeAnalysisManager] ERROR: No MCFileInfo tree found in chain!" << std::endl;
        return;
    }

    auto mcinfo = new NAIA::MCFileInfo();
    mcchain->SetBranchAddress("MCFileInfo", &mcinfo);

    // Loop over all entries of the FileInfo tree
    for (unsigned long long i = 0; i < mcchain->GetEntries(); ++i) {
        mcchain->GetEntry(i);

        unsigned int charge = static_cast<unsigned int>(mcinfo->Charge);
        if (std::find(fCharges.begin(), fCharges.end(), charge) == fCharges.end())
            continue;  // skip if not among studied charges

        TH1D* mc_samp = mcHist[charge][geom]["mc_samp"];
        if (!mc_samp) {
            std::cerr << "[ChargeAnalysisManager] ERROR: missing mc_sample hist for Z=" << charge << std::endl;
            continue;
        }

        double rminDC = mcinfo->GetRMin();
        double rmaxDC = mcinfo->GetRMax();
        double ngen   = mcinfo->GetNGen();

        // Add contributions to the sample histogram
        for (int ibin = 1; ibin <= mc_samp->GetNbinsX(); ++ibin) {
            double rmin = mc_samp->GetBinLowEdge(ibin);
            double rmax = mc_samp->GetBinLowEdge(ibin + 1);

            if (rmin < rminDC) rmin = rminDC;
            if (rmax > rmaxDC) rmax = rmaxDC;
            if (rmin > rmaxDC || rmax < rminDC) continue;

            double weight = ngen * TMath::Log(rmax / rmin) / TMath::Log(rmaxDC / rminDC);
            mc_samp->AddBinContent(ibin, weight);
        }
    }

    delete mcinfo;
}
void ChargeAnalysisManager::Write(const std::map<unsigned int, std::string>& ionPaths,
                                  const std::string& outname,
                                  bool isMC,
                                  unsigned int charge) {

        std::string filename = outname;

        TFile fout(filename.c_str(), "RECREATE");

        auto histIt = histos.find(charge);

        for (const auto& [geom, nameMap] : histIt->second) {
            if (nameMap.empty()) continue;
            fout.mkdir(geom.c_str());
            fout.cd(geom.c_str());
            for (const auto& [name, h] : nameMap) {
                if (h && h->GetEntries() > 0) h->Write();
            }
            fout.cd();
        }
        fout.Close();
}
void ChargeAnalysisManager::WriteDistributions(const std::map<unsigned int, std::string>& ionPaths,
                                               const std::string& outname,
                                               bool isMC) {
    // Map groups of distributions into subsystems
    std::map<std::string, std::vector<std::string>> subsystemMap = {
        {"L1", {"qL1_den", "qL1_den_ltof", "qL1_den_ltof_l9_05", "qL1_den_ltof_l9_1"}},
        {"L1Unb", {"qL1Unb_den", "qL1Unb_den_ltof", "qL1Unb_den_ltof_l9_05", "qL1Unb_den_ltof_l9_1"}},
        {"Tof", {"qtof_den", "qtof_den_ltof", "qtof_den_ltof_l9_05", "qtof_den_ltof_l9_1"}},
        {"TrackCharge", {"qInn_R_den", "qInn_R_den_ltof", "qInn_R_den_ltof_l9_05", "qInn_R_den_ltof_l9_1"}},
        {"L9", {"qUnbL9_R_den", "qUnbL9_R_den_ltof", "qUnbL9_R_den_ltof_l9_05", "qUnbL9_R_den_ltof_l9_1"}}
    };

    for (auto charge : fCharges) {
        std::string ionPath = ionPaths.at(charge);

        // open one file per subsystem
        std::map<std::string, TFile*> files;
        for (auto& [subsys, list] : subsystemMap) {
            std::string filename =
                "/storage/gpfs_ams/ams/users/aubaldi/IonsAnalysis/IonsSelected/" +
                ionPath + "/Distributions/" + subsys + "/" + outname;
            files[subsys] = new TFile(filename.c_str(), "RECREATE");
        }

        // loop histos and write to the right subsystem
        for (const auto& [geom, nameMap] : distributionHistos.at(charge)) {
            for (const auto& [name, h] : nameMap) {
                for (auto& [subsys, list] : subsystemMap) {
                    if (std::find(list.begin(), list.end(), name) != list.end()) {
                        TFile* f = files[subsys];

                        // controlla se la dir esiste già
                        TDirectory* dir = (TDirectory*) f->Get(geom.c_str());
                        if (!dir) {
                            dir = f->mkdir(geom.c_str());
                        }
                        dir->cd();
                        h->Write();
                        f->cd(); // torna alla root del file
                        break;
                    }
                }
            }
        }

        // chiudi i file
        for (auto& [subsys, f] : files) {
            f->Close();
            delete f;
        }
    }
}
void ChargeAnalysisManager::WriteMcTrees(const std::map<unsigned int, std::string>& ionPaths,
                                               const std::string& outname,
                                               bool isMC) {
    for (auto charge : fCharges) {
        std::string ionPath = ionPaths.at(charge);
        std::string filename =
            "/storage/gpfs_ams/ams/users/aubaldi/IonsAnalysis/IonsSelected/" +
            ionPath + "/Passed/" + outname;

        TFile fout(filename.c_str(), "RECREATE");

        for (const auto& [geom, treeMap] : mcTrees.at(charge)) {
            fout.mkdir(geom.c_str());
            fout.cd(geom.c_str());

            for (const auto& [tname, tptr] : treeMap) {
                if (tptr) tptr->Write(tptr->GetName());
            }
            fout.cd();
        }
        for (const auto& [geom, nameMap] : mcHist.at(charge)) {
            fout.cd(geom.c_str());
            for (const auto& [name, h] : nameMap) {
                h->Write();
            }
            fout.cd();
        }
        fout.Close();
    }
}
// Helper: safe access to TemplateSet, returns nullptr if missing
static TemplateSet* _getTemplateSetSafe(std::map<unsigned int, std::map<std::string, TemplateSet>> &mapTS,
                                        unsigned int charge, const std::string &geom) {
    auto it_charge = mapTS.find(charge);
    if (it_charge == mapTS.end()) return nullptr;
    auto it_geom = it_charge->second.find(geom);
    if (it_geom == it_charge->second.end()) return nullptr;
    return &(it_geom->second);
}
void ChargeAnalysisManager::FillTemplateL1(unsigned int charge, const std::string &geom, int ibin, float q, double w) {
    TemplateSet* ts = _getTemplateSetSafe(templateSets, charge, geom);
    if (!ts) return;
    if (ibin < 0 || ibin >= static_cast<int>(ts->L1Templates.size())) return;

    // fill histogram
    TH1D* h = ts->L1Templates[ibin];
    if (h) h->Fill(q, w);

    // fill kernel tree: set storage element and Fill the tree
    ts->L1charge[ibin] = q;
    ts->L1weight[ibin] = w;
    TTree* kt = ts->L1KernelTrees[ibin];
    if (kt) kt->Fill();
}
void ChargeAnalysisManager::FillTemplateL2(unsigned int charge, const std::string &geom, int ibin, float q, double w) {
    TemplateSet* ts = _getTemplateSetSafe(templateSets, charge, geom);
    if (!ts) return;
    if (ibin < 0 || ibin >= static_cast<int>(ts->L2Templates.size())) return;

    TH1D* h = ts->L2Templates[ibin];
    if (h) h->Fill(q, w);

    ts->L2charge[ibin] = q;
    ts->L2weight[ibin] = w;
    TTree* kt = ts->L2KernelTrees[ibin];
    if (kt) kt->Fill();
}
void ChargeAnalysisManager::FillTemplateL1Distribution(unsigned int charge, const std::string &geom, int ibin, float q, double w) {
    TemplateSet* ts = _getTemplateSetSafe(templateSets, charge, geom);
    if (!ts) return;
    if (ibin < 0 || ibin >= static_cast<int>(ts->L1Distributions.size())) return;

    TH1D* h = ts->L1Distributions[ibin];
    if (h) h->Fill(q, w);

    ts->L1Dcharge[ibin] = q;
    ts->L1Dweight[ibin] = w;
    TTree* kt = ts->L1DistKernelTrees[ibin];
    if (kt) kt->Fill();
}
void ChargeAnalysisManager::WriteTemplateDistributions(const std::map<unsigned int, std::string>& ionPaths,
                                                      const std::string& outname,
                                                      bool isMC) {
        for (auto charge : fCharges) {
        auto it = templateSets.find(charge);
        if (it == templateSets.end()) continue;

        const auto &ionPath = ionPaths.at(charge);
        std::string filename =
            "/storage/gpfs_ams/ams/users/aubaldi/IonsAnalysis/Fragmentation/BelowL1/"
            + ionPath + "/" + outname;

        TFile fout(filename.c_str(), "RECREATE", "", ROOT::CompressionSettings(ROOT::kLZ4, 9));
        if (fout.IsZombie()) {
            std::cerr << "Cannot create file " << filename << std::endl;
            continue;
        }

        const auto &geomMap = it->second;
        for (const auto &[geom, ts] : geomMap) {
            TDirectory *geomDir = fout.GetDirectory(geom.c_str());
            if (!geomDir) geomDir = fout.mkdir(geom.c_str());
            geomDir->cd();

            auto writeGroup = [&](const std::vector<TH1D*>& hists,
                                  const std::vector<TTree*>& trees,
                                  const char *dirName) {
                TDirectory *dir = geomDir->GetDirectory(dirName);
                if (!dir) dir = geomDir->mkdir(dirName);
                dir->cd();
                const size_t n = hists.size();
                for (size_t i = 0; i < n; ++i) {
                    if (hists[i]) hists[i]->Write("", TObject::kOverwrite);
                    if (trees[i]) trees[i]->Write("", TObject::kOverwrite);
                }
                geomDir->cd();
            };

            writeGroup(ts.L1Templates, ts.L1KernelTrees, Form("L1TemplateList_%i", charge));
            writeGroup(ts.L2Templates, ts.L2KernelTrees, Form("L2TemplateList_%i", charge));
            writeGroup(ts.L1Distributions, ts.L1DistKernelTrees, Form("L1Distribution_%i", charge));
        }

        fout.Write("", TObject::kOverwrite);
        fout.Close();
    }
}
TH2D* ChargeAnalysisManager::GetHist(unsigned int charge,
                                     const std::string& geom,
                                     const std::string& name) {
    if (histos.count(charge) &&
        histos.at(charge).count(geom) &&
        histos.at(charge).at(geom).count(name)) 
    {
        return histos.at(charge).at(geom).at(name);
    }
    return nullptr;
}
double FillHistos(TH2D &histo, unsigned int charge, double Amass, double weight, NAIA::NAIAChain &chain,  NAIA::Event &event) {
  bool IsMC = chain.IsMC();
  float IGRFCutoff;
  if (!IsMC) {
    const NAIA::RTIInfo &rti_info = chain.GetEventRTIInfo();
    IGRFCutoff = rti_info.MaxIGRFCutoff[1][1];
  } else
    IGRFCutoff = event.mcTruthBase->Primary.GetGenMomentum() / event.mcTruthBase->Primary.Z;
  auto utime = event.header->UTCTime;
  float Beta = event.tofBaseSt->Beta[NAIA::Tof::BetaType::BetaH];
  const float amu = 0.938272075;
  double mass = Amass * amu;
  double BetaRig = mass * Beta / sqrt(1 - Beta * Beta) / charge;

  float EnergyD = 0;
  if (NAIA::ContainsKeys(event.ecalBase->Energy, NAIA::Ecal::EnergyRecoType::EnergyD)) {
    EnergyD = event.ecalBase->Energy[NAIA::Ecal::EnergyRecoType::EnergyD];
  }
  if (BetaRig <= 5) {
    histo.Fill(utime, BetaRig, weight);
	  return BetaRig;
  } else if (IGRFCutoff > 5 && IGRFCutoff <= 30) {
    histo.Fill(utime, IGRFCutoff, weight);
	  return IGRFCutoff;
  } else if (EnergyD > 30 && Efficiency::TrTrackEffSel::IsInsideEcal()(event)) {
    histo.Fill(utime, EnergyD, weight);
	  return EnergyD;
  }
  return EnergyD;
}
