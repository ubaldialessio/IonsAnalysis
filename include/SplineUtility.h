#ifndef SPLINEUTILITY_H
#define SPLINEUTILITY_H
#include "includes.h"
#include "binning.h"
#include "TFitResult.h"

//////spfit///////
#define Nmax 1000
int NODES = 0, NPOINTS = 0;

double SpFuncHermite(double *x, double *par) {
  int N = NODES; // definito globalmente come prima
  double xx = x[0];

  // knots X fissi in par[0..N-1]
  const double* X = par;

  // Y nodali in par[N..2N-1]
  const double* Y = par + N;

  // slopes M in par[2N..3N-1]
  const double* M = par + 2*N;

  // trova il segmento in cui cade xx
  int seg = -1;
  for (int i = 0; i < N-1; ++i) {
    if (xx >= X[i] && xx <= X[i+1]) { seg = i; break; }
  }
  if (seg < 0) {
    // fuori range → extrapolo linearmente
    if (xx < X[0]) return Y[0] + M[0]*(xx - X[0]);
    if (xx > X[N-1]) return Y[N-1] + M[N-1]*(xx - X[N-1]);
  }

  double x0 = X[seg], x1 = X[seg+1];
  double y0 = Y[seg], y1 = Y[seg+1];
  double m0 = M[seg], m1 = M[seg+1];

  double h = x1 - x0;
  double t = (xx - x0) / h;

  // polinomi di Hermite
  double h00 = (1 + 2*t) * pow(1 - t, 2);
  double h10 = t * pow(1 - t, 2);
  double h01 = t*t * (3 - 2*t);
  double h11 = t*t * (t - 1);

  return h00*y0 + h10*h*m0 + h01*y1 + h11*h*m1;
}
void Spline(int n, double *x, double *y, double *b, double *d){
  double u[Nmax];
  d[0] = u[0] = 0;
  if (b) {
	  d[0] = -0.5;
    u[0] = (3/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-b[0]);
  }
  for (int i = 1; i < n-1; i++) {
	  double s = (x[i]-x[i-1])/(x[i+1]-x[i-1]);
    double p = s*d[i-1]+2;
    d[i] = (s-1)/p;
    u[i] = (y[i+1]-y[i])/(x[i+1]-x[i])-(y[i]-y[i-1])/(x[i]-x[i-1]);
    u[i] = (6*u[i]/(x[i+1]-x[i-1])-s*u[i-1])/p;
  }
  double q = 0, v = 0;
  if (b) {
	  q = 0.5;
    v = (3/(x[n-1]-x[n-2]))*(b[1]-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
  }
  d[n-1] = (v-q*u[n-2])/(q*d[n-2]+1);
  for (int i = n-2; i >= 0; i--) d[i] = d[i]*d[i+1]+u[i];
}
double Splint(int n, double *xn, double *yn, double *d, double x){
  int i1 = 0, i2 = n-1;
  while (i2-i1 > 1) {
	int i = (i2+i1)/2;
    if (xn[i] > x) i2 = i;
    else i1 = i;
  }
  double h = xn[i2]-xn[i1];
  if (h == 0) return 0;
  double a  = (xn[i2]-x)/h;
  double b  = (x-xn[i1])/h;
  double a3 = a*a*a;
  double b3 = b*b*b;
  double c  = 0;
  return a*yn[i1]+b*yn[i2]+((a3-a)*d[i1]+(b3-b)*d[i2])*(h*h)/6+c;
}
double SpFunc(double *xp, double *par){
  // Spline function
  // par[0:n-1]   : Nodes x
  // par[n:n*2-1] : Nodes y
  // par[n*2]     : Derivative (dy/dx) of lower boundary
  // par[n*2+1]   : Derivative (dy/dx) of upper boundary
  static double ps[Nmax*2+2] = { 0, 0 };
  static double xn[Nmax];
  static double d2[Nmax];
  int chk = 1;
  for (int i = 0; i < NODES*2+2; i++) if (ps[i] != par[i]) { chk = 0; break; }
  if (!chk) {
	  for (int i = 0; i < NODES*2+2; i++) ps[i] = par[i];
    for (int i = 0; i < NODES; i++) xn[i] = (par[i] > 0) ? TMath::Log10(par[i]) : par[i];
    Spline(NODES, xn, &par[NODES], &par[NODES*2], d2);
  }
  double x = xp[0];
  if (x > 0) x = TMath::Log10(x);
  double y = Splint(NODES, xn, &par[NODES], d2, x);
  if (x < xn[0])       { double dx = x-xn[0];       y = par[NODES]+par[NODES*2]*dx; }
  //if (x > xn[NODES-1]) { double dx = x-xn[NODES-1]; y = par[NODES*2-1]+par[NODES*2+1]*dx;} ORIGINALE
  if (x > xn[NODES-1]) {
    double dx = x - xn[NODES-1];
    y = par[NODES + NODES - 1] + par[NODES * 2 + 1] * dx;
  }
  return y;
}
TH1D *autospline(TH1D *histo, double xmin, double xmax, int llim = 4, int hlim = 10, double xsplit = -1){ //3,7 originale

  int nbins = histo->GetNbinsX()+1;
  double bins[nbins];
  double *ci = new double[nRbins_HighZ], *hx = new double[nRbins_HighZ];
  for(int i=0; i<nbins; ++i) bins[i]=TMath::Log10(histo->GetBinLowEdge(i+1)+1);
  auto histo_log = new TH1D("", "", nbins-1, bins);
  for(int i=1; i<=histo_log->GetNbinsX(); ++i){
    histo_log->SetBinContent(i, histo->GetBinContent(i));
    histo_log->SetBinError(i, histo->GetBinError(i));
  }
  double log_min = TMath::Log10(xmin+1);
  double log_max = TMath::Log10(xmax+1);

  TH1D *final_histo;
  auto final_histo_temp = (TH1D*)hist_rig_highZ->Clone();
  for(int i=0; i<final_histo_temp->GetNbinsX(); i++)
    hx[i]=TMath::Log10(final_histo_temp->GetBinCenter(i+1)+1);
  double bic = 100000, chi2;
  TString title;

  // Special case: constant fit
  if (llim == 1 && hlim == 1) {
    auto constant = new TF1("constant", "[0]", log_min, log_max);
    constant->SetParameter(0, histo_log->GetMean());
    histo_log->Fit(constant, "Q0", "", log_min, log_max);
    chi2 = constant->GetChisquare();
    auto ndf= constant->GetNumberFitPoints()-1;

    const double val = constant->GetParameter(0);
    const double err = constant->GetParError(0);

    final_histo = (TH1D*)final_histo_temp->Clone();
    for (int i = 0; i < final_histo->GetNbinsX(); ++i) {
      final_histo->SetBinContent(i + 1, val);
      final_histo->SetBinError(i + 1, err);
    }

    title = Form("1 knot (constant), Chi2/ndf: %.2f, BIC: N/A", chi2/ndf);
    delete constant;
    final_histo->SetTitle(title);
    final_histo->SetName(title);
    std::cout << "Histogram " << histo->GetName() << " : " << title << std::endl;
    return final_histo;
  }
  
  // ----------------- Hybrid case: constant up to xsplit, spline above -----------------
  if (xsplit > 0) {
      double log_split = TMath::Log10(xsplit+1);

      // --- Fit constant in [xmin, xsplit] ---
      auto constant = new TF1("constant", "[0]", log_min, log_split);
      constant->SetParameter(0, histo_log->GetMean());
      histo_log->Fit(constant, "Q0", "", log_min, log_split);

      const double val = constant->GetParameter(0);
      const double err = constant->GetParError(0);

      // --- Fit spline in [xsplit, xmax] ---
      double best_bic = 1e9;
      TH1D *best_histo = nullptr;
      TString best_title;

      for (NODES = llim; NODES < hlim; NODES++) {
        double knots[NODES];
        for (int i = 0; i < NODES; i++)
            knots[i] = log_split + i * (log_max - log_split) / (NODES - 1);

        auto spline = new TF1("", SpFunc, log_split, log_max, NODES*2+2);

        // fissa i nodi in X
        for (int i = 0; i < NODES; i++) spline->FixParameter(i, knots[i]);

        // fissa il primo nodo in Y = valore costante (continuità di valore)
        spline->FixParameter(NODES, val);

        // Fix lower boundary derivative to 0 so slope matches constant region (C1 cont.)
        // par[2*NODES] is lower derivative in SpFunc convention
        spline->FixParameter(2*NODES, 0);

        // gli altri nodi in Y vengono inizializzati dai contenuti dell’istogramma
        for (int i = 1; i < NODES; i++)
            spline->SetParameter(NODES+i, histo_log->GetBinContent(histo_log->FindBin(knots[i])));

        // keep upper derivative fixed as before (optional)
        spline->FixParameter(NODES*2+1, 0);

        histo_log->Fit(spline, "Q0", "", log_split, log_max);
        NPOINTS = spline->GetNumberFitPoints();
        (TVirtualFitter::GetFitter())->GetConfidenceIntervals(nRbins_HighZ, 1, hx, ci, 0.68);

        auto hybrid = (TH1D*)final_histo_temp->Clone();
        for (int i = 0; i < hybrid->GetNbinsX(); i++) {
            if (hx[i] < log_split) {
                hybrid->SetBinContent(i+1, val);
                hybrid->SetBinError(i+1, err);
            } else {
                hybrid->SetBinContent(i+1, spline->Eval(hx[i]));
                hybrid->SetBinError(i+1, ci[i]);
            }
        }

        // number of free parameters in hybrid:
        // - Y nodes: NODES, but first one is fixed -> NODES-1 free
        // - we fixed both boundary derivatives -> 0 free derivatives
        int k = NODES - 1;
        double chi2 = spline->GetChisquare();
        double bic = chi2 + k * log((double)NPOINTS); // standard BIC

        if (bic < best_bic) {
            best_bic = bic;
            best_histo = (TH1D*)hybrid->Clone();
            best_title = Form("Hybrid (const+%d knots, C1), Chi2/ndf: %f, BIC: %f.",
                              NODES, chi2/ (double)(NPOINTS - k), bic);
        }
        delete spline;
      }

      delete constant;
      if (!best_histo) {
        std::cerr << "autospline: hybrid fit failed for " << histo->GetName() << std::endl;
        return nullptr;
      }
      best_histo->SetTitle(best_title);
      best_histo->SetName(best_title);
      std::cout << "Histogram " << histo->GetName() << " : " << best_title << std::endl;
      return best_histo;
  }


  for(NODES=llim; NODES<hlim; NODES++){
    double knots[NODES];
    for (int i = 0; i < NODES; i++) knots[i] = log_min + i * (log_max-log_min)/(NODES-1);
    auto spline = new TF1("", SpFunc, log_min, log_max, NODES*2+2);
    spline->SetNpx(10000);
    for (int i = 0; i < NODES; i++){
      if (0 & (1<<i)) spline->SetParameter(i, knots[i]); //nodes x
      else            spline->FixParameter(i, knots[i]);
      spline->SetParameter(i+NODES, histo_log->GetBinContent(histo_log->FindBin(knots[i]))); //nodes y
    }
    spline->FixParameter(NODES*2+1, 0);
    spline->FixParameter(2*NODES, 0);
    histo_log->Fit(spline, "Q0", "", log_min, log_max);
    NPOINTS = spline->GetNumberFitPoints();
    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(nRbins_HighZ, 1, hx, ci, 0.68);
    for(int i = 0; i < final_histo_temp->GetNbinsX(); i++){
	  final_histo_temp->SetBinContent(i+1, spline->Eval(hx[i]));
      final_histo_temp->SetBinError(i+1, ci[i]);
    }
    int ndf = (NPOINTS-NODES);
    chi2 = spline->GetChisquare();
    double new_bic = (chi2+NODES*log(NPOINTS))/ndf;
    if(new_bic>0 && new_bic<bic){
	  bic = new_bic;
      final_histo = (TH1D*)final_histo_temp->Clone();
      title = Form("%i knots, Chi2/ndf: %f, BIC: %f.", NODES, chi2/ndf, bic);
    }
    delete spline;
  }

  final_histo->SetTitle(title);
  final_histo->SetName(title);
  std::cout << "Histogram " << histo->GetName() << " : " << title << std::endl;
  return final_histo;
}
TF1 *spfit(TH1D *h, int nodes, double xmin, double xmax){
  NODES = nodes;
  auto func = new TF1("func", SpFunc, xmin, xmax, NODES*2+2);
  for (int i = 0; i < NODES; i++){
    func->FixParameter(i, xmin + i * (xmax-xmin)/(NODES-1)); //nodes x
    func->SetParameter(i+NODES, h->Interpolate(func->GetParameter(i))); //nodes y
  }
  func->FixParameter(NODES*2+1, 0);
  h->Fit(func, "Q0", "", xmin, xmax);
  return h->GetFunction("func");
}
TF1* spfit(TGraphErrors* h,
           int nodes,
           double xmin,
           double xmax,
           const std::vector<double>& knots,
           const TString& knotOpt,
           const std::vector<double>* x_points = nullptr,
           std::vector<double>* y_errors = nullptr) {
    // Decide number of nodes
    if (knotOpt == "lin" || knotOpt == "log")
        NODES = nodes;
    else
        NODES = static_cast<int>(knots.size());

    auto func = new TF1("func", SpFunc, xmin, 3000.0, NODES * 2 + 2);

    // Linear spacing of knots
    if (knotOpt == "lin") {
        for (int i = 0; i < NODES; ++i) {
            double x_knot = xmin + i * (xmax - xmin) / (NODES - 1);
            func->FixParameter(i, x_knot);                     // node x
            func->SetParameter(i + NODES, h->Eval(x_knot));    // node y
        }
    }
    // Logarithmic spacing of knots
    else if (knotOpt == "log") {
        double log_xmin = TMath::Log10(xmin);
        double log_xmax = TMath::Log10(xmax);
        for (int i = 0; i < NODES; ++i) {
            double log_knot = log_xmin + i * (log_xmax - log_xmin) / (NODES - 1);
            double x_knot = TMath::Power(10.0, log_knot);
            func->FixParameter(i, x_knot);
            func->SetParameter(i + NODES, h->Eval(x_knot));
        }
    }
    // Custom knot positions
    else {
        for (std::size_t i = 0; i < knots.size(); ++i) {
            double x_knot = knots[i];
            func->FixParameter(static_cast<int>(i), x_knot);
            func->SetParameter(static_cast<int>(i + knots.size()), h->Eval(x_knot));
        }
    }

    // Optional: fix last slope to 0
    // func->FixParameter(NODES * 2 + 1, 0);

    // Perform fit
    h->Fit(func, "E", "", xmin, xmax);

    // Compute 68% confidence intervals if requested
    if (x_points && y_errors) {
        int n = static_cast<int>(x_points->size());
        std::vector<double> x_arr(n), y_errs(n);

        for (int i = 0; i < n; ++i) {
            x_arr[i] = (*x_points)[i];
        }

        TVirtualFitter::GetFitter()->GetConfidenceIntervals(
            n, 1, x_arr.data(), y_errs.data(), 0.68);

        y_errors->resize(n);
        for (int i = 0; i < n; ++i) {
            (*y_errors)[i] = y_errs[i];
        }
    }

    // Update fit title with chi2/ndf
    double chi2 = func->GetChisquare();
    int ndf = func->GetNDF();
    TString title = Form("#chi^{2}/ndf = %.2f / %d = %.2f", chi2, ndf, (ndf > 0 ? chi2 / ndf : 0.0));

    func->SetName("purityFit");
    func->SetTitle(title);

    return func;
}


//spfit but TH1D*
TH1D* spfitHist(TGraphErrors* h, int nodes, double xmin, double xmax,
            std::vector<double> knots, TString knotOpt,
            const std::vector<double>* x_points = nullptr,
            std::vector<double>* y_errors = nullptr,
            const TH1D* templateHist = nullptr)  // optional: to copy binning
 {
  if (knotOpt=="lin" || knotOpt=="log")
    NODES = nodes;
  else 
    NODES = knots.size();

  auto func = new TF1("func", SpFunc, xmin, xmax, NODES * 2 + 2);

  // --- Fix knot positions and init Y ---
  if (knotOpt=="lin") {
    for (int i = 0; i < NODES; i++) {
      double x_knot = xmin + i * (xmax - xmin) / (NODES - 1);
      func->FixParameter(i, x_knot); 
      func->SetParameter(i + NODES, h->Eval(x_knot)); 
    }
  }
  else if (knotOpt=="log") {
    for (int i = 0; i < NODES; i++) {
      double log_xmin = TMath::Log10(xmin);
      double log_xmax = TMath::Log10(xmax);
      double log_knot = log_xmin + i * (log_xmax - log_xmin) / (NODES - 1);
      double x_knot = TMath::Power(10, log_knot);
      func->FixParameter(i, x_knot);
      func->SetParameter(i + NODES, h->Eval(x_knot));
    }
  }
  else {
    for (int i = 0; i < (int)knots.size(); ++i) {
      double x_knot = knots[i];
      func->FixParameter(i, x_knot);
      func->SetParameter(i + knots.size(), h->Eval(x_knot));
    }
  }

  // Perform fit
  h->Fit(func, "Q0", "", xmin, xmax);

  // ---- Build output histogram ----
  TH1D* out;
  if (templateHist) {
    out = (TH1D*)templateHist->Clone("spfitHist");
    out->Reset();
  } else {
    int nbins = 100; // default binning if no template provided
    out = new TH1D("spfitHist","", nbins, xmin, xmax);
  }

  // Fill contents with function values
  for (int i = 1; i <= out->GetNbinsX(); ++i) {
    double x = out->GetBinCenter(i);
    double val = func->Eval(x);
    out->SetBinContent(i, val);
  }

  // ---- Confidence intervals if requested ----
  if (x_points && y_errors) {
    int n = x_points->size();
    std::vector<Double_t> x_arr(n), y_errs(n);
    for (int i = 0; i < n; ++i)
      x_arr[i] = (*x_points)[i];

    TVirtualFitter::GetFitter()->GetConfidenceIntervals(n, 1, x_arr.data(), y_errs.data(), 0.68);

    y_errors->resize(n);
    for (int i = 0; i < n; ++i)
      (*y_errors)[i] = y_errs[i];
  }

  // ---- Decorate ----
  double chi2 = func->GetChisquare();
  int ndf = func->GetNDF();
  TString title = Form("Spline fit #chi^{2}/ndf = %.2f/%d = %.2f", chi2, ndf, chi2/ndf);
  out->SetTitle(title);
  out->SetName("spfitHist");

  return out;
}

TF1* spfit(TH1D* h, int nodes, double xmin, double xmax,
                std::vector<double> knots, TString knotOpt,
                const std::vector<double>* x_points = nullptr,
                std::vector<double>* y_errors = nullptr) {
  if (!h || h->GetNbinsX() == 0) return nullptr;

  // Convert TH1D to TGraphErrors
  auto g = new TGraphErrors();
  for (int i = 1; i <= h->GetNbinsX(); ++i) {
    double x = h->GetBinCenter(i);
    if (x < xmin || x > xmax) continue;
    double y = h->GetBinContent(i);
    double ey = h->GetBinError(i);
    int p = g->GetN();
    g->SetPoint(p, x, y);
    g->SetPointError(p, 0, ey); // no x-error
  }

  // Call the original spfit with the TGraphErrors
  TF1* result = spfit(g, nodes, xmin, xmax, knots, knotOpt, x_points, y_errors);
  result->SetName("purityFit"); // in case you want to save it later

  delete g;
  return result;
}
TF1* autosplineTF1(TH1D* histo, double xmin, double xmax, int llim = 3, int hlim = 10) {
  int nbins = histo->GetNbinsX() + 1;
  double bins[nbins];
  for (int i = 0; i < nbins; ++i)
    bins[i] = TMath::Log10(histo->GetBinLowEdge(i + 1) + 1);

  auto histo_log = new TH1D("", "", nbins - 1, bins);
  for (int i = 1; i <= histo_log->GetNbinsX(); ++i) {
    histo_log->SetBinContent(i, histo->GetBinContent(i));
    histo_log->SetBinError(i, histo->GetBinError(i));
  }

  double log_min = TMath::Log10(xmin + 1);
  double log_max = TMath::Log10(xmax + 1);

  TF1* best_func = nullptr;
  double best_bic = 1e10;
  TString best_title;

  for (int nodes = llim; nodes <= hlim; ++nodes) {
    NODES = nodes;
    double knots[nodes];
    for (int i = 0; i < nodes; ++i)
      knots[i] = log_min + i * (log_max - log_min) / (nodes - 1);

    auto func = new TF1("spline", SpFunc, log_min, log_max, nodes * 2 + 2);

    for (int i = 0; i < nodes; ++i) {
      func->FixParameter(i, knots[i]); // Fix x positions
      func->SetParameter(i + nodes, histo_log->Interpolate(knots[i])); // Set y guesses
    }

    func->FixParameter(nodes * 2 + 1, 0); // force last param = 0

    histo_log->Fit(func, "Q0", "", log_min, log_max);
    int npoints = func->GetNumberFitPoints();
    int ndf = npoints - nodes;
    double chi2 = func->GetChisquare();
    double bic = (chi2 + nodes * log(npoints)) / ndf;

    if (bic > 0 && bic < best_bic) {
      best_bic = bic;
      best_title = Form("%d knots: #chi^{2}/ndf = %.2f, BIC = %.2f", nodes, chi2 / ndf, bic);
      if (best_func) delete best_func;
      best_func = (TF1*)func->Clone("best_spline");
    }

    delete func;
  }

  delete histo_log;

  if (!best_func) {
    std::cerr << "Spline fit failed for all NODES in range." << std::endl;
    return nullptr;
  }

  best_func->SetTitle(best_title);
  best_func->SetName(best_title);
  std::cout << "Best spline fit: " << best_title << std::endl;
  return best_func;
}
TF1* spfit_autoBIC(TGraphErrors* g, double xmin, double xmax,
                   int llim , int hlim , TString knotOpt,
                   double xsplit = -1,
                   const std::vector<double>* x_points = nullptr,
                   std::vector<double>* y_errors = nullptr) {
  if (!g || xmin >= xmax || llim < 1 || hlim < llim) {
    std::cerr << "[spfit_autoBIC] invalid inputs" << std::endl;
    return nullptr;
  }

  // -------- LOG SPACE -----------
  const double log_xmin = TMath::Log10(xmin+1);
  const double log_xmax = TMath::Log10(xmax+1);

  // punti x per valutare CI
  int npoints = g->GetN();
  std::vector<double> gx(npoints), gy(npoints), gey(npoints);
  for (int i=0; i<npoints; ++i) {
    double xx, yy;
    g->GetPoint(i, xx, yy);
    gx[i] = TMath::Log10(xx+1);
    gy[i] = yy;
    gey[i] = g->GetErrorY(i);
  }

  TF1* bestFit = nullptr;
  double bestBIC = 1e20;
  TString bestTitle;

  // -------------------
  // Caso speciale: costante
  // -------------------
  if (llim == 1 && hlim == 1) {
    auto constant = new TF1("constant", "[0]", log_xmin, log_xmax);
    constant->SetParameter(0, TMath::Mean(npoints, gy.data()));
    g->Fit(constant, "Q0", "", xmin, xmax);
    double chi2 = constant->GetChisquare();
    int ndf   = constant->GetNumberFitPoints()-1;

    bestTitle = Form("Constant fit, #chi^{2}/ndf = %.2f/%d = %.2f",
                     chi2, ndf, chi2/ndf);
    constant->SetTitle(bestTitle);
    std::cout << "[spfit_autoBIC] Selected: " << bestTitle << std::endl;
    return constant;
  }

  // -------------------
  // Caso ibrido: costante sotto xsplit, spline sopra
  // -------------------
  if (xsplit > 0) {
    double log_split = TMath::Log10(xsplit+1);

    // Fit costante su [xmin, xsplit]
    auto constant = new TF1("constant_low", "[0]", log_xmin, log_split);
    constant->SetParameter(0, TMath::Mean(npoints, gy.data()));
    g->Fit(constant, "Q0", "", xmin, xsplit);

    const double val = constant->GetParameter(0);
    const double err = constant->GetParError(0);

    double best_bic = 1e20;
    TF1* bestSpline = nullptr;
    TString best_title;

    for (NODES = llim; NODES <= hlim; ++NODES) {
      std::vector<double> knots(NODES);
      for (int i = 0; i < NODES; ++i)
        knots[i] = log_split + i*(log_xmax - log_split)/(NODES-1);

      auto spline = new TF1("spline", SpFunc, log_split, log_xmax, NODES*2+2);

      // fissa nodi in X
      for (int i = 0; i < NODES; i++)
        spline->FixParameter(i, knots[i]);

      // primo nodo Y = valore costante
      spline->FixParameter(NODES, val);

      // gli altri Y inizializzati da grafico
      for (int i = 1; i < NODES; i++) {
        double yy = g->Eval(TMath::Power(10,knots[i])-1);
        spline->SetParameter(NODES+i, yy);
      }

      // slope finale fissato a 0
      spline->FixParameter(NODES*2+1, 0);

      g->Fit(spline, "Q0", "", xsplit, xmax);

      int NPOINTS = spline->GetNumberFitPoints();
      int ndf = (NPOINTS - (NODES-1));
      double chi2 = spline->GetChisquare();
      double bic = (chi2 + (NODES-1)*std::log(NPOINTS)) / ndf;

      if (bic < best_bic) {
        best_bic = bic;
        delete bestSpline;
        bestSpline = (TF1*)spline->Clone("bestHybrid");
        best_title = Form("Hybrid (const+%d knots), Chi2/ndf: %.2f/%d, BIC=%.2f",
                          NODES, chi2, ndf, bic);
      }
      delete spline;
    }

    delete constant;
    if (!bestSpline) return nullptr;

    bestSpline->SetTitle(best_title);
    std::cout << "[spfit_autoBIC] Selected hybrid: " << best_title << std::endl;
    return bestSpline;
  }

  // -------------------
  // Caso generale: spline pura
  // -------------------
  for (NODES = llim; NODES <= hlim; ++NODES) {
    std::vector<double> knots(NODES);
    for (int i = 0; i < NODES; i++)
      knots[i] = log_xmin + i*(log_xmax - log_xmin)/(NODES-1);

    auto spline = new TF1("spline", SpFunc, log_xmin, log_xmax, NODES*2+2);
    spline->SetNpx(2000);

    // fissa X dei nodi
    for (int i = 0; i < NODES; ++i)
      spline->FixParameter(i, knots[i]);

    // Y inizializzati da grafico
    for (int i = 0; i < NODES; ++i) {
      double yy = g->Eval(TMath::Power(10,knots[i])-1);
      spline->SetParameter(NODES+i, yy);
    }

    // slope finale = 0
    spline->FixParameter(NODES*2+1, 0);

    g->Fit(spline, "Q0", "", xmin, xmax);
    int NPOINTS = spline->GetNumberFitPoints();
    int ndf = (NPOINTS - NODES);
    double chi2 = spline->GetChisquare();
    double bic  = (chi2 + NODES*std::log(NPOINTS)) / ndf;

    if (bic > 0 && bic < bestBIC) {
      bestBIC = bic;
      delete bestFit;
      bestFit = (TF1*)spline->Clone("bestSpline");
      bestTitle = Form("%d knots, Chi2/ndf: %.2f/%d, BIC=%.2f",
                       NODES, chi2, ndf, bic);
    }

    delete spline;
  }

  if (!bestFit) {
    std::cerr << "[spfit_autoBIC] No valid fits found! Fallback constant." << std::endl;
    auto constant = new TF1("constant_fallback", "[0]", xmin, xmax);
    constant->SetParameter(0, g->GetMean(2));
    g->Fit(constant, "Q0", "", xmin, xmax);
    return constant;
  }

  bestFit->SetTitle(bestTitle);
  std::cout << "[spfit_autoBIC] Selected best: " << bestTitle << std::endl;
  return bestFit;
}

namespace SplineUtility {
    enum Efficiency : int { TofEff, TrackEff, TriggerEff, L1Eff, L1UnbEff, TrackChEff, Acc, L9Eff, DaqEff, L1ChargeCut};
    std::pair<double,double> SetFitLimits(int Eff, unsigned int Charge);
    std::pair<int,int> SetKnots(int Eff, unsigned int Charge);
}
std::pair<double, double> SplineUtility::SetFitLimits(int Eff, unsigned int Charge) {
    double xmin = 2.15;
    double xmax = 1000; // reasonable default

    switch (static_cast<SplineUtility::Efficiency>(Eff)) {
        case TofEff:
            // All charges use default xmax = 1000
            break;
        case TrackEff:
            switch (Charge) {
                case 14: xmax = 30; break;
                case 15:
                case 16: xmax = 30; break;
                default: xmax = 30; break;
            }
            break;
        case TriggerEff:
            switch (Charge) {
                case 14: xmax = 60; break;
                case 16: xmax = 90; break;
                case 26: xmax = 1000; break;
                default: xmax = 1000; break;
            }
            break;
        case L1Eff:
            switch (Charge) {
                case 14: xmax = 120; break; 
                case 15:
                case 16: 
                case 18:
                case 20:
                default: xmax = 100; break;
            }
            break;
        case L1UnbEff:
            xmax = 100;
            break;
        case TrackChEff:
            switch (Charge) {
                case 14: xmax = 100; break;
                case 16: xmax = 1000; break;
                default: xmax = 1000; break;
            }
            break;
        case L1ChargeCut:
            switch (Charge) {
                case 14: xmin = 0.8; xmax = 100; break;
                case 15: xmin = 0.8; xmax = 1000; break;
                case 16: xmin = 0.8; xmax = 1000; break;
                default: xmin = 0.8; xmax = 1000; break;
            }
        case Acc:
            //xmin = 2.15;
            xmax = 2000;
            break;
        case L9Eff:
            xmax = 500;
            break;
        case DaqEff:
            xmin = 0.8;
            xmax = 2000;
            break;
        default:
            // Optional: log warning or assert
            break;
    }
    return {xmin, xmax};
}
std::pair<int, int> SplineUtility::SetKnots(int Eff, unsigned int Charge) {
    int kmin = 2, kmax = 3; // default values
    switch (static_cast<SplineUtility::Efficiency>(Eff)) {
        case TofEff:
            switch (Charge) {
                case 14: kmin = 2; kmax = 3; break;
                case 16: kmin = 2; kmax = 3; break;
                default: kmin = 3; kmax = 5; break;
            }
            break;
        case TrackEff:
            switch (Charge) {
                case 14: kmin = 2; kmax = 5; break;
                case 16: kmin = 2; kmax = 5; break;
                default: kmin = 2; kmax = 5; break;
            }
            break;
        case TriggerEff:
            switch (Charge) {
                case 14: kmin = 2; kmax = 5; break;
                case 16: kmin = 2; kmax = 5; break;
                default: kmin = 2; kmax = 5; break;
            }
            break;
        case L1Eff:
            switch (Charge) {
                case 14: kmin = 2; kmax = 5; break;
                case 16: kmin = 2; kmax = 5; break;
                default: kmin = 2; kmax = 3; break;
            }
            break;
        case L1UnbEff:
              switch (Charge) {
                case 14: kmin = 2; kmax = 5; break;
                case 16: kmin = 2; kmax = 5; break;
                default: kmin = 2; kmax = 3; break;
            }
            break;
        case L1ChargeCut:
            switch (Charge) {
                case 14: kmin = 1; kmax = 1; break;
                case 15: kmin = 2; kmax = 3; break;
                case 16: kmin = 1; kmax = 1; break;
                default: kmin = 2; kmax = 3; break;
            }
            break;
        case TrackChEff:
            switch (Charge) {
                case 14: kmin = 4 ; kmax = 6; break;
                case 16: kmin = 1; kmax = 1; break;
                default: kmin = 3; kmax = 6; break;
            }
            break;
        case Acc:
            kmin = 2;
            kmax = 60;
            break;
        case L9Eff:
            kmin = 2;
            kmax = 4;
            break;
        case DaqEff:
            switch (Charge) {
                case 14: kmin = 3; kmax = 4; break;
                case 16: kmin = 4; kmax = 5; break;
                case 18: kmin = 3; kmax = 7; break;
                case 20: kmin = 3; kmax = 4; break;
                case 26: kmin = 2; kmax = 4; break;
                default: kmin = 3; kmax = 4; break;
            }
            break;
        default:
            // Optional: add logging or error
            break;
    }
    return {kmin, kmax};
}
#endif
