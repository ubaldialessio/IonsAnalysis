void plot_unfolded_flux() {
    TFile* file = TFile::Open("Phosphorus_unfolded.root");
    TTree* tree = (TTree*)file->Get("unfolded_flux");

    double x, y, y_err;
    std::string* graph_name = nullptr;
    std::string* graph_title = nullptr;
    std::string* x_title = nullptr;
    std::string* y_title = nullptr;

    tree->SetBranchAddress("x", &x);
    tree->SetBranchAddress("y", &y);
    tree->SetBranchAddress("y_err", &y_err);
    tree->SetBranchAddress("graph_name", &graph_name);
    tree->SetBranchAddress("graph_title", &graph_title);
    tree->SetBranchAddress("x_title", &x_title);
    tree->SetBranchAddress("y_title", &y_title);

    Long64_t nentries = tree->GetEntries();

    std::vector<double> vx, vy, vex, vey;
    for (Long64_t i = 0; i < nentries; ++i) {
        tree->GetEntry(i);
        vx.push_back(x);
        vy.push_back(y);
        vex.push_back(0);     
        vey.push_back(y_err);
    }

    TGraphErrors* gr = new TGraphErrors(nentries, vx.data(), vy.data(), vex.data(), vey.data());


    TCanvas* c1 = new TCanvas("c1", "Unfolded Flux", 800, 600);
    gr->SetMarkerStyle(20);
    gr->Draw("AP");
    c1->SaveAs("unfolded_flux.pdf");

    file->Close();
}
