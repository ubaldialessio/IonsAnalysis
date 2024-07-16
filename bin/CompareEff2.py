# -*- coding: utf-8 -*-
import ROOT
import matplotlib.pyplot as plt
import numpy as np

def getIonName(charge):
    ion_names = {
        1: "Pr",
        2: "He",
        3: "Li",
        4: "Be",
        5: "B",
        6: "C",
        7: "N",
        8: "O",
        9: "F",
        10: "Ne",
        11: "Na",
        12: "Mg",
        13: "Al",
        14: "Si",
        15: "P",
        16: "S"
    }
    return ion_names.get(charge, "")

def read_efficiencies(filename):
    data_efficiencies = {}
    file = ROOT.TFile(filename)
    if not file or file.IsZombie():
        print("Error opening file: {}".format(filename))
        return None

    for key in file.GetListOfKeys():
        hist = key.ReadObj()
        if isinstance(hist, ROOT.TH1D) and ("dataEff" in hist.GetName() or "mcEff" in hist.GetName()):
            data_efficiencies[hist.GetName()] = hist

    file.Close()
    return data_efficiencies

def plot_comparison(charge_list, option):
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']

    for charge in charge_list:
        ion_name = getIonName(charge)
        if not ion_name:
            print("Charge {} not recognized.".format(charge))
            continue

        filename = "/storage/gpfs_ams/ams/users/aubaldi/IonsAnalysis/bin/IonsSelected/{}/{}_output.root".format(ion_name, ion_name)
        data_efficiencies = read_efficiencies(filename)

        if not data_efficiencies or len(data_efficiencies) != 8:
            print("Error reading efficiencies from file: {}".format(filename))
            continue

        # Extract data and mc efficiencies
        data_eff_l1 = data_efficiencies["dataEff_l1"]
        data_eff_tf = data_efficiencies["dataEff_tf"]
        data_eff_tr = data_efficiencies["dataEff_tr"]
        data_eff_tk = data_efficiencies["dataEff_tk"]

        mc_eff_l1 = data_efficiencies["mcEff_l1"]
        mc_eff_tf = data_efficiencies["mcEff_tf"]
        mc_eff_tr = data_efficiencies["mcEff_tr"]
        mc_eff_tk = data_efficiencies["mcEff_tk"]

        data_efficiencies = [data_eff_l1, data_eff_tf, data_eff_tr, data_eff_tk]
        mc_efficiencies = [mc_eff_l1, mc_eff_tf, mc_eff_tr, mc_eff_tk]

        # Plotting with matplotlib
        fig, axs = plt.subplots(2, 2, figsize=(10, 8))
        fig.suptitle('Efficiency Comparison - Ion {}'.format(ion_name), fontsize=16)

        for i in range(4):
            row = i // 2
            col = i % 2
            data_values = [data_efficiencies[i].GetBinContent(j) for j in range(1, data_efficiencies[i].GetNbinsX() + 1)]
            mc_values = [mc_efficiencies[i].GetBinContent(j) for j in range(1, mc_efficiencies[i].GetNbinsX() + 1)]
            x_labels = ['Bin {}'.format(j) for j in range(1, len(data_values) + 1)]
            axs[row, col].bar(x_labels, data_values, color=colors[i], label='Data', alpha=0.7)
            axs[row, col].bar(x_labels, mc_values, color='red', label='MC', alpha=0.5)
            axs[row, col].set_title('Efficiency Type {}'.format(i+1))
            axs[row, col].set_ylabel('Efficiency')
            axs[row, col].legend()
            axs[row, col].set_ylim(0, 1)

        plt.tight_layout()
        plt.savefig('Comparison_{}.pdf'.format(ion_name))
        plt.close()

if __name__ == "__main__":
    import sys
    if len(sys.argv) < 3:
        print("Usage: python CompareEff.py <charge1> <charge2> ... <charge n> <data/mc/ratio>")
        sys.exit(1)

    charges = [int(charge) for charge in sys.argv[1:-1]]
    option = sys.argv[-1]

    plot_comparison(charges, option)
