from MyPyLib_v4 import *

def plot_multifit( data_exp, save_as, datas_fit=None, labels_fit=None, E_of_C=None, C_of_E=None, plot_range_x=None, plot_range_y=None ):


    sample_format_dict_1 = {
        "label"      : "Messdaten",          
        "fmt"        : 'o', 
        "color"      : "grey",                               
        "markersize" : 1.0, 
        "linewidth"  : 1,
        "capsize"    : 0,
        "alpha"      : 0.8                                   
    }

    sample_format_dict_gauss = {
        "label"      : None,          
        "fmt"        : '-', 
        "color"      : None,                               
        "markersize" : 0.7, 
        "linewidth"  : 2,
        "capsize"    : 0,
        "alpha"      : 1                                   
    }

    writtings = {
        "title"       : None,
        "x_ax_label"  : r"Kanal-Index $K$ / 1",
        "y_ax_label"  : r"Ereignisse $N$ / 1"
    }
    
    general_format_dict = standard_format_dict.copy()
    zoom_params         = no_zooming.copy()
    colorbar_params     = no_colorbar.copy()
    extra_label         = no_extra_label.copy()

    general_format_dict["custom_x_range"] = plot_range_x
    general_format_dict["custom_y_range"] = plot_range_y

    if ((E_of_C == None) or (C_of_E == None)):
        extra_xaxis = no_extra_xaxis.copy()
    else:
        extra_xaxis         = {
            "do_axis"   :   True,
            "f_new(old)":   E_of_C,
            "f_old(new)":   C_of_E,
            "label"     :   r"Energie $E$ / keV"
        }
    
    all_data                = [ data_exp ]                               
    all_sample_format_dicts = [ sample_format_dict_1 ]

    if( isinstance(datas_fit, list) ): 
        for i in range(len(datas_fit)):
            colors_fit                          = sns.color_palette("magma", len(datas_fit))
            sample_format_dict_gauss["label"]   = labels_fit[i]
            sample_format_dict_gauss["color"]   = colors_fit[i]
            all_sample_format_dicts.append( sample_format_dict_gauss.copy() )
            all_data.append(datas_fit[i])

    save_plot = True, "../Figures/"+save_as+".jpg"                                     
    ultimate_plot_advanced (all_data, writtings, zoom_params, colorbar_params, extra_label, extra_xaxis, save_plot, all_sample_format_dicts, general_format_dict)

def fitfunc_gauss_bg( x, mu, FWHM, area, a, b ):
    sigma = FWHM / ( 2*np.sqrt(2*np.log(2)) )
    return fitfunc_gauss_area(x, mu, sigma, area) + fitfunc_linear(x, a, b)

def fitfunc_gauss_double_bg( x, mu1, FWHM1, area1, mu2, FWHM2, area2, a, b ):
    sigma1 = FWHM1 / ( 2*np.sqrt(2*np.log(2)) )
    sigma2 = FWHM2 / ( 2*np.sqrt(2*np.log(2)) )
    return fitfunc_gauss_area(x, mu1, sigma1, area1) + fitfunc_gauss_area(x, mu2, sigma2, area2) + fitfunc_linear(x, a, b)

def gauss_guess( dataset, ROI ):
    x_data      = dataset[0]
    y_data      = dataset[2]
    mask_ROI    = ((x_data < ROI[1])  &  (x_data > ROI[0]))
    x_data      = x_data[mask_ROI]
    y_data      = y_data[mask_ROI]

    amplitude   = max(y_data)
    mu          = np.sum( x_data*y_data ) / np.sum(y_data)
    sigma       = np.sqrt( np.abs( np.sum( x_data**2 *y_data ) / np.sum(y_data) - mu**2 ) ) 
    area        = amplitude * np.sqrt(2*np.pi)*sigma
    FWHM        = 2*np.sqrt(2*np.log(2)) * sigma

    return [mu, FWHM, area]



# --------------------------- Import Data + Raw Plots ---------------------------

data_names_HPGe = [
    "HPGe_BG_300s.txt",
    "HPGe_Co_300s_20mm.txt",
    "HPGe_Cs_300s_100mm.txt",
    "HPGe_Eu_300s_150mm.txt"
]
prefix = "HPGe"

channel     = np.loadtxt("../Data/"+data_names_HPGe[0])[:-1 , 0]
bg          = np.loadtxt("../Data/"+data_names_HPGe[0])[:-1 , 1]
counts_Co   = np.loadtxt("../Data/"+data_names_HPGe[1])[:-1 , 1]
counts_Cs   = np.loadtxt("../Data/"+data_names_HPGe[2])[:-1 , 1]
counts_Eu   = np.loadtxt("../Data/"+data_names_HPGe[3])[:-1 , 1]

counts_Co_clean, counts_Co_clean_err     = counts_Co - bg, np.sqrt(counts_Co+bg+1)
counts_Cs_clean, counts_Cs_clean_err     = counts_Cs - bg, np.sqrt(counts_Cs+bg+1)
counts_Eu_clean, counts_Eu_clean_err     = counts_Eu - bg, np.sqrt(counts_Eu+bg+1)

data_Co     = [ channel, None, counts_Co_clean, counts_Co_clean_err ]
data_Cs     = [ channel, None, counts_Cs_clean, counts_Cs_clean_err ]
data_Eu     = [ channel, None, counts_Eu_clean, counts_Eu_clean_err ]

data_bg_noerr   = [ channel, None, bg, None ]
data_Co_noerr   = [ channel, None, counts_Co - bg, None ]
data_Cs_noerr   = [ channel, None, counts_Cs - bg, None ]
data_Eu_noerr   = [ channel, None, counts_Eu - bg, None ]

plot_multifit(data_Co_noerr, prefix+"_Co")
plot_multifit(data_bg_noerr, prefix+"_bg")
#plot_multifit(data_Cs_noerr, prefix+"_Cs")
#plot_multifit(data_Eu_noerr, prefix+"_Eu")


# --------------------------- Fit Cobalt Data ---------------------------


#plot_multifit(data_Co_noerr, prefix+"_Co_1173", plot_range_x=[11850, 11950])
#plot_multifit(data_Co_noerr, prefix+"_Co_1332", plot_range_x=[13460, 13560])

range_Co_1173 = [11850, 11950]
fitbook_Co_1173 = {
    "data_set"          : data_Co,
    "region_of_interest": range_Co_1173,
    "fit_function"      : fitfunc_gauss_bg,
    "params_guess"      : [*gauss_guess(data_Co, range_Co_1173), 0, 0],
    "region_of_fit"     : None,
    "fit_density"       : 400
}
range_Co_1332 = [13460, 13560]
fitbook_Co_1332 = {
    "data_set"          : data_Co,
    "region_of_interest": range_Co_1332,
    "fit_function"      : fitfunc_gauss_bg,
    "params_guess"      : [*gauss_guess(data_Co, range_Co_1332), 0, 0],
    "region_of_fit"     : None,
    "fit_density"       : 400
}

fitdata_Co_1, params_Co_1, params_Co_1_err, chi_Co_1 = ultimate_fit( fitbook_Co_1173 )
fitdata_Co_2, params_Co_2, params_Co_2_err, chi_Co_2 = ultimate_fit( fitbook_Co_1332 )

fitdata_Co     = [fitdata_Co_1, fitdata_Co_2]
fitlabels_Co   = [ 
    f"1173keV-Peak $(\\chi^2 = {chi_Co_1:.1f})$",
    f"1332keV-Peak $(\\chi^2 = {chi_Co_2:.1f})$"
]

#plot_multifit(data_Co_noerr, prefix+"_Co", fitdata_Co, fitlabels_Co )
#plot_multifit(data_Co_noerr, prefix+"_Co_1173", fitdata_Co, fitlabels_Co, plot_range_x=range_Co_1173 )
#plot_multifit(data_Co_noerr, prefix+"_Co_1332", fitdata_Co, fitlabels_Co, plot_range_x=range_Co_1332 )



# --------------------------- Fit Caesium Data ---------------------------


#plot_multifit(data_Cs_noerr, prefix+"_Cs_661", plot_range_x=[6660, 6760])

range_Cs_661 = [6660, 6760]
fitbook_Cs_661 = {
    "data_set"          : data_Cs,
    "region_of_interest": range_Cs_661,
    "fit_function"      : fitfunc_gauss_bg,
    "params_guess"      : [*gauss_guess(data_Cs, range_Cs_661), 0, 0],
    "region_of_fit"     : None,
    "fit_density"       : 400
}

fitdata_Cs, params_Cs, params_Cs_err, chi_Cs = ultimate_fit( fitbook_Cs_661 )
fitlabel_Cs = [ 
    f"661keV-Peak $(\\chi^2 = {chi_Cs:.1f})$",
]

#plot_multifit(data_Cs_noerr, prefix+"_Cs", [fitdata_Cs], fitlabel_Cs )
#plot_multifit(data_Cs_noerr, prefix+"_Cs_661", [fitdata_Cs], fitlabel_Cs, plot_range_x=range_Cs_661 )


# --------------------------- First Energy Kalibration ---------------------------


known_energies      = np.array([ 661.7, 1173.2, 1332.5 ]) 
peak_positions      = np.array([params_Cs, params_Co_1, params_Co_2 ])[: , 0]
peak_positions_err  = np.array([params_Cs_err, params_Co_1_err, params_Co_2_err ])[: , 0]
dataset_energy      = [ known_energies, None, peak_positions, peak_positions_err ]

fitbook_energy = { 
    "data_set"          : dataset_energy,
    "region_of_interest": None,
    "fit_function"      : fitfunc_linear,
    "params_guess"      : None,
    "region_of_fit"     : None,
    "fit_density"       : 400
}
fitdata_energy, params_energy, params_energy_err, chi_energy = ultimate_fit(fitbook_energy)


def plot_energie_kalib(exp_data, fit_data, chi, id):
    
    sample_format_dict_1 = {
        "label"      : r"Messdaten",          
        "fmt"        : 'o', 
        "color"      : sns.color_palette("dark")[0],                               
        "markersize" : 4, 
        "linewidth"  : 1,
        "capsize"    : 0,
        "alpha"      : 1                                   
    }
    sample_format_dict_2 = {
        "label"      : f"Linearer Fit $(\\chi^2 = {chi:.1f})$",                       
        "fmt"        : '--', 
        "color"      : sns.color_palette("dark")[9],        
        "markersize" : 4, 
        "linewidth"  : 1,
        "capsize"    : 0,
        "alpha"      : 1
    }  

    writtings = {
        "title"       : None,
        "y_ax_label"  : r"Kanal-Index $K$ / 1",
        "x_ax_label"  : r"Energie $E$ / keV"
    }
    
    general_format_dict = standard_format_dict.copy()
    zoom_params         = no_zooming.copy()
    colorbar_params     = no_colorbar.copy()
    extra_label         = no_extra_label.copy()
    extra_xaxis         = no_extra_xaxis.copy()
 
    
    all_data                = [ exp_data, fit_data ]                               
    all_sample_format_dicts = [ sample_format_dict_1, sample_format_dict_2 ]

    save_plot = True, "../Figures/"+prefix+"_Kalibration_E_"+id+".jpg"                                     
    ultimate_plot_advanced (all_data, writtings, zoom_params, colorbar_params, extra_label, extra_xaxis, save_plot, all_sample_format_dicts, general_format_dict)

plot_energie_kalib( dataset_energy, fitdata_energy, chi_energy, "1" )

a, a_err = params_energy[0], params_energy_err[0]
b, b_err = params_energy[1], params_energy_err[1]

print(prefix+" Energie-Kalibration 1: K = a*E + b")
print(f"a = ({a:.4f} +- {a_err:.4f}) / keV")
print(f"b = ({b:.2f} +- {b_err:.2f}) / 1")

def C_of_E( E ):
    return a*E + b
def E_of_C( C ):
    return (C-b) / a


# --------------------------- Fit Europium Data ---------------------------


#plot_multifit(data_Eu_noerr, prefix+"_Eu_5000", E_of_C=E_of_C, C_of_E=C_of_E, plot_range_x=[0, 5000])
#plot_multifit(data_Eu_noerr, prefix+"_Eu_10000", E_of_C=E_of_C, C_of_E=C_of_E, plot_range_x=[7500, 10000], plot_range_y=[-50,500])
#plot_multifit(data_Eu_noerr, prefix+"_Eu_15000", E_of_C=E_of_C, C_of_E=C_of_E, plot_range_x=[10000, 15000], plot_range_y=[-50,500])


range_Eu_121    = [1200 , 1250 ]
fitbook_Eu_121  = {
    "data_set"          : data_Eu,
    "region_of_interest": range_Eu_121,
    "fit_function"      : fitfunc_gauss_bg,
    "params_guess"      : [*gauss_guess(data_Eu, range_Eu_121), 0, 0],
    "region_of_fit"     : None,
    "fit_density"       : 400
}
range_Eu_244    = [2450 , 2550]
fitbook_Eu_244  = {
    "data_set"          : data_Eu,
    "region_of_interest": range_Eu_244,
    "fit_function"      : fitfunc_gauss_bg,
    "params_guess"      : [*gauss_guess(data_Eu, range_Eu_244), 0, 0],
    "region_of_fit"     : None,
    "fit_density"       : 400
}
range_Eu_344    = [3450 , 3550]
fitbook_Eu_344  = {
    "data_set"          : data_Eu,
    "region_of_interest": range_Eu_344,
    "fit_function"      : fitfunc_gauss_bg,
    "params_guess"      : [*gauss_guess(data_Eu, range_Eu_344), 0, 0],
    "region_of_fit"     : None,
    "fit_density"       : 400
}
range_Eu_411    = [4120 , 4200]
fitbook_Eu_411  = {
    "data_set"          : data_Eu,
    "region_of_interest": range_Eu_411,
    "fit_function"      : fitfunc_gauss_bg,
    "params_guess"      : [*gauss_guess(data_Eu, range_Eu_411), 0, 0],
    "region_of_fit"     : None,
    "fit_density"       : 400
}
range_Eu_443    = [4450 , 4550]
fitbook_Eu_443  = {
    "data_set"          : data_Eu,
    "region_of_interest": range_Eu_443,
    "fit_function"      : fitfunc_gauss_bg,
    "params_guess"      : [*gauss_guess(data_Eu, range_Eu_443), 0, 0],
    "region_of_fit"     : None,
    "fit_density"       : 400
}
range_Eu_778    = [7850 , 7950]
fitbook_Eu_778  = {
    "data_set"          : data_Eu,
    "region_of_interest": range_Eu_778,
    "fit_function"      : fitfunc_gauss_bg,
    "params_guess"      : [*gauss_guess(data_Eu, range_Eu_778), 0, 0],
    "region_of_fit"     : None,
    "fit_density"       : 400
}
range_Eu_867    = [8700 , 8900]
fitbook_Eu_867  = {
    "data_set"          : data_Eu,
    "region_of_interest": range_Eu_867,
    "fit_function"      : fitfunc_gauss_bg,
    "params_guess"      : [*gauss_guess(data_Eu, range_Eu_867), 0, 0],
    "region_of_fit"     : None,
    "fit_density"       : 400
}
range_Eu_964    = [9700 , 9850]
fitbook_Eu_964  = {
    "data_set"          : data_Eu,
    "region_of_interest": range_Eu_964,
    "fit_function"      : fitfunc_gauss_bg,
    "params_guess"      : [*gauss_guess(data_Eu, range_Eu_964), 0, 0],
    "region_of_fit"     : None,
    "fit_density"       : 400
}
range_Eu_1085   = [10960 , 11050]
fitbook_Eu_1085 = {
    "data_set"          : data_Eu,
    "region_of_interest": range_Eu_1085,
    "fit_function"      : fitfunc_gauss_bg,
    "params_guess"      : [*gauss_guess(data_Eu, range_Eu_1085), 0, 0],
    "region_of_fit"     : None,
    "fit_density"       : 400
}
range_Eu_1112   = [11230 , 11320]
fitbook_Eu_1112 = {
    "data_set"          : data_Eu,
    "region_of_interest": range_Eu_1112,
    "fit_function"      : fitfunc_gauss_bg,
    "params_guess"      : [*gauss_guess(data_Eu, range_Eu_1112), 0, 0],
    "region_of_fit"     : None,
    "fit_density"       : 400
}
range_Eu_1408   = [14200 , 14400]
fitbook_Eu_1408 = {
    "data_set"          : data_Eu,
    "region_of_interest": range_Eu_1408,
    "fit_function"      : fitfunc_gauss_bg,
    "params_guess"      : [*gauss_guess(data_Eu, range_Eu_1408), 0, 0],
    "region_of_fit"     : None,
    "fit_density"       : 400
}


fitdata_Eu_121, params_Eu_121, params_Eu_121_err, chi_Eu_121 = ultimate_fit( fitbook_Eu_121 )
fitdata_Eu_244, params_Eu_244, params_Eu_244_err, chi_Eu_244 = ultimate_fit( fitbook_Eu_244 )
fitdata_Eu_344, params_Eu_344, params_Eu_344_err, chi_Eu_344 = ultimate_fit( fitbook_Eu_344 )
fitdata_Eu_411, params_Eu_411, params_Eu_411_err, chi_Eu_411 = ultimate_fit( fitbook_Eu_411 )
fitdata_Eu_443, params_Eu_443, params_Eu_443_err, chi_Eu_443 = ultimate_fit( fitbook_Eu_443 )
fitdata_Eu_778, params_Eu_778, params_Eu_778_err, chi_Eu_778 = ultimate_fit( fitbook_Eu_778 )
fitdata_Eu_867, params_Eu_867, params_Eu_867_err, chi_Eu_867 = ultimate_fit( fitbook_Eu_867 )
fitdata_Eu_964, params_Eu_964, params_Eu_964_err, chi_Eu_964 = ultimate_fit( fitbook_Eu_964 )
fitdata_Eu_1085, params_Eu_1085, params_Eu_1085_err, chi_Eu_1085 = ultimate_fit( fitbook_Eu_1085 )
fitdata_Eu_1112, params_Eu_1112, params_Eu_1112_err, chi_Eu_1112 = ultimate_fit( fitbook_Eu_1112 )
fitdata_Eu_1408, params_Eu_1408, params_Eu_1408_err, chi_Eu_1408 = ultimate_fit( fitbook_Eu_1408 )

fitdata_Eu_5000     = [fitdata_Eu_121, fitdata_Eu_244, fitdata_Eu_344, fitdata_Eu_411, fitdata_Eu_443]
fitdata_Eu_10000    = [fitdata_Eu_778, fitdata_Eu_867, fitdata_Eu_964] 
fitdata_Eu_15000    = [fitdata_Eu_1085, fitdata_Eu_1112, fitdata_Eu_1408]
fitdata_Eu          = fitdata_Eu_5000 + fitdata_Eu_10000 + fitdata_Eu_15000

fitlabels_Eu_5000   = [ 
    f"121keV-Peak $(\\chi^2 = {chi_Eu_121:.1f})$",
    f"244keV-Peak $(\\chi^2 = {chi_Eu_244:.1f})$",
    f"344keV-Peak $(\\chi^2 = {chi_Eu_344:.1f})$",
    f"411keV-Peak $(\\chi^2 = {chi_Eu_411:.1f})$",
    f"443keV-Peak $(\\chi^2 = {chi_Eu_443:.1f})$"
]
fitlabels_Eu_10000  = [ 
    f"778keV-Peak $(\\chi^2 = {chi_Eu_778:.1f})$",
    f"867keV-Peak $(\\chi^2 = {chi_Eu_867:.1f})$",
    f"964keV-Peak $(\\chi^2 = {chi_Eu_964:.1f})$"
]
fitlabels_Eu_15000  = [ 
    f"1085keV-Peak $(\\chi^2 = {chi_Eu_1085:.1f})$",
    f"1112keV-Peak $(\\chi^2 = {chi_Eu_1112:.1f})$",
    f"1408keV-Peak $(\\chi^2 = {chi_Eu_1408:.1f})$"
]
fitlabels_Eu        = [
    f"121keV-Peak",
    f"244keV-Peak",
    f"344keV-Peak",
    f"411keV-Peak",
    f"443keV-Peak",
    f"778keV-Peak",
    f"867keV-Peak",
    f"964keV-Peak",
    f"1085keV-Peak",
    f"1112keV-Peak",
    f"1408keV-Peak" 
]


#plot_multifit(data_Eu_noerr, prefix+"_Eu", fitdata_Eu, fitlabels_Eu, E_of_C=E_of_C, C_of_E=C_of_E)
#plot_multifit(data_Eu_noerr, prefix+"_Eu_5000", fitdata_Eu_5000, fitlabels_Eu_5000, E_of_C=E_of_C, C_of_E=C_of_E, plot_range_x=[0, 5000])
#plot_multifit(data_Eu_noerr, prefix+"_Eu_10000", fitdata_Eu_10000, fitlabels_Eu_10000, E_of_C=E_of_C, C_of_E=C_of_E, plot_range_x=[7500, 10000], plot_range_y=[-50,500])
#plot_multifit(data_Eu_noerr, prefix+"_Eu_15000", fitdata_Eu_15000, fitlabels_Eu_15000, E_of_C=E_of_C, C_of_E=C_of_E, plot_range_x=[10000, 15000], plot_range_y=[-50,500])


# --------------------------- Second Energy Kalibration ---------------------------


params_Eu       = [params_Eu_121, params_Eu_244, params_Eu_344, params_Eu_411, params_Eu_443, params_Eu_778, params_Eu_867, params_Eu_964, params_Eu_1085, params_Eu_1112, params_Eu_1408]
params_Eu_err   = [params_Eu_121_err, params_Eu_244_err, params_Eu_344_err, params_Eu_411_err, params_Eu_443_err, params_Eu_778_err, params_Eu_867_err, params_Eu_964_err, params_Eu_1085_err, params_Eu_1112_err, params_Eu_1408_err]
energies_Eu     = [121.7817, 244.6974, 344.2785, 411.1165, 443.9606, 778.9045, 867.380, 964.057, 1085.837, 1112.076, 1408.013]
energies_Eu_err = [0.0003,   0.0008,   0.0012,   0.0012,   0.0016,   0.0024,   0.003,   0.005,   0.010,    0.003,    0.003]

known_energies      = np.array([ 661.7, 1173.2, 1332.5, *energies_Eu ]) 
peak_positions      = np.array([params_Cs, params_Co_1, params_Co_2, *params_Eu ])[: , 0]
peak_positions_err  = np.array([params_Cs_err, params_Co_1_err, params_Co_2_err, *params_Eu_err ])[: , 0]
dataset_energy      = [ known_energies, None, peak_positions, peak_positions_err ]

fitbook_energy = { 
    "data_set"          : dataset_energy,
    "region_of_interest": None,
    "fit_function"      : fitfunc_linear,
    "params_guess"      : None,
    "region_of_fit"     : None,
    "fit_density"       : 2
}
fitdata_energy, params_energy, params_energy_err, chi_energy = ultimate_fit(fitbook_energy)

plot_energie_kalib( dataset_energy, fitdata_energy, chi_energy, "2" )

a, a_err = params_energy[0], params_energy_err[0]
b, b_err = params_energy[1], params_energy_err[1]

print(prefix+" Energie-Kalibration 2: K = a*E + b")
print(f"a = ({a:.6f} +- {a_err:.6f}) / keV")
print(f"b = ({b:.4f} +- {b_err:.4f}) / 1")

def C_of_E( E ):
    return a*E + b
def E_of_C( C ):
    return (C-b) / a

def E_of_C_err( C, C_err ):
    E       = (C-b) / a
    E_err   = np.sqrt( (C_err/a)**2 + (b_err/a)**2 + (E*a_err/a)**2 )
    return E, E_err

# --------------------------- Final Plots ---------------------------


plot_multifit(data_Co_noerr, prefix+"_Co", fitdata_Co, fitlabels_Co, E_of_C=E_of_C, C_of_E=C_of_E)
plot_multifit(data_Co_noerr, prefix+"_Co_1173", fitdata_Co, fitlabels_Co, E_of_C=E_of_C, C_of_E=C_of_E, plot_range_x=range_Co_1173 )
plot_multifit(data_Co_noerr, prefix+"_Co_1332", fitdata_Co, fitlabels_Co, E_of_C=E_of_C, C_of_E=C_of_E, plot_range_x=range_Co_1332 )

plot_multifit(data_Cs_noerr, prefix+"_Cs", [fitdata_Cs], fitlabel_Cs, E_of_C=E_of_C, C_of_E=C_of_E )
plot_multifit(data_Cs_noerr, prefix+"_Cs_661", [fitdata_Cs], fitlabel_Cs, E_of_C=E_of_C, C_of_E=C_of_E, plot_range_x=range_Cs_661 )

plot_multifit(data_Eu_noerr, prefix+"_Eu", fitdata_Eu, fitlabels_Eu, E_of_C=E_of_C, C_of_E=C_of_E)
plot_multifit(data_Eu_noerr, prefix+"_Eu_5000", fitdata_Eu_5000, fitlabels_Eu_5000, E_of_C=E_of_C, C_of_E=C_of_E, plot_range_x=[0, 5000])
plot_multifit(data_Eu_noerr, prefix+"_Eu_10000", fitdata_Eu_10000, fitlabels_Eu_10000, E_of_C=E_of_C, C_of_E=C_of_E, plot_range_x=[7500, 10000], plot_range_y=[-50,500])
plot_multifit(data_Eu_noerr, prefix+"_Eu_15000", fitdata_Eu_15000, fitlabels_Eu_15000, E_of_C=E_of_C, C_of_E=C_of_E, plot_range_x=[10000, 15000], plot_range_y=[-50,500])


# --------------------------- HPGe Intrinsische Halbwertsbreite  ---------------------------


FWHM_Eu         = np.array(params_Eu)[: , 1]
FWHM_Eu_err     = np.array(params_Eu_err)[: , 1]
E_FWHM_Eu, E_FWHM_Eu_err    = E_of_C_err( FWHM_Eu, FWHM_Eu_err )

mu_Eu           = np.array(params_Eu)[: , 0]
mu_Eu_err       = np.array(params_Eu_err)[: , 0]
E_gamma_Eu, E_gamma_Eu_err  = E_of_C_err( mu_Eu, mu_Eu_err )

data_intrins = [E_gamma_Eu, E_gamma_Eu_err, E_FWHM_Eu, E_FWHM_Eu_err]

def fitfunc_intrins( E_gamma, alpha, E_e ):
    return np.sqrt( alpha**2 * E_gamma + E_e**2 )

fitbook_intrins = {
    "data_set"          : data_intrins,
    "region_of_interest": None,
    "fit_function"      : fitfunc_intrins,
    "params_guess"      : None,
    "region_of_fit"     : None,
    "fit_density"       : 400
}

fitdata_intrins, params_intrins, params_intrins_err, chi_intrins = ultimate_fit( fitbook_intrins )

print(prefix+" Intrinsische Breite: FWHM = (a^2 *E + b^2 )^(1/2)")
print(f"a = ({params_intrins[0]:.4f} +- {params_intrins_err[0]:.4f}) / 1")
print(f"b = ({params_intrins[1]:.4f} +- {params_intrins_err[1]:.4f}) * keV")

def plot_FWHM_intrins(exp_data, fit_data, chi):
    
    sample_format_dict_1 = {
        "label"      : r"Messdaten",          
        "fmt"        : 'o', 
        "color"      : sns.color_palette("dark")[0],                               
        "markersize" : 4, 
        "linewidth"  : 1,
        "capsize"    : 2,
        "alpha"      : 1                                   
    }
    sample_format_dict_2 = {
        "label"      : f"Fitfunktion $(\\chi^2 = {chi:.1f})$",                       
        "fmt"        : '--', 
        "color"      : sns.color_palette("dark")[9],        
        "markersize" : 4, 
        "linewidth"  : 1,
        "capsize"    : 0,
        "alpha"      : 1
    }  

    writtings = {
        "title"       : None,
        "y_ax_label"  : r"Halbwertsbreite $\Delta E$ / keV",
        "x_ax_label"  : r"Energie des Photo-Peaks $E_\gamma$ / keV"
    }
    
    general_format_dict = standard_format_dict.copy()
    zoom_params         = no_zooming.copy()
    colorbar_params     = no_colorbar.copy()
    extra_label         = no_extra_label.copy()
    extra_xaxis         = no_extra_xaxis.copy()
 
    
    all_data                = [ exp_data, fit_data ]                               
    all_sample_format_dicts = [ sample_format_dict_1, sample_format_dict_2 ]

    save_plot = True, "../Figures/"+prefix+"_FWHM_intrins.jpg"                                     
    ultimate_plot_advanced (all_data, writtings, zoom_params, colorbar_params, extra_label, extra_xaxis, save_plot, all_sample_format_dicts, general_format_dict)

plot_FWHM_intrins(data_intrins, fitdata_intrins, chi_intrins)


# --------------------------- Peak-To-Total ---------------------------


    # korrekt, dass wir hier das bg-Spektrum abziehen?
area_Cs, area_Cs_err    = params_Cs[2], params_Cs_err[2]
total_Cs, total_Cs_err  = np.sum(counts_Cs_clean), np.sqrt(np.sum(counts_Cs_clean))     # hier illegale Fehlerrechnung!!
area_Co, area_Co_err    = params_Co_1[2]+params_Co_2[2], np.sqrt(params_Co_1_err[2]**2 + params_Co_2_err[2]**2)
total_Co, total_Co_err  = np.sum(counts_Co_clean), np.sqrt(np.sum(counts_Co_clean))

PeakTotal_Cs            = area_Cs/total_Cs
PeakTotal_Cs_err        = PeakTotal_Cs * np.sqrt( (area_Cs_err/area_Cs)**2 + (total_Cs_err/total_Cs)**2 )

PeakTotal_Co            = area_Co/total_Co
PeakTotal_Co_err        = PeakTotal_Co * np.sqrt( (area_Co_err/area_Co)**2 + (total_Co_err/total_Co)**2 )

print(prefix + f" PeakTotal Cs: ({PeakTotal_Cs*100:.2f} +- {PeakTotal_Cs_err*100:.2f})"+r" %")
print(prefix + f" PeakTotal Co: ({PeakTotal_Co*100:.2f} +- {PeakTotal_Co_err*100:.2f})"+r" %")


# --------------------------- Detektor-Effizienz ---------------------------


dist_Cs, dist_Cs_err    = 100, 5        # in mm
dist_Eu, dist_Eu_err    = 150, 5        # in mm
radius_HPGe             = 55.7 / 2      # in mm

coverage_Cs             = (np.pi * radius_HPGe**2) / (4*np.pi * dist_Cs**2)
coverage_Cs_err         = coverage_Cs * 2*dist_Cs_err/dist_Cs

coverage_Eu             = (np.pi * radius_HPGe**2) / (4*np.pi * dist_Eu**2)
coverage_Eu_err         = coverage_Eu * 2*dist_Eu_err/dist_Eu

print(prefix + f" Coverage Cs:  ({coverage_Cs*100:.2f} +- {coverage_Cs_err*100:.2f}) " + r"%")
print(prefix + f" Coverage Eu:  ({coverage_Eu*100:.2f} +- {coverage_Eu_err*100:.2f}) " + r"%")

sources_age             = (4*365 + 31 + 30 + 3)/365     # in Jahren 

halftime_Cs             = 30.007                        # in Jahren (Nundat)
halftime_Eu             = 13.517                        # in Jahren (Nundat)

init_activ_Cs           = 405 *10**3                    # in Bq
init_activ_Eu           = 709 *10**3                    # in Bq

acitv_Cs                = init_activ_Cs * 2**(-sources_age/halftime_Cs)
acitv_Eu                = init_activ_Eu * 2**(-sources_age/halftime_Eu)

T                       = 300                           # in s

exp_counts_Cs           = acitv_Cs * coverage_Cs 
exp_counts_Cs_err       = acitv_Cs * coverage_Cs_err
meas_counts_Cs          = area_Cs / T
meas_counts_Cs_err      = area_Cs_err / T

rel_itens_Eu            = np.array([28.53, 7.55, 26.59, 2.237, 2.827, 12.93, 4.23, 14.51, 10.11, 13.67, 20.87]) /100
rel_itens_Eu_err        = np.array([0.16,  0.04, 0.20,  0.013, 0.014, 0.08,  0.03, 0.07,  0.05,  0.08,  0.09]) /100
exp_counts_Eu           = acitv_Eu * coverage_Eu * rel_itens_Eu
exp_counts_Eu_err       = exp_counts_Eu * np.sqrt( (coverage_Eu_err/coverage_Eu)**2 + (rel_itens_Eu_err/rel_itens_Eu)**2 )
meas_counts_Eu          = np.array(params_Eu)[: , 2] / T
meas_counts_Eu_err      = np.array(params_Eu_err)[: , 2] / T

abs_eff_Cs              = meas_counts_Cs/exp_counts_Cs
abs_eff_Cs_err          = abs_eff_Cs * np.sqrt( (meas_counts_Cs_err/meas_counts_Cs)**2 + (exp_counts_Cs_err/exp_counts_Cs)**2 )

abs_eff_Eu              = meas_counts_Eu/exp_counts_Eu
abs_eff_Eu_err          = abs_eff_Eu * np.sqrt( (meas_counts_Eu_err/meas_counts_Eu)**2 + (exp_counts_Eu_err/exp_counts_Eu)**2 )

print(prefix + f" Abs. Eff. Cs 661keV-Peak: ({abs_eff_Cs*100:.2f} +- {abs_eff_Cs_err*100:.2f})" + r"%")

data_abs_eff                = [E_gamma_Eu, E_gamma_Eu_err, abs_eff_Eu*100, abs_eff_Eu_err*100]
E_gamma_Cs, E_gamma_Cs_err  = E_of_C_err(params_Cs[0], params_Cs_err[0])
data_abs_eff_Cs             = [ np.array([E_gamma_Cs]), np.array([E_gamma_Cs_err]), np.array([abs_eff_Cs*100]), np.array([abs_eff_Cs_err*100]) ]


    # Leo ----> https://www.sciencedirect.com/science/article/pii/016890028990394X?via=ihub S. 299
def fitfunc_eff_3( E, a_1, a_2):
    x = np.log(1.022/E)
    return  np.exp(a_1*x + a_2*x**2)

fitbook_abs_eff = {
    "data_set"          : data_abs_eff,
    "region_of_interest": None,
    "fit_function"      : fitfunc_eff_3,
    "params_guess"      : None,
    "region_of_fit"     : None,
    "fit_density"       : 400
}
fitdata_eff, params_eff, params_eff_err, chi_eff = ultimate_fit( fitbook_abs_eff )

print(prefix + " Effizienzkurve nach Leo: (4)")
print(f"a = ({params_eff[0]:.4f} +- {params_eff_err[0]:.4f}) / 1")
print(f"b = ({params_eff[1]:.4f} +- {params_eff_err[1]:.4f}) / 1")

def plot_abs_eff():
    
    sample_format_dict_Eu = {
        "label"      : r"Messdaten $^{152}$Eu",          
        "fmt"        : 'o', 
        "color"      : sns.color_palette("dark")[0],                               
        "markersize" : 4, 
        "linewidth"  : 1,
        "capsize"    : 2,
        "alpha"      : 1                                   
    }
    sample_format_dict_Cs = {
        "label"      : r"Messdaten $^{137}$Cs",          
        "fmt"        : 'o', 
        "color"      : sns.color_palette("dark")[6],                               
        "markersize" : 4, 
        "linewidth"  : 1,
        "capsize"    : 2,
        "alpha"      : 1                                   
    }
    sample_format_dict_2 = {
        "label"      : f"Fitfunktion $(\\chi^2 = {chi_eff:.1f})$",                       
        "fmt"        : '--', 
        "color"      : sns.color_palette("dark")[9],        
        "markersize" : 4, 
        "linewidth"  : 1,
        "capsize"    : 0,
        "alpha"      : 1
    }  

    writtings = {
        "title"       : None,
        "y_ax_label"  : r"Absolute Effizienz $\epsilon(E)$ / "+r"\%",
        "x_ax_label"  : r"Energie des Photo-Peaks $E_\gamma$ / keV"
    }
    
    general_format_dict = standard_format_dict.copy()
    zoom_params         = no_zooming.copy()
    colorbar_params     = no_colorbar.copy()
    extra_label         = no_extra_label.copy()
    extra_xaxis         = no_extra_xaxis.copy()
 
    
    all_data                = [ data_abs_eff, data_abs_eff_Cs, fitdata_eff ]                               
    all_sample_format_dicts = [ sample_format_dict_Eu, sample_format_dict_Cs, sample_format_dict_2]

    save_plot = True, "../Figures/"+prefix+"_abs_Effizienz.jpg"                                     
    ultimate_plot_advanced (all_data, writtings, zoom_params, colorbar_params, extra_label, extra_xaxis, save_plot, all_sample_format_dicts, general_format_dict)

plot_abs_eff()


# --------------------------- Cobalt Aktivität ---------------------------


    # Will Aktivität der Cobald Quelle aus bekannte Energie-Effizienz-Kurve bestimmen
    # Aktivität <----- Raumwinkel + korregierter Fläche der Cs-Gaußkurve <------- Fläche Gaußkurve + Effizienzkurve

mu_Co1, FWHM_Co1, area_Co1, trash1, trash2 = params_Co_1
mu_Co2, FWHM_Co2, area_Co2, trash1, trash2 = params_Co_2

Energy_axis, Energy_axis_err    = E_of_C_err(channel, np.full(len(channel),0))
weight_axis         = fitfunc_eff_3(Energy_axis, *params_eff)
weight_axis_err     = None   

area_Co1_corr       = np.sum( fitfunc_gauss_bg( channel, mu_Co1, FWHM_Co1, area_Co1, 0, 0 ) / weight_axis )
area_Co2_corr       = np.sum( fitfunc_gauss_bg( channel, mu_Co2, FWHM_Co2, area_Co2, 0, 0 ) / weight_axis )


halftime_Co             = 5.271                         # in Jahren (Nundat)
init_activ_Co           = 67 *10**3                     # in Bq
acitv_Co                = init_activ_Co * 2**(-sources_age/halftime_Co)

dist_Co, dist_Co_err    = 20, 5        # in mm
coverage_Co             = (np.pi * radius_HPGe**2) / (4*np.pi * dist_Co**2)
coverage_Co_err         = coverage_Co * 2*dist_Co_err/dist_Co

activ_Co_meas           = (area_Co1_corr+area_Co2_corr) / coverage_Co
activ_Co_meas_err       = activ_Co_meas * coverage_Co_err/coverage_Co

print(prefix + f" Activity Co lit:  ({acitv_Co*10**(-3):.2f}) kBq")
print(prefix + f" Activity Co meas: ({activ_Co_meas*10**(-3):.2f} +-  {activ_Co_meas_err*10**(-3):.2f}) kBq")



# --------------------------- Tabellen-Erzeugung ---------------------------


Gauss_params        = np.array([params_Co_1, params_Co_2, params_Cs, *params_Eu]).T
Gauss_params_err    = np.array([params_Co_1_err, params_Co_2_err, params_Cs_err, *params_Eu_err]).T

Gauss_params        = np.concatenate(( [known_energies], Gauss_params ))
Gauss_params_err    = np.concatenate(( [np.zeros(len(known_energies))], Gauss_params ))

Gauss_header        = [
    r"$E_{lit}$ / keV",
    r"$\mu$ / Kanal",
    r"$\sigma_{FWHM}$ / Kanal",
    r"$A$ / 1",
    r"$a$ / Kanal$^{-1}$",
    r"$b$ / 1"
]
array_to_latex_2( "../Data/Params_"+prefix+"_Gaussfits.txt", Gauss_params, Gauss_params_err, Gauss_header, ".4f")



Eu_lit_params       = np.array([energies_Eu, rel_itens_Eu*100])
Eu_lit_params_err   = np.array([energies_Eu_err, rel_itens_Eu_err*100])
Eu_lit_header       = [
    r"$E_{lit}$ / keV",
    r"rel. Intens $I$ / $\%$"
]

array_to_latex_2( "../Data/Params_Eu_lit.txt", Eu_lit_params, Eu_lit_params_err, Eu_lit_header, ".4f")


print(f"Aktivität Co: {acitv_Co/1000:.2f} kBq")
print(f"Aktivität Cs: {acitv_Cs/1000:.2f} kBq")
print(f"Aktivität Eu: {acitv_Eu/1000:.2f} kBq")