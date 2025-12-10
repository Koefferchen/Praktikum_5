from MyPyLib_v3 import *

def plot_multifit( data_exp, save_as, datas_fit=None, labels_fit=None, E_of_C=None, C_of_E=None ):


    sample_format_dict_1 = {
        "label"      : "Messdaten",          
        "fmt"        : 'o', 
        "color"      : "grey",                               
        "markersize" : 0.5, 
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

def fitfunc_gauss_bg( x, mu, sigma, amp, a, b ):
    return fitfunc_gauss_amplitude(x, mu, sigma, amp) + fitfunc_linear(x, a, b)

def fitfunc_gauss_double_bg( x, mu1, sigma1, amp1, mu2, sigma2, amp2, a, b ):
    return fitfunc_gauss_amplitude(x, mu1, sigma1, amp1) + fitfunc_gauss_amplitude(x, mu2, sigma2, amp2) + fitfunc_linear(x, a, b)


# --------------------------- Plot Raw Data ---------------------------


data_Voraufgabe = np.loadtxt("../Data/spectrum_Voraufgabe.txt")
channel         = data_Voraufgabe[:, 0]
counts          = data_Voraufgabe[:, 1]

data_exp_noerr = [channel, None, counts, None]
data_exp       = [channel, None, counts, np.sqrt(counts)]

plot_multifit(data_exp, "Voraufgabe.jpg")
plot_multifit(data_exp, "Voraufgabe.jpg")


# --------------------------- Fit Data ---------------------------


fitbook_1 = {
    "data_set"          : data_exp,
    "region_of_interest": [250, 410],
    "fit_function"      : fitfunc_gauss_bg,
    "params_guess"      : [350, 50, 500, 0, 0],
    "region_of_fit"     : None,
    "fit_density"       : 400
}
fitbook_2 = {
    "data_set"          : data_exp,
    "region_of_interest": [700, 1100],
    "fit_function"      : fitfunc_gauss_bg,
    "params_guess"      : [900, 60, 300, 0, 0],
    "region_of_fit"     : None,
    "fit_density"       : 400
}
fitbook_5 = {
    "data_set"          : data_exp,
    "region_of_interest": [500, 700],
    "fit_function"      : fitfunc_gauss_bg,
    "params_guess"      : [600, 20, 60, 0, 0],
    "region_of_fit"     : None,
    "fit_density"       : 400
}
fitbook_3_4 = {
    "data_set"          : data_exp,
    "region_of_interest": [3000, 4500],
    "fit_function"      : fitfunc_gauss_double_bg,
    "params_guess"      : [3350, 60, 50, 3850, 100, 100, 0, 0],
    "region_of_fit"     : None,
    "fit_density"       : 400
}

fitdata_1, params_1, params_err_1, chi_1 = ultimate_fit( fitbook_1 )
fitdata_2, params_2, params_err_2, chi_2 = ultimate_fit( fitbook_2 )
fitdata_5, params_5, params_err_5, chi_5 = ultimate_fit( fitbook_5 )
fitdata_3_4, params_3_4, params_err_3_4, chi_3_4 = ultimate_fit( fitbook_3_4 )

fitdata_all     = [fitdata_1,fitdata_5, fitdata_2, fitdata_3_4]
fitlabels_all   = [ f"Fit 1 $(\\chi^2 = {chi_1:.1f} )$",
                    f"Fit 2 $(\\chi^2 = {chi_2:.1f})$",
                    f"Fit 5 $(\\chi^2 = {chi_5:.1f})$",
                    f"Fit 3+4 $(\\chi^2 = {chi_3_4:.1f})$"
]

plot_multifit(data_exp, "Voraufgabe_Fit.jpg", fitdata_all, fitlabels_all)
