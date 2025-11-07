
from MyPyLib_v2 import *


# -------------------------- Plotting of raw energy spectrums --------------------------

dateinamen = [  "01_Na_Spektrum_Links", 
                "02_Na_Spektrum_Kalibration_Links",
                "03_Ba_Spektrum_Kalibration_Links",
                "04_Na_Spektrum_Rechts",
                "05_Na_Spektrum_Kalibration_Rechts",
                "06_Ba_Spektrum_Kalibration_Rechts"  ]

   

def plot_raw_data( n ):
    
    data = np.loadtxt("../Data/"+dateinamen[n]+".txt")
    channels    = data[ : , 0]
    counts      = data[ : , 1]
    if ((n == 2) or (n == 5)):
        color = sns.color_palette("bright")[0]
    else:
        color = sns.color_palette("bright")[0]



    sample_format_dict_1 = {
        "label"      : None,          
        "fmt"        : 'o', 
        "color"      : color,                               
        "markersize" : 0.7, 
        "linewidth"  : 1,
        "capsize"    : 0,
        "alpha"      : 0.6                                   
    }

    writtings = {
        "title"       : None,
        "x_ax_label"  : r"Kanal-Index $K$ [1]",
        "y_ax_label"  : r"Ereignisse $N$ [1]"
    }
    
    general_format_dict = standard_format_dict.copy()
    zoom_params         = no_zooming.copy()
    colorbar_params     = no_colorbar.copy()
    extra_label         = no_extra_label.copy()
    extra_xaxis         = no_extra_xaxis.copy()

    data_set_1  = channels, None, counts, None
    
    all_data                = [ data_set_1 ]                               
    all_sample_format_dicts = [ sample_format_dict_1 ]

    save_plot = True, "../Figures/"+dateinamen[n]+".jpg"                                     
    ultimate_plot_advanced (all_data, writtings, zoom_params, colorbar_params, extra_label, extra_xaxis, save_plot, all_sample_format_dicts, general_format_dict)
plot_raw_data(1)
plot_raw_data(1)
#for j in range(6):
#    plot_raw_data(j)


# -------------------------- Gauss-Plotting Na  --------------------------


def plot_fit_data_Na(exp_data, fit_data, chi, parity, id, E_of_C=None, C_of_E=None ):

    sample_format_dict_1 = {
        "label"      : "Messdaten",          
        "fmt"        : 'o', 
        "color"      : sns.color_palette("bright")[9],                               
        "markersize" : 0.7, 
        "linewidth"  : 1,
        "capsize"    : 0,
        "alpha"      : 0.6                                   
    }
    sample_format_dict_gauss = {
        "label"      : f"Fit an $511$keV-Linie $(\\chi^2 = {chi:.1f})$ ",          
        "fmt"        : '-', 
        "color"      : sns.color_palette("magma", 2)[0],                               
        "markersize" : 0.7, 
        "linewidth"  : 2,
        "capsize"    : 0,
        "alpha"      : 1                                   
    }

    writtings = {
        "title"       : None,
        "x_ax_label"  : r"Kanal-Index $K$ [1]",
        "y_ax_label"  : r"Ereignisse $N$ [1]"
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
            "label"     :   r"Energie $E$ [keV]"
        }
    
    all_data                = [ exp_data, fit_data ]                               
    all_sample_format_dicts = [ sample_format_dict_1, sample_format_dict_gauss ]

    save_plot = True, "../Figures/0"+id+"_Na_"+parity+"_511keV_Fit.jpg"                                     
    ultimate_plot_advanced (all_data, writtings, zoom_params, colorbar_params, extra_label, extra_xaxis, save_plot, all_sample_format_dicts, general_format_dict)


# -------------------------- Gauss-Fit Na (Left) --------------------------


data_left = np.loadtxt("../Data/"+dateinamen[1]+".txt")
Na_channels_left    = data_left[ : , 0]
Na_counts_left      = data_left[ : , 1]
Na_counts_err_left  = np.sqrt(Na_counts_left)+1
exp_data_Na_left    = [Na_channels_left, None, Na_counts_left, Na_counts_err_left ]
exp_noerr_Na_left   = [Na_channels_left, None, Na_counts_left, None ]

fitbook_Na_left = {
    "data_set"          : exp_data_Na_left,
    "region_of_interest": [7000, 7800],
    "fit_function"      : gauss_fit,
    "params_guess"      : [7400, 700, 175],
    "region_of_fit"     : None,
    "fit_density"       : 400
}
fit_data_Na_left, params_Na_left_511, params_err_Na_left_511, chi_Na_left = ultimate_fit( fitbook_Na_left )

plot_fit_data_Na( exp_noerr_Na_left, fit_data_Na_left, chi_Na_left, "Links", "2" )


# -------------------------- Gauss-Fit Na (Right) --------------------------


data_right = np.loadtxt("../Data/"+dateinamen[4]+".txt")
Na_channels_right    = data_right[ : , 0]
Na_counts_right      = data_right[ : , 1]
Na_counts_err_right  = np.sqrt(Na_counts_right)+1
exp_data_Na_right    = [Na_channels_right, None, Na_counts_right, Na_counts_err_right ]
exp_noerr_Na_right    = [Na_channels_right, None, Na_counts_right, None ]

fitbook_Na_right = {
    "data_set"          : exp_data_Na_right,
    "region_of_interest": [6700, 7700],
    "fit_function"      : gauss_fit,
    "params_guess"      : [7300, 700, 120],
    "region_of_fit"     : None,
    "fit_density"       : 400
}
fit_data_Na_right, params_Na_right_511, params_err_Na_right_511, chi_Na_right = ultimate_fit( fitbook_Na_right )

plot_fit_data_Na( exp_noerr_Na_right, fit_data_Na_right, chi_Na_right, "Rechts", "5" )



# -------------------------- Gauss-Plotting Ba --------------------------


def plot_fit_data_Ba( exp_data, fit_datas, chis, parity, id, E_of_C=None, C_of_E=None ):

    fit_colors  = sns.color_palette("magma", 6)
    fit_lables  = [ f"Fit an 31keV-Linie $(\\chi^2 = {chis[0]:.1f})$",
                    f"Fit an 81keV-Linie $(\\chi^2 = {chis[1]:.1f})$",
                    f"Fit an 302keV-Linie $(\\chi^2 = {chis[2]:.1f})$",
                    f"Fit an 356keV-Linie $(\\chi^2 = {chis[3]:.1f})$", ]

    sample_format_dict_1 = {
        "label"      : "Messdaten",          
        "fmt"        : 'o', 
        "color"      : sns.color_palette("bright")[9],                               
        "markersize" : 0.7, 
        "linewidth"  : 1,
        "capsize"    : 0,
        "alpha"      : 0.6                                   
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
        "x_ax_label"  : r"Kanal-Index $K$ [1]",
        "y_ax_label"  : r"Ereignisse $N$ [1]"
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
            "label"     :   r"Energie $E$ [keV]"
        }

    
    all_data                = [ exp_data, *fit_datas ]                               
    all_sample_format_dicts = [ sample_format_dict_1 ]

    for i in range(len(fit_datas)):
        sample_format_dict_gauss["label"] = fit_lables[i]
        sample_format_dict_gauss["color"] = fit_colors[i]
        all_sample_format_dicts.append( sample_format_dict_gauss.copy() )

    save_plot = True, "../Figures/0"+id+"_Ba_"+parity+"_Fits.jpg"                                     
    ultimate_plot_advanced (all_data, writtings, zoom_params, colorbar_params, extra_label, extra_xaxis, save_plot, all_sample_format_dicts, general_format_dict)


# -------------------------- Gauss-Fits Ba (Left) --------------------------


data = np.loadtxt("../Data/"+dateinamen[2]+".txt")
Ba_channels_left    = data[ : , 0]
Ba_counts_left      = data[ : , 1]
Ba_counts_err_left  = np.sqrt(Ba_counts_left)+1

exp_data_Ba_left = [Ba_channels_left, None, Ba_counts_left, Ba_counts_err_left ]
exp_noerr_Ba_left = [Ba_channels_left, None, Ba_counts_left, None ]

fitbook_Ba_left_31keV = {
    "data_set"          : exp_data_Ba_left,
    "region_of_interest": [320, 480],
    "fit_function"      : gauss_fit,
    "params_guess"      : [400, 100, 6500],
    "region_of_fit"     : None,
    "fit_density"       : 400
}
fitbook_Ba_left_81keV = {
    "data_set"          : exp_data_Ba_left,
    "region_of_interest": [1080, 1290],
    "fit_function"      : gauss_fit,
    "params_guess"      : [1250, 80, 1600],
    "region_of_fit"     : None,
    "fit_density"       : 400
}
fitbook_Ba_left_302keV = {
    "data_set"          : exp_data_Ba_left,
    "region_of_interest": [4210, 4510],
    "fit_function"      : gauss_fit,
    "params_guess"      : [4350, 250, 180],
    "region_of_fit"     : None,
    "fit_density"       : 400
}
fitbook_Ba_left_356keV = {
    "data_set"          : exp_data_Ba_left,
    "region_of_interest": [4750, 5270],
    "fit_function"      : gauss_fit,
    "params_guess"      : [5150, 350, 400],
    "region_of_fit"     : None,
    "fit_density"       : 400
}

fit_data_Ba_left_31keV, params_Ba_left_31keV, params_err_Ba_left_31keV, chi_Ba_left_31keV = ultimate_fit( fitbook_Ba_left_31keV )
fit_data_Ba_left_81keV, params_Ba_left_81keV, params_err_Ba_left_81keV, chi_Ba_left_81keV = ultimate_fit( fitbook_Ba_left_81keV )
fit_data_Ba_left_302keV, params_Ba_left_302keV, params_err_Ba_left_302keV, chi_Ba_left_302keV = ultimate_fit( fitbook_Ba_left_302keV )
fit_data_Ba_left_356keV, params_Ba_left_356keV, params_err_Ba_left_356keV, chi_Ba_left_356keV = ultimate_fit( fitbook_Ba_left_356keV )

fit_datas_Ba_left   = [fit_data_Ba_left_31keV, fit_data_Ba_left_81keV, fit_data_Ba_left_302keV, fit_data_Ba_left_356keV]
chis_Ba_left        = [chi_Ba_left_31keV, chi_Ba_left_81keV, chi_Ba_left_302keV, chi_Ba_left_356keV]

plot_fit_data_Ba( exp_noerr_Ba_left, fit_datas_Ba_left, chis_Ba_left, "Links", "3")

Kalibration_params_left = np.array([ params_Ba_left_31keV, params_Ba_left_81keV, params_Ba_left_302keV, params_Ba_left_356keV, params_Na_left_511 ])
Kalibration_params_left_err = np.array([ params_err_Ba_left_31keV, params_err_Ba_left_81keV, params_err_Ba_left_302keV, params_err_Ba_left_356keV, params_err_Na_left_511])

array_to_latex( "../Data/Params_Kalibration_Links.txt", Kalibration_params_left.T, Kalibration_params_left_err.T)


# -------------------------- Gauss-Fits Ba (Right) --------------------------


data = np.loadtxt("../Data/"+dateinamen[5]+".txt")
Ba_channels_right   = data[ : , 0]
Ba_counts_right     = data[ : , 1]
Ba_counts_err_right = np.sqrt(Ba_counts_right)+1

exp_data_Ba_right   = [Ba_channels_right, None, Ba_counts_right, Ba_counts_err_right ]
exp_noerr_Ba_right  = [Ba_channels_right, None, Ba_counts_right, None ]

fitbook_Ba_right_31keV = {
    "data_set"          : exp_data_Ba_right,
    "region_of_interest": [290, 440],
    "fit_function"      : gauss_fit,
    "params_guess"      : [400, 100, 5500],
    "region_of_fit"     : None,
    "fit_density"       : 400
}
fitbook_Ba_right_81keV = {
    "data_set"          : exp_data_Ba_right,
    "region_of_interest": [1060, 1250],
    "fit_function"      : gauss_fit,
    "params_guess"      : [1250, 80, 1600],
    "region_of_fit"     : None,
    "fit_density"       : 400
}
fitbook_Ba_right_302keV = {
    "data_set"          : exp_data_Ba_right,
    "region_of_interest": [4210, 4510],
    "fit_function"      : gauss_fit,
    "params_guess"      : [4350, 250, 180],
    "region_of_fit"     : None,
    "fit_density"       : 400
}
fitbook_Ba_right_356keV = {
    "data_set"          : exp_data_Ba_right,
    "region_of_interest": [4750, 5270],
    "fit_function"      : gauss_fit,
    "params_guess"      : [5150, 350, 400],
    "region_of_fit"     : None,
    "fit_density"       : 400
}

fit_data_Ba_right_31keV, params_Ba_right_31keV, params_err_Ba_right_31keV, chi_Ba_right_31keV = ultimate_fit( fitbook_Ba_right_31keV )
fit_data_Ba_right_81keV, params_Ba_right_81keV, params_err_Ba_right_81keV, chi_Ba_right_81keV = ultimate_fit( fitbook_Ba_right_81keV )
fit_data_Ba_right_302keV, params_Ba_right_302keV, params_err_Ba_right_302keV, chi_Ba_right_302keV = ultimate_fit( fitbook_Ba_right_302keV )
fit_data_Ba_right_356keV, params_Ba_right_356keV, params_err_Ba_right_356keV, chi_Ba_right_356keV = ultimate_fit( fitbook_Ba_right_356keV )

fit_datas_Ba_right   = [fit_data_Ba_right_31keV, fit_data_Ba_right_81keV, fit_data_Ba_right_302keV, fit_data_Ba_right_356keV]
chis_Ba_right        = [chi_Ba_right_31keV, chi_Ba_right_81keV, chi_Ba_right_302keV, chi_Ba_right_356keV]

plot_fit_data_Ba( exp_noerr_Ba_right, fit_datas_Ba_right, chis_Ba_right, "Rechts", "6")

Kalibration_params_right        = np.array([ params_Ba_right_31keV, params_Ba_right_81keV, params_Ba_right_302keV, params_Ba_right_356keV, params_Na_right_511 ])
Kalibration_params_right_err    = np.array([ params_err_Ba_right_31keV, params_err_Ba_right_81keV, params_err_Ba_right_302keV, params_err_Ba_right_356keV, params_err_Na_right_511])

array_to_latex( "../Data/Params_Kalibration_Rechts.txt", Kalibration_params_right.T, Kalibration_params_right_err.T)



# -------------------------- Energie Kalibraiton Plotting --------------------------


def plot_energie_kalib(exp_data, fit_data, chi, parity):
    
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
        "y_ax_label"  : r"Kanal-Index $K$ [1]",
        "x_ax_label"  : r"Energie $E$ [keV]"
    }
    
    general_format_dict = standard_format_dict.copy()
    zoom_params         = no_zooming.copy()
    colorbar_params     = no_colorbar.copy()
    extra_label         = no_extra_label.copy()
    extra_xaxis         = no_extra_xaxis.copy()
 
    
    all_data                = [ exp_data, fit_data ]                               
    all_sample_format_dicts = [ sample_format_dict_1, sample_format_dict_2 ]

    save_plot = True, "../Figures/Energy_Kalibration_Fit_"+parity+".jpg"                                     
    ultimate_plot_advanced (all_data, writtings, zoom_params, colorbar_params, extra_label, extra_xaxis, save_plot, all_sample_format_dicts, general_format_dict)

def energy_of_channel( channel, channel_err, params_kalib, params_kalib_err ):
    a           = params_kalib[0]
    b           = params_kalib[1]
    a_err       = params_kalib_err[0]
    b_err       = params_kalib_err[1]
    energy      = (channel - b)/a 
    energy_err  = np.sqrt( (channel_err/a)**2 + (b_err/a)**2 + (energy*a_err/a)**2 )
    return energy, energy_err


    # Source: https://www.nndc.bnl.gov/ensnds/133/Cs/ec_decay_10.551_y.pdf:  80.9979 (11),  302.8508 (5), 356.0129 (7)
    # Source: https://www.ezag.com/wp-content/uploads/2023/08/Ba-133.pdf:    30.85
    # Source: https://www-nds.iaea.org/xgamma_standards/genergies1.htm?utm_source=chatgpt.com: 511

energy_array        = np.array([ 30.85, 80.9979, 302.8508, 356.0129, 511.0 ])
energy_err_array    = np.array([ 0.0, 0.0011, 0.0005, 0.0007, 0.0])

# -------------------------- Energie Kalibraiton Left --------------------------


channel_array       = Kalibration_params_left[ : , 0 ]
channel_array_err   = Kalibration_params_left_err[ : , 0]

fitbook_kalib_left = {
    "data_set"          : [energy_array, None, channel_array, channel_array_err],
    "region_of_interest": None,
    "fit_function"      : linear_fit,
    "params_guess"      : None,
    "region_of_fit"     : None,
    "fit_density"       : 400
}
fitdata_kalib_left, params_kalib_left, params_err_kalib_left, chi_kalib_left = ultimate_fit(fitbook_kalib_left)
expdata_kalib_left  = energy_array, None, channel_array, channel_array_err      

plot_energie_kalib(expdata_kalib_left, fitdata_kalib_left, chi_kalib_left, "Links")

E_peaks_left, E_peaks_err_left = energy_of_channel(channel_array, channel_array_err, params_kalib_left, params_err_kalib_left)
print("E     links: ", E_peaks_left)
print("E_err links: ", E_peaks_err_left)


# -------------------------- Energie Kalibraiton Right --------------------------


channel_array       = Kalibration_params_right[ : , 0 ]
channel_array_err   = Kalibration_params_right_err[ : , 0]

fitbook_kalib_right = {
    "data_set"          : [energy_array, None, channel_array, channel_array_err],
    "region_of_interest": None,
    "fit_function"      : linear_fit,
    "params_guess"      : None,
    "region_of_fit"     : None,
    "fit_density"       : 400
}
fitdata_kalib_right, params_kalib_right, params_err_kalib_right, chi_kalib_right = ultimate_fit(fitbook_kalib_right)
expdata_kalib_right  = energy_array, None, channel_array, channel_array_err      

plot_energie_kalib(expdata_kalib_right, fitdata_kalib_right, chi_kalib_right, "Rechts")

E_peaks_right, E_peaks_err_right = energy_of_channel(channel_array, channel_array_err, params_kalib_right, params_err_kalib_right)
print("E     rechts: ", E_peaks_right)
print("E_err rechts: ", E_peaks_err_right)


array_to_latex("../Data/Params_E_kalib.txt", np.array([params_kalib_left, params_kalib_right]).T, np.array([params_err_kalib_left, params_err_kalib_right]).T, ".4f" )
array_to_latex("../Data/Results_E.txt", np.array([energy_array, E_peaks_left, E_peaks_right ]), np.array([energy_err_array, E_peaks_err_left, E_peaks_err_right]), ".4f")



# -------------------------- Plotting mit Energie-Achse --------------------------


a_L          = params_kalib_left[0]
b_L          = params_kalib_left[1]
a_R          = params_kalib_right[0]
b_R          = params_kalib_right[1]

def C_L(E_L):
    return a_L*E_L + b_L 
def C_R(E_R):
    return a_R*E_R + b_R
def E_L(C_L):
    return (C_L - b_L)/a_L 
def E_R(C_R):
    return (C_R - b_R)/a_R 


def plot_raw_kalibrated( dateiname, E_of_C, C_of_E ):
    
    data = np.loadtxt("../Data/"+dateiname+".txt")
    channels    = data[ : , 0]
    counts      = data[ : , 1]

    color = sns.color_palette("bright")[0]


    sample_format_dict_1 = {
        "label"      : None,          
        "fmt"        : 'o', 
        "color"      : color,                               
        "markersize" : 0.7, 
        "linewidth"  : 1,
        "capsize"    : 0,
        "alpha"      : 0.6                                   
    }

    writtings = {
        "title"       : None,
        "x_ax_label"  : r"Kanal-Index $K$ [1]",
        "y_ax_label"  : r"Ereignisse $N$ [1]"
    }
    
    general_format_dict = standard_format_dict.copy()
    zoom_params         = no_zooming.copy()
    colorbar_params     = no_colorbar.copy()
    extra_label         = no_extra_label.copy()
    extra_xaxis         = {
        "do_axis"   :   True,
        "f_new(old)":   E_of_C,
        "f_old(new)":   C_of_E,
        "label"     :   r"Energie $E$ [keV]"
    }

    data_set_1  = channels, None, counts, None
    
    all_data                = [ data_set_1 ]                               
    all_sample_format_dicts = [ sample_format_dict_1 ]

    save_plot = True, "../Figures/"+dateiname+".jpg"                                     
    ultimate_plot_advanced (all_data, writtings, zoom_params, colorbar_params, extra_label, extra_xaxis, save_plot, all_sample_format_dicts, general_format_dict)


plot_raw_kalibrated( "07_Na_Spektrum_Fenster_Links", E_L, C_L )
plot_raw_kalibrated( "08_Na_Spektrum_Fenster_Rechts", E_R, C_R )


    # überschreibe frühere Gaußfitplots, da nun Energieskalierung bekannt:
plot_fit_data_Na( exp_noerr_Na_left, fit_data_Na_left, chi_Na_left, "Links", "2", E_L, C_L )
plot_fit_data_Na( exp_noerr_Na_right, fit_data_Na_right, chi_Na_right, "Rechts", "5", E_R, C_R )
plot_fit_data_Ba( exp_noerr_Ba_left, fit_datas_Ba_left, chis_Ba_left, "Links", "3", E_L, C_L)
plot_fit_data_Ba( exp_noerr_Ba_right, fit_datas_Ba_right, chis_Ba_right, "Rechts", "6", E_R, C_R)


# -------------------------- Ba mit mit CFD-Schwelle auf 0 --------------------------

plot_raw_kalibrated( "09_Ba_Spektrum_CFD_Links", E_L, C_L )
plot_raw_kalibrated( "10_Ba_Spektrum_CFD_Rechts", E_R, C_R )
plot_raw_kalibrated( "12_Ba_Links_356_Fenster", E_L, C_L )
plot_raw_kalibrated( "13_Ba_Rechts_81_Fenster", E_R, C_R )


