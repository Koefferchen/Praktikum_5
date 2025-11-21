
from MyPyLib_v2 import *
from scipy.special import erf

# ------------------------ Promptkurve plotten ------------------------

def plot_promptkurven( exp_data, fit_datas=None, chis=None, T_of_C=None, C_of_T=None ):

    sample_format_dict_1 = {
        "label"      : "Messdaten",          
        "fmt"        : 'o', 
        "color"      : sns.color_palette("bright")[9],                               
        "markersize" : 1.0, 
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

    general_format_dict["custom_x_range"] = [True, 0, 10000]
    
    if ((T_of_C == None) or (C_of_T == None)):
        extra_xaxis = no_extra_xaxis.copy()
    else:
        extra_xaxis         = {
            "do_axis"   :   True,
            "f_new(old)":   T_of_C,
            "f_old(new)":   C_of_T,
            "label"     :   r"Zeit $t$ [ns]"
        }


    all_sample_format_dicts = [ sample_format_dict_1 ]

    if(fit_datas == None):
        save_plot = True, "../Figures/Na_Promptkurven.jpg"
        all_data = [exp_data]
    else:
        save_plot = True, "../Figures/Na_Promptkurven_T.jpg"
        all_data = [ exp_data, *fit_datas ] 
        fit_colors  = sns.color_palette("magma", len(fit_datas))                           
        for i in range(len(fit_datas)):
            sample_format_dict_gauss["label"] = f"$\\chi^2 = {chis[i]:.1f}$"
            sample_format_dict_gauss["color"] = fit_colors[i]
            all_sample_format_dicts.append( sample_format_dict_gauss.copy() )

                     
    ultimate_plot_advanced (all_data, writtings, zoom_params, colorbar_params, extra_label, extra_xaxis, save_plot, all_sample_format_dicts, general_format_dict)



pro_data_txt    = np.loadtxt("../Data/11_Na_Promptkurve_16ns.txt")
pro_channel     = pro_data_txt[ : , 0]
pro_count       = pro_data_txt[ : , 1]
pro_count_err   = np.sqrt(pro_count)+1
pro_data        = [pro_channel, None, pro_count, pro_count_err]
pro_data_noerr  = [pro_channel, None, pro_count, None]

    # Plot der Rohdaten
plot_promptkurven(pro_data_noerr)


# ------------------------ Promptkurve fitten ------------------------



density = 400

fitbook_1 = {
    "data_set"          : pro_data,
    "region_of_interest": [500, 1200],
    "fit_function"      : gauss_fit,
    "params_guess"      : [800, 200, 40],
    "region_of_fit"     : None,
    "fit_density"       : density
}
fitbook_2 = {
    "data_set"          : pro_data,
    "region_of_interest": [1800, 2400],
    "fit_function"      : gauss_fit,
    "params_guess"      : [2100, 200, 60],
    "region_of_fit"     : None,
    "fit_density"       : density
}
fitbook_3 = {
    "data_set"          : pro_data,
    "region_of_interest": [3050, 3650],
    "fit_function"      : gauss_fit,
    "params_guess"      : [3300, 200, 40],
    "region_of_fit"     : None,
    "fit_density"       : density
}
fitbook_4 = {
    "data_set"          : pro_data,
    "region_of_interest": [4200, 5000],
    "fit_function"      : gauss_fit,
    "params_guess"      : [4700, 200, 40],
    "region_of_fit"     : None,
    "fit_density"       : density
}
fitbook_5 = {
    "data_set"          : pro_data,
    "region_of_interest": [5500, 6200],
    "fit_function"      : gauss_fit,
    "params_guess"      : [5800, 200, 40],
    "region_of_fit"     : None,
    "fit_density"       : density
}
fitbook_6 = {
    "data_set"          : pro_data,
    "region_of_interest": [6900, 7500],
    "fit_function"      : gauss_fit,
    "params_guess"      : [7200, 200, 40],
    "region_of_fit"     : None,
    "fit_density"       : density
}

fitbooks        = [fitbook_1, fitbook_2, fitbook_3, fitbook_4, fitbook_5, fitbook_6]
pro_fit_data_all= []
pro_params      = np.zeros((len(fitbooks), 3))
pro_params_err  = np.zeros((len(fitbooks), 3))
pro_chis        = np.zeros(len(fitbooks))

for i, book in enumerate(fitbooks):
    pro_fit_data_i, pro_params[i], pro_params_err[i], pro_chis[i] = ultimate_fit(book)
    pro_fit_data_all.append(pro_fit_data_i)


    # überarbeiteter Plot, jetzt mit Fit
plot_promptkurven(pro_data_noerr, pro_fit_data_all, pro_chis )


array_to_latex("../Data/Params_Promtkurve_Fits.txt", pro_params.T, pro_params_err.T, ".4f")


# ------------------------ Zeit_Achse fitten ------------------------


pro_peaks       = pro_params[ : , 0]
pro_peaks_err   = pro_params_err[ : , 0]
time_array      = np.array([0, 16, 32, 48, 64, 80])
peak_data       = [time_array, None, pro_peaks, pro_peaks_err]

fitbook_time = {
    "data_set"          : peak_data,
    "region_of_interest": None,
    "fit_function"      : linear_fit,
    "params_guess"      : None,
    "region_of_fit"     : None,
    "fit_density"       : 2
}
time_fit, time_params, time_params_err, time_chi = ultimate_fit( fitbook_time )


def plot_energie_kalib(exp_data, fit_data, chi):
    
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
        "x_ax_label"  : r"Zeit $t$ [ns]"
    }
    
    general_format_dict = standard_format_dict.copy()
    zoom_params         = no_zooming.copy()
    colorbar_params     = no_colorbar.copy()
    extra_label         = no_extra_label.copy()
    extra_xaxis         = no_extra_xaxis.copy()
 
    
    all_data                = [ exp_data, fit_data ]                               
    all_sample_format_dicts = [ sample_format_dict_1, sample_format_dict_2 ]

    save_plot = True, "../Figures/Zeit-Kalibration.jpg"                                     
    ultimate_plot_advanced (all_data, writtings, zoom_params, colorbar_params, extra_label, extra_xaxis, save_plot, all_sample_format_dicts, general_format_dict)
plot_energie_kalib(peak_data, time_fit, time_chi)


a       = time_params[0]
b       = time_params[1]
a_err   = time_params_err[0]
b_err   = time_params_err[1]

print("Zeitkalibration: Channel = a*Zeit + b")
print(f"a = ({a:.2f}+- {a_err:.2f})/ns")
print(f"b = ({b:.2f}+- {b_err:.2f})")

def C_of_T(time):
    return a*time + b
def T_of_C(channel):
    return (channel - b)/a


    # überarbeiteter Plot, jetzt mit Zeitachse
plot_promptkurven(pro_data_noerr, pro_fit_data_all, pro_chis, T_of_C, C_of_T )



# ------------------------ Lebenskurve plotten ------------------------


decay_data_txt  = np.loadtxt("../Data/14_Ba_Lebenskurve.txt")
decay_channel   = decay_data_txt[:, 0]
decay_counts    = decay_data_txt[:, 1]
decay_counts_err= np.sqrt(decay_counts)+1
decay_data_time = [T_of_C(decay_channel), None, decay_counts, decay_counts_err]
decay_data_noerr= [decay_channel, None, decay_counts, None]


def plot_lebenskurve( exp_data, fit_data=None, chi=None, T_of_C=None, C_of_T=None ):

    if(chi == None):
        chi = 0

    sample_format_dict_1 = {
        "label"      : "Messdaten",          
        "fmt"        : 'o', 
        "color"      : sns.color_palette("bright")[9],                               
        "markersize" : 1.0, 
        "linewidth"  : 1,
        "capsize"    : 0,
        "alpha"      : 0.6                                   
    }

    sample_format_dict_fit = {
        "label"      : f"Fit ($\\chi^2 = {chi:.1f}$)",          
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
    
    if ((T_of_C == None) or (C_of_T == None)):
        extra_xaxis = no_extra_xaxis.copy()
    else:
        extra_xaxis         = {
            "do_axis"   :   True,
            "f_new(old)":   T_of_C,
            "f_old(new)":   C_of_T,
            "label"     :   r"Zeit $t$ [ns]"
        }

    if(fit_data == None):
        save_plot               = True, "../Figures/Ba_Lebenskurve.jpg"
        all_data                = [exp_data]
        all_sample_format_dicts = [ sample_format_dict_1 ]

    else:
        save_plot               = True, "../Figures/Ba_Lebenskurve_T.jpg"
        all_data                = [ exp_data, fit_data ] 
        all_sample_format_dicts = [ sample_format_dict_1, sample_format_dict_fit ]
                     
    ultimate_plot_advanced (all_data, writtings, zoom_params, colorbar_params, extra_label, extra_xaxis, save_plot, all_sample_format_dicts, general_format_dict)


    # Plotten der Rohdaten mit Zeit-Achse 
plot_lebenskurve(decay_data_noerr, None, None, T_of_C, C_of_T)


# ------------------------ Lebenskurve fitten ------------------------


    # Fitfunktion nach Anleitung. Eigentlich müsst vor erf() noch "np.sqrt(np.pi)/2" stehen.
def faltung_fit( t, tau, t_0, sigma, I):
    return I/(2*tau) * np.exp( (sigma**2 - 2*tau*(t-t_0))/(2*tau**2) ) * ( 1 + erf( 2**(-1/2) *((t-t_0)/sigma - sigma/tau ) ) )

fitbook_leben = {
    "data_set"          : decay_data_time,
    "region_of_interest": None,
    "fit_function"      : faltung_fit,
    "params_guess"      : [9, 8, 5, 100/(2*9)],
    "region_of_fit"     : None,
    "fit_density"       : 400
}

decay_fit, decay_params, decay_params_err, decay_chi = ultimate_fit( fitbook_leben )
decay_fit[0] = C_of_T(decay_fit[0])

plot_lebenskurve(decay_data_noerr, decay_fit, decay_chi, T_of_C, C_of_T)


tau, t_0, sigma, I                  = decay_params
tau_err, t_0_err, sigma_err, I_err  = decay_params_err
t_halb, t_halb_err                  = tau * np.log(2), tau_err * np.log(2)

print("Lebenskurve Fitparameter:")
print(f"tau   = ({tau:.2f} +- {tau_err:.2f})ns ")
print(f"t_halb= ({t_halb:.2f} +- {t_halb_err:.2f})ns ")
print(f"t_0   = ({t_0:.2f} +- {t_0_err:.2f})ns")
print(f"sigma = ({sigma:.2f} +- {sigma_err:.2f})ns")
print(f"I     = ({I:.2f} +- {I_err:.2f})")


# ------------------------ Zeitauflösung bestimmen ------------------------


def T_of_C_mit_err(channel, channel_err):
    T = (channel - b)/a
    T_err = np.sqrt( (channel_err/a)**2 + (b_err/a)**2 + (T*a_err/a)**2 )
    return T, T_err

sigmas      = pro_params[ : , 1 ]
sigmas_err  = pro_params_err[ : , 1 ]

FWHMs       = 2 * sigmas *np.sqrt( 2 * np.log(2))
FWHMs_err   = 2 * sigmas_err *np.sqrt( 2 * np.log(2))

FWHM_av, FWHM_av_err = stichproben_varianz( FWHMs )

print(FWHMs)
print(FWHMs_err)
print(f"Zeitauflösung: FWHM = ({FWHM_av:.2f} +- {FWHM_av_err:.2f}) Kanal")

FWHM_av_T, FWHM_av_T_err = T_of_C_mit_err(FWHM_av, FWHM_av_err)
print(f"Zeitauflösung: FWHM = ({FWHM_av_T:.2f} +- {FWHM_av_T_err:.2f}) ns")

