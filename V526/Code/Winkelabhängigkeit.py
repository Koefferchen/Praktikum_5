from MyPyLib_v3 import *



# ---------------------------- Energie Kalibration Raw Plot ----------------------------


data_kalib      = np.loadtxt("../Data/09_Eu_Kalibration_90grad.txt")
data_kalib_bg   = np.loadtxt("../Data/08_Eu_Hintergrund_90grad.txt")

channel_kalib       = data_kalib[ : , 0]
counts_kalib        = data_kalib[ : , 1] - data_kalib_bg[ : , 1]
counts_kalib_err    = np.sqrt( (np.sqrt(data_kalib[ : , 1]))**2 + (np.sqrt(data_kalib_bg[ : , 1]))**2 ) +1 

exp_data_kalib_noerr  = [channel_kalib, None, counts_kalib, None]
exp_data_kalib        = [channel_kalib, None, counts_kalib, counts_kalib_err]



def plot_kalib( exp_data, fit_datas=None, chis=None, E_of_C=None, C_of_E=None ):

    fit_colors  = sns.color_palette("magma", 7)

    if( isinstance(chis,  np.ndarray) ):
        fit_lables  = [ f"keV-Linie $(\\chi^2 = {chis[0]:.1f})$",
                        f"keV-Linie $(\\chi^2 = {chis[1]:.1f})$",
                        f"121keV-Linie $(\\chi^2 = {chis[2]:.1f})$",
                        f"244keV-Linie $(\\chi^2 = {chis[3]:.1f})$",
                        f"344keV-Linie $(\\chi^2 = {chis[3]:.1f})$",
                        f"778keV-Linie $(\\chi^2 = {chis[3]:.1f})$",
                        f"1112keV-Linie $(\\chi^2 = {chis[3]:.1f})$" ]

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

    all_sample_format_dicts = [ sample_format_dict_1 ]

    if( isinstance(fit_datas, list) ):
        all_data            = [ exp_data, *fit_datas ] 
        for i in range(len(fit_datas)):
            sample_format_dict_gauss["label"] = fit_lables[i]
            sample_format_dict_gauss["color"] = fit_colors[i]
            all_sample_format_dicts.append( sample_format_dict_gauss.copy() )

    else:
        all_data            = [ exp_data ] 

    save_plot = True, "../Figures/Eu_Kalibration.jpg"                                     
    ultimate_plot_advanced (all_data, writtings, zoom_params, colorbar_params, extra_label, extra_xaxis, save_plot, all_sample_format_dicts, general_format_dict)


plot_kalib(exp_data_kalib_noerr)


# ---------------------------- Energie Kalibration Fits Europium ----------------------------

def area( amplitude, sigma ):
    return amplitude * np.sqrt(2*np.pi) * sigma


fitboook_peak_1 = { # K_alpha ?
    "data_set"          : exp_data_kalib,
    "region_of_interest": [240, 430],
    "fit_function"      : fitfunc_gauss_area,
    "params_guess"      : [350, 50, area( 300, 50)],
    "region_of_fit"     : None,
    "fit_density"       : 400
}
fitboook_peak_2 = { # ???
    "data_set"          : exp_data_kalib,
    "region_of_interest": [500, 700],
    "fit_function"      : fitfunc_gauss_area,
    "params_guess"      : [600, 50, area( 50, 50)],
    "region_of_fit"     : None,
    "fit_density"       : 400
}
fitboook_peak_3 = { # 121.7817 keV (Anleitung)
    "data_set"          : exp_data_kalib,
    "region_of_interest": [850, 1050],
    "fit_function"      : fitfunc_gauss_area,
    "params_guess"      : [900, 50, area( 150, 50)],
    "region_of_fit"     : None,
    "fit_density"       : 400
}
fitboook_peak_4 = { # 244.6974 keV
    "data_set"          : exp_data_kalib,
    "region_of_interest": [1550, 2100],
    "fit_function"      : fitfunc_gauss_area,
    "params_guess"      : [1850, 50, area( 50, 50)],
    "region_of_fit"     : None,
    "fit_density"       : 400
}
fitboook_peak_5 = { # 344.2785 keV
    "data_set"          : exp_data_kalib,
    "region_of_interest": [2300, 2800],
    "fit_function"      : fitfunc_gauss_area,
    "params_guess"      : [2550, 60, area( 60, 60)],
    "region_of_fit"     : None,
    "fit_density"       : 400
}
fitboook_peak_6 = { # 778.9045 keV
    "data_set"          : exp_data_kalib,
    "region_of_interest": [5200, 5900],
    "fit_function"      : fitfunc_gauss_area,
    "params_guess"      : [5550, 100, area( 15, 100)],
    "region_of_fit"     : None,
    "fit_density"       : 400
}
fitboook_peak_7 = { # 964.057 keV
    "data_set"          : exp_data_kalib,
    "region_of_interest": [7400, 8100],
    "fit_function"      : fitfunc_gauss_area,
    "params_guess"      : [7700, 100, area( 15, 100)],
    "region_of_fit"     : None,
    "fit_density"       : 400
}
fitbooks    = [fitboook_peak_1, fitboook_peak_2, fitboook_peak_3, fitboook_peak_4, fitboook_peak_5, fitboook_peak_6, fitboook_peak_7]
params_len  = 3 

fitdata_all     = []
params_i        = np.zeros((len(fitbooks), params_len))
params_err_i    = np.zeros((len(fitbooks), params_len))
chis_i          = np.zeros(len(fitbooks))

for i, book in enumerate(fitbooks):
    fitdata_i, params_i[i], params_err_i[i], chis_i[i] = ultimate_fit(book)
    fitdata_all.append(fitdata_i)


plot_kalib(exp_data_kalib_noerr, fitdata_all, chis_i)

array_to_latex( "../Data/Params_E_Kalibration_Gauss.txt", params_i.T, params_err_i.T, ".4f")

# ---------------------------- Energie Kalibration Fits Caesium ----------------------------


data_Cs = np.loadtxt("../Data/01_Caesium_0mm.txt") 

fitbook_Cs = { # 661.7 keV
    "data_set"          : [ data_Cs[:,0], None, data_Cs[:,1], np.sqrt(data_Cs[:,1])+1 ],
    "region_of_interest": [4500, 5500],
    "fit_function"      : fitfunc_gauss_area,
    "params_guess"      : [4850, 170, 400*400 ],
    "region_of_fit"     : None,
    "fit_density"       : 400
}

dummy1, params_Cs, params_Cs_err, chi_Cs = ultimate_fit(fitbook_Cs)
print(f"Cs 0mm Fit: chi^2 = {chi_Cs:.1f}")



# ---------------------------- Energie Kalibration Bestimmung ----------------------------

    # energies(uncercainties) according to Anleitung
known_energies      = np.array([ 661.7, 121.7818, 244.6974, 344.2785, 778.9045, 1112.076])  #964.057 ])
known_energies_err  = np.array([ 0, 0.0003, 0.0008, 0.0012, 0.0024, 0.003 ])
peak_positions      = np.concatenate( ( [params_Cs[0]], params_i[ 2: , 0] ) )
peak_positions_err  = np.concatenate( ( [params_Cs_err[0]], params_err_i[ 2: , 0] ) )
dataset_energy      = [ known_energies, known_energies_err, peak_positions, peak_positions_err ]

fitbook_energy = { 
    "data_set"          : dataset_energy,
    "region_of_interest": None,
    "fit_function"      : fitfunc_linear,
    "params_guess"      : None,
    "region_of_fit"     : None,
    "fit_density"       : 400
}
fitdata_energy, params_energy, params_energy_err, chi_energy = ultimate_fit(fitbook_energy)


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

    save_plot = True, "../Figures/Energie_Kalibration.jpg"                                     
    ultimate_plot_advanced (all_data, writtings, zoom_params, colorbar_params, extra_label, extra_xaxis, save_plot, all_sample_format_dicts, general_format_dict)


plot_energie_kalib( dataset_energy, fitdata_energy, chi_energy )

a, a_err = params_energy[0], params_energy_err[0]
b, b_err = params_energy[1], params_energy_err[1]

print("Energie-Kalibration: K = a*E + b")
print(f"a = ({a:.4f} +- {a_err:.4f}) / keV")
print(f"b = ({b:.2f} +- {b_err:.2f}) / 1")

def C_of_E( E ):
    return a*E + b
def E_of_C( C ):
    return (C-b) / a


plot_kalib(exp_data_kalib_noerr, fitdata_all, chis_i, E_of_C, C_of_E)


# ---------------------------- Energie Effizenz Kalibration ----------------------------


known_rel_intens        = np.array([28.53, 7.55, 26.59, 12.93, 13.68])
known_rel_intens_err    = np.array([0.16, 0.04, 0.20, 0.08, 0.08 ])

meas_intens             = params_i[ 2: , 2]
meas_intens_err         = params_err_i[ 2: , 2]

efficiency              = meas_intens / known_rel_intens
efficiency_err          = efficiency * np.sqrt( (meas_intens_err/meas_intens)**2 + (known_rel_intens_err/known_rel_intens)**2 )

mask_max                = (efficiency == max(efficiency))
norm_efficiency         = efficiency / max(efficiency)
norm_efficiency_err     = np.sqrt( (efficiency_err/max(efficiency))**2 + (efficiency_err[mask_max] * efficiency/max(efficiency)**2)**2 )

meas_channel            = params_i[ 2: , 0]
meas_channel_err        = params_err_i[ 2: , 0]

data_eff                = [ meas_channel, meas_channel_err, norm_efficiency, norm_efficiency_err ]

def effizienz_kurve( x, mu, sigma, area, a, b):
    return fitfunc_gauss_area(x, np.abs(mu), sigma, np.abs(area) ) - (x-b)/np.abs(a) 

def effizienz_kurve_2( x, a, b, c, d):
    return a + b*x + c*x**2 + d*x**3 

fitboook_eff = { 
    "data_set"          : data_eff,
    "region_of_interest": None,
    "fit_function"      : effizienz_kurve,
    "params_guess"      : [1800, 1000, area(1, 1200), 12000, 1000],
    "region_of_fit"     : None,
    "fit_density"       : 400
}

fitboook_eff_2 = { 
    "data_set"          : data_eff,
    "region_of_interest": None,
    "fit_function"      : effizienz_kurve_2,
    "params_guess"      : None,
    "region_of_fit"     : None,
    "fit_density"       : 400
}

fitdata_eff, params_eff, params_eff_err, chi_eff = ultimate_fit(fitboook_eff)

def plot_energie_effizienz():
    
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
        "label"      : f"Empirischer Fit",                       
        "fmt"        : '--', 
        "color"      : sns.color_palette("dark")[9],        
        "markersize" : 4, 
        "linewidth"  : 1,
        "capsize"    : 0,
        "alpha"      : 1
    }  

    writtings = {
        "title"       : None,
        "x_ax_label"  : r"Kanal-Index $K$ / 1",
        "y_ax_label"  : r"Effizienz $\epsilon$ / 1"
    }
    
    general_format_dict = standard_format_dict.copy()
    zoom_params         = no_zooming.copy()
    colorbar_params     = no_colorbar.copy()
    extra_label         = no_extra_label.copy()
    extra_xaxis         = {
        "do_axis"   :   True,
        "f_new(old)":   E_of_C,
        "f_old(new)":   C_of_E,
        "label"     :   r"Energie $E$ / keV" }
 
    
    all_data                = [ data_eff, fitdata_eff ]                               
    all_sample_format_dicts = [ sample_format_dict_1, sample_format_dict_2 ]

    save_plot = True, "../Figures/Energie_Effizienz.jpg"                                     
    ultimate_plot_advanced (all_data, writtings, zoom_params, colorbar_params, extra_label, extra_xaxis, save_plot, all_sample_format_dicts, general_format_dict)

plot_energie_effizienz()

print("mu, sigma, area, a, b = \n", params_eff)
print("uncertainsties: \n", params_eff_err)

array_to_latex( "../Data/Params_Eff_Fit.txt", np.array([ params_eff]), np.array([ params_eff_err]), ".4f") 

