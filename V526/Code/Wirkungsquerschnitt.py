
from MyPyLib_v3 import *


# ---------------------- Plotte Rohdaten ----------------------


dateinamen = [
    "01_Caesium_0mm",
    "02_Caesium_1mm",
    "03_Caesium_5mm",
    "07_Caesium_5mm_verdreht_45grad",
    "04_Caesium_10mm",
    "05_Caesium_20mm",
    "06_Caesium_30mm"
]
dateinamen = dateinamen[::-1]


data_len = len(dateinamen)

def plot_raw_data( n ):
    
    data = np.loadtxt("../Data/"+dateinamen[n]+".txt")
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
    extra_xaxis         = no_extra_xaxis.copy()

    data_set_1  = channels, None, counts, None
    all_data                = [ data_set_1 ]                               
    all_sample_format_dicts = [ sample_format_dict_1 ]

    save_plot = True, "../Figures/"+dateinamen[n]+".jpg"                                     
    ultimate_plot_advanced (all_data, writtings, zoom_params, colorbar_params, extra_label, extra_xaxis, save_plot, all_sample_format_dicts, general_format_dict)

plot_raw_data(0)

    # plotting raw data
for i in range(data_len):
    plot_raw_data(i)


# ---------------------- Fitte Rohdaten ----------------------


params_len = 3

fitbook_vorlage = {
    "data_set"          : None,
    "region_of_interest": [4500, 5500],
    "fit_function"      : fitfunc_gauss_area,
    "params_guess"      : [4850, 170, 400*400 ],
    "region_of_fit"     : None,
    "fit_density"       : 400
}
fitbooks = []

for i in range(data_len):
    data = np.loadtxt("../Data/"+dateinamen[i]+".txt")
    channels    = data[ : , 0]
    counts      = data[ : , 1]
    fitbook_vorlage["data_set"] = [channels, None, counts, np.sqrt(counts)+1]
    fitbooks.append( fitbook_vorlage.copy() )

fitdata_all     = []
params_i        = np.zeros((len(fitbooks), params_len))
params_err_i    = np.zeros((len(fitbooks), params_len))
chis_i          = np.zeros(len(fitbooks))

for i, book in enumerate(fitbooks):
    fitdata_i, params_i[i], params_err_i[i], chis_i[i] = ultimate_fit(book)
    fitdata_all.append(fitdata_i)


# ---------------------- Plotte Rohdaten + Fit alle zusammen ----------------------


def plot_fitted_data(E_of_C = None, C_of_E = None):

    labels = [
        "0mm",
        "1mm",
        "5mm",
        r"$5\sqrt{2}$mm",
        "10mm",
        "20mm",
        "30mm"
    ]
    labels = labels[::-1]
    
    sample_format_dict_data = {
        "label"      : None,          
        "fmt"        : 'o', 
        "color"      : None,                               
        "markersize" : 0.7, 
        "linewidth"  : 1,
        "capsize"    : 0,
        "alpha"      : 0.3                                   
    }
    sample_format_dict_fit = {
        "label"      : None,          
        "fmt"        : '-', 
        "color"      : None,                               
        "markersize" : 0.7, 
        "linewidth"  : 1,
        "capsize"    : 0,
        "alpha"      : 1.0                                   
    }
    
    general_format_dict = standard_format_dict.copy()
    zoom_params         = no_zooming.copy()
    colorbar_params     = no_colorbar.copy()
    extra_label         = no_extra_label.copy()
    general_format_dict["fig_side_lengh"] = [9,5]

    if ((E_of_C == None) or (C_of_E == None)):
            extra_xaxis = no_extra_xaxis.copy()
    else:
        extra_xaxis         = {
            "do_axis"   :   True,
            "f_new(old)":   E_of_C,
            "f_old(new)":   C_of_E,
            "label"     :   r"Energie $E$ / keV"
        }

    shift       = 100
    all_data    = []
    all_sample_format_dicts = []
    colors_bright   = sns.color_palette("bright", data_len)
    colors_dark     = sns.color_palette("dark", data_len)

    for i in range(data_len):
        data = np.loadtxt("../Data/"+dateinamen[i]+".txt")
        channels    = data[ : , 0]
        counts      = data[ : , 1]
        x_fit, x_fit_err, y_fit, y_fit_err = fitdata_all[i]

        all_data.append( [channels, None, counts +i*shift, None] )
        all_data.append( [x_fit, x_fit_err, y_fit + i * shift, y_fit_err] )

        sample_format_dict_data["color"] = colors_bright[i]
        sample_format_dict_fit["color"]  = colors_dark[i]
        sample_format_dict_data["label"] = labels[i] + f" $(\\chi^2 = {chis_i[i]:.1f})$"

        all_sample_format_dicts.append( sample_format_dict_data.copy() )
        all_sample_format_dicts.append( sample_format_dict_fit.copy() )

    writtings = {
        "title"       : None,
        "x_ax_label"  : r"Kanal-Index $K$ / 1",
        "y_ax_label"  : f"Ereignisse $(N + {shift:.0f}n)$ / 1"
    }

    save_plot = True, "../Figures/00_Caesium_ALLmm.jpg"                                     
    ultimate_plot_advanced (all_data, writtings, zoom_params, colorbar_params, extra_label, extra_xaxis, save_plot, all_sample_format_dicts, general_format_dict)

plot_fitted_data()


# ---------------------- Exportiere Fit-Parameter nach Latex ----------------------


array_to_latex( "../Data/Params_TotWirkQuerschnitt_Gauss.txt", params_i.T, params_err_i.T, ".4f")


# ---------------------- Bestimme mu und sigma_tot ----------------------

Area_0      = params_i[ data_len-1, 2 ]
Area_0_err  = params_err_i[ data_len-1, 2 ]
Areas       = params_i[ :(data_len-1) , 2]
Areas_err   = params_err_i[ :(data_len-1) , 2]

Intens      = Areas / Area_0
Intens_err  = Intens * np.sqrt( (Areas_err/Areas)**2 + (Area_0_err/Area_0)**2 )

log_Intens      = - np.log(Intens)
log_Intens_err  = Intens_err / Intens

thickness       = np.array([30, 20.16, 10.12, 5*np.sqrt(2), 5, 1]) 
thickness_err   = np.array([0.02, 0.02, 0.02, 0.02*np.sqrt(2), 0.02, 0.02])

fitbook_lin = {
    "data_set"          : [thickness, thickness_err, log_Intens, log_Intens_err],
    "region_of_interest": [0,25],
    "fit_function"      : fitfunc_linear,
    "params_guess"      : None,
    "region_of_fit"     : None,
    "fit_density"       : 2
}

fitdata_lin, params_lin, params_lin_err, chi_lin = ultimate_fit( fitbook_lin )


def plot_lambert_beer():
    
    sample_format_dict_1 = {
        "label"      : r"Messdaten",          
        "fmt"        : 'o', 
        "color"      : sns.color_palette("dark")[0],                               
        "markersize" : 1, 
        "linewidth"  : 1,
        "capsize"    : 2,
        "alpha"      : 1                                   
    }
    sample_format_dict_2 = {
        "label"      : f"Linearer Fit $(\\chi^2 = {chi_lin:.1f})$",                       
        "fmt"        : '--', 
        "color"      : sns.color_palette("dark")[9],        
        "markersize" : 4, 
        "linewidth"  : 1,
        "capsize"    : 0,
        "alpha"      : 1
    }  

    writtings = {
        "title"       : None,
        "x_ax_label"  : r"Materialdicke $d$ / mm",
        "y_ax_label"  : r"Log. Intensität $\log{ \left(\frac{I_0}{I(d)} \right)}$ / 1"
    }
    
    general_format_dict = standard_format_dict.copy()
    zoom_params         = no_zooming.copy()
    colorbar_params     = no_colorbar.copy()
    extra_label         = no_extra_label.copy()
    extra_xaxis         = no_extra_xaxis.copy()
 
    exp_data            = [thickness, thickness_err, log_Intens, log_Intens_err]
    all_data                = [ exp_data, fitdata_lin ]                               
    all_sample_format_dicts = [ sample_format_dict_1, sample_format_dict_2 ]

    save_plot = True, "../Figures/Lambert-Beer_Fit.jpg"                                     
    ultimate_plot_advanced (all_data, writtings, zoom_params, colorbar_params, extra_label, extra_xaxis, save_plot, all_sample_format_dicts, general_format_dict)

plot_lambert_beer()


mu, mu_err  = params_lin[0], params_lin_err[0]
b, b_err    = params_lin[1], params_lin_err[1]

print("Linear Fit Params: log(I_0/I) = mu*x + b")
print(f"mu = ({mu:.4f}+-{mu_err:.4f}) / mm")
print(f"b = ({b:.6f}+-{b_err:.6f}) / 1")


# https://periodic-table.rsc.org/element/13/aluminium: 
#       density     = 2.70 g / cm^3 
#       atomic mass = 26,982  u    
#       ordnungszahl= 13
# (wiki) 1u         = 1.66053906892(52)e-27 kg = 1.66053906892(52)e-24 g

# https://www.spektrum.de/lexikon/physik/aluminium/416_
#       density     =  2,699  g / cm^3 
#       atomic mass = 26,9815  u    


density     = 2.699 
atomic_mass = 26.9815 *  1.66053907 * 10**(-24)
ordn_zahl   = 13
n_Al        = density / atomic_mass *10**6
n_e         = n_Al * ordn_zahl
print(f"Teilchendichte:   n_Al = {n_Al:.3e} / m^3")
print(f"Elektronendichte: n_e  = {n_e:.3e} / m^3")

sigma_tot       = 1000*mu / n_e    # mm ---> m 
sigma_tot_err   = 1000*mu_err/ n_e
   
print(f"sigma_tot = ({sigma_tot*10**28:.3e} +- {sigma_tot_err*10**28:.3e})  barn") # 1 barn = 10**(-28) m^2



alpha_EM    = 7.2973525643 * 10**(-3)
m_e         = 510.99895069 
E_Cs        = 661.7
h_bar_c     = 197.3 * 10**(3) * 10**(-15) # 1 in keV * m
def sigma_KN( E ):
    x = E/m_e
    f_1 = 2*x*(2+x*(1+x)*(8+x)) / (1+2*x)**2
    f_2 = ( (x-2)*x -2 )*np.log(1+2*x)
    return np.pi * alpha_EM**2 / (m_e**2 * x**3) * (f_1 + f_2) * h_bar_c**2

print(sigma_KN(E_Cs)*10**(28))