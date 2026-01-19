from MyPyLib_v4 import *


data = np.loadtxt("../Data/Protokoll_Driftstrom.txt")

R           = 10**6

U           = data[:, 0]                    # alles in V
U_err       = 0.01 * U
U_ohne      = data[:, 1] *10**(-3)
U_ohne_err  = data[:, 2] *10**(-3)
U_mit       = data[:, 3] *10**(-3)
U_mit_err   = data[:, 4] *10**(-3)

I_ohne      = U_ohne/R      *10**9          # alles in nA
I_ohne_err  = U_ohne_err/R  *10**9          
I_mit       = U_mit/R       *10**9
I_mit_err   = U_mit_err/R   *10**9

dataset_ohne    = [U, U_err, I_ohne, I_ohne_err ]
dataset_mit     = [U[:-2], U_err[:-2], I_mit[:-2], I_mit_err[:-2] ]

def fitfunc_exp(x, a, b, c):
    return a * np.exp(x/b) + c

fitbook = {
    "data_set"          : dataset_mit,
    "region_of_interest": None,
    "fit_function"      : fitfunc_exp,
    "params_guess"      : [1, 2000, 200],
    "region_of_fit"     : None,
    "fit_density"       : 400
}
fitdata, params, params_err, chi = ultimate_fit(fitbook)


def ultimate_plot_ohne( ):
    
    sample_format_dict_1 = {
        "label"      : r"Messwerte",          
        "fmt"        : '.', 
        "color"      : "black",                               
        "markersize" : 1, 
        "linewidth"  : 1,
        "capsize"    : 1.5,
        "alpha"      : 1                                   
    }
 

    writtings = {
        "title"       : None,
        "x_ax_label"  : r"Hochspannung $U_{dk}$ / V",
        "y_ax_label"  : r"Driftstrom $I$ / nA"
    }
    
    general_format_dict = standard_format_dict.copy()
    zoom_params         = no_zooming.copy()
    colorbar_params     = no_colorbar.copy()
    extra_label         = no_extra_label.copy()
    extra_xaxis         = no_extra_xaxis.copy()
    

    
    all_data                = [ dataset_ohne ]                               
    all_sample_format_dicts = [ sample_format_dict_1 ]

    save_plot = True, "../Figures/Driftstrom_ohne.jpg"                                     
    ultimate_plot_advanced (all_data, writtings, zoom_params, colorbar_params, extra_label, extra_xaxis, save_plot, all_sample_format_dicts, general_format_dict)

init_ultimate_plotting()
ultimate_plot_ohne()



def ultimate_plot_mit( ):
    
    sample_format_dict_1 = {
        "label"      : r"Messwerte",          
        "fmt"        : '.', 
        "color"      : "black",                               
        "markersize" : 1, 
        "linewidth"  : 1,
        "capsize"    : 1.5,
        "alpha"      : 1                                   
    }
    sample_format_dict_2 = {
        "label"      : f"Exp. Fit ($\\chi^2 = {chi:.1f}$)",          
        "fmt"        : '-', 
        "color"      : sns.color_palette("dark")[2],                               
        "markersize" : 4, 
        "linewidth"  : 4,
        "capsize"    : 0,
        "alpha"      : 0.5                                   
    }

    writtings = {
        "title"       : None,
        "x_ax_label"  : r"Hochspannung $U_{dk}$ / V",
        "y_ax_label"  : r"Driftstrom $I$ / nA"
    }
    
    general_format_dict = standard_format_dict.copy()
    zoom_params         = no_zooming.copy()
    colorbar_params     = no_colorbar.copy()
    extra_xaxis         = no_extra_xaxis.copy()

    t1 = f"a = ({params[0]*10**6:.1f} $\\pm$ {params_err[0]*10**6:.1f}) mA \n"
    t2 = f"b = ({params[1]:.0f} $\\pm$ {params_err[1]:.0f}) V \n"
    t3 = f"c = ({params[2]:.1f} $\\pm$ {params_err[2]:.1f}) nA "
    
    extra_label         = {
        "do_label"  :   True,
        "position"  :   [0.03, 0.8],
        "font_size" :   10,
        "content"   :   t1+t2+t3
    }
    

    
    all_data                = [ dataset_mit, fitdata ]                               
    all_sample_format_dicts = [ sample_format_dict_1, sample_format_dict_2 ]

    save_plot = True, "../Figures/Driftstrom_mit.jpg"                                     
    ultimate_plot_advanced (all_data, writtings, zoom_params, colorbar_params, extra_label, extra_xaxis, save_plot, all_sample_format_dicts, general_format_dict)

init_ultimate_plotting()
ultimate_plot_mit()

print("Exp. Fit:    f(x) = a * exp(x/b) + c ")
print(f"a = ({params[0]*10**6:.4f} +- {params_err[0]*10**6:.4f}) mA ")
print(f"b = ({params[1]:.4f} +- {params_err[1]:.4f}) V ")
print(f"c = ({params[2]:.4f} +- {params_err[2]:.4f}) nA ")





# ------------------------ Myonenfluss expected ------------------------

fluss       = 100

x, x_err    = 36.0 /100, 0.5 /100
y, y_err    = 4.8 /100, 0.2 /100

area, area_err  = x * y, np.sqrt( (x_err*y)**2 + (y_err*x)**2 )

n, n_err    = fluss * area, fluss * area_err

print(f"\nn_exp = ({n:.4f} +- {n_err:.4f}) 1/s ")
