from MyPyLib_v2 import *



# ---------------------------- Plot Counts - Bg for 35 degree ----------------------------


data_deg    = np.loadtxt("../Data/10_Caesium_35grad.txt", skiprows=0)
data_bg     = np.loadtxt("../Data/11_Caesium_35grad_Hintergrund.txt", skiprows=0)
channel     = data_deg[ : , 0 ] 
counts      = data_deg[ : , 1 ] - data_bg[ : , 1]
counts_err  = np.sqrt(   np.sqrt(data_deg[ : , 1 ])**2  +  np.sqrt(data_bg[ : , 1 ])**2 ) + 1   # !!!!!!!!!!!!




def ultimate_plot():
    
    sample_format_dict_1 = {
        "label"      : r"Messwerte 1",          
        "fmt"        : '.', 
        "color"      : sns.color_palette("dark")[1],                               
        "markersize" : 4, 
        "linewidth"  : 1,
        "capsize"    : 0,
        "alpha"      : 1                                   
    }

    writtings = {
        "title"       : r"Title",
        "x_ax_label"  : r"Time $t$ [$s$]",
        "y_ax_label"  : r"Angle $\theta$ [$rad$]"
    }
    
    general_format_dict = standard_format_dict.copy()
    zoom_params         = no_zooming.copy()
    colorbar_params     = no_colorbar.copy()
    extra_label         = no_extra_label.copy()
    extra_xaxis         = no_extra_xaxis.copy()
    
    data_set_1  = channel,     None,   counts,   None       
    
    all_data                = [ data_set_1 ]                               
    all_sample_format_dicts = [ sample_format_dict_1 ]

    save_plot = True, "../Figures/text_30grad.jpg"                                     
    ultimate_plot_advanced (all_data, writtings, zoom_params, colorbar_params, extra_label, extra_xaxis, save_plot, all_sample_format_dicts, general_format_dict)

ultimate_plot()
ultimate_plot()