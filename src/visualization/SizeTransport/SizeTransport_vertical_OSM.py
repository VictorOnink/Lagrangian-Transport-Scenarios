import settings
import utils
import visualization.visualization_utils as vUtils
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from advection_scenarios import advection_files
import numpy as np
import string
import cmocean.cm as cmo
from Analysis.CMEMS_mediterranean_mean_MLD import CMEMS_mediterranean_mean_MLD


class SizeTransport_vertical_OSM:
    def __init__(self, scenario, figure_direc, time_selection, season, rho_list=[920], with_mld=True,
                 shore='all', fixed_resus=False, resus_time=7):
        # Figure Parameters
        self.fig_size = (16, 10)
        self.fig_shape = (1, 1)
        self.x_label = 'Normalized Concentration'
        self.y_label = 'Depth (m)'
        self.ax_ticklabel_size = 16
        self.ax_label_size = 18
        self.legend_size = 16
        self.title_size = 24
        self.line_width = 3
        self.xmin, self.xmax = 1e-7, 1e0 + 0.1
        self.ymin, self.ymax = 3e3, 1e0
        self.ax_range = self.xmax, self.xmin, self.ymax, self.ymin
        self.number_of_plots = 1
        self.rho_line_dict = {30: 'dashed', 920: '-', 980: 'dashed', 1020: '-'}
        self.season_list = ['JFM', 'AMJ', 'JAS', 'OND']
        self.with_mld = with_mld
        # Data parameters
        self.output_direc = figure_direc + 'vertical_profile/'
        self.data_direc = settings.DATA_OUTPUT_DIREC + 'concentrations/SizeTransport/'
        utils.check_direc_exist(self.output_direc)
        self.prefix = 'vertical_concentration'
        # Simulation parameters
        self.scenario = scenario
        self.size_list = np.array([5000, 2500, 1250, 625, 313, 156, 78, 39, 20, 10, 5, 2]) * settings.SIZE_FACTOR
        self.time_selection = time_selection
        self.season = season
        self.year = 2010 + self.time_selection
        self.rho_list = rho_list
        self.tau = 0.0
        self.shore = shore
        self.fixed_resus = fixed_resus
        self.resus_time = resus_time

    def plot(self):
        # Loading the data
        output_dict = {}
        for rho in self.rho_list:
            output_dict[rho] = {}
            for size in self.size_list:
                data_dict = vUtils.SizeTransport_load_data(scenario=self.scenario, prefix=self.prefix,
                                                           data_direc=self.data_direc, fixed_resus=self.fixed_resus,
                                                           size=size, rho=rho, tau=self.tau, resus_time=self.resus_time)
                output_dict[rho][size] = data_dict[utils.analysis_simulation_year_key(self.time_selection)]
        depth_bins = data_dict['depth']

        # Averaging by season
        conc_type = {'all': 'concentration', 'offshore': 'concentration_offshore',
                     'nearshore': 'concentration_nearshore'}[self.shore]
        for rho in self.rho_list:
            for size in self.size_list:
                for month in np.arange(0, 12, 3):
                    month_stack = np.vstack([output_dict[rho][size][month][conc_type],
                                             output_dict[rho][size][month + 1][conc_type],
                                             output_dict[rho][size][month + 2][conc_type]])
                    output_dict[rho][size][month][conc_type] = np.nanmean(month_stack, axis=0)

        # Normalizing the profiles
        for rho in self.rho_list:
            for size in self.size_list:
                for month in range(0, 12):
                    output_dict[rho][size][month][conc_type] /= np.sum(output_dict[rho][size][month][conc_type])

        # setting the zero values to nan to clear up the plots
        for rho in self.rho_list:
            for size in self.size_list:
                for month in range(0, 12):
                    selection = output_dict[rho][size][month][conc_type] == 0
                    output_dict[rho][size][month][conc_type][selection] = np.nan

        # Get the mean MLD depth data
        if self.with_mld:
            MLD_data = CMEMS_mediterranean_mean_MLD().calculate_seasonal_mean()[self.year]
            MLD_mean = dict.fromkeys(self.season_list)
            for season in self.season_list:
                MLD_mean[season] = np.nanmean(MLD_data[season])

        # Creating the figure
        ax = vUtils.base_figure(fig_size=self.fig_size, ax_range=self.ax_range, x_label=self.x_label,
                                y_label=self.y_label, ax_ticklabel_size=self.ax_ticklabel_size,
                                ax_label_size=self.ax_label_size, shape=self.fig_shape, plot_num=self.number_of_plots,
                                log_yscale=True, log_xscale=True, all_x_labels=True, all_y_labels=True,
                                legend_axis=True, width_ratios=[1, 0.3])

        # Labelling the subfigures
        for index_ax in range(self.number_of_plots):
            ax[index_ax].set_title(self.season, fontsize=self.title_size, fontweight='bold')

        # Adding in a legend
        cmap_list, label_list = [], []
        for index_size, size in enumerate(self.size_list):
            cmap_list.append(vUtils.discrete_color_from_cmap(index_size, subdivisions=self.size_list.__len__()))
            label_list.append(self.legend_label(size))
        size_colors = [plt.plot([], [], c=cmap_list[i], label=label_list[i], linestyle='-')[0] for i in
                       range(cmap_list.__len__())]
        rho_lines = [plt.plot([], [], c='k', label=r'$\rho=$' + str(rho) + r' kg m$^{-3}$', linestyle=self.rho_line_dict[rho])[0]
                     for rho in self.rho_list]
        if self.with_mld:
            mld_line = [plt.plot([], [], c='r', label=r'Mean MLD', linestyle='-')[0]]
            ax[-1].legend(handles=rho_lines + mld_line + size_colors, fontsize=self.legend_size)
        else:
            ax[-1].legend(handles=rho_lines + size_colors, fontsize=self.legend_size)
        ax[-1].axis('off')

        # The actual plotting
        if self.season == 'Winter':
            selection_name, selection_ind = 'JFM', 0
        elif self.season == 'Summer':
            selection_name, selection_ind = 'JAS', 5
        for rho in self.rho_list:
            for index_size, size in enumerate(self.size_list):
                c = cmap_list[index_size]
                if np.sum(~np.isnan(output_dict[rho][size][selection_ind][conc_type])) < 2:
                    linestyle, markerstyle = None, 'o'
                else:
                    linestyle, markerstyle = self.rho_line_dict[rho], None
                ax[0].plot(output_dict[rho][size][selection_ind][conc_type], depth_bins, linestyle=linestyle,
                           c=c, marker=markerstyle, linewidth=self.line_width)
                # Adding in a horizontal line for the MLD
                if self.with_mld:
                    ax[0].axhline(y=MLD_mean[selection_name], color='r', linestyle='-', linewidth=self.line_width)

        plt.savefig(self.file_name(), bbox_inches='tight')
        plt.close('all')

    def file_name(self, file_type='.png'):
        base = 'OSM_season_{}_year_{}'.format(self.season, self.time_selection)
        return self.output_direc + base + file_type

    @staticmethod
    def legend_label(size):
        return r'd = {:.3f} mm'.format(size * 1e3)
