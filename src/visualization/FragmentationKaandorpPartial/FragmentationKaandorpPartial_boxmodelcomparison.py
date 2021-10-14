import settings
import utils
import visualization.visualization_utils as vUtils
import matplotlib.pyplot as plt
import string
from Analysis.FragmentationKaandorp_boxmodel import FragmentationKaandorp_box_model
from copy import deepcopy


class FragmentationKaandorpPartial_boxmodelcomparison:
    def __init__(self, figure_direc, scenario, shore_time, lambda_frag, rho, sink=True, sim_length=2, month_step=4):
        # Data parameters
        self.output_direc = figure_direc + 'size_distribution/'
        self.data_direc = utils.get_output_directory(server=settings.SERVER) + 'size_distribution/FragmentationKaandorpPartial/'
        utils.check_direc_exist(self.output_direc)
        self.prefix = 'size_distribution'
        # Simulation parameters
        self.scenario = scenario
        self.shore_time = shore_time
        self.lambda_frag = lambda_frag
        self.rho = rho
        self.class_num = settings.SIZE_CLASS_NUMBER
        self.sink = sink
        if not self.sink:
            self.count, self.mass = 'particle_number', 'particle_mass'
        else:
            self.count, self.mass = 'particle_number_sink', 'particle_mass_sink'
        self.sim_length = sim_length
        self.month_step = month_step
        self.reservoir_list = ['total', 'beach', 'adrift']
        # Figure parameters
        self.fig_size = (14, 10)
        self.fig_shape = (3, 2)
        self.x_label = 'Size (mm)'
        self.y_label = r'Particle Number (n)'
        self.twiny_label = r'Particle Mass (g)'
        self.ax_ticklabel_size = 12
        self.ax_label_size = 14
        self.legend_size = 12
        self.xmin, self.xmax = 1e-3, 2e2
        self.ymin, self.ymax = 1e1, 1e5
        self.ax_range = self.xmax, self.xmin, self.ymax, self.ymin
        self.twin_ymin, self.twin_ymax = 1e-2, 1e5
        self.twin_ax_range = self.xmax, self.xmin, self.twin_ymax, self.twin_ymin
        self.number_of_plots = self.fig_shape[0] * self.fig_shape[1]
        self.cmap = 'viridis'

    def plot(self):
        # Getting the sizes of the size classes, and we convert from meters to mm
        size_classes = utils.size_range(size_class_number=self.class_num, units='mm')

        # Loading the parcels data
        base_dict = {self.count: {}, self.mass: {}}
        data_dict = {'total': deepcopy(base_dict), 'adrift': deepcopy(base_dict), 'beach': deepcopy(base_dict)}
        data = vUtils.FragmentationKaandorpPartial_load_data(scenario=self.scenario, prefix=self.prefix,
                                                             data_direc=self.data_direc, shore_time=self.shore_time,
                                                             lambda_frag=self.lambda_frag, rho=self.rho,
                                                             postprocess=True)
        time_indices = data['beach'][self.mass].keys()
        for time in time_indices:
            data_dict['total'][self.count][time] = data['adrift'][self.count][time] + data['beach'][self.count][time]
            data_dict['total'][self.mass][time] = data['adrift'][self.mass][time] + data['beach'][self.mass][time]
            for variable in ['adrift', 'beach']:
                data_dict[variable][self.count][time] = data[variable][self.count][time]
                data_dict[variable][self.mass][time] = data[variable][self.mass][time]

        # Loading the box model data
        box_model_data = FragmentationKaandorp_box_model(sim_length=self.sim_length, lambda_f=388).load_box_model(rerun=True)
        box_mass, box_number = box_model_data['mass'], box_model_data['number']
        box_time = box_model_data['mass'].keys()

        # Creating the figure
        ax, twin_ax = vUtils.base_figure(fig_size=self.fig_size, ax_range=self.ax_range, x_label=self.x_label,
                                         y_label=self.y_label, ax_ticklabel_size=self.ax_ticklabel_size,
                                         ax_label_size=self.ax_label_size, shape=self.fig_shape,
                                         plot_num=self.number_of_plots, log_yscale=True, log_xscale=True,
                                         all_x_labels=True, all_y_labels=True, add_twinx=True,
                                         twinx_y_label=self.twiny_label, twinx_ax_range=self.twin_ax_range,
                                         log_twinxscale=True, legend_axis=True, width_ratios=[1, 1, 0.5])

        # Labelling the subfigures
        for index_ax in range(self.number_of_plots):
            ax[index_ax].set_title(subfigure_title(index_ax), fontsize=self.ax_label_size)

        # Plotting the model distributions from the parcels analysis
        # subdivision_number = int(12 * self.sim_length // self.month_step)
        # for index_time, time in enumerate(time_indices):
        #     if index_time % self.month_step == 0 and index_time <= 12 * self.sim_length:
        #         c = vUtils.discrete_color_from_cmap(index=index_time // self.month_step, subdivisions=subdivision_number, cmap=self.cmap)
        #         ax[0].plot(size_classes, data_dict['total'][self.count][time], linestyle='-', c=c, label='Month {}'.format(index_time))
        #         twin_ax[1].plot(size_classes, data_dict['total'][self.mass][time], linestyle='-', c=c)
        #
        # for index_time, time in enumerate(time_indices):
        #     if index_time % self.month_step == 0 and index_time <= 12 * self.sim_length:
        #         c = vUtils.discrete_color_from_cmap(index=index_time // self.month_step, subdivisions=subdivision_number, cmap=self.cmap)
        #         ax[2].plot(size_classes, data_dict['beach'][self.count][time], linestyle='-', c=c)
        #         twin_ax[3].plot(size_classes, data_dict['beach'][self.mass][time], linestyle='-', c=c)
        #
        # for index_time, time in enumerate(time_indices):
        #     if index_time % self.month_step == 0 and index_time <= 12 * self.sim_length:
        #         c = vUtils.discrete_color_from_cmap(index=index_time // self.month_step, subdivisions=subdivision_number, cmap=self.cmap)
        #         ax[4].plot(size_classes, data_dict['adrift'][self.count][time], linestyle='-', c=c)
        #         twin_ax[5].plot(size_classes, data_dict['adrift'][self.mass][time], linestyle='-', c=c)

        lines, labels = ax[0].get_legend_handles_labels()

        # Plotting the model distributions from the box model
        step_div = 4 * self.month_step
        subdivision_number = int(52 * self.sim_length // step_div)
        for index_time, time in enumerate(box_time):
            if index_time % (4 * self.month_step) == 0 and index_time <= 52 * self.sim_length and type(time) == int:
                c = vUtils.discrete_color_from_cmap(index=index_time // step_div, subdivisions=subdivision_number, cmap=self.cmap)
                ax[0].plot(size_classes, box_number[time]['total'], linestyle='--', c=c)
                twin_ax[1].plot(size_classes, box_mass[time]['total'], linestyle='--', c=c)

        for index_time, time in enumerate(box_time):
            if index_time % 4 * self.month_step == 0 and index_time <= 52 * self.sim_length and type(time) == int:
                c = vUtils.discrete_color_from_cmap(index=index_time // step_div, subdivisions=subdivision_number, cmap=self.cmap)
                ax[2].plot(size_classes, box_number[time]['beach'], linestyle='--', c=c)
                twin_ax[3].plot(size_classes, box_mass[time]['beach'], linestyle='--', c=c)

        for index_time, time in enumerate(box_time):
            if index_time % 4 * self.month_step == 0 and index_time <= 52 * self.sim_length and type(time) == int:
                c = vUtils.discrete_color_from_cmap(index=index_time // step_div, subdivisions=subdivision_number, cmap=self.cmap)
                ax[4].plot(size_classes, box_number[time]['ocean'] + box_number[time]['coastal'], linestyle='--', c=c)
                twin_ax[5].plot(size_classes, box_mass[time]['ocean'] + box_mass[time]['coastal'], linestyle='--', c=c)

        # Adding a legend
        subdivision_number = int(12 * self.sim_length // self.month_step)
        line_colors = [plt.plot([], [], c=vUtils.discrete_color_from_cmap(ind, subdivisions=subdivision_number, cmap=self.cmap),
                                label=label, linestyle='-')[0] for ind, label in enumerate(labels)]
        ax[-1].legend(handles=line_colors, fontsize=self.legend_size, loc='upper right')

        # Saving the figure
        str_format = self.shore_time, self.rho, self.lambda_frag, self.sink
        file_name = self.output_direc + 'boxmodel_comparison-ST={}-rho={}-lambda_f={}_sink={}.png'.format(*str_format)
        plt.savefig(file_name, bbox_inches='tight', dpi=300)


def subfigure_title(index):
    alphabet = string.ascii_lowercase
    subtitle_list = ['Counts', 'Mass']
    reservoir_list = ['Total', 'Beach', 'Adrift']
    return '({}) {} - {}'.format(alphabet[index], reservoir_list[index // 2], subtitle_list[index % 2])




