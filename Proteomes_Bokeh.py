# %% IMPORTING PACKAGES =======================================================
# =============================================================================

# Importing supporting packages
import os
import sys
import math

# Importing main packages
import pandas as pd
import numpy as np
import seaborn as sns

# Importing Bokeh packages
from bokeh.plotting import figure
from bokeh.io import curdoc
from bokeh.models import ColumnDataSource, CategoricalColorMapper, HoverTool
from bokeh.models import Div
from bokeh.models.widgets import Spinner, Toggle, CheckboxButtonGroup, Button
from bokeh.models.widgets import Panel, Tabs
from bokeh.layouts import row, column

# Importing pyProtamine UDFs
import Proteomes_engine as sperm

# Importing Scikit Learn packages for incomming DBSCAN clustering
from sklearn.preprocessing import StandardScaler

# %% IMPORTING INPUT DATA =====================================================
# =============================================================================

# Importing Project, Experimental Design dictionaries
p = np.load(f'{os.getcwd()}\\Output\\Project_dict.npy', allow_pickle='TRUE').item()
e = np.load(f'{os.getcwd()}\\Output\\Experiment_dict.npy', allow_pickle='TRUE').item()

# Import the pre-DBSCAN clustering DataFrame
df = pd.read_excel(f'{os.getcwd()}\\Output\\Table S1.xlsx', engine='openpyxl')

# Adding "empty" 'Cluster' column to the pre-DBSCAN DataFrame
df['Cluster'] = 'Unclustered'

# Getting the scaling constant (to start the plot  with a nice circle radius)
scaler = StandardScaler()
scaler.fit_transform(df[['m (Da)']])
scale = scaler.scale_[0]

# %% PREPARING TOOLTIP ========================================================
# =============================================================================

# Preparing the tooltips in a list of tuples
tooltips_list = [
                 ('m (Da)',     '@{m (Da)}{int}'),
                 ('RT (min)',   '@{RT (min)}{0.0}'),
                 ('z log2 Int', '@{z log2 Int}{0.0}'),
                 ('Cluster',    '@{Cluster}'),
                 ('Node',     '@{Node}'),
                 ]

# Creating a hover tool legend
ht = HoverTool(tooltips=tooltips_list)

del tooltips_list

# %% PREPARING PALETTE AND CCM ================================================
# =============================================================================

# Listing unique items from 'Cluster' to get sorted list of factors (hue)
cluster_list = sorted(list(df['Cluster'].unique()))

# Defining our color palette based on the number of clusters
pal = sns.color_palette('hls', len(cluster_list))  # 'husl'

# HEX-converting the palette to make it suitable for Bokeh
pal = list(pal.as_hex())

# Creating the Categorical Color Mapper for 'Cluster'
ccmapper = CategoricalColorMapper(factors=cluster_list, palette=pal)

del pal, cluster_list

# %% INITIATING COLUMN DATA SOURCE AND BOKEH APP ==============================
# =============================================================================

# Creating our Column Data Source from the input DatFrame
cds = ColumnDataSource(df)

# Cleaning current document to avoid app concatenation within the *.html file
curdoc().clear()

# %% PREPARING BOKEH FIGURES ==================================================
# =============================================================================

# Initiating empty dictionary to store unpcomming Bokeh figures
figures_dict = {}


# Defining supporting UDF to shrink code
def bokeh_figure(fig_name, x, y, option):

    # TODO
    cs = len(df.loc[df["Cluster"] != "Unclustered", "Cluster"].unique())
    t = len(df["Cluster"])
    c = sum(df["Cluster"] != "Unclustered")
    uc = sum(df["Cluster"] == "Unclustered")

    # Initiating Bokeh figure
    fig = figure(x_axis_label=x,
                 y_axis_label=y,
                 tools='pan,box_zoom,wheel_zoom,box_select,reset,save',
                 plot_width=1900, plot_height=800,
                 title=f'Clusters: {cs} - PrSMs: {t} - Clustered: {c} ({100*c/t:.1f}%) - Unclustered: {uc} ({100*uc/t:.1f}%)',
                 )

    # Adding the hover legend tool to the Bokeh figure
    fig.add_tools(ht)

    # Using radius= as circle "size" if the input parameter option='radius'
    if option == 'radius':

        # Adding circle glyph to the Bokeh figure (using radius=)
        circle_glyph = fig.circle(
                                  x=x, y=y, source=cds, alpha=0.9,
                                  radius=(p['ε'] * scale) / 2,
                                  color={'field': 'Cluster',
                                         'transform': ccmapper},
                                  )

    # Using size= as circle "size" if the input parameter option='size'
    if option == 'size':

        # Adding circle glyph to the Bokeh figure (using size=)
        circle_glyph = fig.circle(
                                  x=x, y=y, source=cds, alpha=0.9,
                                  size=2,
                                  color={'field': 'Cluster',
                                         'transform': ccmapper},
                                  )

    # NOTICE: Usually we just do fig.circle, but since we need to acces the
    # circle glyphs below, we must store it as circle_glyph

    # Storing the Bokeh figure and its circle glyphs as a tuple
    figures_dict[fig_name] = (fig, circle_glyph)


# Generating our Bokeh figures batch
bokeh_figure(fig_name='f1', x='m (Da)', y='z log2 Int', option='radius')
bokeh_figure(fig_name='f2', x='m (Da)', y='RT (min)', option='radius')
bokeh_figure(fig_name='f3', x='RT (min)', y='z log2 Int', option='size')

# %% PREPARING PANELS AND TABS ================================================
# =============================================================================

# Seetting one Panel with a single Figure in each Tab
tab1 = Panel(child=figures_dict['f1'][0], title='Intensity vs. Mass')
tab2 = Panel(child=figures_dict['f2'][0], title='Retention Time vs. Mass')
tab3 = Panel(child=figures_dict['f3'][0], title='Intensity vs. Retention Time')

# %% INITIATING BOKEH WIDGETS =================================================
# =============================================================================

# Getting list of unique Nodes
nodes = list(df['Node'].unique())

# Creating Check Box Button Group Widget for Nodes
node_cbbg = CheckboxButtonGroup(labels=nodes,
                                active=[nodes.index(i) for i in nodes])

# Creating Toggle Button Widget for show/hide unclustered PrSMs
unclustered_toggle = Toggle(label='Hide unclustered')

# Creating Spinner Widget for the DBSCAN eps hyperparameter
ε_spinner = Spinner(title='ε', width=100, value=p['ε'],
                    low=0.0005, high=0.0050, step=0.0001)

# Creating Spinner Widget for the DBSCAN min_samples hyperparameter
n_spinner = Spinner(title='n_min', width=100, value=p['n_min'],
                    low=1, high=50, step=1)

# Getting RT floor & ceil for low= and high= Spinner parameters
t_floor = math.floor(min(df['RT (min)']))
t_ceil = math.ceil(max(df['RT (min)']))

# Creating Spinner Widget for low and high RT filtering
tmin_spinner = Spinner(title='Min. RT', width=100,
                       value=p['rt_lower_lim'],
                       low=t_floor, high=t_ceil, step=1)
tmax_spinner = Spinner(title='Max. RT', width=100,
                       value=p['rt_upper_lim'],
                       low=t_floor, high=t_ceil, step=1)

# Getting abundance floor & ceil for low= and high= Spinner parameters
i_floor = math.floor(min(df['z log2 Int']))
i_ceil = math.ceil(max(df['z log2 Int']))

# Creating Spinner Widget for low and high abundance filtering
imin_spinner = Spinner(title='Min. I', width=100,
                       value=p['z-i_lower_lim'],
                       low=i_floor, high=i_ceil, step=0.1)
imax_spinner = Spinner(title='Max. I', width=100,
                       value=p['z-i_upper_lim'],
                       low=i_floor, high=i_ceil, step=0.1)

# Button to stop the Bokeh server
stop_toggle = Toggle(label="Online", button_type="success")

del t_floor, t_ceil, i_floor, i_ceil

# %% INITIATING BOKEH CALLBACKS ===============================================
# =============================================================================


# Initiating main Bokeh Call Back
def cb(attr, old, new):

    # Getting values of each Widget
    ε, n = ε_spinner.value, n_spinner.value
    t_min, t_max, = tmin_spinner.value, tmax_spinner.value
    i_min, i_max, = imin_spinner.value, imax_spinner.value

    # Getting active nodes within the Check Box Button Group Widget
    active_list = [node_cbbg.labels[i] for i in node_cbbg.active]

    # Duplicating the input df in order to get the "Callback" df
    df_cb = df.copy()

    # TODO
    m0 = (df_cb['RT (min)'] >= t_min) & (df_cb['RT (min)'] <= t_max)
    m2 = (df_cb['z log2 Int'] >= i_min) & (df_cb['z log2 Int'] <= i_max)

    # Masking the DataFrame looking for "active" features
    m3 = df_cb['Node'].isin(active_list)

    # Filtering DataFrame
    df_cb = df_cb[m0 & m2 & m3].copy()

    # Calling DBSCAN_clustering UDF on the "Callback" DataFrame
    sperm.DBSCAN_clustering(data=df_cb, eps=ε, min_samples=n)

    # If the 'Hide unclustered' toggle button is active...
    if unclustered_toggle.active:

        # ... change the toggle button tag
        unclustered_toggle.label = "Show unclustered"

        # Masking DataFrame looking for selected window
        m0 = df_cb['Cluster'] != 'Unclustered'

        # Filtering DataFrame
        df_cb = df_cb[m0].copy()

    # If the 'Hide unclustered' toggle button is not active...
    if not unclustered_toggle.active:

        # ... do not change the toggle button tag
        unclustered_toggle.label = "Hide unclustered"
    # Listing unique items from 'Cluster' to get sorted list of factors (hue)
    factors_cb = sorted(list(df_cb['Cluster'].unique()))

    # Defining our color palette based on the number of clusters
    pal_cb = sns.color_palette('hls', len(factors_cb))  # 'husl'

    # HEX-converting the palette to make it suitable for Bokeh
    pal_cb = list(pal_cb.as_hex())

    # Creating the Categorical Color Mapper for the filtered-in features
    ccm_cb = CategoricalColorMapper(palette=pal_cb, factors=factors_cb)

    # Updating .palette and .factors attributes of the Categorical Color Mapper
    ccmapper.palette = ccm_cb.palette
    ccmapper.factors = ccm_cb.factors

    # Creating the (new) column data source
    cds_cb = ColumnDataSource(df_cb)

    # Updating the (old) column data source
    cds.data.update(cds_cb.data)

    # TODO
    cs = len(df_cb.loc[df_cb["Cluster"] != "Unclustered", "Cluster"].unique())
    t = len(df_cb["Cluster"])
    c = sum(df_cb["Cluster"] != "Unclustered")
    uc = sum(df_cb["Cluster"] == "Unclustered")

    # Creating a nice title reporting the number of clusters
    figures_dict['f1'][0].title.text = f'Clusters: {cs} - PrSMs: {t} - Clustered: {c} ({100*c/t:.1f}%) - Unclustered: {uc} ({100*uc/t:.1f}%)'
    figures_dict['f2'][0].title.text = f'Clusters: {cs} - PrSMs: {t} - Clustered: {c} ({100*c/t:.1f}%) - Unclustered: {uc} ({100*uc/t:.1f}%)'
    figures_dict['f3'][0].title.text = f'Clusters: {cs} - PrSMs: {t} - Clustered: {c} ({100*c/t:.1f}%) - Unclustered: {uc} ({100*uc/t:.1f}%)'

    # Updating circle glyphs radius according the new ε hyperparameter value
    figures_dict['f1'][1].glyph.radius = (ε * scaler.scale_[0]) / 2
    figures_dict['f2'][1].glyph.radius = (ε * scaler.scale_[0]) / 2
    # figures_dict['f3'][1].glyph.radius = (ε * scaler.scale_[0]) / 2


# TODO
def stop_cb(attr, old, new):

    # If the 'Hide unclustered' toggle button is active...
    if stop_toggle.active:

        # ... change the toggle button tag
        stop_toggle.label = "Click again to stop Bokeh Server"
        stop_toggle.button_type = "warning"

    # If the 'Hide unclustered' toggle button is not active...
    if not stop_toggle.active:

        # ... do not change the toggle button tag
        stop_toggle.label = "Online"
        stop_toggle.button_type = "success"

        sys.exit()  # Stop the server


# %% ATTACHING BOKEH CALLBACKS ================================================
# =============================================================================

# Attaching CallBack to the Toggle Button a Checkbox Button Group widgets
for wg in [stop_toggle]:
    # wg.on_change('active', lambda attr, old, new: select_movies())
    wg.on_change('active', stop_cb)
del wg

# Attaching CallBack to the Toggle Button a Checkbox Button Group widgets
for wg in [unclustered_toggle, node_cbbg]:
    # wg.on_change('active', lambda attr, old, new: select_movies())
    wg.on_change('active', cb)
del wg

# Attaching CallBack to the Spinner widgets
for wg in [ε_spinner, n_spinner, tmin_spinner, tmax_spinner, imin_spinner, imax_spinner]:
    wg.on_change('value', cb)
del wg

# %% CREATING BOKEH LAYOUT ====================================================
# =============================================================================

# Creating a layout combining the slider and the figure
layout = column(
                row(unclustered_toggle, Div(text='Node: '), node_cbbg, Div(text='Stop Bokeh Server: '), stop_toggle),
                row(Div(text='Filter by retention time: '), tmin_spinner, tmax_spinner, Div(text='Filter by intensity: '), imin_spinner, imax_spinner, Div(text='DBSCAN parameters: '), ε_spinner, n_spinner),
                row(Tabs(tabs=[tab1, tab2, tab3]))
                )

# %% CONCLUDING ===============================================================
# =============================================================================

# Adding our layout and a title to our current document
curdoc().add_root(layout)
curdoc().title = 'Proteomes Journal'
