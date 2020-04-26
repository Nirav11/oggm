import os
import geopandas as gpd
import oggm
from oggm import cfg, utils, tasks, graphics, workflow, tasks
from oggm.utils import get_demo_file
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm as colormap
import xarray as xr


try:
    import salemf
except ImportError:
    pass


# Read the default parameter file and make them available to other OGGM tools via cfg.PARAMS
cfg.initialize(logging_level='WORKFLOW')
cfg.PARAMS['prcp_scaling_factor'], cfg.PARAMS['ice_density'], cfg.PARAMS['continue_on_error']
cfg.PATHS['working_dir'] = utils.gettempdir(dirname='/home/nirab/OGGM/OGGM-GettingStarted', reset=True)

rgi_ids = ['RGI60-11.01328', 'RGI60-11.00897']

#gdirs = workflow.init_glacier_directories(rgi_ids, from_prepro_level=2, prepro_border=80)
#graphics.plot_domain(gdirs, figsize=(8, 7))
#gdir = rgi_ids[1]
gdir= []
# Plotting the glacier over the DEM file
for rgi in rgi_ids:

    f, axs = plt.subplots(2, 2, figsize=(8, 6))
    for ax, border in zip(np.array(axs).flatten(), [10, 80, 160, 250]):
        gdirs = workflow.init_glacier_directories(rgi,
                                                    from_prepro_level=1,prepro_border=border)

        graphics.plot_domain(gdirs, ax=ax, title='Border: {}'.format(border),
                                add_colorbar=False,
                                lonlat_contours_kwargs={'add_tick_labels':False})
    gdir = gdir + gdirs                                           
    print("dir==",gdir)
plt.show()
# run the glacier_masks task on all gdirs
#workflow.execute_entity_task(tasks.glacier_masks, gdirs)

#gdir = gdirs[1]
#print('Path to the DEM:', gdir.get_filepath('dem'))


#    print('Path to the masks:', gdir.get_filepath('gridded_data'))  
# Path to the masks: /home/nirab/OGGM/OGGM-GettingStarted/per_glacier/RGI60-11/RGI60-11.00/RGI60-11.00897/gridded_data.nc

# Apply several tasks sequentially (i.e. one after an other) on our glacier list
""" list_talks = [
         tasks.compute_centerlines, #Compute the centerlines following Kienholz et al., (2014).
         tasks.initialize_flowlines, #Computes more physical Inversion Flowlines from geometrical Centerlines
         tasks.compute_downstream_line, #Computes the Flowline along the unglaciated downstream topography
         tasks.catchment_area, #Compute the catchment areas of each tributary line.
         tasks.catchment_width_geom, #Compute geometrical catchment widths for each point of the flowlines
         tasks.catchment_width_correction, #Corrects for NaNs and inconsistencies in the geometrical widths.
         tasks.compute_downstream_bedshape #The bedshape obtained by fitting a parabola to the line's normals and downstream altitude
         ]
 
for task in list_talks:
        workflow.execute_entity_task(task, gdirs)

for agdir in gdirs:
    graphics.plot_centerlines(agdir, figsize=(8, 7), use_flowlines=True, add_downstream=True)
    graphics.plot_catchment_areas(agdir, figsize=(8, 7))
    graphics.plot_catchment_width(agdir, corrected=True, figsize=(8, 7))

# Location of Monthly Climate Data for the Glacier
#fpath = gdir.get_filepath('climate_monthly')
#print(fpath)

#ds = xr.open_dataset(fpath)
#print(ds)
# Data is in hydrological years
# -> let's just ignore the first and last calendar years
#ds.temp.resample(time='AS').mean()[1:-1].plot()

plt.show() 

#workflow.execute_entity_task(tasks.local_t_star, gdirs);
#workflow.execute_entity_task(tasks.mu_star_calibration, gdirs); 

#Each flowline now knows what area will contribute to its surface mass-balance and ice flow. 
#gdir = gdirs[0]
# Accordingly, it is possible to compute each glacier cross-section's width, and correct it
#  so that the total glacier area and elevation distribution is conserved


 """