#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plotting routines for 2D and 3D look directions of the instrument onboard 
spacecraft.

.. toctree::
  :maxdepth: 2   
  attitude_2d_altitude
  attitude_2d_on_map
  attitude_3d_connect_to_subpoint
  attitude_3d_connect_to_target
  attitude_3d_ground_quiver
  display_observation_geometry
  earth_radius_at_latitude
  fov_plotter
  plot_attitude_accuracy
  plot_slew_rri
  trajectory_plotter_2d_map

@author: ceren
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as dtformat
import cartopy.crs as ccrs
import cartopy
from matplotlib.collections import PolyCollection
from cartopy.mpl.patch import geos_to_path
from matplotlib.ticker import MultipleLocator
from mpl_toolkits.mplot3d import axes3d
import matplotlib.patches as mpatches
import matplotlib.ticker as ticker
import matplotlib.cm as cm
import itertools
from matplotlib import rcParams
import matplotlib.gridspec as gridspec
from matplotlib.offsetbox import AnchoredText
from matplotlib.cm import ScalarMappable
from matplotlib.ticker import FixedLocator, FixedFormatter
from matplotlib.legend_handler import HandlerBase

import cavsiopy.miscellaneous as misc
import cavsiopy.ephemeris_importer as ei

class MyAxes3D(axes3d.Axes3D):
    """ Class to draw 3D grids for the counter-clockwise angles of 3D plots."""
    def __init__(self, baseObject, sides_to_draw):
        self.__class__ = type(baseObject.__class__.__name__,
                              (self.__class__, baseObject.__class__),
                              {})
        self.__dict__ = baseObject.__dict__
        self.sides_to_draw = list(sides_to_draw)
        self.mouse_init()

    def set_some_features_visibility(self, visible):
        for t in self.w_zaxis.get_ticklines() + self.w_zaxis.get_ticklabels():
            t.set_visible(visible)
        self.w_zaxis.line.set_visible(visible)
        self.w_zaxis.pane.set_visible(visible)
        self.w_zaxis.label.set_visible(visible)

    def draw(self, renderer):
        # set visibility of some features False
        self.set_some_features_visibility(False)
        # draw the axes
        super(MyAxes3D, self).draw(renderer)
        # set visibility of some features True.
        # This could be adapted to set your features to desired visibility,
        # e.g. storing the previous values and restoring the values
        self.set_some_features_visibility(True)

        zaxis = self.zaxis
        draw_grid_old = zaxis.axes._draw_grid
        # disable draw grid
        zaxis.axes._draw_grid = False

        tmp_planes = zaxis._PLANES

        if 'l' in self.sides_to_draw :
            # draw zaxis on the left side
            zaxis._PLANES = (tmp_planes[2], tmp_planes[3],
                             tmp_planes[0], tmp_planes[1],
                             tmp_planes[4], tmp_planes[5])
            zaxis.draw(renderer)
        if 'r' in self.sides_to_draw :
            # draw zaxis on the right side
            zaxis._PLANES = (tmp_planes[3], tmp_planes[2],
                             tmp_planes[1], tmp_planes[0],
                             tmp_planes[4], tmp_planes[5])
            zaxis.draw(renderer)

        zaxis._PLANES = tmp_planes

        # disable draw grid
        zaxis.axes._draw_grid = draw_grid_old
        
def display_observation_geometry(title, time_array, Px, Py, Pz, 
                                 x, y, z, target_name, *V, step = 60, 
                                 vc=['red', '#12e193', 'red', '#12e193', \
                                     'black', 'blue'], 
                                 ls=['solid', 'solid', 'solid', 'solid', \
                                     'solid', 'solid'], 
                                 lw = [1.5, 1.5, 1.5, 1.5, 1.5, 1.5], 
                                 arrlen = [.5, .5, .5, .5, .75, 1],
                                 labels=['$Dipole_{1}$','$Dipole_{2}$',\
                                         'Boresight','$Ray_{LOS}$'], 
                                 markers=['_','_',  r'$\longrightarrow$',\
                                          r'$\longrightarrow$'], 
                                 colors=['red', '#12e193', 'black', 'blue'], 
                                 edgecolors=['red', '#12e193', 'black','blue'],
                                 arrowhead = [0.01, 0.01, 0.01, 0.01, \
                                              0.45, 0.45], 
                                 sct_kwargs= {'alpha': 1,'edgecolor': 'black',\
                                              'c': 'lightgrey', 'marker': '*',\
                                                  's': 180}, 
                                 loc = 'upper center'):
    """
    Plotting function to display observation geometry.
    The defaults for the keywords are based on RRI plots.

    Parameters
    ----------
    title : str
        plot title.
    time_array : datetime.datetime
        experiment time interval as datetime array.
    Px : float
        geodetic longitude of the target.
    Py : float
        geodetic latitude of the target.
    Pz : float
        altitude of the target.
    x : numpy.ndarray[float]
        spacecraft longitude(degrees).
    y : numpy.ndarray[float]
        spacecraft latitude (degrees).
    z : numpy.ndarray[float]
        spacecraft altitude (km).
    target_name : str
        target name.
    *V : numpy.ndarray[float]
        vectors to be plotted (need to be in ENU coordinate system).
        vec_args = M1_enu, M3_enu, M2_enu, M4_enu, RRI_enu, los_enu_arr
    step : int
        time between the vectors. The default is 60. (60 seconds)
    vc : list, optional
        vector colors. The default is ['red', '#12e193', 'red', '#12e193',\
                                       'black', 'blue'].
    ls : list, optional
        linestyles for vectors. The default is ['solid', 'solid', 'solid',\
                                                'solid', 'solid', 'solid'].
    lw : list, optional
        linewidth of vectors. The default is [1.5, 1.5, 1.5, 1.5, 1.5, 1.5].
    arrlen : list, optional
        length for vectors. The default is [.5, .5, .5, .5, .75, 1].
    labels : list, optional
        vector labels for legend. The default is ['$Dipole_{1}$',\
                                                  '$Dipole_{2}$', 'Boresight',\
                                                      '$Ray_{LOS}$'].
    markers : list, optional
        legend markers. The default is ['_','_',  r'$\longrightarrow$',\
                                        r'$\longrightarrow$'].
    colors : list, optional
        colors for labels in the legend. The default is ['red', '#12e193',\
                                                         'black', 'blue'].
    edgecolors : list, optional
        edgecolors for labels in the legend. The default is ['red', '#12e193',\
                                                             'black','blue'].
    arrowhead : list, optional
        how large is the arrow head. The default is [0.01, 0.01, 0.01, 0.01, \
                                                     0.45, 0.45].
    sct_kwargs : dict, optional
        target marker specifics.The default is {'alpha': 1,\
                                                'edgecolor':'black',\
                                                'c': 'lightgrey', \
                                                'marker': '*', 's': 180}.
    loc : str
        legend location. uses same keywords as matplotlib legend. \
            The default is 'upper right'.

    Returns
    -------
    fig : figure.Figure
        Figure object of matplotlib.figure module.
    ax : axes
        Axes object of matplotlib.
        
    Examples
    --------
    vec_args= M1_enu, M3_enu, M2_enu, M4_enu, RRI_enu, los_enu_arr
    connected_plot = op.display_observation_geometry(title_vec, \
                                                     time_array_1sec,\
                                                     pLon, pLat, OH, \
                                                     Lon, Lat, Alt,\
                                                     'Ottawa', *vec_args, n=45)

    """

    a=0; 
    Pz = Pz/100 ; z = z/100 # for better visualization, divide altitude by 100
    fig = plt.figure(figsize=(12, 12))
    ax=plt.axes(projection='3d')
    for vec in V:
        if len(x)>1:
            a+=1
            i=a-1
            ax.quiver(x[::step], y[::step], z[::step],
                          vec[:,0][::step], vec[:,1][::step], vec[:,2][::step],
                          length = arrlen[i], color = vc[i], 
                          linestyle = ls[i], lw = lw[i],
                          arrow_length_ratio = arrowhead[i])
            
            Px_vec = [Px]*len(x)
            Py_vec = [Py]*len(y)
            Pz_vec = [0.07]*len(z)
            for i in range(0, len(x), step):
                lon = [Px_vec[i], x[i]]
                lat = [Py_vec[i], y[i]]
                alt = [Pz_vec[i], z[i]]
                ax.plot(lon,lat,alt,linestyle = 'dashed', color = '#95d0fc', 
                        alpha = .7, lw=1)
                
            xmin, xmax = misc.set_3Dplot_limits(Px, x, 5)
            ymin, ymax = misc.set_3Dplot_limits(Py, y, 5)
            zmin, zmax = misc.set_3Dplot_limits(Pz, z, 5)
            ax.set_xlim(xmin, xmax)
            ax.set_ylim(ymin, ymax)
            ax.set_zlim(0, np.round(max(z)+0.5))

            plt.xticks(fontsize= 18)
            plt.yticks(fontsize= 18)

        else:
            a+=1
            i=a-1
            ax.quiver(x, y, z,
                  vec[0], vec[1], vec[2], 
                  length=arrlen[i], color=vc[i])
            coord_info="".join(['Lat: ', str(y), '\N{DEGREE SIGN}\n',
                              'Lon: ', str(x), '\N{DEGREE SIGN}\n',
                              'Alt: ', str(z), 
                              '(\N{DEGREE SIGN} x 100 km/\N{DEGREE SIGN})'])
            ax.text2D(0.05, 0.95, coord_info, fontsize=14, va='top', ha='left',
                    rotation_mode='anchor',
                    transform=ax.transAxes)
            Pz= [4]
            lon = [Px, x]
            lat = [Py, y]
            alt = [Pz, z]
   
    #  insert the legend
    misc.put_legend_fnt(ax, 4, loc, 0.2, 0.525, 0.90, 13,
                      labels=labels, linestyles= ls, 
                      markers=markers, colors=colors, 
                      edgecolors=edgecolors)  
      
    time1 = time_array[0].strftime('%H:%M:%S')
    time2 = time_array[-1].strftime('%H:%M:%S')
    
    #  insert the start and end times
    ax.text(x[0]-0.6, y[0]-0.95, z[0], time1, size=14)   
    ax.text(x[-1]-0.5, y[-1]-0.2, z[-1], time2, size = 14)
    
    # target
    ax.scatter3D(Px, Py, 0.01, marker = 'o', edgecolor="black", 
              facecolor='lightgrey', s= 180)
    if x[0] < Px:
        ax.text(Px-.5, Py, 0.5, target_name, c='k', size=14, zorder=20)
    else:
        ax.text(Px+.5, Py, 0.5, target_name, c='k', size=14, zorder=20)
    
    proj_ax = plt.figure().add_subplot(111, projection=ccrs.PlateCarree())
        
    proj_ax.set_xlim(ax.get_xlim())
    proj_ax.set_ylim(ax.get_ylim())
    # a useful function
    concat = lambda iterable: list(itertools.chain.from_iterable(iterable))
    # geometry
    target_projection = proj_ax.projection
    feature = cartopy.feature.NaturalEarthFeature('physical', 'land', '10m')
    geoms = feature.geometries()
    # specify boundaries
    boundary = proj_ax._get_extent_geom()
    # Project the geometry into target projection geometry
    geoms = [target_projection.project_geometry(geom, feature.crs) \
             for geom in geoms]
    
    # do not append invalid geometry
    geoms2 = []
    for i in range(len(geoms)) :
        if geoms[i].is_valid :
            geoms2.append(geoms[i])
    
    geoms = geoms2
    
    # intersection
    geoms = [boundary.intersection(geom) for geom in geoms]
    # geometry to path to polygon to collection
    paths = concat(geos_to_path(geom) for geom in geoms)
    polys = concat(path.to_polygons() for path in paths)
    lc = PolyCollection(polys, edgecolor='black', facecolor='green', 
                        closed=True, alpha=0.1)
    # arrange the view
    # nadir
    ax.view_init(azim=-90, elev=10)
    plt.tight_layout()   
    plt.close(proj_ax.figure)
    from matplotlib import rcParams
    rcParams['axes.labelpad'] = 15  
    plt.axis('off')
    plt.grid(b=None)
    
    plt.tight_layout()

    return fig, ax

def attitude_3d_ground_quiver(title, time_array, Px, Py, Pz, x, y, z, dir_vec, 
                                target_name, *V, step = 60, 
                                vc=['red', '#12e193', 'red', '#12e193', \
                                  'black', '#0165fc'], 
                                arrlen = [0.25, 0.25, 0.25, 0.25, 0.75, 0.75],
                                ls = ['solid', 'solid', 'solid', 'solid', \
                                  'solid',  'solid'],
                                labels=['$Dipole_{1}$','$Dipole_{2}$',\
                                       'Boresight','$Ray_{LOS}$'], 
                                markers=['_','$--$',r'$\longrightarrow$', \
                                     r'$\longrightarrow$','H','H'], 
                                legend_colors=['red', '#12e193','black',\
                                           '#0165fc' ,'None', 'black'], 
                                legend_edgecolors=['red', '#12e193','black',\
                                               '#0165fc' ,'black', 'black'], 
                                loc = 'upper center'):
    """
    Plotting function to show projected look direction of RRI on the ground.
    The defaults for the keywords are based on RRI plots.

    Parameters
    ----------
    title : str
        plot title.
    time_array : datetime.datetime
        experiment time interval as datetime array.
    Px : float
        geodetic longitude of the target.
    Py : float
        geodetic latitude of the target.
    Pz : float
        altitude of the target.
    x : numpy.ndarray[float]
        spacecraft longitude.
    y : numpy.ndarray[float]
        spacecraft latitude.
    z : numpy.ndarray[float]
        spacecraft altitude.
    dir_vec : numpy.ndarray[float]
        look direction vector in ENU.
    target_name : str
        target name.
    *V : numpy.ndarray
        vectors to be plotted (need to be in ENU coordinate system).
        vec_args = M1_enu, M3_enu, M2_enu, M4_enu, RRI_enu, los_enu_arr
    step : float, optional
        Step in seconds to plot vectors. The default is 60.
    vc : list, optional
        vector colors. 
        The default is ['red','#12e193', 'red', '#12e193', 'black', '#0165fc'].
    arrlen : TYPE, optional
        DESCRIPTION. 
        The default is [0.25, 0.25, 0.25, 0.25, 0.75, 0.75].
    ls : list, optional
        linestyles for vectors.
        The default is ['solid', 'solid', 'solid', 'solid','solid',  'solid'].
    labels : list, optional
        vector labels for legend.
        The default is ['$Dipole_{1}$','$Dipole_{2}$','Boresight''$Ray_{LOS}$'].
    markers : list, optional
        legend markers.
        The default is ['_','$--$',r'$\longrightarrow$',\
                        r'$\longrightarrow$','H','H'].
    legend_colors : list, optional
        colors for labels in the legend. 
        The default is ['red', '#12e193','black','#0165fc' ,'None', 'black'].
    legend_edgecolors : list, optional
        edgecolors for labels in the legend.
        The default is ['red', '#12e193','black','#0165fc' ,'black', 'black'].
    loc : str
        legend location. uses same keywords as matplotlib legend.
        The default is 'upper center'.

    Returns
    -------
    fig : figure.Figure
        Figure object of matplotlib.figure module.
    ax : axes
        Axes object of matplotlib.
        
    Examples
    --------
    V = M1_enu, M3_enu, M2_enu, M4_enu, RRI_enu, los_enu_arr
    ground_quiver_3d = op.attitude_3d_ground_quiver(title, time_array,\
                                                  pLon, pLat, OH, \
                                                  Lon, Lat, Alt, RRI_enu, 
                                                  'Ottawa', *V, step = 60) 
    """

    a=0 ; start_time = time_array[0]
    Pz = Pz/100 ; z = z/100
    
    fig = plt.figure(figsize=(12.06, 10.79))
    fig.suptitle(title)
    ax=plt.axes(projection='3d')
    for vec in V:
        if len(V)>1:
            lim1 = len(x)
            lim2 = len(vec)
            start_time=time_array[0]
            a+=1
            i=a-1
            # 3D quiver on the trajectory
            if np.less(lim1,lim2)==False:
                ax.quiver(x[0:lim2:step], y[0:lim2:step], z[0:lim2:step],
                          vec[:,0][::step], vec[:,1][::step], vec[:,2][::step],
                          length=arrlen[i], color=vc[i], 
                          linestyle=ls[i])
            else:
                ax.quiver(x[::step], y[::step], z[::step],
                          vec[:,0][::step], vec[:,1][::step], vec[:,2][::step],
                          length=arrlen[i], color=vc[i], 
                          linestyle=ls[i])               

        else:
            a+=1
            i=a-1
            ax.quiver(x, y, z,
                  vec[0], vec[1], vec[2], 
                  length=arrlen[i], color=vc[i])
            coord_info="".join(['Time (UTC): ', 
                                time_array.strftime('%H:%M:%S'),'\n',
                                'Lat: ', str(y), '\N{DEGREE SIGN}\n',
                                'Lon: ', str(x), '\N{DEGREE SIGN}\n',
                                'Alt: ', str(z), '( x 100 km)'])
            ax.text(0.05, 0.95, 0.2, coord_info, fontsize=14, va='top', 
                    ha='left', rotation_mode='anchor', transform=ax.transAxes)

    # 2D quiver on the ground
    for i in range(0, len(x), step):    
        ax.quiver(x[i], y[i], 0, dir_vec[:,0][i], dir_vec[:,1][i], 0, 
                  length=2, color='black', linewidths=0.7, 
                  arrow_length_ratio=0.2, zorder= 50)
        # connect to the target
        lon = [Px, x[i]]
        lat = [Py, y[i]]
        alt = [0, 0]
        ax.plot(lon, lat, alt, linestyle = 'dashed', color = '#95d0fc', lw=1.5)

    # determine plot limits
    xmin, xmax = misc.set_3Dplot_limits(Px, x, 1)
    ymin, ymax = misc.set_3Dplot_limits(Py, y, 1)
    zmin, zmax = misc.set_3Dplot_limits(Pz, z, 1)
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_zlim(0, np.round(max(z)+0.5))

    # target
    ax.scatter3D(Px, Py, 0.01, marker = 'o', edgecolor="black", 
              facecolor='lightgrey', s= 180)
    ax.text(Px, Py, 0.02, target_name, c='k', size=14, zorder=20)

    for i in range(0, len(x), step):
        if y[i] < 45.4215:
            ax.scatter(x[i], y[i], 0.3,
            marker='H', edgecolor='black', color='black', s=60)
        else:
            ax.scatter(x[i], y[i], 0.3,
            marker='H', edgecolor='black', color='None', s=60)
    
    time1=start_time.strftime('%H:%M:%S')
    time2=time_array[-1].strftime('%H:%M:%S')
    
    #  insert the start and end times
    ax.text(x[0]-0.6, y[0]-0.95, 1, time1, size=14, zorder = 55)
    ax.text(x[-1]-0.5, y[-1]-0.2, 1, time2, size = 14, zorder = 55)

    #  insert the legend
    misc.put_legend_fnt(ax, 6, loc, 0.2, 0.525, 0.90, 13,
                     labels = labels, linestyles = ls, 
                     markers = markers, colors = legend_colors, 
                     edgecolors = legend_edgecolors)  
    
    proj_ax = plt.figure().add_subplot(111, projection=ccrs.PlateCarree()) 
    proj_ax.set_xlim(ax.get_xlim())
    proj_ax.set_ylim(ax.get_ylim())

    concat = lambda iterable: list(itertools.chain.from_iterable(iterable))
    # geometry
    target_projection = proj_ax.projection
    feature = cartopy.feature.NaturalEarthFeature('physical', 'land', '10m')
    geoms = feature.geometries()
    # specify boundaries
    boundary = proj_ax._get_extent_geom()
    # Project the geometry into target projection geometry
    geoms = [target_projection.project_geometry(geom, feature.crs) \
             for geom in geoms]
    
    # do not append invalid geometry
    geoms2 = []
    for i in range(len(geoms)) :
        if geoms[i].is_valid :
            geoms2.append(geoms[i])  
    geoms = geoms2
    
    # intersection
    geoms = [boundary.intersection(geom) for geom in geoms]
    # geometry to path to polygon to collection
    paths = concat(geos_to_path(geom) for geom in geoms)
    polys = concat(path.to_polygons() for path in paths)
    lc = PolyCollection(polys, edgecolor='black', facecolor='lightgreen', 
                        closed=True, alpha=0.1)
    
    # overlay subplots
    ax.add_collection3d(lc, zs=0)
    # arrange the view
    ax.view_init(azim=-35, elev=15)
    
    plt.close(proj_ax.figure)

    tmp_planes = ax.zaxis._PLANES 
    ax.zaxis._PLANES = (tmp_planes[2], tmp_planes[3],
                             tmp_planes[0], tmp_planes[1],
                             tmp_planes[4], tmp_planes[5])
    
    rcParams['axes.labelpad'] = 30 
    ax.set_xlabel('Geographic Longitude (\N{DEGREE SIGN})', fontsize=18)
    ax.set_ylabel('Geographic Latitude (\N{DEGREE SIGN})', fontsize=18)
    ax.zaxis.set_rotate_label(False)  # disable automatic rotation
    ax.set_zlabel('Altitude ( x 100 km)', rotation=90, fontsize=18)

    
    ax.tick_params(axis='x', which='major', labelsize=18)
    ax.tick_params(axis='y', which='major', labelsize=18)
    ax.tick_params(axis='z', which='major', labelsize=18)

    
    plt.tight_layout()
    plt.show()

    return fig, ax

def attitude_3d_connect_to_target(title, time_array,                         
                             Px, Py, Pz, x, y, z, target_name, *V, step = 60,
                             vc=['red', '#12e193', 'red', '#12e193', 'black', \
                                           'blue'],
                            ls=['solid', 'solid', 'solid', 'solid',\
                                 'solid','solid'], 
                            lw = [2.5, 2.5, 2.5, 2.5, 2, 2],
                            arrlen = [0.5, 0.5, 0.5, 0.5, 1, 1],
                            labels=['$Dipole_{1}$','$Dipole_{2}$', \
                                    'Boresight','$Ray_{LOS}$','N','S'], 
                            markers=['_','_',  r'$\longrightarrow$', \
                                r'$\longrightarrow$', 'H','H'], 
                            colors=['red', '#12e193', 'black','blue','None', \
                                    'black'], 
                            edgecolors=['red', '#12e193', 'black','blue',\
                                        'black', 'black'],
                             arrowhead = [0.01, 0.01, 0.01, 0.01,  0.25, 0.25],
                             sct_kwargs = {'alpha': 1, 'edgecolor': 'black', \
                                           'c': 'lightgrey', 'marker': '*',
                                           's': 180}, 
                            loc = 'upper center' ):
    
    """
    Plotting function in 3D to connect the spacecraft location with the target.
    The defaults for the keywords are based on RRI plots.

    Parameters
    ----------
    title : str
        plot title.
    time_array : datetime.datetime
        experiment time interval as datetime array.
    Px : float
        geodetic longitude of the target.
    Py : float
        geodetic latitude of the target.
    Pz : float
        altitude of the target (km).
    x : numpy.ndarray
        spacecraft longitude.
    y : numpy.ndarray
        spacecraft latitude.
    z : numpy.ndarray
        spacecraft altitude.
    target_name : str
        target name.
    *V : numpy.ndarray
        vectors to be plotted (need to be in ENU coordinate system).
        vec_args = M1_enu, M3_enu, M2_enu, M4_enu, RRI_enu, los_enu_arr
    step : int
        time between the vectors. The default is 60. (60 seconds)
    vc : list, optional
        vector colors. The default is ['red', '#12e193', 'red', '#12e193',\
                                       'black', 'blue'].
    ls : list, optional
        linestyles for vectors. The default is ['solid', 'solid', 'solid',\
                                                'solid', 'solid', 'solid'].
    lw : list, optional
        linewidth of vectors. The default is [2.5, 2.5, 2.5, 2.5, 2, 2].
    arrlen : list, optional
        length for vectors. The default is [0.5, 0.5, 0.5, 0.5, 1, 1].
    labels : list, optional
        vector labels for legend. The default is ['$Dipole_{1}$',\
                                                  '$Dipole_{2}$', 'Boresight',\
                                                      '$Ray_{LOS}$'].
    markers : list, optional
        legend markers. The default is ['_','_',  r'$\longrightarrow$',\
                                        r'$\longrightarrow$'].
    colors : list, optional
        colors for labels in the legend. The default is ['red', '#12e193',\
                                                         'black', 'blue'].
    edgecolors : list, optional
        edgecolors for labels in the legend. The default is ['red', '#12e193',\
                                                             'black','blue'].
    arrowhead : list, optional
        how large is the arrow head. The default is [0.01, 0.01, 0.01, 0.01, \
                                                     0.25, 0.25].
    sct_kwargs : dict, optional
        target marker specifics.The default is {'alpha': 1,\
                                                'edgecolor':'black',\
                                                'c': 'lightgrey', \
                                                'marker': '*', 's': 180}.
    loc : str, optional
        Location of the legend. The default is 'upper center'.

    Returns
    -------
    fig : figure.Figure
        Figure object of matplotlib.figure module.
        
    ax : axes
        Axes object of matplotlib.
        
    Examples
    --------
    vec_args= M1_enu, M3_enu, M2_enu, M4_enu, RRI_enu, los_enu_arr
    connected_plot, ax = ap.attitude_3d_connect_to_target(title_vec, \
                                                     time_array_1sec,\
                                                     pLon, pLat, pH, \
                                                     Lon, Lat, Alt,\
                                                     target_name, \
                                                    *vec_args, step = 45)
    """
    a=0
    fig = plt.figure(figsize=(12, 12))
    ax=plt.axes(projection='3d')
    for vec in V:
        if len(x)>1:
            a+=1
            i=a-1
            ax.quiver(x[::step], y[::step], z[::step],
                          vec[:,0][::step], vec[:,1][::step], vec[:,2][::step],
                          length = arrlen[i], color = vc[i], 
                          linestyle = ls[i], linewidth = lw[i],
                          arrow_length_ratio = arrowhead[i])
            
            Px_vec = [Px]*len(x)
            Py_vec = [Py]*len(y)
            Pz_vec = [0.07]*len(z)
            for i in range(0, len(x), step):
                lon = [Px_vec[i], x[i]]
                lat = [Py_vec[i], y[i]]
                alt = [Pz_vec[i], z[i]]
                ax.plot(lon,lat,alt,linestyle = 'dashed', color = '#95d0fc', 
                        alpha = .7, linewidth=1)
                
            xmin, xmax = misc.set_3Dplot_limits(Px, x, 5)
            ymin, ymax = misc.set_3Dplot_limits(Py, y, 5)
            zmin, zmax = misc.set_3Dplot_limits(Pz, z, 5)
            ax.set_xlim(xmin, xmax)
            ax.set_ylim(ymin, ymax)
            ax.set_zlim(0, np.round(max(z)+0.5))
            plt.xticks(fontsize= 18)
            plt.yticks(fontsize= 18)

        else:
            a+=1
            i=a-1
            ax.quiver(x, y, z,
                  vec[0], vec[1], vec[2], 
                  length=arrlen[i], color=vc[i])
            coord_info="".join(['Lat: ', str(y), '\N{DEGREE SIGN}\n',
                              'Lon: ', str(x), '\N{DEGREE SIGN}\n',
                              'Alt: ', str(z), 
                              '(\N{DEGREE SIGN} x 100 km/\N{DEGREE SIGN})'])
            ax.text2D(0.05, 0.95, coord_info,  fontsize=14, va='top', ha='left',
                    rotation_mode='anchor',
                    transform=ax.transAxes)
            plt.xticks(fontsize= 18)
            plt.yticks(fontsize= 18)
   
    misc.put_legend_fnt(ax, 5, loc, 1.5, 0.5, 1, 12,
                     labels=labels, linestyles= ls, 
                     markers=markers, colors=colors, 
                     edgecolors=edgecolors)
    
    for i in range(0, len(x), step):
        if y[i] < Py:
            ax.scatter(x[i], y[i], 0.1,
            marker='H', edgecolor='black', color='black', s=60)

        else:
            ax.scatter(x[i], y[i], 0.1,
            marker='H', edgecolor='black', color='None', s=60)

    ax.plot3D(x, y, 0.01, c = 'black', linewidth= 1.25)

    #  insert the legend
    misc.put_legend_fnt(ax, 6, loc, 0.2, 0.525, 0.90, 13,
                     labels=labels, linestyles= ls, 
                     markers=markers, colors=colors, 
                     edgecolors=edgecolors)  
      
    time1 = time_array[0].strftime('%H:%M:%S')
    time2 = time_array[-1].strftime('%H:%M:%S')
    
    #  insert the start and end times
    ax.text(x[0]-0.6, y[0]-0.95, 0.3, time1, size=14)   
    ax.text(x[-1]-0.5, y[-1]-0.2, 0.3, time2, size = 14)
    
    # target
    ax.scatter3D(Px, Py, 0.01, marker = 'o', edgecolor="black", 
              facecolor='lightgrey', s= 180)
    if x[0] < Px:
        ax.text(Px-.5, Py, 0.5, target_name, c='k', size=14, zorder=20)
    else:
        ax.text(Px+.5, Py, 0.5, target_name, c='k', size=14, zorder=20)
    
   
    proj_ax = plt.figure().add_subplot(111, projection=ccrs.PlateCarree())
        
    proj_ax.set_xlim(ax.get_xlim())
    proj_ax.set_ylim(ax.get_ylim())

    concat = lambda iterable: list(itertools.chain.from_iterable(iterable))

    # geometry
    target_projection = proj_ax.projection
    feature = cartopy.feature.NaturalEarthFeature('physical', 'land', '10m')
    geoms = feature.geometries()
    # specify boundaries
    boundary = proj_ax._get_extent_geom()
    # Project the geometry into target projection geometry
    geoms = [target_projection.project_geometry(geom, feature.crs) \
             for geom in geoms]
    
    # do not append invalid geometry
    geoms2 = []
    for i in range(len(geoms)) :
        if geoms[i].is_valid :
            geoms2.append(geoms[i])
    
    geoms = geoms2
    
    # intersection
    geoms = [boundary.intersection(geom) for geom in geoms]
    # geometry to path to polygon to collection
    paths = concat(geos_to_path(geom) for geom in geoms)
    polys = concat(path.to_polygons() for path in paths)
    lc = PolyCollection(polys, edgecolor='black', facecolor='green', 
                        closed=True, alpha=0.1)
    
    # overlay subplots
    ax.add_collection3d(lc, zs=0)
    # arrange the view
    # nadir
    ax.view_init(azim=-128, elev=22)
    plt.tight_layout()   
    plt.close(proj_ax.figure)

    rcParams['axes.labelpad'] = 15  
    ax.set_xlabel('Geographic Longitude (\N{DEGREE SIGN})', fontsize=18)
    ax.set_ylabel('Geographic Latitude (\N{DEGREE SIGN})', fontsize=18)
    ax.zaxis.set_rotate_label(False)  # disable automatic rotation
    ax.set_zlabel('Altitude ( x 100 km)', rotation=90, fontsize=18)

    
    ax.tick_params(axis='x', which='major', labelsize = 18)
    ax.tick_params(axis='y', which='major', labelsize=18)
    ax.tick_params(axis='z', which='major', labelsize=18)

    plt.tight_layout()

    return fig, ax

def attitude_3d_connect_to_subpoint(title, time_array, Px, Py, Pz, x, y, z, 
                                   target_name, *V, step = 60, 
                                vc=['red', '#12e193', 'red', '#12e193', \
                                  'black', '#0165fc'], 
                                arrlen = [0.25, 0.25, 0.25, 0.25, 0.75, 0.75],
                                ls = ['solid', 'solid', 'solid', 'solid', \
                                  'solid',  'solid'],
                                labels=['$Dipole_{1}$','$Dipole_{2}$',\
                                       'Boresight','$Ray_{LOS}$'], 
                                markers=['_','$--$',r'$\longrightarrow$', \
                                     r'$\longrightarrow$','H','H'], 
                                legend_colors=['red', '#12e193','black',\
                                           '#0165fc' ,'None', 'black'], 
                                legend_edgecolors=['red', '#12e193','black',\
                                               '#0165fc' ,'black', 'black'], 
                                loc = 'upper center'):
    """
    Plotting function to show projected look direction of RRI on the ground.
    The defaults for the keywords are based on RRI plots.

    Parameters
    ----------
    title : str
        plot title.
    time_array : datetime.datetime
        experiment time interval as datetime array.
    Px : float
        geodetic longitude of the target.
    Py : float
        geodetic latitude of the target.
    Pz : float
        altitude of the target.
    x : numpy.ndarray[float]
        spacecraft longitude.
    y : numpy.ndarray[float]
        spacecraft latitude.
    z : numpy.ndarray[float]
        spacecraft altitude.
    target_name : str
        target name.
    *V : numpy.ndarray
        vectors to be plotted (need to be in ENU coordinate system).
        vec_args = M1_enu, M3_enu, M2_enu, M4_enu, RRI_enu, los_enu_arr
    step : float, optional
        Step in seconds to plot vectors. The default is 60.
    vc : list, optional
        vector colors. 
        The default is ['red','#12e193', 'red', '#12e193', 'black', '#0165fc'].
    arrlen : TYPE, optional
        DESCRIPTION. 
        The default is [0.25, 0.25, 0.25, 0.25, 0.75, 0.75].
    ls : list, optional
        linestyles for vectors.
        The default is ['solid', 'solid', 'solid', 'solid','solid',  'solid'].
    labels : list, optional
        vector labels for legend.
        The default is ['$Dipole_{1}$','$Dipole_{2}$','Boresight''$Ray_{LOS}$'].
    markers : list, optional
        legend markers.
        The default is ['_','$--$',r'$\longrightarrow$',\
                        r'$\longrightarrow$','H','H'].
    legend_colors : list, optional
        colors for labels in the legend. 
        The default is ['red', '#12e193','black','#0165fc' ,'None', 'black'].
    legend_edgecolors : list, optional
        edgecolors for labels in the legend.
        The default is ['red', '#12e193','black','#0165fc' ,'black', 'black'].
    loc : str
        legend location. uses same keywords as matplotlib legend.
        The default is 'upper center'.

    Returns
    -------
    fig : figure.Figure
        Figure object of matplotlib.figure module.
        
    Examples
    --------
    V = M1_enu, M3_enu, M2_enu, M4_enu, RRI_enu, los_enu_arr
    connect_to_subpoint= ap.attitude_3d_connect_to_subpoint(title, time_array,\
                                                  pLon, pLat, OH, \
                                                  Lon, Lat, Alt,
                                                  'Ottawa', *V, step = 60) 
    """

    a=0 ; start_time = time_array[0]
    Pz = Pz/100 ; z = z/100
    
    fig = plt.figure(figsize=(12.06, 10.79))
    fig.suptitle(title)
    ax=plt.axes(projection='3d')
    for vec in V:
        if len(V)>1:
            lim1 = len(x)
            lim2 = len(vec)
            start_time=time_array[0]
            a+=1
            i=a-1
            # 3D quiver on the trajectory
            if np.less(lim1,lim2)==False:
                ax.quiver(x[0:lim2:step], y[0:lim2:step], z[0:lim2:step],
                          vec[:,0][::step], vec[:,1][::step], vec[:,2][::step],
                          length=arrlen[i], color=vc[i], 
                          linestyle=ls[i])
            else:
                ax.quiver(x[::step], y[::step], z[::step],
                          vec[:,0][::step], vec[:,1][::step], vec[:,2][::step],
                          length=arrlen[i], color=vc[i], 
                          linestyle=ls[i])               

        else:
            a+=1
            i=a-1
            ax.quiver(x, y, z,
                  vec[0], vec[1], vec[2], 
                  length=arrlen[i], color=vc[i])
            coord_info="".join(['Time (UTC): ', 
                                time_array.strftime('%H:%M:%S'),'\n',
                                'Lat: ', str(y), '\N{DEGREE SIGN}\n',
                                'Lon: ', str(x), '\N{DEGREE SIGN}\n',
                                'Alt: ', str(z), '( x 100 km)'])
            ax.text(0.05, 0.95, 0.2, coord_info, fontsize=14, va='top', 
                    ha='left', rotation_mode='anchor', transform=ax.transAxes)


    # determine plot limits          # 
    xmin, xmax = misc.set_3Dplot_limits(Px, x, 1)
    ymin, ymax = misc.set_3Dplot_limits(Py, y, 1)
    zmin, zmax = misc.set_3Dplot_limits(Pz, z, 1)
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_zlim(0, np.round(max(z)+0.5))
 
    for i in range(0, len(x), step):
        lon = [x[i], x[i]]
        lat = [y[i], y[i]]
        alt = [0, z[i]]
        ax.plot(lon,lat,alt,linestyle = 'dashed', color = 'k', 
                alpha = .7, linewidth=1)

    # target
    ax.scatter3D(Px, Py, 0.01, marker = 'o', edgecolor="black", 
              facecolor='lightgrey', s= 180)
    ax.text(Px, Py, 0.02, target_name, c='k', size=14, zorder=20)

    for i in range(0, len(x), step):
        if y[i] < 45.4215:
            ax.scatter(x[i], y[i], 0.3,
            marker='H', edgecolor='black', color='black', s=60)
        else:
            ax.scatter(x[i], y[i], 0.3,
            marker='H', edgecolor='black', color='None', s=60)
    
    time1=start_time.strftime('%H:%M:%S')
    time2=time_array[-1].strftime('%H:%M:%S')
    
    #  insert the start and end times
    ax.text(x[0]-0.6, y[0]-0.95, 1, time1, size=14, zorder = 55)
    ax.text(x[-1]-0.5, y[-1]-0.2, 1, time2, size = 14, zorder = 55)

    #  insert the legend
    misc.put_legend_fnt(ax, 6, loc, 0.2, 0.525, 0.90, 13,
                     labels = labels, linestyles = ls, 
                     markers = markers, colors = legend_colors, 
                     edgecolors = legend_edgecolors)  
    
    proj_ax = plt.figure().add_subplot(111, projection=ccrs.PlateCarree()) 
    proj_ax.set_xlim(ax.get_xlim())
    proj_ax.set_ylim(ax.get_ylim())

    concat = lambda iterable: list(itertools.chain.from_iterable(iterable))
    # geometry
    target_projection = proj_ax.projection
    feature = cartopy.feature.NaturalEarthFeature('physical', 'land', '10m')
    geoms = feature.geometries()
    # specify boundaries
    boundary = proj_ax._get_extent_geom()
    # Project the geometry into target projection geometry
    geoms = [target_projection.project_geometry(geom, feature.crs) \
             for geom in geoms]
    
    # do not append invalid geometry
    geoms2 = []
    for i in range(len(geoms)) :
        if geoms[i].is_valid :
            geoms2.append(geoms[i])  
    geoms = geoms2
    
    # intersection
    geoms = [boundary.intersection(geom) for geom in geoms]
    # geometry to path to polygon to collection
    paths = concat(geos_to_path(geom) for geom in geoms)
    polys = concat(path.to_polygons() for path in paths)
    lc = PolyCollection(polys, edgecolor='black', facecolor='lightgreen', 
                        closed=True, alpha=0.1)
    
    # overlay subplots
    ax.add_collection3d(lc, zs=0)
    # arrange the view
    ax.view_init(azim=-35, elev=15)
    
    plt.close(proj_ax.figure)

    tmp_planes = ax.zaxis._PLANES 
    ax.zaxis._PLANES = (tmp_planes[2], tmp_planes[3],
                             tmp_planes[0], tmp_planes[1],
                             tmp_planes[4], tmp_planes[5])
    
    rcParams['axes.labelpad'] = 30 
    ax.set_xlabel('Geographic Longitude (\N{DEGREE SIGN})', fontsize=18)
    ax.set_ylabel('Geographic Latitude (\N{DEGREE SIGN})', fontsize=18)
    ax.zaxis.set_rotate_label(False)  # disable automatic rotation
    ax.set_zlabel('Altitude ( x 100 km)', rotation=90, fontsize=18)

    
    ax.tick_params(axis='x', which='major', labelsize=18)
    ax.tick_params(axis='y', which='major', labelsize=18)
    ax.tick_params(axis='z', which='major', labelsize=18)

    
    plt.tight_layout()
    plt.show()

    return fig, ax


def earth_radius_at_latitude(latitude):
    """
    Function to calculate Earth radius at latitude to correctly scale 
    the ellipse height for FOV plotter.

    Parameters
    ----------
    latitude : float
        Geodetic latitude (degrees).

    Returns
    -------
    R : float
        Earth's radius in meters.

    """
    # Earth's radius at equator and poles (in meters)
    a = 6378137.0
    b = 6356752.3

    # Latitude in radians
    latitude_rad = np.deg2rad(latitude)

    # Earth's radius at the given latitude
    R = np.sqrt((((a**2 * np.cos(latitude_rad))**2 + (b**2 * np.sin(latitude_rad))**2) /
                 ((a * np.cos(latitude_rad))**2 + (b * np.sin(latitude_rad))**2)))
    return R


def plot_attitude_accuracy(file_RRI):
    """
    Plots the attitude accuracy using the accuracy data in RRI data files.
    0=Dropout, 1=Rough, 2=Coarse, 3=Moderate, 4=Fine, 9=NaN.
    
    Parameters
    ----------
    file_RRI : str
        Filename of the RRI data file including the path.

    Returns
    -------
    fig : figure.Figure
        Figure object of matplotlib.figure module.
    ax : axes
        Axes object of matplotlib.  

    """
    
    dict_rri = ei.rri_ephemeris(file_RRI)
    acc = dict_rri['accuracy'] 
    time_array = dict_rri['time_array']
    
    myFmt = dtformat.DateFormatter("%H:%M:%S")
    seconds=dtformat.MicrosecondLocator(interval=10000000)
    
    fig, ax = plt.subplots()
    ax.plot(time_array[::5], acc[::5], color='grey', marker='.', \
            markersize=10)   
    ax.xaxis.set_major_formatter(myFmt)
    ax.xaxis.set_minor_locator(seconds)
    ax.set_ylim(-0.1,5.3)
    ax.grid(True, linestyle='--')
    
    ax.yaxis.set_major_locator(MultipleLocator(1))
    ax.set_ylabel('Accuracy', color="k", fontsize=14)
    at = AnchoredText(" 0=Dropout, 1=Rough, 2=Coarse, 3=Moderate, 4=Fine, 9=NaN", 
                      prop=dict(size=12), frameon=False, loc='upper left') 
    ax.add_artist(at)

    ax.tick_params(axis='both', which='major', labelsize=14)
    ax.get_yaxis().set_label_coords(-0.075,0.5)
    
    
    return fig, ax

def trajectory_plotter_2d_map(time_array, extent, x, y, z, Px, Py, XV, \
                              target_name, index_ca):
    """
    Plots the spacecraft trajectory on map.
    Prints the trajectory information on screen: altitude increasing/decreasing
    spacecraft going towards: NE, SW, NW, SE.

    Parameters
    ----------
    time_array : datetime.datetime
        experiment time interval as datetime array.
    extent : TYPE
        DESCRIPTION.
    x : numpy.ndarray[float]
        spacecraft longitude(degrees).
    y : numpy.ndarray[float]
        spacecraft latitude (degrees).
    z : numpy.ndarray[float]
        spacecraft altitude (km).
    Px : float
        geodetic longitude of the target.
    Py : float
        geodetic latitude of the target.
    target_name : str
        Target name for labels.
    index_ca : int
        index of the point of closest approach. can be found by using 
        miscellaneous.find_index

    Returns
    -------
    fig : figure.Figure
        Figure object of matplotlib.figure module.
    ax : axes
        Axes object of matplotlib.        

    """
                
    w = 0.004; OH = min(z)-0.01;
    start_time = time_array[0]
    central_lon, central_lat = misc.coverage(extent)
    
    fig = plt.figure(figsize=(12, 15))
    proj = ccrs.Orthographic(central_longitude=central_lon, 
                                              central_latitude=central_lat)
    ax = plt.subplot(projection=proj)
        
    ax.set_extent(extent, ccrs.PlateCarree())
    ax.add_feature(cartopy.feature.OCEAN, color='white', alpha=0.2, zorder=0)
    ax.add_feature(cartopy.feature.LAND, edgecolor='grey', 
                   color='silver', alpha=0.3, zorder=10)
    ax.add_feature(cartopy.feature.LAKES, color='white', alpha=0.6, zorder=15)
    
    ax.plot(x, y, color="orange", linestyle = 'dashed', linewidth = '2.5', 
                zorder = 20, transform=ccrs.PlateCarree())
    
    ax.text(x[0]-2.5, y[0]-0.5, start_time.strftime("%H:%M:%S"), 
             transform=ccrs.PlateCarree(), fontsize= 12) # trajectory start time
    
    # mark the point of closest approach
    ca_time = time_array[index_ca][0]
    ax.text(x[index_ca]-2.5, y[index_ca]-0.5, ca_time.strftime("%H:%M:%S"), 
             transform=ccrs.PlateCarree(), fontsize= 12, zorder=25) # trajectory start time
    
    ax.text(x[-1]-2.5, y[-1]+0.5, time_array[-1].strftime("%H:%M:%S"),  
             transform=ccrs.PlateCarree(), fontsize= 12) # trajectory end time
    
    ax.scatter(Px, Py, marker = '*', edgecolor="black", \
               facecolor= 'red',  s= 150, transform=ccrs.PlateCarree(), \
                   label= target_name)
    
    
    ax.scatter(x[index_ca], y[index_ca], marker = 'x', color="black",  s= 200, 
                transform=ccrs.PlateCarree(), zorder = 25,
                label= 'Closest approach') 
    
    # insert legend
    list_edgecolor  = ["black", "orange", "black", "black"]
    list_color  = ["black", "orange", "red",  "black"]
    list_mak   = [r'$\longrightarrow$', '_','*',  'x']
    list_lab    = ['$V_{S/C}$', '$Trajectory_{S/C}$', 
                   target_name, 'Closest approach']
    
    class MarkerHandler(HandlerBase):
        def create_artists(self, legend, tup, xdescent, ydescent,
                            width, height, fontsize, trans):
            return [ax.scatter([width/2], [height/2.],
                           marker=tup[2], color=tup[1], edgecolor=tup[0], 
                           transform=trans)]

    ax.legend(list(zip(list_edgecolor, list_color, list_mak)), list_lab, 
              handler_map={tuple:MarkerHandler()}, fontsize=10, 
              loc='upper right').set_zorder(102) 
    
    ax.text(-0.1, 0.5, 'Geographic Latitude (\N{DEGREE SIGN})', 
             va='bottom', ha='center',
            rotation='vertical', rotation_mode='anchor',
            transform=ax.transAxes, fontsize= 14)
    ax.text(0.5, -0.12, 'Geographic Longitude (\N{DEGREE SIGN})', 
             va='bottom', ha='center',
            rotation='horizontal', rotation_mode='anchor',
            transform=ax.transAxes, fontsize= 14)
    
    misc.sc_direction_plotter(ax, x, y, z)
    
    ax.set_aspect('auto', adjustable=None)  
    gl = ax.gridlines(draw_labels=True, dms=True, x_inline=False, \
                      y_inline=False, linestyle=':', zorder=15)
    gl.top_labels = gl.right_labels = False
    gl.xlabel_style = {'size': 14, 'color': 'black'}
    gl.ylabel_style = {'size': 14, 'color': 'black'}

    return fig, ax


def attitude_2d_on_map(time_array, extent, x, y, z, Px, Py, V, \
                              target_name, inst_name, index_ca, step = 60):
    """
    Plots the spacecraft trajectory and instrument pointing direction vectors 
    on map.

    Parameters
    ----------
    time_array : datetime.datetime
        experiment time interval as datetime array.
    extent : list
        extent = [Lonmin, Lonmax, Latmin, Latmax].
    x : numpy.ndarray[float]
        spacecraft longitude(degrees).
    y : numpy.ndarray[float]
        spacecraft latitude (degrees).
    z : numpy.ndarray[float]
        spacecraft altitude (km).
    Px : float
        geodetic longitude of the target.
    Py : float
        geodetic latitude of the target.
    V : numpy.ndarray[float]
        vectors to be plotted (need to be in ENU coordinate system).
    inst_name : str
        Instrument name for labels.
    target_name : str
        Target name for labels.
    index_ca : int
        index of the point of closest approach. can be found by using 
        miscellaneous.find_index
    step : int
        Time in seconds between pointing direction vectors. The default is 60.
        

    Returns
    -------
    fig : figure.Figure
        Figure object of matplotlib.figure module.
    ax : axes
        Axes object of matplotlib.        

    """
                
    w = 0.004
    start_time = time_array[0]
    central_lon, central_lat = misc.coverage(extent)
    
    fig = plt.figure(figsize=(12, 15))
    proj = ccrs.PlateCarree()
    ax = plt.subplot(projection=proj)
        
    ax.set_extent(extent, ccrs.PlateCarree())
    ax.add_feature(cartopy.feature.OCEAN, color='white', alpha=0.2, zorder=0)
    ax.add_feature(cartopy.feature.LAND, edgecolor='grey', 
                   color='silver', alpha=0.3, zorder=10)
    ax.add_feature(cartopy.feature.LAKES, color='white', alpha=0.6, zorder=15)
    
    # plot the RRI pointing direction.
    ax.quiver(x[::step], y[::step], 
              V[::step,0], V[::step,1],  angles='xy', headaxislength = 5,
              transform=ccrs.PlateCarree(), label = 'RRI boresight',
              color='black', scale = 10, width=w, axes=ax, zorder=20)
    
    # spacecraft trajectory
    ax.plot(x, y, color="orange", linestyle = 'dashed', linewidth = '2.5', 
                zorder = 20, transform=ccrs.PlateCarree())
    
    ax.text(x[0]-2.5, y[0]-0.5, start_time.strftime("%H:%M:%S"), 
             transform=ccrs.PlateCarree(), fontsize= 12) # trajectory start time
    
    
    ax.text(x[-1]-2.5, y[-1]+0.5, time_array[-1].strftime("%H:%M:%S"),  
             transform=ccrs.PlateCarree(), fontsize= 12) # trajectory end time
    
    # target location
    ax.scatter(Px, Py, marker = '*', edgecolor="black", \
               facecolor= 'red', s= 200, transform=ccrs.PlateCarree(), \
                   label= target_name)

    # beginning of the pass
    ax.scatter(x[0], y[0], marker = 'x', color='navy', s=150, linewidth=2.0,
               transform=ccrs.PlateCarree(), label ='Beginning of the pass', 
               zorder = 25) 
    
    list_edgecolor  = ["black", "navy", "orange", "black"]
    list_color  = ["black","navy", "orange", "red"]
    list_mak   = [r'$\longrightarrow$', 'x', '_', '*']
    list_lab    = ['RRI boresight', 'Start', '$Trajectory_{S/C}$', target_name]
    
    class MarkerHandler(HandlerBase):
        def create_artists(self, legend, tup, xdescent, ydescent,
                            width, height, fontsize, trans):
            return [ax.scatter([width/2], [height/2.],
                           marker=tup[2], color=tup[1], edgecolor=tup[0], 
                           transform=trans)]

    ax.legend(list(zip(list_edgecolor, list_color, list_mak)), list_lab, 
              handler_map={tuple:MarkerHandler()}, fontsize=10, 
              loc='upper right').set_zorder(102) 
    
    ax.text(-0.1, 0.5, 'Geographic Latitude (\N{DEGREE SIGN})', 
             va='bottom', ha='center',
            rotation='vertical', rotation_mode='anchor',
            transform=ax.transAxes, fontsize= 14)
    ax.text(0.5, -0.12, 'Geographic Longitude (\N{DEGREE SIGN})', 
             va='bottom', ha='center',
            rotation='horizontal', rotation_mode='anchor',
            transform=ax.transAxes, fontsize= 14)
        
    ax.set_aspect('auto', adjustable=None)  
    gl = ax.gridlines(draw_labels=True, dms=True, x_inline=False, \
                      y_inline=False, linestyle=':', zorder=15)
    gl.top_labels = gl.right_labels = False
    gl.xlabel_style = {'size': 14, 'color': 'black'}
    gl.ylabel_style = {'size': 14, 'color': 'black'}

    return fig, ax

def attitude_2d_altitude(time_array, extent, x, y, z, Px, V, \
                         inst_name, target_name, x_axis = 'lon', step = 60):
    """
    

    Parameters
    ----------
    time_array : datetime.datetime
        experiment time interval as datetime array.
    extent : list
        extent = [Lonmin, Lonmax, Altmin, Altmax].
    x : numpy.ndarray[float]
        spacecraft longitude(degrees).
    y : numpy.ndarray[float]
        spacecraft latitude (degrees).
    z : numpy.ndarray[float]
        spacecraft altitude (km).
    Px : float
        Location of the target on the x_axis (lat or lon) (degrees).
    V : numpy.ndarray[float]
        vectors to be plotted (need to be in ENU coordinate system).
    inst_name : str
        Instrument name for labels.
    target_name : str
        Target name for labels.
    x_axis : str, optional
        Parameter to be plotted on the x-axis (lat or lon). 
        The default is 'lon'.
    step : int
        Time in seconds between pointing direction vectors. The default is 60.

    Returns
    -------
    fig : figure.Figure
        Figure object of matplotlib.figure module.
    ax : axes
        Axes object of matplotlib.    

    """
    
   
    w = 0.004; OH = min(z)-0.01;
    start_time = time_array[0]
    
    fig = plt.figure(figsize=(12, 15))
    ax = plt.subplot()
    ax.set_xlim(extent[0:2])
    ax.set_ylim(extent[2:4])
    
    if x_axis =='lon':
        ax.quiver(x[::step], z[::step],  V[::step,0], V[::step,2],  
                  label = 'RRI boresight', color='black', scale=10, \
                      width=w, axes=ax, zorder=20)  
        ax.set_xlabel('Geographic Longitude (\N{DEGREE SIGN})', fontsize=14)
        ax.set_ylabel('Spacecraft Altitude (km)', fontsize=14) 
        # spacecraft trajectory
        ax.plot(x, z, color='orange', linewidth = 2.5, linestyle = 'dashed')    
        # target location
        ax.scatter(Px, OH, marker='*', edgecolor= 'black', color='red',
                    s=200, label = target_name)
        # beginning of the trajectory
        ax.scatter(x[0], z[0], marker = 'x', color='navy', s=150, 
                   label ='Beginning of the pass', zorder = 25) 
        # trajectory start time
        ax.text(x[0]-2.5, z[0]-0.5, start_time.strftime("%H:%M:%S"), 
                  fontsize= 12) 
        # trajectory end time
        ax.text(x[-1]-2.5, z[-1]+0.5, time_array[-1].strftime("%H:%M:%S"),  
                 fontsize= 12) 
    elif x_axis =='lat':
        ax.quiver(y[::step], z[::step],  V[::step,1], V[::step,2],  
                  label = 'RRI boresight', color='black', scale=10, \
                      width=w, axes=ax, zorder=20)  
        ax.set_xlabel('Geographic Latitude (\N{DEGREE SIGN})', fontsize=14)
        ax.set_ylabel('Spacecraft Altitude (km)', fontsize=14)     
        # target location        
        ax.scatter(Px, OH, marker='*', edgecolor= 'black', color='red',
                    s=200, label = target_name)
        # beginning of the trajectory        
        ax.scatter(y[0], z[0], marker = 'x', color='navy', s=150, 
                   label ='Beginning of the pass', zorder = 25)  
        # trajectory start time
        ax.text(y[0]-2.5, z[0]-0.5, start_time.strftime("%H:%M:%S"), 
                 fontsize= 12) 
        # trajectory end time
        ax.text(y[-1]-2.5, z[-1]+0.5, time_array[-1].strftime("%H:%M:%S"),  
                 fontsize= 12) 

    list_edgecolor  = ["black", "navy", "orange", "black"]
    list_color  = ["black","navy", "orange", "red"]
    list_mak   = [r'$\longrightarrow$', 'x', '_', '*']
    list_lab    = ['RRI boresight', 'Start', '$Trajectory_{S/C}$', target_name]
    
    class MarkerHandler(HandlerBase):
        def create_artists(self, legend, tup, xdescent, ydescent,
                            width, height, fontsize, trans):
            return [ax.scatter([width/2], [height/2.],
                           marker=tup[2], color=tup[1], edgecolor=tup[0], 
                           transform=trans)]
    
    ax.legend(list(zip(list_edgecolor, list_color, list_mak)), list_lab, 
              handler_map={tuple:MarkerHandler()}, fontsize=10, loc='upper right') 
    ax.set_aspect('auto', adjustable=None)  

    ax.grid(True, linestyle='--')
    plt.xticks(fontsize= 14)
    plt.yticks(fontsize= 14)
    
    return fig, ax

def fov_plotter(extent, time_array, x, y, z, fov_deg, 
                px, py, step= 90, inst_name = 'FAI', target_name = 'ICEBEAR'):
    
    """
    Plot the FOV of the instrument for Nadir look directions.
    Caution: This code ONLY works for NADIR pointing instruments.
    
    Parameters
    ----------
    extent : list
        [Lonmin, Lonmax, Latmin, Latmax].    
    time_array : numpy.ndarray[datetime]
        Experiment time interval.
    x : numpy.ndarray[float]
        Spacecraft longitude (degrees).    
    y : numpy.ndarray[float]
        Spacecraft latitude (degrees).
    z : numpy.ndarray[float]
        Spacecraft altitude (km).
    fov_deg : float
        Field-of-view angle (degrees).
    px : float
        Target longitude (degrees).
    py : float
        Target latitude (degrees).
    step: float, optional
        Time interval in seconds to plot the inst. vector.The default is 90.
    inst_name : string
        Name of the instrument. The default is 'FAI'.
    target_name : string
        Target name. The default is 'ICEBEAR'.
    Returns
    -------
    fig: figure.Figure
        Figure object of matplotlib.figure module.
    ax : TYPE
        DESCRIPTION.
    """
    
    central_lon, central_lat = misc.coverage(extent)
    crs = ccrs.Orthographic(central_longitude=central_lon, 
                            central_latitude=central_lat)
    title = '-'.join([inst_name, target_name + ' Coordinated experiment\n' +\
                      time_array[0].strftime('%Y-%m-%d')])
    fig = plt.figure(figsize=(12,8))
    ax = plt.axes(projection=crs)  
    plt.suptitle(title, fontsize = 28)
    ax.set_extent([extent[0]-20, extent[1]+40, extent[2]-20, 90])
    
    # ax.set_global()
    # ax.stock_img()
    
    ax.add_feature(cartopy.feature.OCEAN, color='white', alpha=1, zorder=0)
    ax.add_feature(cartopy.feature.LAND, edgecolor='white', 
                    color='silver', alpha=0.3, zorder=10)
    ax.add_feature(cartopy.feature.LAKES, color='white', alpha=1, zorder=0)

      
    num_colors = len(x) + 1
    cmap = cm.get_cmap('Set1', num_colors)
    
    i , j = 0, 0
    
    for i in range(0, len(x), step):
    
        # Generate a random color index
        i += 1
        # color_index = np.random.randint(0, num_colors)
        color_index = i
        # Get the color from the colormap based on the index
        color = cmap(color_index)
        facecolor = (*color[:3], 0.2)  # Add alpha to the base color
        edgecolor = (*color[:3], 1)   # Fully opaque edge     
        
        pos = [z[i], y[i], x[i]] # alt (m), lat, long
        radius = pos[0] * np.sin(np.deg2rad(fov_deg/2)) 
        
        ax.tissot(rad_km=radius, lons=pos[2], lats=pos[1], n_samples=64, \
                     facecolor=facecolor, edgecolor=edgecolor, 
                     linewidth=2, alpha = 0.3, zorder = 15)
    
        ax.text(x[i] + 5, y[i], (time_array[i].strftime("%H:%M:%S")), c=color, 
                size=18, weight = 'bold', transform=ccrs.PlateCarree(), 
                zorder=20)
        
                
    ax.scatter(px, py, marker = '*', s=400, c= 'magenta', zorder=20, 
               transform = ccrs.PlateCarree(), label = target_name)    
    
  
    ax.plot(x, y, linestyle='dashed', c='r', zorder=30,
            transform = ccrs.PlateCarree(), label = 'Trajectory')
    
    
    ax.scatter(x[::step], y[::step], marker = 'x', s = 100, color = 'k', lw = 2,
               label = 'Nadir Pointing', zorder = 20,
               transform = ccrs.PlateCarree())
    
    ax.legend(fontsize = 18, loc = 'upper right')
    ax.coastlines(linewidth=0.15) 
    gl = ax.gridlines(draw_labels=True, linewidth=1, alpha=0.3, linestyle='--', 
                 crs = ccrs.PlateCarree())
    
    gl.xlocator = ticker.FixedLocator(np.arange(extent[0]-20,
                                                 extent[1]+40,
                                                 10))
    gl.ylocator = ticker.FixedLocator(np.arange(extent[2]-20, 90,10))
    gl.xlabel_style = {'size': 24}
    gl.ylabel_style = {'size': 24}

    return fig, ax


def plot_slew_rri(ax, ylim_min, ylim_max, panel_number, time_array, slew, \
                  slew_angle, cb_axis = 'no', time='no'):
    """
    Function to plot the temporal change of slew accuracy for the criterion
    set by slew_angle.

    Parameters
    ----------
    ax : axes
        Axes object of matplotlib.  
    ylim_min : float
        DESCRIPTION.
    ylim_max : float
        DESCRIPTION.
    panel_number : int
        panel number in slew criteria plot.
    time_array : datetime.datetime
        experiment time interval as datetime array.
    slew : numpy.ndarray[float]
        slew parameter. Can be [-2, -1, -0.5, 0, 0.5, 1, 2]
    slew_angle : float
        criteria for slew (degrees).
    cb_axis : str, optional
        color bar axes display ('yes'/'no'). The default is 'no'.
    time : str, optional
        x-axis time display ('yes'/'no'). The default is 'no'.

    Returns
    -------
    None.

    """
   
    # ticks for the color map
    levels = [-2, -1, -0.1, 0, 0.1, 1, 2]
    # Create a continuous colormap
    cmap = plt.cm.get_cmap('PiYG_r')
    norm = plt.Normalize(min(levels), max(levels))
    sm = ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    
    range_max = len(slew)-1
    
    # patch rectangles according to the slew value
    y = ylim_min + 2.5 * panel_number
    height = ( ylim_max - ylim_min )  * 0.05
    for i, value in enumerate(slew):
        if i != range_max:
            x = time_array[i]       
            width = time_array[i+1]-time_array[i]
            rect = mpatches.Rectangle((x, y), width, height,
                                     edgecolor=None,
                                     facecolor=sm.to_rgba(value),
                                     fill=True, lw=0)
            c1 = ax.add_patch(rect)
            
    # dummy plot to properly insert time axis
    dummy_plot = ax.scatter(time_array, slew,
                                c=slew, cmap=cmap,
                                norm=norm, 
                                s=0)
    ax.annotate(slew_angle + '\N{DEGREE SIGN}', (time_array[-15], y+.25), \
                   fontsize=12)
    # beautify the axes
    ax.set_yticks([])
    if time =='yes':
        myFmt = dtformat.DateFormatter("%H:%M:%S")
        seconds = dtformat.MicrosecondLocator(interval=10000000)
        ax.xaxis.set_major_formatter(myFmt)
        ax.xaxis.set_minor_locator(seconds)
        ax.set_xticks([], minor=True)
    if cb_axis =='yes':
        cbar = plt.colorbar(sm, ticks=levels, boundaries=levels,
                            orientation = 'horizontal')
        labels = ['Back\nSlew', '1 Dipole\nBack', 
                  '', 'None', ' ',
                  '1 Dipole\nFront', 'Slew']
        new_pos = [-1.5, -0.5, -0.1, 0, 0.1, 0.5, 1.5]
        cbar.ax.xaxis.set_major_locator(FixedLocator(new_pos))
        cbar.ax.xaxis.set_major_formatter(FixedFormatter(labels))
        cbar.ax.tick_params(labelsize=18) 
   
    ax.patch.set_alpha(0.001)
    
    return    
 
