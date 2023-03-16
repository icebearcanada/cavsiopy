#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 17 14:53:01 2021
Plotting codes
@author: ceren
"""
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy
from matplotlib.collections import PolyCollection
from cartopy.mpl.patch import geos_to_path
from matplotlib.ticker import MultipleLocator
from mpl_toolkits.mplot3d import axes3d
import itertools
from matplotlib import rcParams
import cavsiopy.misc as misc

def indices_and_intervals(start_time, time_data, intval):
    """ Generates date for usage with plot titles and output filenames,
    determines the index of the last used data row and finds the number of
    intervals for cases we do not use the complete data set.

    Parameters
    ----------
    start_time: The beginning of the passage
    time_data : Time array for the whole passage
    intval : how many seconds do we want for intervals

    Returns
    -------
    date : date as a string
    last_data_index: index of the last used data row
    number_of_intervals: number of intervals for cases
    seconds: tick locator for seconds
    """
    SizeArr=len(time_data)
    last_data_index=SizeArr-(SizeArr%intval)
    number_of_intervals=int((SizeArr-(SizeArr%intval)) / intval)

    return last_data_index, number_of_intervals

def coverage(extent):
    """ Generates the plot parameters

    Parameters
    ----------
    extent: list
        [Lon_min, Lon_max, Lat_min, Lat_max]

    Returns
    -------
    central_lon, central lat: float
    pre-defined: float
            left (.01), width (.5), bottom (.01), height (.5)
    """

    central_lon = np.mean(extent[:2])
    central_lat = np.mean(extent[2:])
    left, width = .01, .5
    bottom, height = .01, .5

    return central_lon, central_lat, left, width, bottom, height

def set_3Dplot_limits(P, coord, how_far):
    """ Generates the limit of axes wrt the selected point on the ground for
    3D plots.

    Parameters
    ----------
        P: float
            coordinate of the ground point at x, y or z
        coord: float
            lat, lon or alt
        how_far: float
            how far is the axis limit from the point on the ground

    Returns
    -------
        coord_min: float
        coord_max: float

    Examples
    --------
    xmin, xmax = set_3Dplot_limits(Px, x)
    """

    if np.size(coord)>1:
        if min(coord)<P:
            coord_min = min(coord)-how_far
        else:
            coord_min = P - how_far

        if max(coord)<P:
            coord_max = P + how_far
        else:
            coord_max = max(coord) + how_far
    else:
        if coord<P:
            coord_min = coord-how_far
        else:
            coord_min = P - how_far

        if coord<P:
            coord_max = P + how_far
        else:
            coord_max = coord + how_far

    return coord_min, coord_max

def put_legend_fnt(ax, n, location, labelspace, anchorx, anchory, fontsize,
                labels={}, linestyles={},markers={}, colors={}, edgecolors={}):

    """
    Inserts legend to figure for cases, which need customized symbols.

    Parameters
    ----------
    ax: str
        axis
    n: float
        how many columns should the legend have
    location: str
        location of the legend
    labelspace: float
        space between the labels (values can be between 0 and 1).
    anchorx, anchory: float
        coordinate of the lower left corner of legend box
        (values can be between 0 and 1)

    The following are for the customization of the symbols we want to display:

    labels: str
        labels to display in the legend
    linestyles: str
        linestyles used in plotting the parameters
    markers: str
        markers used in plotting the parameters
    colors: str
        colors used in plotting the parameters.
        can be also given in RGB
    edgecolors: str
        edgecolors used in plotting the parameters

    Returns
    -------
        None
    """
    from matplotlib.legend_handler import HandlerBase
    class MarkerHandler(HandlerBase):
        def create_artists(self, legend, tup, xdescent, ydescent,
                            width, height, fontsize, trans):
            return [ax.scatter([width/2], [height/2],
                           marker=tup[2], color=tup[1], edgecolor=tup[0],
                           s=80, transform=trans)]

    # just use the default location if there is no anchorx and anchory
    if anchorx == None and anchory == None:
        ax.legend(list(zip(edgecolors, colors, markers)), labels,
                  handler_map={tuple:MarkerHandler()}, fontsize=fontsize,
                  labelspacing=labelspace, handletextpad=0, handlelength=2,
                  ncol=n, loc=location, frameon = False)

    else:
        ax.legend(list(zip(edgecolors, colors, markers)), labels,
                  handler_map={tuple:MarkerHandler()}, fontsize=fontsize,
                  labelspacing=labelspace, handlelength=2,
                  handletextpad=0, ncol=n, loc=location,
                  bbox_to_anchor=(anchorx,anchory),
                  frameon = False)

    return

def plotter_with_map_on_latlon_grid(fig, spec, extent, x, y, XV, n, m,
                 time_array, Ox, Oy, end_ind, d, dec,
                 x_name, y_name, title, kword, kword2, spec_col):

    """ Plots 2D maps for instrument pointing direction

    Parameters
    ----------
    fig: int
        figure number
    spec: class
        figure grid properties
        ex: spec = gridspec.GridSpec(ncols=5, nrows=1, figure = fig)
    x: float
       longitude
    y: float
        latitude
    XV: float
        instrument pointing vector to be plotted
    n: int
        column number of the East (X) component of the inst. pointing vector
        ex: for enu : n = 0
    m: int
        column number of the North (X) component of the inst. pointing vector
        ex: for enu : m = 1
    time_array: datetime
        time array of spacecraft pass
    Ox: float
        ground target longitude
    Oy: float
        ground target latitude
    end_ind: int
        index of the last data point
    d: int
        how many pointing vectors would you like to see
    dec: float
        offset of time info from the trajectory (in degrees lon)
    x_name: str
        x axis label
    y_name: str
        y axis label
    kword: str
        hide or display x and y axis labels
        can be 'left_on', 'bottom_on', 'right_on', else (everything is on)
    kword2: str
        'legend_on' : plots the legend to the upper right of the fig.
    spec_col: int
        column number to plot

    Returns
    -------
        None
    """

    central_lon, central_lat, left, width, bottom, height = coverage(extent)
    ax = fig.add_subplot(spec[0, spec_col], projection=\
                         ccrs.Orthographic(central_longitude=central_lon,
                                       central_latitude=central_lat))


    ax.set_extent(extent,ccrs.PlateCarree())
    ax.add_feature(cartopy.feature.OCEAN, color='white', alpha=1, zorder=0)
    ax.add_feature(cartopy.feature.LAND, edgecolor='white',
    color='silver', alpha=0.3, zorder=10)
    ax.add_feature(cartopy.feature.LAKES, color='white', alpha=1, zorder=0)

    # plot the RRI pointing direction.
    ax.quiver(x[::d], y[::d],
              XV[::d,n], XV[::d,m],  headaxislength = 7.5,
              transform=ccrs.PlateCarree(), label = 'RRI',
              color='black', scale=6, width=.006, axes=ax, zorder=20)

    ax.plot(x,y, transform=ccrs.PlateCarree(),
                color='#40a368', linestyle = 'solid')
    ax.scatter(Ox, Oy, transform=ccrs.PlateCarree(),
               marker='*',  edgecolor='k', color= '#fdaa48', s=200,
               label='Ottawa')
    ax.scatter(x[0], y[0], transform=ccrs.PlateCarree(),
               marker = 'x', color='magenta', s=150,
               label ='Beginning of the pass')

    for i in range(0, end_ind, 45):
        ax.text(x[i]+0.5, y[i], (time_array[i].strftime("%H:%M:%S")),
            transform=ccrs.PlateCarree(),
              c='b',size=10, zorder=20)
    if kword2 =='legend_on':
        ax.legend(loc='upper right')


    gl = ax.gridlines(draw_labels=True, dms=True,
    x_inline=False, y_inline=False,
    linestyle=':', zorder=15)
    if kword == 'left_on':
        gl.top_labels = gl.right_labels = False
    elif kword == 'bottom_on':
        gl.top_labels = gl.right_labels = gl.left_labels = False
    elif kword == 'right_on':
        gl.top_labels = gl.left_labels = False
    else:
        gl.top_labels = gl.right_labels = True
    gl.xlabel_style = {'size': 14, 'color': 'black'}
    gl.ylabel_style = {'size': 14, 'color': 'black'}

    plt.xticks(fontsize= 14)
    plt.yticks(fontsize= 14)
    plt.title(title)

    return

class MyAxes3D(axes3d.Axes3D):
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

def vector_direction_plotter_ground_trajectory_3D(title, time_array, end_ind, ind1, 
                             Px, Py, Pz, x, y, z, arrlen, loc, target, *V, 
                             vector_colors={}, linestyles={}, linewidth = {},
                             labels={}, markers={}, colors={}, edgecolors={},
                             arrowhead = {},
                             sct_kwargs = {}):
    
    """v: vx, vy, vz """

    a=0
    fig = plt.figure(figsize=(12, 12))
    ax=plt.axes(projection='3d')
    for vec in V:
        if np.size(x)>1:
            # plt.suptitle(title, fontsize=20, y=0.95)
            a+=1
            i=a-1
            ax.quiver(x[::ind1], y[::ind1], z[::ind1],
                          vec[:,0][::ind1], vec[:,1][::ind1], vec[:,2][::ind1],
                          length = arrlen[i], color = vector_colors[i], 
                          linestyle = linestyles[i], linewidth = linewidth[i],
                          arrow_length_ratio = arrowhead[i])
            
            # linewidth = 2.5, arrow_length_ratio = 0.35
            Px_vec = [Px]*len(x)
            Py_vec = [Py]*len(y)
            Pz_vec = [0.07]*len(z)
            for i in range(0, len(x), ind1):
                lon = [Px_vec[i], x[i]]
                lat = [Py_vec[i], y[i]]
                alt = [Pz_vec[i], z[i]]
                ax.plot(lon,lat,alt,linestyle = 'dashed', color = '#95d0fc', 
                        alpha = .7, linewidth=1)
                
            xmin, xmax = set_3Dplot_limits(Px, x, 5)
            ymin, ymax = set_3Dplot_limits(Py, y, 5)
            zmin, zmax = set_3Dplot_limits(Pz, z, 5)
            ax.set_xlim(xmin, xmax)
            ax.set_ylim(ymin, ymax)
            ax.set_zlim(0, np.round(max(z)+0.5))
            # ax.scatter3D(Px, Py, 0.07, **sct_kwargs)
            plt.xticks(fontsize= 18)
            plt.yticks(fontsize= 18)

        else:
            # plt.suptitle(title, fontsize=14, y=0.97)
            a+=1
            i=a-1
            ax.quiver(x, y, z,
                  vec[0], vec[1], vec[2], 
                  length=arrlen[i], color=vector_colors[i])
            coord_info="".join(['Lat: ', str(y), '\N{DEGREE SIGN}\n',
                              'Lon: ', str(x), '\N{DEGREE SIGN}\n',
                              'Alt: ', str(z), 
                              '(\N{DEGREE SIGN} x 100 km/\N{DEGREE SIGN})'])
            ax.text2D(0.05, 0.95, coord_info,  fontsize=14, va='top', ha='left',
                    rotation_mode='anchor',
                    transform=ax.transAxes)
            Pz= [4]
            lon = [Px, x]
            lat = [Py, y]
            alt = [Pz, z]
            # ax.plot(lon,lat,alt,'k--',alpha=1, linewidth=0.5)
            # ax.scatter3D(Px, Py, 0.07, **sct_kwargs)
            plt.xticks(fontsize= 18)
            plt.yticks(fontsize= 18)
   
    misc.put_legend_fnt(ax, 5, loc, 1.5, 0.5, 1, 12,
                     labels=labels, linestyles= linestyles, 
                     markers=markers, colors=colors, 
                     edgecolors=edgecolors)
    
    for i in range(0, len(x), ind1):
        # lon = [x[i], x[i]]
        # lat = [y[i], y[i]]
        # alt = [0, z[i]]
        if y[i] < Py:
            ax.scatter(x[i], y[i], 0.1,
            marker='H', edgecolor='black', color='black', s=60)

        else:
            ax.scatter(x[i], y[i], 0.1,
            marker='H', edgecolor='black', color='None', s=60)

    ax.plot3D(x, y, 0.01, c = 'black', linewidth= 1.25)

    #  insert the legend
    misc.put_legend_fnt(ax, 6, loc, 0.2, 0.525, 0.90, 13,
                     labels=labels, linestyles= linestyles, 
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
        ax.text(Px-.5, Py, 0.5, target, c='k', size=14, zorder=20)
    else:
        ax.text(Px+.5, Py, 0.5, target, c='k', size=14, zorder=20)
    
    import itertools
    
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
    geoms = [target_projection.project_geometry(geom, feature.crs) for geom in geoms]
    
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
    from matplotlib import rcParams
    rcParams['axes.labelpad'] = 15  
    ax.set_xlabel('Geographic Longitude (\N{DEGREE SIGN})', fontsize=18)
    ax.set_ylabel('Geographic Latitude (\N{DEGREE SIGN})', fontsize=18)
    ax.zaxis.set_rotate_label(False)  # disable automatic rotation
    ax.set_zlabel('Altitude ( x 100 km)', rotation=90, fontsize=18)

    
    ax.tick_params(axis='x', which='major', labelsize = 18)
    ax.tick_params(axis='y', which='major', labelsize=18)
    ax.tick_params(axis='z', which='major', labelsize=18)
    # fig.tight_layout(rect=[0, 0, 0.9,1])
    plt.tight_layout()
    # fig.subplots_adjust(right=0.65)
    f= title +'.png'
    # plt.savefig(f)
    # plt.show()
    # plt.close(fig)

    return fig

def vector_direction_plotter_connect_to_target(title, time_array, end_ind, ind1, 
                             Px, Py, Pz, x, y, z, arrlen, loc, target, *V, 
                             vector_colors={}, linestyles={}, linewidth = {},
                             labels={}, markers={}, colors={}, edgecolors={},
                             arrowhead = {},
                             sct_kwargs = {}):
    
    """v: vx, vy, vz """

    a=0
    fig = plt.figure(figsize=(12, 12))
    ax=plt.axes(projection='3d')
    for vec in V:
        if np.size(x)>1:
            # plt.suptitle(title, fontsize=20, y=0.95)
            a+=1
            i=a-1
            ax.quiver(x[::ind1], y[::ind1], z[::ind1],
                          vec[:,0][::ind1], vec[:,1][::ind1], vec[:,2][::ind1],
                          length = arrlen[i], color = vector_colors[i], 
                          linestyle = linestyles[i], linewidth = linewidth[i],
                          arrow_length_ratio = arrowhead[i])
            
            # linewidth = 2.5, arrow_length_ratio = 0.35
            Px_vec = [Px]*len(x)
            Py_vec = [Py]*len(y)
            Pz_vec = [0.07]*len(z)
            for i in range(0, len(x), ind1):
                lon = [Px_vec[i], x[i]]
                lat = [Py_vec[i], y[i]]
                alt = [Pz_vec[i], z[i]]
                ax.plot(lon,lat,alt,linestyle = 'dashed', color = '#95d0fc', 
                        alpha = .7, linewidth=1)
                
            xmin, xmax = set_3Dplot_limits(Px, x, 5)
            ymin, ymax = set_3Dplot_limits(Py, y, 5)
            zmin, zmax = set_3Dplot_limits(Pz, z, 5)
            ax.set_xlim(xmin, xmax)
            ax.set_ylim(ymin, ymax)
            ax.set_zlim(0, np.round(max(z)+0.5))
            # ax.scatter3D(Px, Py, 0.07, **sct_kwargs)
            plt.xticks(fontsize= 18)
            plt.yticks(fontsize= 18)

        else:
            # plt.suptitle(title, fontsize=14, y=0.97)
            a+=1
            i=a-1
            ax.quiver(x, y, z,
                  vec[0], vec[1], vec[2], 
                  length=arrlen[i], color=vector_colors[i])
            coord_info="".join(['Lat: ', str(y), '\N{DEGREE SIGN}\n',
                              'Lon: ', str(x), '\N{DEGREE SIGN}\n',
                              'Alt: ', str(z), 
                              '(\N{DEGREE SIGN} x 100 km/\N{DEGREE SIGN})'])
            ax.text2D(0.05, 0.95, coord_info,  fontsize=14, va='top', ha='left',
                    rotation_mode='anchor',
                    transform=ax.transAxes)
            Pz= [4]
            lon = [Px, x]
            lat = [Py, y]
            alt = [Pz, z]
            # ax.plot(lon,lat,alt,'k--',alpha=1, linewidth=0.5)
            # ax.scatter3D(Px, Py, 0.07, **sct_kwargs)
            plt.xticks(fontsize= 18)
            plt.yticks(fontsize= 18)
   
    # misc.put_legend_fnt(ax, 5, loc, 1.5, 0.5, 1, 12,
    #                  labels=labels, linestyles= linestyles, 
    #                  markers=markers, colors=colors, 
    #                  edgecolors=edgecolors)
    
    # if y[0] < Py:
    #     ax.scatter(x[0], y[0], 0.1,
    #     marker='H', edgecolor='black', color='black', s=30)

    # else:
    #     ax.scatter(x[0], y[0], 0.1,
    #     marker='H', edgecolor='black', color='None', s=30)
        
    # if y[-1] < Py:
    #     ax.scatter(x[-1], y[-1], 0.1,
    #     marker='H', edgecolor='black', color='black', s=30)

    # else:
    #     ax.scatter(x[-1], y[-1], 0.1,
    #     marker='H', edgecolor='black', color='None', s=30)

    # ax.plot3D(x, y, 0.01, c = 'black', linewidth= 1.25)

    # #  insert the legend
    # misc.put_legend_fnt(ax, 6, loc, 0.2, 0.525, 0.90, 13,
    #                  labels=labels, linestyles= linestyles, 
    #                  markers=markers, colors=colors, 
    #                  edgecolors=edgecolors)  
      
    # time1 = time_array[0].strftime('%H:%M:%S')
    # time2 = time_array[-1].strftime('%H:%M:%S')
    
    # #  insert the start and end times
    # ax.text(x[0]-0.6, y[0]-0.95, 0.3, time1, size=14)   
    # ax.text(x[-1]-0.5, y[-1]-0.2, 0.3, time2, size = 14)
    
    # target
    ax.scatter3D(Px, Py, 0.01, marker = 'o', edgecolor="black", 
              facecolor='lightgrey', s= 180)
    if x[0] < Px:
        ax.text(Px-.5, Py, 0.5, target, c='k', size=14, zorder=20)
    else:
        ax.text(Px+.5, Py, 0.5, target, c='k', size=14, zorder=20)
    
    import itertools
    
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
    geoms = [target_projection.project_geometry(geom, feature.crs) for geom in geoms]
    
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
    # ax.add_collection3d(lc, zs=0)
    
    # arrange the view
    # nadir
    ax.view_init(azim=-90, elev=10)
    plt.tight_layout()   
    plt.close(proj_ax.figure)
    from matplotlib import rcParams
    rcParams['axes.labelpad'] = 15  
    ax.set_xlabel('Geographic Longitude (\N{DEGREE SIGN})', fontsize=18)
    ax.set_ylabel('Geographic Latitude (\N{DEGREE SIGN})', fontsize=18)
    ax.zaxis.set_rotate_label(False)  # disable automatic rotation
    ax.set_zlabel('Altitude ( x 100 km)', rotation=90, fontsize=18)

    plt.axis('off')
    plt.grid(b=None)
    
    # ax.tick_params(axis='x', which='major', labelsize = 18)
    # ax.tick_params(axis='y', which='major', labelsize=18)
    # ax.tick_params(axis='z', which='major', labelsize=18)
    # fig.tight_layout(rect=[0, 0, 0.9,1])
    plt.tight_layout()
    # fig.subplots_adjust(right=0.65)
    f= title +'.png'
    # plt.savefig(f)
    # plt.show()
    # plt.close(fig)

    return fig

def pointing_through_3D_trajectory(start_time,input_time, end_ind, step,
                 Px, Py, Pz, x, y, z, pointing, arrlen, x_avg, y_avg, z_avg,
                 loc, *V,vector_colors={}, labels={}, linestyles={}, markers={},
                 colors={}, edgecolors={}):

    a=0
    fig = plt.figure(figsize=(12.06, 10.79))
    ax=plt.axes(projection='3d')
    for vec in V:
        if np.size(V)>1:
            lim1 = len(x)
            lim2 = len(vec)
            start_time=input_time[0]
            # title=" ".join(['RRI Cross-Dipoles and Antenna Boresight on',
            #                 start_time.strftime('%d-%b-%Y')])
            # plt.suptitle(title, fontsize=16, y=0.95)
            a+=1
            i=a-1
# =============================================================================
# 3D quiver on the trajectory
# =============================================================================
            if np.less(lim1,lim2)==False:
                ax.quiver(x[0:lim2:step], y[0:lim2:step], z[0:lim2:step],
                          vec[:,0][::step], vec[:,1][::step], vec[:,2][::step],
                          length=arrlen[i], color=vector_colors[i],
                          linestyle=linestyles[i])
            else:
                ax.quiver(x[::step], y[::step], z[::step],
                          vec[:,0][::step], vec[:,1][::step], vec[:,2][::step],
                          length=arrlen[i], color=vector_colors[i],
                          linestyle=linestyles[i])

        else:
            a+=1
            i=a-1
            ax.quiver(x, y, z,
                  vec[0], vec[1], vec[2],
                  length=arrlen[i], color=vector_colors[i])
            coord_info="".join(['Time (UTC): ',
                                input_time.strftime('%H:%M:%S'),'\n',
                                'Lat: ', str(y), '\N{DEGREE SIGN}\n',
                                'Lon: ', str(x), '\N{DEGREE SIGN}\n',
                                'Alt: ', str(z), '( x 111 km)'])
            ax.text(0.05, 0.95, 0.2, coord_info, fontsize=14, va='top', ha='left',
                    rotation_mode='anchor', transform=ax.transAxes)


    cond1 = np.isclose(0, pointing[0,2], atol=0.05)
    cond2 = np.isclose(0, pointing[0,0], atol=0.05) and \
            np.isclose(0, pointing[0,1], atol=0.05)
    cond3 = pointing[0,2] < 0
# =============================================================================
# 2D quiver on the ground
# =============================================================================
    for i in range(0, len(x), step):
         if cond1==True or cond3 ==True:
             ax.quiver(x[::step], y[::step], 0,
                       pointing[:,0][::step], pointing[:,1][::step], 0,
                       length=1, color='black',
                       linewidths=0.7, arrow_length_ratio=0.1, zorder= 50)
         elif cond2==True:
             ax.quiver(x[::step], y[::step], z[::step],
                   pointing[:,0][::step], pointing[:,1][::step],
                   pointing[:,2][::step],
                   length=0.75, color='black',
                   linewidths=1, arrow_length_ratio=0.2, zorder= 50)

    # determine plot limits
    xmin, xmax = set_3Dplot_limits(Px, x, 1)
    ymin, ymax = set_3Dplot_limits(Py, y, 1)
    zmin, zmax = set_3Dplot_limits(Pz, z, 1)
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_zlim(0, np.round(max(z)+0.5))

    # ground target
    ax.scatter3D(Px, Py, 0.01, marker = '*', edgecolor="black",
              facecolor='lightgrey', s= 180)
    ax.text(Px, Py, 0.02, 'Ottawa', c='k', size=14, zorder=20)

    for i in range(0, len(x), step):
        lon = [x[i], x[i]]
        lat = [y[i], y[i]]
        alt = [0, z[i]]
        if y[i] < 45.4215:
            ax.scatter(x[i], y[i], 0.1,
            marker='H', edgecolor='black', color='black', s=60)
            ax.plot(lon,lat,alt,'k--',alpha=1, linewidth=1)
        else:
            ax.scatter(x[i], y[i], 0.1,
            marker='H', edgecolor='black', color='None', s=60)
            ax.plot(lon,lat,alt,'k--',alpha=1, linewidth=1)

    time1=start_time.strftime('%H:%M:%S')
    time2=input_time[-1].strftime('%H:%M:%S')

    #  insert the start and end times
    ax.text(x[0]-0.6, y[0]-0.95, 0.3, time1, size=14)
    ax.text(x[-1]-0.5, y[-1]-0.2, 0.3, time2, size = 14)

    #  insert the legend
    misc.put_legend_fnt(ax, 5, loc, 0.2, 0.525, 0.90, 13,
                     labels=labels, linestyles= linestyles,
                     markers=markers, colors=colors,
                     edgecolors=edgecolors)

    import itertools

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
    geoms = [target_projection.project_geometry(geom, feature.crs) for geom in geoms]

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

    rcParams['axes.labelpad'] = 15
    ax.set_xlabel('Geographic Longitude (\N{DEGREE SIGN})', fontsize=14)
    ax.set_ylabel('Geographic Latitude (\N{DEGREE SIGN})', fontsize=14)
    ax.zaxis.set_rotate_label(False)  # disable automatic rotation
    ax.set_zlabel('Altitude ( x 100 km)', rotation=90, fontsize=14)

    ax.tick_params(axis='x', which='major', labelsize=14)
    ax.tick_params(axis='y', which='major', labelsize=14)
    ax.tick_params(axis='z', which='major', labelsize=14)

    plt.tight_layout()
    plt.show()

    return

def plotter_with_map(fig, spec, extent, x, y, XV, n, m,
                 time_array, Ox, Oy, end_ind, d, dec,
                 x_name, y_name, title, kword, kword2, spec_col):

    """ Plots 2D maps for RRI pointing direction

    Parameters
    ----------
    fig: figure number
    spec: figure grid properties
        ex: spec = gridspec.GridSpec(ncols=5, nrows=1, figure = fig)
    x, y = Lon, Lat: longitude and latitude
    XV: RRI pointing vector to be plotted
    n: column number of the East (X) component of the RRI pointing vector
        ex: for enu : n = 0
    m: column number of the North (X) component of the RRI pointing vector
        ex: for enu : m = 1
    time_array: time array of RRI pass
    Ox: ground target longitude
    Oy: ground target latitude
    end_ind: index of the last data point
    d: how many RRI vectors would you like to see
    dec: offset of time info from the trajectory (in degrees lon)
    x_name: x axis label
    y_name: y axis label
    kword: hide or display x and y axis labels
        'left_on', 'bottom_on', 'right_on', else (everything is on)
    kword2: 'legend_on' : plots the legend to the upper right of the fig.
    spec_col: column number to plot

    Returns
    -------
    figure as a subplot specified by gridspec"""

    central_lon, central_lat, left, width, bottom, height = coverage(extent)
    ax = fig.add_subplot(spec[0, spec_col], projection=\
                         ccrs.Orthographic(central_longitude=central_lon,
                                       central_latitude=central_lat))


    ax.set_extent(extent,ccrs.PlateCarree())
    ax.add_feature(cartopy.feature.OCEAN, color='white', alpha=1, zorder=0)
    ax.add_feature(cartopy.feature.LAND, edgecolor='white',
    color='silver', alpha=0.3, zorder=10)
    ax.add_feature(cartopy.feature.LAKES, color='white', alpha=1, zorder=0)

    # plot the RRI pointing direction.
    ax.quiver(x[::d], y[::d],
              XV[::d,n], XV[::d,m],  headaxislength = 5,
              transform=ccrs.PlateCarree(), label = 'RRI',
              color='black', scale=8.75, width=.009, axes=ax, zorder=20)

    ax.plot(x,y, transform=ccrs.PlateCarree(),
                color='#40a368', linestyle = 'solid')
    ax.scatter(Ox, Oy, transform=ccrs.PlateCarree(),
               marker='*',  edgecolor='k', color= '#fdaa48', s=200,
               label='Ottawa')
    ax.scatter(x[0], y[0], transform=ccrs.PlateCarree(),
               marker = 'x', color='magenta', s=150,
               label ='Beginning of the pass', zorder = 20)

    for i in range(0, end_ind, 90):
        ax.text(x[i]+0.5, y[i], (time_array[i].strftime("%H:%M:%S")),
            transform=ccrs.PlateCarree(),
              c='b',size=10, zorder=20)
    if kword2 =='legend_on':
        ax.legend(loc='upper right')


    gl = ax.gridlines(draw_labels=True, dms=True,
    x_inline=False, y_inline=False,
    linestyle=':', zorder=15)
    if kword == 'left_on':
        gl.top_labels = gl.right_labels = False
    elif kword == 'bottom_on':
        gl.top_labels = gl.right_labels = gl.left_labels = False
    elif kword == 'right_on':
        gl.top_labels = gl.left_labels = False
    else:
        gl.top_labels = gl.right_labels = True
    gl.xlabel_style = {'size': 14, 'color': 'black'}
    gl.ylabel_style = {'size': 14, 'color': 'black'}

    plt.xticks(fontsize= 14)
    plt.yticks(fontsize= 14)
    plt.title(title, fontsize = 16, y = 1)

    return ax
