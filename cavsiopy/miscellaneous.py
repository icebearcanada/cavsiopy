#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This module contains functions for: 
1- legend insertion to figures (put_legend, put_legend_fnt)
2- combining images (combine_horizontal, combine_vertical)
3- marking specific times on map(mark_on_map, mark_beginning_on_map) 
    and altitude plots(mark_beginning, mark_altitude_plots)
4- plotting the direction of the spacecraft
5- finding index of desired values in parameters: find_index

.. toctree::
  :maxdepth: 2   
  combine_horizontal
  combine_vertical
  coverage
  find_index
  indices_and_intervals
  mark_altitude_plots
  mark_beginning
  mark_beginning_on_map
  mark_on_map
  put_legend_fnt
  sc_direction_plotter
  set_3Dplot_limits
@author: ceren

"""
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.tri as tri

def indices_and_intervals(start_time, time_data, intval):
    """ Generates date for usage with plot titles and output filenames,
    determines the index of the last used data row and finds the number of
    intervals for cases we do not use the complete data set.

    Parameters
    ----------
    start_time: datetime.datetime
        The beginning of the passage
    time_data : datetime.datetime 
        Time array for the whole passage
    intval : int
        how many seconds do we want for intervals

    Returns
    -------
    date : str
        date as a string.
    last_data_index: int
        index of the last used data row
    number_of_intervals: int
        number of intervals for cases
    seconds: int
        tick locator for seconds
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
    central_lon : float
        Central longitude (degrees)
    central lat : float
        Central latitude (degrees)
    pre-defined : float
            left (.01), width (.5), bottom (.01), height (.5)
    """

    central_lon = np.mean(extent[:2])
    central_lat = np.mean(extent[2:])

    return central_lon, central_lat

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

# =============================================================================
# combine images
# =============================================================================
import PIL
def combine_horizontal(f, fnew):
    """
    Horizontally combines images with arbitrary image size as a new image and
    saves the new image with the given file name.
    
    Parameters:
    -----------
    f: Name of image files, list
    fnew: Output file name
    
    Returns:
    -----------
    Saves horizontally combined image as fnew to the directory
    Reference: https://stackoverflow.com/a/30228789/15165141"""
    
    list_im = f
    # open all images in the list
    imgs    = [ PIL.Image.open(i) for i in list_im ]
    # find the minimum size in given images
    min_shape = sorted( [(np.sum(i.size), i.size ) for i in imgs])[0][1]
    # resize all images to match the image with minimum size
    # and stack horizontally
    imgs_comb = np.hstack([np.asarray(i.resize(min_shape)) for i in imgs])
    
    # Create an image memory from an object 
    imgs_comb = PIL.Image.fromarray(imgs_comb)
    # save new image
    imgs_comb.save(fnew)
    
def combine_vertical(f, fnew):
    """
    Vertically combines images with arbitrary image size as a new image and
    saves the new image with the given file name.
    
    Parameters:
    -----------
    f: Name of image files, list
    fnew: Output file name
    
    Returns:
    -----------
    Saves vertically combined image as fnew to the directory
    Reference: https://stackoverflow.com/a/30228789/15165141 """
    
    list_im = f
    # open all images in the list
    imgs    = [ PIL.Image.open(i) for i in list_im ]
    # find the minimum size in given images    
    min_shape = sorted( [(np.sum(i.size), i.size ) for i in imgs])[0][1]
    # resize all images according to the minimum size
    imgs_comb = np.vstack([np.asarray(i.resize(min_shape)) for i in imgs])
    
    # Create an image memory from an object 
    imgs_comb = PIL.Image.fromarray(imgs_comb)
    # save new image
    imgs_comb.save(fnew)

# =============================================================================
# find the indices of close values for specific angle values
# =============================================================================
def find_index(data, closeto, tol, verbose = 'off'):
    
    """ 
    Finds the index of numbers that are close to the value we are searching for
    with a tolerance 'tol'.
 
    Parameters:
    -----------
    data: np.ndarray[float]
        input data for inspection, array
    closeto: float
        The value we are searching for
    tol: float
        Amount of tolerance
    verbose : string, optional
        Verbose. The default is 'off'.
    
    Returns:
    --------
    ind_close: list
        Index of the close values
    """
    
    # initialize the index array of close points
    ind_close=[]
    
    # are there any close numbers within a tolerance of tol?
    boolean = np.isclose(closeto, data, atol=tol)
    # no close values
    if np.all(boolean==0):
        if verbose != 'off':
            print("".join(['No values found close to', str(closeto)]))
            ind_close=np.empty(0)
    
    # if there are close values
    else:
        # find how many close values
        how_many_close_vals=np.count_nonzero(boolean)
        # or np.sum(boolean!=0)
        # append the index of close values in ind_close
        ind_close=np.where(boolean==1)
        if verbose != 'off':
            print(' '.join(['found', str(how_many_close_vals),\
                            'values close to', str(closeto), 'within',\
                            str(tol), 'tolerance.\n']))
        # col gives the number of close values   
        row, col = np.shape(ind_close)
       
    return ind_close

# =============================================================================
# mark the beginning of the pass and insert the location of the selected city 
# to the plot
# =============================================================================

def mark_on_map(ax, Px, Py, color='cyan', edgecolor='black', marker='*', 
                transform=ccrs.PlateCarree(), **kwargs):
     
    """
    Function to mark locations on maps.
    
    Parameters
    ----------
    ax : Axes
        axes.
    Px : float
        Longitude of the point to be marked (degrees).
    Py : float
        Latitude of the point to be marked (degrees).
    color : str, optional
        Color of the point to be marked. Default is 'cyan'.
    edgecolor : str, optional
        Edge color of the point to be marked. Default is 'black'.
    transform : object, optional
        cartopy projection. The default is ccrs.PlateCarree().
    **kwargs : dict
        Additional keyword arguments for scatter plot.

    Returns
    -------
    None.
    """

    # add desired ground location to plot
    ax.scatter(Px, Py, color=color, edgecolor = edgecolor, marker=marker, 
               s=150, transform=transform, zorder= 20)
                  
    return

def mark_beginning_on_map(ax, x, y, z, transform=ccrs.PlateCarree()):
    
    """
    Function to mark the beginning of the pass and tp show the direction of the
    spacecraft with an arrow on map projections.

    Parameters
    ----------
    ax : Axes
        axes.
    x : float
        Longitude of the beginning of the trajectory (degrees).
    y : float
        Latitude of the beginning of the trajectory (degrees).
    z : float, optional
        Latitude of the beginning of the trajectory (km).
    transform : object, optional
        cartopy projection. The default is ccrs.PlateCarree().

    Returns
    -------
    None.

    """
    
    sc_direction_plotter(ax, x, y, z, transform = transform, zorder= 20)
    
    return

def mark_altitude_plots(ax, Px, Py, color='cyan', edgecolor = 'k', marker='*', 
                        **kwargs):
    """
    Function to mark locations on altitude plots.
    
    Parameters
    ----------
    ax : Axes
        axes.
    Px : float
        Longitude of the point to be marked (degrees).
    Py : float
        Altitude of the point to be marked (km).
    color : str, optional
        Color of the point to be marked. Default is 'cyan'.
    edgecolor : str, optional
        Edge color of the point to be marked. Default is 'black'.
    marker : str, optional
        Marker symbol. The default is start '*'.
    **kwargs : dict
        Additional keyword arguments for scatter plot.

    Returns
    -------
    None.
    """
     
     
    # add desired ground location to plot
    ax.scatter(Px, Py, color = color, edgecolor = edgecolor, marker=marker, 
               s=150, zorder= 20)
           
    return

def mark_beginning(ax, x, y, z):
    """
    Plot to mark the beginning of the pass and show the direction of the
    spacecraft with an arrow for altitude plots.

    Parameters
    ----------
    ax : Axes
        axes.
    x : float
        Longitude of the beginning of the trajectory (degrees).
    y : float
        Latitude of the beginning of the trajectory (degrees).
    z : float
        Altitude of the beginning of the trajectory (km).

    Returns
    -------
    None.

    """
    
    sc_direction_plotter(ax, x, y, z, zorder= 20)
    
    return
 

def sc_direction_plotter(ax, Lon, Lat, Alt):
    """
    Plots an arrow to depict the direction of the spacecraft velocity.
    Prints out trajectory information: altitude increasing/decreasing, 
    going towards SE, NW, SW, NE.

    Parameters
    ----------
    ax : axes
        Axes object of matplotlib. 
    Lon : numpy.ndarray[float]
        Longitude of the beginning of the trajectory (degrees).
    Lat : numpy.ndarray[float]
        Latitude of the beginning of the trajectory (degrees).
    Alt : numpy.ndarray[float]
        Altitude of the beginning of the trajectory (km).

    Returns
    -------
    None.

    """

    L=np.size(Alt)
    # calculate the differences to determine where the spacecraft is heading to
    # diff_alt  > 0 = descending in altitude
    # diff_alt  < 0 = ascending in altitude
    diff_alt= Alt[0]-Alt[5]
    # diff_lon  > 0 = heading west (as we are westward from Greenwhich)
    # diff_lon  < 0 = heading east 
    diff_lon= Lon[0]-Lon[5]
    # diff_lat  > 0 = heading south (as we are northward from equator)
    # diff_lat  < 0 = heading north
    diff_lat= Lat[0]-Lat[5]
    # determine lat_min and max to put the arrows to the start of the passage
    Lat_min, Lat_max=np.min(Lat), np.max(Lat)
    # specify arrrow length and head size using constants
    # c1: arrow length to denote the altitudinal direction
    # c2: arrow length to denote the lat-long plane direction
    # c3: head size for altitudinal direction
    # c4: head size for lat-long plane direction
    c1=1 ; c2= -1.5; c3=0.3; c4=0.3;
    
    # arrows to denote the changes in altitude, , color: blue
    if (diff_alt>0 and diff_lon>0 and diff_lat>0):
        # ax.arrow(Lon[0],Lat_max, 0, -c1,
        # head_width=c3, head_length=c3, fc='b', ec='b', zorder=25)
        print("Altitude decreasing-Going towards: SW")
    elif (diff_alt>0 and diff_lon<0 and diff_lat>0):
        # ax.arrow(Lon[0],Lat_max, 0, -c1,
        # head_width=c3, head_length=c3, fc='b', ec='b', zorder=25)
        print("Altitude decreasing-Going towards: SE")
    elif (diff_alt>0 and diff_lon<0 and diff_lat<0):
        # ax.arrow(Lon[L-1],Lat_max, 0, c1,
        # head_width=c3, head_length=c3, fc='b', ec='b', zorder=25)
        print("Altitude decreasing-Going towards: NE")
    elif (diff_alt>0 and diff_lon>0 and diff_lat<0):
        # ax.arrow(Lon[L-1],Lat_max, 0, c1,
        # head_width=c3, head_length=c3, fc='b', ec='b', zorder=25)
        print("Altitude decreasing-Going towards: NW")
    elif (diff_alt<0 and diff_lon>0 and diff_lat>0):
        # ax.arrow(Lon[0],Lat_max, 0, -c1,
        #      head_width=c3, head_length=c3, fc='b', ec='b', zorder=25)
        print("Altitude increasing-Going towards: SW")
    elif (diff_alt<0 and diff_lon<0 and diff_lat>0):
        # ax.arrow(Lon[0],Lat_max, 0, -c1,
        # head_width=c3, head_length=c3, fc='b', ec='b', zorder=25)
        print("Altitude increasing-Going towards: SE")
    elif (diff_alt<0 and diff_lon<0 and diff_lat<0):
        # ax.arrow(Lon[L-1],Lat_max, 0, c1,
        # head_width=c3, head_length=c3, fc='b', ec='b', zorder=25)
        print("Altitude increasing-Going towards: NE")
    elif (diff_alt<0 and diff_lon>0 and diff_lat<0):
        # ax.arrow(Lon[L-1],Lat_max, 0, c1,
        # head_width=c3, head_length=c3, fc='b', ec='b', zorder=25)
        print("Altitude increasing-Going towards: NW")
    
    # arrows for denoting the direction, color: black
    
    if (diff_alt>0 and diff_lon>0 and diff_lat>0):
        ax.arrow(Lon[0],Lat[0]-.2, diff_lon*c2, diff_lat*c2,
             head_width=c4, head_length=c4, fc='k', ec='k', 
             transform=ccrs.PlateCarree(), zorder=25)
        
    elif (diff_alt>0 and diff_lon<0 and diff_lat>0):
        ax.arrow(Lon[0],Lat[0]-.2, diff_lon*c2, diff_lat*c2,
        head_width=c4, head_length=c4, fc='k', ec='k', 
        transform=ccrs.PlateCarree(), zorder=25)

    elif (diff_alt>0 and diff_lon<0 and diff_lat<0):
        ax.arrow(Lon[0],Lat[0]-.2, diff_lon*c2, diff_lat*c2,
        head_width=c4, head_length=c4, fc='k', ec='k', 
        transform=ccrs.PlateCarree(), zorder=25)
        
    elif (diff_alt>0 and diff_lon>0 and diff_lat<0):
        ax.arrow(Lon[0],Lat[0]-.2, diff_lon*c2, diff_lat*c2,
        head_width=c4, head_length=c4, fc='k', ec='k', 
        transform=ccrs.PlateCarree(), zorder=25)
        
    elif (diff_alt<0 and diff_lon>0 and diff_lat>0):
        ax.arrow(Lon[0],Lat[0]-.2, diff_lon*c2, diff_lat*c2,
             head_width=c4, head_length=c4, fc='k', ec='k', 
             transform=ccrs.PlateCarree(), zorder=25)
        
    elif (diff_alt<0 and diff_lon<0 and diff_lat>0):
        ax.arrow(Lon[0],Lat[0]-.2, diff_lon*c2, diff_lat*c2,
        head_width=c4, head_length=c4, fc='k', ec='k', 
        transform=ccrs.PlateCarree(), zorder=25)
       
    elif (diff_alt<0 and diff_lon<0 and diff_lat<0):
        ax.arrow(Lon[0],Lat[0]-.2, diff_lon*c2, diff_lat*c2,
        head_width=c4, head_length=c4, fc='k', ec='k', 
        transform=ccrs.PlateCarree(), zorder=25)
        
    elif (diff_alt<0 and diff_lon>0 and diff_lat<0):
        ax.arrow(Lon[0],Lat[0]+.1, diff_lon*c2, diff_lat*c2,
        head_width=c4, head_length=c4, fc='k', ec='k', 
        transform=ccrs.PlateCarree(), zorder=25)
        
        
    return
