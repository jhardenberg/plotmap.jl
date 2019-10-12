using PyPlot, PyCall
using NetCDF
	
"""
    plotmap(fname, var; ...)

Easy plotting of gridded datasets on a global map.
Plots variable `var` from a netcdf file `fname`.
Optional arguments are available to control the details of the plot.

# Arguments
- 'fname::String': netcdf filename containing the data to plot
- 'var'::String  : name of the variable to plot

Optional:
- 'lon::String': name of the lon variable ("lon")
- 'lat::String': name of the lat variable ("lat")
- 'lonb::String': name of the lon bounds variable ("lon_bnds")
- 'latb::String': name of the lat boundsvariable ("lat_bnds")
- 'titles':      title string ("")
- 'cstep': 	 divisions of the colorbar axis ([])
- 'cmap':        colormap ("RdBu_r")
- 'proj':        projection. One of ["platecarree", "robinson", "mollweide"]. Defaults to "platecarree".
- 'cpad':        padding (shift) of the colorbar (0.08)
- 'sub':         matplotlib subplot option (e.g. "221" for the first panel of 2x2 subplots) ("111")
- 'clabel':      label of the colorbar (defaults to the units string read from the netcdf file)
- 'cdir':        direction of the colorbar. One of ["horizontal", "vertical"]. ("horizontal")
- 'cscale':      scaling of the colorbar (0.65)
- 'cfs':         colorbar ticks font size (12)
- 'lfs':         colorbar label font size (12)
- 'tfs':         title font size (14)
- 'tpad':        padding (shift) of the title string
- 'tweight':     weight of title font. One of [ 'normal" "bold" "heavy" "light" "ultrabold" "ultralight"]. ("normal")
- 'grid':        grid spacing (defaults to [60,30]). set to empty [] to remove gridlines.
- 'region::NTuple{4,Int64}':       region to plot in format (lon1, lon2, lat1, lat2). Defaults to global.
- 'style':       one of ["pcolormesh" "contourf"]. Defaults to "pcolormesh".
- 'levels':      contour plot levels. Can be an array or a number of levels (auto)
- 'extend':      plot levels outside `levels` range. One of ["neither", "both", "min", "max"]. Default: "neither".

"""
function plotmap(fname::String, var::String; lon="lon", lat="lat", lonb="lon_bnds", latb="lat_bnds", titles="", cstep=[], cmap="RdBu_r", proj="", cpad=0.08, tpad=24, sub=111, clabel="NONE", cdir="horizontal", cscale=0.65, tfs=14, cfs=12, lfs=12, tweight="normal", grid=[60,30], region=(), style="pcolormesh", levels=0, extend="neither")

# pcolormesh needs cell boundaries
if style=="pcolormesh"
    try 
        lonb=ncread(fname, lonb);
    catch
        println("Did not find lon bounds in file, reconstructing")
        lonv=ncread(fname, lon);
        lonb=zeros(2,length(lonv))
        lonb[1,2:end]=0.5*(lonv[2:end]+lonv[1:(end-1)])
        lonb[1,1]=lonv[1]-(lonv[2]-lonv[1])*0.5
        lonb[2,end]=lonv[end]+(lonv[end]-lonv[end-1])*0.5
    end
    try
        latb=ncread(fname, latb);
    catch
        println("Did not find lat bounds in file, reconstructing")
        latv=ncread(fname, lat);
        latb=zeros(2,length(latv))
        latb[1,2:end]=0.5*(latv[2:end]+latv[1:(end-1)])
        latb[1,1]=latv[1]-(latv[2]-latv[1])*0.5
        latb[2,end]=latv[end]+(latv[end]-latv[end-1])*0.5
        if latb[1,1]>89; latb[1,1]=90 ; end 
        if latb[1,1]<-89; latb[1,1]=-90 ; end 
        if latb[2,end]>89; latb[2,end]=90 ; end 
        if latb[2,end]<-89; latb[2,end]=-90 ; end 
    end
    lonv=vcat(lonb[1,:],lonb[2,end])
    latv=vcat(latb[1,:],latb[2,end])
else
    lonv=ncread(fname, lon);
    latv=ncread(fname, lat);
end

data=ncread(fname, var);
units=ncgetatt(fname, var, "units")
if clabel=="NONE" clabel=units end

plotmap(lonv, latv, data; titles=titles, cstep=cstep, cmap=cmap, proj=proj, cpad=cpad, tpad=tpad, sub=sub, clabel=clabel, cdir=cdir, cscale=cscale, tfs=tfs, cfs=cfs, lfs=lfs, tweight=tweight, grid=grid, region=region, style=style, levels=levels, extend=extend)

end

"""
    plotmap(lon, lat, data; ...)

Easy plotting of gridded datasets on a global map.
Plots data in 2D array `data` with longitudes `lon` and latitudes `lat`.
Optional arguments are available to control the details of the plot.

# Arguments
- 'data::Array{Float32,2}': data to plot. If a 3D array is passed, only the first frame is plotted: `data[:,:,1]`
- 'lon::Array{Float64,1}' : longitudes
- 'lat::Array{Float64,1}' : latitudes

Optional:
- 'titles':      title string ("")
- 'cstep':       divisions of the colorbar axis ([])
- 'cmap':        colormap ("RdBu_r")
- 'proj':        projection. One of ["platecarree", "robinson", "mollweide"]. Defaults to "platecarree".
- 'cpad':        padding (shift) of the colorbar (0.08)
- 'sub':         matplotlib subplot option (e.g. "221" for the first panel of 2x2 subplots) ("111")
- 'clabel':      label of the colorbar (defaults to the units string read from the netcdf file)
- 'cdir':        direction of the colorbar. One of ["horizontal", "vertical"]. ("horizontal")
- 'cscale':      scaling of the colorbar (0.65)
- 'cfs':         colorbar ticks font size (12)
- 'lfs':         colorbar label font size (12)
- 'tfs':         title font size (14)
- 'tpad':        padding (shift) of the title string
- 'tweight':     weight of title font. One of [ 'normal" "bold" "heavy" "light" "ultrabold" "ultralight"]. ("normal")
- 'grid':        grid spacing (defaults to [60,30]). set to empty [] to remove gridlines.
- 'region::NTuple{4,Int64}':       region to plot in format (lon1, lon2, lat1, lat2). Defaults to global.
- 'style':       one of ["pcolormesh" "contourf"]. Defaults to "pcolormesh".
- 'levels':      contour plot levels. Can be an array or a number of levels (auto)
- 'extend':      plot levels outside `levels` range. One of ["neither", "both", "min", "max"]. Default: "neither".

"""
function plotmap(lon, lat, data; titles="", cstep=[], cmap="RdBu_r", proj="", cpad=0.08, tpad=24, sub=111, clabel="", cdir="horizontal", cscale=0.65, tfs=14, cfs=12, lfs=12, tweight="normal", grid=[60,30], region=(), style="pcolormesh", levels=0, extend="neither")

#cf = pyimport("cartopy.feature")

dd = size(data)

if length(dd)==3 data=data[:,:,1] end
# Needed below for correct plotting
if style=="pcolormesh"
    if length(lon)==(dd[1]+1) data=data' end
else
    if length(lon)==dd[1] data=data' end
end

ccrs = pyimport("cartopy.crs")
cutil = pyimport("cartopy.util")

if proj=="robinson"
    proj=ccrs.Robinson()
    dlabels=false
elseif proj == "mollweide"
    proj=ccrs.Mollweide()
    dlabels=false
else
    proj=ccrs.PlateCarree()
    dlabels=true
end

ax = subplot(sub, projection=proj)
if length(region)>0 ax.set_extent(region, crs=ccrs.PlateCarree()) end
ax.coastlines()
xlocvec=vcat(-vcat(grid[1]:grid[1]:180)[end:-1:1], vcat(0:grid[1]:180))
ylocvec=vcat(-vcat(grid[2]:grid[2]:90)[end:-1:1], vcat(0:grid[2]:90))

if dlabels
    ax.gridlines(linewidth=1, color="gray", alpha=0.5, linestyle="--", draw_labels=true, xlocs=xlocvec, ylocs=ylocvec)
else
    ax.gridlines(linewidth=1, color="gray", alpha=0.5, linestyle="--", xlocs=xlocvec, ylocs=ylocvec)
end

if style=="contourf"
    data_cyc, lon_cyc = cutil.add_cyclic_point(data, coord=lon)
else
    data_cyc=data
    lon_cyc=lon
end

if style=="contourf"
    if levels==0 
        contourf(lon_cyc, lat, data_cyc, transform=ccrs.PlateCarree(), cmap=cmap, extend=extend)
    else
        contourf(lon_cyc, lat, data_cyc, transform=ccrs.PlateCarree(), cmap=cmap, levels=levels, extend=extend)
    end
else
    pcolormesh(lon_cyc, lat, data_cyc, transform=ccrs.PlateCarree(), cmap=cmap)
end

if length(titles)>1 title(titles, pad=tpad, fontsize=tfs, weight=tweight) end
cbar=colorbar(orientation=cdir, extend="both", pad=cpad, label=clabel, shrink=cscale)
cbar.set_label(label=clabel,size=lfs)
cbar.ax.tick_params(labelsize=cfs) 
if length(cstep)>0 cbar.set_ticks(cstep) end
tight_layout()

end
                  
