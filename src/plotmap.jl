using PyPlot, PyCall
using NetCDF
pyimport_conda("cartopy.crs", "cartopy", "conda-forge")
ssl = pyimport_conda("ssl","ssl")
ssl._create_default_https_context = ssl._create_unverified_context
	
"""
    plotmap(fname, var; ...)

Easy plotting of gridded datasets on a global map.
Plots variable `var` from a netcdf file `fname`.
Optional arguments are available to control the details of the plot.

# Arguments
* data::Array{Float32,2} : data to plot. If a 3D array is passed, only the first frame is plotted: `data[:,:,1]`
* lon::Array{Float64,1} : longitudes
* lat::Array{Float64,1} : latitudes

Optional:
* title :       title string ("")
* cstep :       divisions of the colorbar axis ([])
* cmap :        colormap ("RdBu_r")
* proj :        projection. One of ["platecarree", "robinson", "mollweide"]. Defaults to "platecarree".
* cpad :        padding (shift) of the colorbar (0.08)
* sub :         matplotlib subplot option (e.g. "221" for the first panel of 2x2 subplots) ("111")
* clabel :      label of the colorbar (defaults to the units string read from the netcdf file)
* cdir :        direction of the colorbar. One of ["horizontal", "vertical"]. ("horizontal")
* cscale :      scaling of the colorbar (0.65)
* cfs :         colorbar ticks font size (12)
* lfs :         colorbar label font size (12)
* tfs :         title font size (14)
* borders :     national borders (false)
* tpad :        padding (shift) of the title string
* tweight :     weight of title font. One of [ 'normal" "bold" "heavy" "light" "ultrabold" "ultralight"]. ("normal")
* grid :        grid spacing (defaults to [60,30]). set to empty [] to remove gridlines.
* region::NTuple{4,Int64}:       region to plot in format (lon1, lon2, lat1, lat2). Defaults to global.
* style :       one of ["pcolormesh" "contourf"]. Defaults to "pcolormesh".
* levels :      contour plot levels. Can be an array or a number of levels (auto)
* extend :      plot levels outside `levels` range. One of ["neither", "both", "min", "max"]. Default: "neither".

"""
function plotmap(fname::String, var::String; lon="lon", lat="lat", lonb="lon_bnds", latb="lat_bnds", titles="", title="", cstep=[], cmap="RdBu_r", proj="", cpad=0.08, tpad=24, sub=111, clabel="NONE", cdir="horizontal", cscale=0.65, tfs=14, cfs=12, lfs=12, tweight="normal", grid=[60,30], region=(), style="pcolormesh", levels=0, extend="neither", borders=false)

# pcolormesh needs cell boundaries
if style=="pcolormesh"
    try 
        lonb=ncread(fname, lonb);
        lonv=vcat(lonb[1,:],lonb[2,end])
    catch
        lonv=ncread(fname, lon);
    end
    try
        latb=ncread(fname, latb);
        latv=vcat(latb[1,:],latb[2,end])
    catch
        latv=ncread(fname, lat);
    end
else
    lonv=ncread(fname, lon);
    latv=ncread(fname, lat);
end

data=ncread(fname, var);
units=ncgetatt(fname, var, "units")
if clabel=="NONE" clabel=units end

plotmap(lonv, latv, data; title=title, titles=titles, cstep=cstep, cmap=cmap, proj=proj, cpad=cpad, tpad=tpad, sub=sub, clabel=clabel, cdir=cdir, cscale=cscale, tfs=tfs, cfs=cfs, lfs=lfs, tweight=tweight, grid=grid, region=region, style=style, levels=levels, extend=extend, borders=borders)

end

"""
    plotmap(lon, lat, data; ...)

Easy plotting of gridded datasets on a global map.
Plots data in 2D array `data` with longitudes `lon` and latitudes `lat`.
Optional arguments are available to control the details of the plot.

# Arguments
* data::Array{Float32,2} : data to plot. If a 3D array is passed, only the first frame is plotted: `data[:,:,1]`
* lon::Array{Float64,1} : longitudes
* lat::Array{Float64,1} : latitudes

Optional:
* title :       title string ("")
* cstep :       divisions of the colorbar axis ([])
* cmap :        colormap ("RdBu_r")
* proj :        projection. One of ["platecarree", "robinson", "mollweide"]. Defaults to "platecarree".
* cpad :        padding (shift) of the colorbar (0.08)
* sub :         matplotlib subplot option (e.g. "221" for the first panel of 2x2 subplots) ("111")
* clabel :      label of the colorbar (defaults to the units string read from the netcdf file)
* cdir :        direction of the colorbar. One of ["horizontal", "vertical"]. ("horizontal")
* cscale :      scaling of the colorbar (0.65)
* cfs :         colorbar ticks font size (12)
* lfs :         colorbar label font size (12)
* tfs :         title font size (14)
* borders :     national borders (false)
* tpad :        padding (shift) of the title string
* tweight :     weight of title font. One of [ 'normal" "bold" "heavy" "light" "ultrabold" "ultralight"]. ("normal")
* grid :        grid spacing (defaults to [60,30]). set to empty [] to remove gridlines.
* region::NTuple{4,Int64}:       region to plot in format (lon1, lon2, lat1, lat2). Defaults to global.
* style :       one of ["pcolormesh" "contourf"]. Defaults to "pcolormesh".
* levels :      contour plot levels. Can be an array or a number of levels (auto)
* extend :      plot levels outside `levels` range. One of ["neither", "both", "min", "max"]. Default: "neither".

"""
function plotmap(lon, lat, data; title="", titles="", cstep=[], cmap="RdBu_r", proj="", cpad=0.08, tpad=24, sub=111, clabel="", cdir="horizontal", cscale=0.65, tfs=14, cfs=12, lfs=12, tweight="normal", grid=[60,30], region=(), style="pcolormesh", levels=0, extend="neither", borders=false)

#cf = pyimport("cartopy.feature")

dd = size(data)

if length(dd)==3 data=data[:,:,1] end
if style=="pcolormesh"
    if length(lon) in dd
        #println("pcolormesh needs cell boundaries, reconstructing lon")
        lonb=zeros(2,length(lon))
        lonb[1,2:end]=0.5*(lon[2:end]+lon[1:(end-1)])
        lonb[1,1]=lon[1]-(lon[2]-lon[1])*0.5
        lonb[2,end]=lon[end]+(lon[end]-lon[end-1])*0.5
        lon=vcat(lonb[1,:],lonb[2,end])
    end
    if length(lat) in dd
        #println("pcolormesh needs cell boundaries, reconstructing lat")
        latb=zeros(2,length(lat))
        latb[1,2:end]=0.5*(lat[2:end]+lat[1:(end-1)])
        latb[1,1]=lat[1]-(lat[2]-lat[1])*0.5
        latb[2,end]=lat[end]+(lat[end]-lat[end-1])*0.5
        if latb[1,1]>89; latb[1,1]=90 ; end
        if latb[1,1]<-89; latb[1,1]=-90 ; end
        if latb[2,end]>89; latb[2,end]=90 ; end
        if latb[2,end]<-89; latb[2,end]=-90 ; end
        lat=vcat(latb[1,:],latb[2,end])
    end
    if length(lon)==(dd[1]+1) data=data' end
else
    if length(lon)==dd[1] data=data' end
end

ccrs = pyimport("cartopy.crs")
cutil = pyimport("cartopy.util")
cfeat = pyimport("cartopy.feature")                  

cm=0
if length(region)>0
    cm=(region[1]+region[2])/2
end
if proj=="robinson"
    proj=ccrs.Robinson()
    dlabels=false
elseif proj == "mollweide"
    proj=ccrs.Mollweide()
    dlabels=false
elseif proj == "polar_north"
    proj=ccrs.NorthPolarStereo()
    region = (-180, 180, 45, 90)
    dlabels=false
elseif proj == "polar_south"
    proj=ccrs.SouthPolarStereo()
    region = (-180, 180, -90, -45)
    dlabels=false
else
    proj=ccrs.PlateCarree(central_longitude=cm)
    dlabels=true
end

ax = subplot(sub, projection=proj)
if length(region)>0 ax.set_extent(region, crs=ccrs.PlateCarree()) end
ax.coastlines()
if borders ax.add_feature(cfeat.BORDERS) end
xlocvec=vcat(-vcat(grid[1]:grid[1]:180)[end:-1:1], vcat(0:grid[1]:180))
ylocvec=vcat(-vcat(grid[2]:grid[2]:90)[end:-1:1], vcat(0:grid[2]:90))

if dlabels
    ax.gridlines(linewidth=1, color="gray", alpha=0.5, linestyle="--", draw_labels=true, xlocs=xlocvec, ylocs=ylocvec)
else
    ax.gridlines(linewidth=1, color="gray", alpha=0.5, linestyle="--", xlocs=xlocvec, ylocs=ylocvec)
end

if style=="contourf"
    data_cyc, lon_cyc = cutil.add_cyclic_point(data, coord=lon)
    if levels==0 
        contourf(lon_cyc, lat, data_cyc, transform=ccrs.PlateCarree(), cmap=cmap, extend=extend)
    else
        contourf(lon_cyc, lat, data_cyc, transform=ccrs.PlateCarree(), cmap=cmap, levels=levels, extend=extend)
    end
else
    pcolormesh(lon, lat, data, transform=ccrs.PlateCarree(), cmap=cmap)
end

if length(cstep)>0 clim(cstep[1],cstep[end]); end
if length(title)>1 PyPlot.title(title, pad=tpad, fontsize=tfs, weight=tweight) end
if length(titles)>1 PyPlot.title(titles, pad=tpad, fontsize=tfs, weight=tweight) end
cbar=colorbar(orientation=cdir, extend="both", pad=cpad, label=clabel, shrink=cscale)
cbar.set_label(label=clabel,size=lfs)
cbar.ax.tick_params(labelsize=cfs) 
if length(cstep)>0 cbar.set_ticks(cstep) end
tight_layout()

end
