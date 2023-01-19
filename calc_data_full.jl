using NetCDF
using FFTW
#using Plots
#gr()

filename = "/scratch/nr2489/riceWRF_3km/wrfout_d01_2016-08-02_00:15:00.nc"

## LENGTHSCALES
Lhi = 700000        # filter wavelength in km
Llo = 100000        # filter wavelength in km
dx = 3000           # X-gridscale for derivative
dy = 3000           # Y-gridscale for derivative
	            # Z lengthscale is calculated exactly from level height

zb = 10 # z batch size
zmax = size(ncread(filename,"U",start=[1,1,1,1],count=[1,1,-1,1]))[3]

###############DEBUGGING CODE
stz = 170
zmax -= (stz-1)
nz = zmax
ln = zb         # length of batch

nx = size(ncread(filename,"V",start=[1,1,1,1],count=[-1,1,1,1]))[1]
ny = size(ncread(filename,"W",start=[1,1,1,1],count=[1,-1,1,1]))[2]
xx = reshape((0.5/nx):1/nx:(1-0.5/nx), (nx,1))
xs = 0:1/nx:1
yy = reshape((0.5/ny):1/ny:(1-0.5/ny), (1,ny))
ys = 0:1/ny:1
ll = (0.5/nz):1/nz:(1-0.5/nz)
ls = 0:1/nz:1
function fw(x, L)
    if x < 0
        return 0.0
    elseif x > L
        return 1.0
    else
        return exp(-(L/x)^2*exp(-L/(L-x)))
    end
end
Lw = 50/nx      # 50 grid cells ~ 150km
wx = fw.(xx,Lw) .* fw.(1 .- xx,Lw)
wy = fw.(yy,Lw) .* fw.(1 .- yy,Lw)

bin = (16,16,1)     # Condense array to 100x82x179, for hor. res. of ~50km
mx = ceil(Int,nx/bin[1])
my = ceil(Int,ny/bin[2])
mz = ceil(Int,nz/bin[3])

UU = zeros((mx,my,mz))
UV = zeros((mx,my,mz))
UW = zeros((mx,my,mz))
VV = zeros((mx,my,mz))
VW = zeros((mx,my,mz))
WW = zeros((mx,my,mz))

function coarsen(f,bin)
    bx,by,bz = bin
    nx,ny,nz = size(f)
    mx = ceil(Int,nx/bx)
    my = ceil(Int,ny/by)
    mz = ceil(Int,nz/bz)
    out = zeros(mx, my, mz)
    for i = 1:mx
        ex = min(i*bx, nx)
        for j = 1:my
            ey = min(j*by, ny)
            for k = 1:mz
                ez = min(k*bz, nz)
                tmp = view(f,1+bx*(i-1):ex, 1+by*(j-1):ey, 1+bz*(k-1):ez)
                out[i,j,k] = sum(tmp)/length(tmp)
            end
        end
    end
    return out
end

function low_pass!(f, Lx, Ly, Llo)
    nx,ny,nz = size(f)
    ff = rfft(f,(1,2))

    for x = 1:div(nx,2)+1 # Because rfft only stores half of the coefficients
        kx = (x-1)
        for y = 1:ny
            ky = (mod(y-1+ny/2,ny)-ny/2)
            # Gaussian filter
            ff[x,y,:] .*= exp(-0.25*((Llo/Lx*kx)^2+(Llo/Ly*ky)^2))
            # Spectral boxcar
            #ff[x,y,:] .*= (hypot(kx/Lx, ky/Ly) > pi/Llo ? 0.0 : 1.0)
        end
    end

    f .= irfft(ff,nx,(1,2))
    return f
end
function ddx(f, dx)
    # Spectrally differentiate f in terms of x
    nx,ny,nz = size(f)
    n2x = div(nx,2)+1
    fh = rfft(f,(1,))
    kk = 1im*(0:n2x-1)	# rfft collapses the first dimension so it's only half as long
    fh .*= kk
    fh ./= dx
    return irfft(fh,nx,(1,))
end
function ddy(f, dy)
    # Spectrally differentiate f in terms of y
    nx,ny,nz = size(f)
    n2y = div(ny,2)+1
    fh = rfft(f,(2,))
    kk = 1im*(0:n2y-1)'
    fh .*= kk
    fh ./= dy
    return irfft(fh,ny,(2,))
end
function ddz(f, z)
    # Finite-difference approximation to the derivative df/dz, expressed in terms of a third indexing parameter
    nx,ny,nz = size(f)
    fh = zeros((nx,ny,nz))
    println(size(fh),size(z))
    for i = 2:nz-1
        fh[:,:,i-1] .= (f[:,:,i+1] .- f[:,:,i-1])./(z[:,:,i+1] .- z[:,:,i-1])
    end
    return fh
end

fout = "dataslice.nc4"
isfile(fout) && rm(fout)
UUatts = Dict("longname" => "UU Reynolds flux", "units" => "m^2 s^-2")
UVatts = Dict("longname" => "UV Reynolds flux", "units" => "m^2 s^-2")
UWatts = Dict("longname" => "UW Reynolds flux", "units" => "m^2 s^-2")
VVatts = Dict("longname" => "VV Reynolds flux", "units" => "m^2 s^-2")
VWatts = Dict("longname" => "VW Reynolds flux", "units" => "m^2 s^-2")
WWatts = Dict("longname" => "WW Reynolds flux", "units" => "m^2 s^-2")
Dxatts = Dict("longname" => "gw drag in x", "units" => "m s^-2")
Dyatts = Dict("longname" => "gw drag in y", "units" => "m s^-2")
Dzatts = Dict("longname" => "gw drag in z", "units" => "m s^-2")
nccreate(fout,"UU","X",mx,"Y",my,"L",mz,atts=UUatts)
nccreate(fout,"UV","X",mx,"Y",my,"L",mz,atts=UVatts)
nccreate(fout,"UW","X",mx,"Y",my,"L",mz,atts=UWatts)
nccreate(fout,"VV","X",mx,"Y",my,"L",mz,atts=VVatts)
nccreate(fout,"VW","X",mx,"Y",my,"L",mz,atts=VWatts)
nccreate(fout,"WW","X",mx,"Y",my,"L",mz,atts=WWatts)
nccreate(fout,"Dx","X",mx,"Y",my,"L",mz-2,atts=Dxatts)
nccreate(fout,"Dy","X",mx,"Y",my,"L",mz-2,atts=Dyatts)
nccreate(fout,"Dz","X",mx,"Y",my,"L",mz-2,atts=Dzatts)

########## DEBUG FILE ##########
fdebug = "debug.nc4"
isfile(fdebug) && rm(fdebug)
U_raw_atts = Dict("longname" => "U raw", "units" => "m/s")
U_cent_atts = Dict("longname" => "U on centered grid", "units" => "m/s")
U_wind_atts = Dict("longname" => "U windowed", "units" => "m/s")
U_lowpass_atts = Dict("longname" => "U low-pass filtered", "units" => "m/s")
UU_raw_atts = Dict("longname" => "UU raw", "units" => "m/s")
UU_lowpass_atts = Dict("longname" => "UU low-pass filtered", "units" => "m/s")
UU_coarse_atts = Dict("longname" => "UU coarse-grained (final)", "units" => "m/s")
nccreate(fdebug,"U_raw","Xl",nx+1,"Y",ny,atts=U_raw_atts)
nccreate(fdebug,"U_cent","X",nx,"Y",ny,atts=U_cent_atts)
nccreate(fdebug,"U_wind","X",nx,"Y",ny,atts=U_wind_atts)
nccreate(fdebug,"U_lowpass","X",nx,"Y",ny,atts=U_lowpass_atts)
nccreate(fdebug,"UU_raw","X",nx,"Y",ny,atts=UU_raw_atts)
nccreate(fdebug,"UU_lowpass","X",nx,"Y",ny,atts=UU_lowpass_atts)
nccreate(fdebug,"UU_coarse","Xs",mx,"Ys",my,atts=UU_coarse_atts)

U_raw = zeros(nx+1,ny)
U_cent = zeros(nx,ny)
U_wind = zeros(nx,ny)
U_lowpass = zeros(nx,ny)
UU_raw = zeros(nx,ny)
UU_lowpass = zeros(nx,ny)
UU_coarse = zeros(mx,my)

for zm = 1:zb:zmax
    println(zm)

    st = zm         # start -- no batch overlapping at the moment
    lr = min(ln,zmax-st)	# actual batch size, respecting array edges
    nd = st+lr-1    # end of batch

    # Read in data
    u_raw = ncread(filename, "U", start=[1,1,st,1], count=[-1,-1,lr,1])[:,:,:,1]
    v_raw = ncread(filename, "V", start=[1,1,st,1], count=[-1,-1,lr,1])[:,:,:,1]
    w_raw = ncread(filename, "W", start=[1,1,st,1], count=[-1,-1,lr+1,1])[:,:,:,1]
   
    ###
    U_raw .= u_raw[:,:,1] 
    #println("Finished reading in data")

    # Linearly interpolate onto shared non-staggered grid
    U = 0.5*(u_raw[2:nx+1,:,:]+u_raw[1:nx,:,:])
    u_raw = nothing
    V = 0.5*(v_raw[:,2:ny+1,:]+v_raw[:,1:ny,:])
    v_raw = nothing
    W = 0.5*(w_raw[:,:,2:lr+1]+w_raw[:,:,1:lr])
    w_raw = nothing

    ###
    U_cent .= U[:,:,1]
    #println("Finished interpolating data")

    # Windowing
    U .*= wx.*wy
    V .*= wx.*wy
    W .*= wx.*wy

    ###
    U_wind .= U[:,:,1]

    #println("Finished windowing")

    # Large filter to isolate gravity waves
    Lx = nx*3000        # width of domain in x
    Ly = ny*3000        # width of domain in y
    uf = copy(U); vf = copy(V); wf = copy(W)
    low_pass!(uf,Lx,Ly,Lhi)
    low_pass!(vf,Lx,Ly,Lhi)
    low_pass!(wf,Lx,Ly,Lhi)

    ###
    U_lowpass .= uf[:,:,1]
    println(size(U_lowpass))
    #println("Finished high-pass filtering")

    # Calculate Reynolds stresses

    uub = U.*U .- uf.*uf
    uvb = U.*V .- uf.*vf
    uwb = U.*W .- uf.*wf
    vvb = V.*V .- vf.*vf
    vwb = V.*W .- vf.*wf
    wwb = W.*W .- wf.*wf

    ###
    UU_raw .= uub[:,:,1]

    low_pass!(uub,Lx,Ly,Llo)
    low_pass!(uvb,Lx,Ly,Llo)
    low_pass!(uwb,Lx,Ly,Llo)
    low_pass!(vvb,Lx,Ly,Llo)
    low_pass!(vwb,Lx,Ly,Llo)
    low_pass!(wwb,Lx,Ly,Llo)

    ###
    UU_lowpass .= uub[:,:,1]

    UU[:,:,st:nd] = coarsen(uub,bin)
    UV[:,:,st:nd] = coarsen(uvb,bin)
    UW[:,:,st:nd] = coarsen(uwb,bin)
    VV[:,:,st:nd] = coarsen(vvb,bin)
    VW[:,:,st:nd] = coarsen(vwb,bin)
    WW[:,:,st:nd] = coarsen(wwb,bin)
    #println(UU[1:10,1:10,1])

    ###
    UU_coarse .= UU[:,:,1]
    #println("Finished calculating Reynolds stresses")

end

# Calculating divergence of Reynolds stress tensor
# First and last levels in Z are rubbish so we get rid of them
z_raw = ncread(filename, "PHB", start=[1,1,stz,1], count=[-1,-1,-1,1])[:,:,:,1]/9.81 # Geopotential = gz
Z = 0.5*(z_raw[:,:,2:nz+1]+z_raw[:,:,1:nz])
Z = coarsen(Z,bin)
#z_raw = nothing

Dx = ddx(UU,dx) .+ ddy(UV,dy) .+ ddz(UW,Z)
Dy = ddx(UV,dx) .+ ddy(VV,dy) .+ ddz(VW,Z)
Dz = ddx(UW,dx) .+ ddy(VW,dy) .+ ddz(WW,Z)
Dx = Dx[:,:,2:zmax-1]
Dy = Dy[:,:,2:zmax-1]
Dz = Dz[:,:,2:zmax-1]

#println("Finished calculating divergence of Reynolds stress tensor")

# Write to file

ncwrite(UU,fout,"UU")
ncwrite(UV,fout,"UV")
ncwrite(UW,fout,"UW")
ncwrite(VV,fout,"VV")
ncwrite(VW,fout,"VW")
#ncwrite(vf,fout,"V")
#ncwrite(vf,fout,"W")
ncwrite(Dx,fout,"Dx")
ncwrite(Dy,fout,"Dy")
ncwrite(Dz,fout,"Dz")
ncclose(fout)

################# UNBUG ###########
println(size(U_raw))
ncwrite(U_raw,fdebug,"U_raw")
ncwrite(U_cent,fdebug,"U_cent")
ncwrite(U_wind,fdebug,"U_wind")
ncwrite(U_lowpass,fdebug,"U_lowpass")
ncwrite(UU_raw,fdebug,"UU_raw")
ncwrite(UU_lowpass,fdebug,"UU_lowpass")
ncwrite(UU_coarse,fdebug,"UU_coarse")
