using NetCDF
using FFTW
#using Plots
#gr()

filename = "/scratch/nr2489/riceWRF_3km/wrfout_d01_2016-08-02_00:15:00.nc"

zb = 10 # z batch size
zmax = size(ncread(filename,"U",start=[1,1,1,1],count=[1,1,-1,1]))[3]

st = zm         # start -- no batch overlapping at the moment
ln = zb         # length of batch
nd = st+ln      # end of batch
nx = size(ncread(filename,"V",start=[1,1,1,1],count=[-1,1,1,1]))[1]
ny = size(ncread(filename,"W",start=[1,1,1,1],count=[1,-1,1,1]))[2]
n2x = div(nx,2)+1
n2y = div(ny,2)+1
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
mx = ceil(nx/bin[1])
my = ceil(ny/bin[2])
mz = ceil(nz/bin[3])
UU = zeros(mx,my,mz)
UV = zeros(mx,my,mz)
UW = zeros(mx,my,mz)
VV = zeros(mx,my,mz)
VW = zeros(mx,my,mz)
WW = zeros(mx,my,mz)

function high_pass!(f, Lx, Ly, Lhi)
    nx,ny,nz = size(f)
    ff = rfft(f,(1,2))
    for x = 1:div(nx,2) # Because rfft only stores half of the coefficients
        for y = 1:ny
            for l = 1:nz
                kx = x/Lx
                ky = (mod(y+ny/2,ny)-ny/2)/Ly
                if sqrt(kx^2 + ky^2) < pi/Lhi
                    ff[x,y,l] = 0.0
                end
            end
        end
    end
    f .= irfft(ff,nx,(1,2))
end

function low_pass!(f, Lx, Ly, Llo, bin)
    bx,by,bz = bin
    nx,ny,nz = size(f)
    ff = rfft(f,(1,2))
    for x = 1:div(nx,2) # Because rfft only stores half of the coefficients
        for y = 1:ny
            for l = 1:nz
                kx = x/Lx
                ky = (mod(y+ny/2,ny)-ny/2)/Ly
                if sqrt(kx^2 + ky^2) > pi/Llo
                    ff[x,y,l] = 0.0
                end
            end
        end
    end
    f .= irfft(ff,nx,(1,2))
    mx = ceil(nx/bx)
    my = ceil(ny/by)
    mz = ceil(nz/bz)
    out = zeros(mx, my, mz)
    for i = 1:mx
        ex = min((bx+1)*i, nx)
        for j = 1:my
            ey = min((by+1)*j, ny)
            for k = 1:mz
                # Coarsening resolution now that high frequencies have been removed
                ez = min((bz+1)*k, nz)
                tmp = view(f,(bx*i:ex, by*j:ey, bz*k:ez))
                out[i,j,k] = sum(tmp)/length(tmp)
            end
        end
        return out
    end
end
function ddx(f, dx)
    # Spectrally differentiate f in terms of x
    nx,ny,nz = size(f)
    fh = rfft(f,(1,))
    kk = [1im*mod(k-n2x,nx)+n2x for k=1:nx]
    fh .*= kk
    fh ./= dx
    return irfft(fh,nx,(1,))
end
function ddy(f, dy)
    # Spectrally differentiate f in terms of y
    nx,ny,nz = size(f)
    fh = rfft(f,(2,))
    kk = [1im*mod(k-n2y,ny)+n2y for k=1:ny]
    fh .*= kk
    fh ./= dy
    return irfft(fh,ny,(2,))
end
function ddz(f, z)
    # Finite-difference approximation to the derivative df/dz, expressed in terms of a third indexing parameter
    nx,ny,nz = size(f)
    fh = zeros((nx,ny,nz))
    for i = 2:nz-1
        fh[:,:,i-1] .= (f[:,:,i+1] .- f[:,:,i-1])./(z[:,:,i+1] .- z[:,:,i-1])
    end
end

fout = "dataslice.nc4"
isfile(fout) && rm(fout)
UUatts = Dict("longname" => "UU Reynolds flux", "units" => "m^2 s^-2")
UVatts = Dict("longname" => "UV Reynolds flux", "units" => "m^2 s^-2")
UWatts = Dict("longname" => "UW Reynolds flux", "units" => "m^2 s^-2")
VVatts = Dict("longname" => "VV Reynolds flux", "units" => "m^2 s^-2")
VWatts = Dict("longname" => "VW Reynolds flux", "units" => "m^2 s^-2")
WWatts = Dict("longname" => "WW Reynolds flux", "units" => "m^2 s^-2")
Uatts  = Dict("longname" => "U gravity wave field", "units" => "m/s")
Vatts  = Dict("longname" => "V gravity wave field", "units" => "m/s")
Watts  = Dict("longname" => "W gravity wave field", "units" => "m/s")
Dxatts = Dict("longname" => "gw drag in x", "units" => "m s^-2")
Dyatts = Dict("longname" => "gw drag in y", "units" => "m s^-2")
Dzatts = Dict("longname" => "gw drag in z", "units" => "m s^-2")
nccreate(fout,"UU","X",nx,"Y",ny,"L",zmax,atts=UUatts)
nccreate(fout,"UV","X",nx,"Y",ny,"L",zmax,atts=UVatts)
nccreate(fout,"UW","X",nx,"Y",ny,"L",zmax,atts=UWatts)
nccreate(fout,"VV","X",nx,"Y",ny,"L",zmax,atts=VVatts)
nccreate(fout,"VW","X",nx,"Y",ny,"L",zmax,atts=VWatts)
nccreate(fout,"WW","X",nx,"Y",ny,"L",zmax,atts=WWatts)
nccreate(fout,"U","X",nx,"Y",ny,"L",zmax,atts=Uatts)
nccreate(fout,"V","X",nx,"Y",ny,"L",zmax,atts=Vatts)
nccreate(fout,"W","X",nx,"Y",ny,"L",zmax,atts=Watts)
nccreate(fout,"Dx","X",nx,"Y",ny,"L",zmax-2,atts=Dxatts)
nccreate(fout,"Dy","X",nx,"Y",ny,"L",zmax-2,atts=Dyatts)
nccreate(fout,"Dz","X",nx,"Y",ny,"L",zmax-2,atts=Dzatts)

for zm = 1:zb:zmax

    # Read in data
    u_raw = ncread(filename, "U", start=[1,1,st,1], count=[-1,-1,ln,1])[:,:,:,1]
    v_raw = ncread(filename, "V", start=[1,1,st,1], count=[-1,-1,ln,1])[:,:,:,1]
    w_raw = ncread(filename, "W", start=[1,1,st,1], count=[-1,-1,ln+1,1])[:,:,:,1]
    z_raw = ncread(filename, "PHB", start=[1,1,st,1], count=[-1,-1,ln+1,1])[:,:,:,1]

    #println("Finished reading in data")

    # Linearly interpolate onto shared non-staggered grid
    U = 0.5*(u_raw[2:nx+1,:,:]+u_raw[1:nx,:,:])
    u_raw = nothing
    V = 0.5*(v_raw[:,2:ny+1,:]+v_raw[:,1:ny,:])
    v_raw = nothing
    W = 0.5*(w_raw[:,:,2:ln+1]+w_raw[:,:,1:ln])
    w_raw = nothing
    Z = 0.5*(z_raw[:,:,2:ln+1]+z_raw[:,:,1:ln])
    z_raw = nothing

    #println("Finished interpolating data")

    # Windowing
    U .*= wx.*wy
    V .*= wx.*wy
    W .*= wx.*wy

    #println("Finished windowing")

    # High-pass filter
    Lhi = 2000          # filter wavelength in km
    Lx = nx*3           # width of domain in x
    Ly = ny*3           # width of domain in y
    uf = U; vf = V; wf = W
    high_pass!(uf,Lx,Ly,Lhi)
    high_pass!(vf,Lx,Ly,Lhi)
    high_pass!(wf,Lx,Ly,Lhi)

    #println("Finished high-pass filtering")

    # Calculate Reynolds stresses

    uub = uf.*uf
    uvb = uf.*vf
    uwb = uf.*wf
    vvb = vf.*vf
    vwb = vf.*wf
    wwb = wf.*wf
    Llo = 500           # filter wavelength in km
    UU[:,:,st:nd] = low_pass!(uub,Lx,Ly,Llo,bin)
    UV[:,:,st:nd] = low_pass!(uvb,Lx,Ly,Llo,bin)
    UW[:,:,st:nd] = low_pass!(uwb,Lx,Ly,Llo,bin)
    VV[:,:,st:nd] = low_pass!(vvb,Lx,Ly,Llo,bin)
    VW[:,:,st:nd] = low_pass!(vwb,Lx,Ly,Llo,bin)
    WW[:,:,st:nd] = low_pass!(wwb,Lx,Ly,Llo,bin)

    #println("Finished calculating Reynolds stresses")

end

# Extra! Extra! Calculating divergence of Reynolds stress tensor
# First and last levels in Z are rubbish so we get rid of them
Dx = ddx(UU,dx) .+ ddy(UV,dy) .+ ddz(UW,Z)
Dy = ddx(UV,dx) .+ ddy(VV,dy) .+ ddz(VW,Z)
Dz = ddx(UW,dx) .+ ddy(VW,dy) .+ ddz(WW,Z)
Dx = Dx[:,:,2:zmax-1]
Dy = Dy[:,:,2:zmax-1]
Dz = Dz[:,:,2:zmax-1]

#println("Finished calculating divergence of Reynolds stress tensor")

# Write to file
println(zm)
ncwrite(UU,fout,"UU")
ncwrite(UV,fout,"UV")
ncwrite(UW,fout,"UW")
ncwrite(VV,fout,"VV")
ncwrite(VW,fout,"VW")
ncwrite(WW,fout,"WW")
ncwrite(uf,fout,"U")
ncwrite(vf,fout,"V")
ncwrite(vf,fout,"W")
ncwrite(Dx,fout,"Dx")
ncwrite(Dy,fout,"Dy")
ncwrite(Dz,fout,"Dz")
ncclose(fout)