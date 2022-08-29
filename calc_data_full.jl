using NetCDF
using FFTW
#using Plots
#gr()

filename = "/scratch/nr2489/riceWRF_3km/wrfout_d01_2016-08-02_00:15:00.nc"

zb = 10 # z batch size
zmax = size(ncread(filename,"U",start=[1,1,1,1],count=[1,1,-1,1]),3)

nx = size(v_raw)[1]
ny = size(w_raw)[2]
nz = size(u_raw)[3]
n2x = div(nx,2)+1
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
wx = fw.(xx,Lw) .* fw.(1 .- xx,Lw)
wy = fw.(yy,Lw) .* fw.(1 .- yy,Lw)

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

function low_pass!(f, Lx, Ly, Llo)
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
nccreate(fout,"UU","X",nx,"Y",ny,"L",nz,atts=UUatts)
nccreate(fout,"UV","X",nx,"Y",ny,"L",nz,atts=UVatts)
nccreate(fout,"UW","X",nx,"Y",ny,"L",nz,atts=UWatts)
nccreate(fout,"VV","X",nx,"Y",ny,"L",nz,atts=VVatts)
nccreate(fout,"VW","X",nx,"Y",ny,"L",nz,atts=VWatts)
nccreate(fout,"WW","X",nx,"Y",ny,"L",nz,atts=WWatts)
nccreate(fout,"U","X",nx,"Y",ny,"L",nz,atts=Uatts)
nccreate(fout,"V","X",nx,"Y",ny,"L",nz,atts=Vatts)
nccreate(fout,"W","X",nx,"Y",ny,"L",nz,atts=Watts)

for zm = 1:zb:zmax

    # Read in data
    u_raw = ncread(filename, "U", start=[1,1,zm,1], count=[-1,-1,zb,1])[:,:,:,1]
    v_raw = ncread(filename, "V", start=[1,1,zm,1], count=[-1,-1,zb,1])[:,:,:,1]
    w_raw = ncread(filename, "W", start=[1,1,zm,1], count=[-1,-1,zb+1,1])[:,:,:,1]

    #println("Finished reading in data")

    # Linearly interpolate onto shared non-staggered grid
    U = 0.5*(u_raw[2:nx+1,:,:]+u_raw[1:nx,:,:])
    u_raw = nothing
    V = 0.5*(v_raw[:,2:ny+1,:]+v_raw[:,1:ny,:])
    v_raw = nothing
    W = 0.5*(w_raw[:,:,2:nz+1]+w_raw[:,:,1:nz])
    w_raw = nothing

    #println("Finished interpolating data")

    # Windowing
    Lw = 50/nx      # 50 grid cells ~ 150km
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

    UU = uf.*uf
    UV = uf.*vf
    UW = uf.*wf
    VV = vf.*vf
    VW = vf.*wf
    WW = wf.*wf
    Llo = 500           # filter wavelength in km
    low_pass!(UU,Lx,Ly,Llo)
    low_pass!(UV,Lx,Ly,Llo)
    low_pass!(UW,Lx,Ly,Llo)
    low_pass!(VV,Lx,Ly,Llo)
    low_pass!(VW,Lx,Ly,Llo)
    low_pass!(WW,Lx,Ly,Llo)

    #println("Finished calculating Reynolds stresses")

    # Extra! Extra! Calculating divergence of Reynolds stress tensor
    T = ncread(filename, "T", start=[1,1,zm,1], count=[-1,-1,zb,1])[:,:,:,1]
    P = ncread(filename, "PB", start=[1,1,zm,1], count=[-1,-1,zb,1])[:,:,:,1]
    Φ_raw = ncread(filename, "PHB", start=[1,1,zm,1], count=[-1,-1,zb+1,1])[:,:,:,1]
    Φ = 0.5*(Φ_raw[:,:,2:nz].+Φ_raw[:,:,1:nz-1])
    Φ_raw = nothing
    g = 9.81            # m/s^2
    mm = 29970          # kg/mol
    R = 8.315           # kg m^2/(s^2 mol K)
    k = R/(g*mm)
    Z = zeros((nx,ny,nz))
    # dz = k T/P dΦ
    for i = 2:nz
        # this is for debugging mostly
        @. Z[:,:,i] = Z[:,:,i-1] + k*(T[:,:,i]+T[:,:,i])/(P[:,:,i]+P[:,:,i-1])*(Φ[:,:,i]-Φ[:,:,i-1])
    end
    # d/dz = dl/dz d/dl = P/(k T) d/dΦ
    function dZ(f, T, P, Φ, k)
        nz,ny,nx = shape(f)
        out = zeros((nz,ny,nx))
        for i = 2:nz-1
        @. out[:,:,i] = 1/k*P[:,:,i]/T[:,:,i]*(f[:,:,i+1]-f[:,:,i-1])/(Φ[:,:,i+1]-Φ[:,:,i-1])
        end
        out[:,:,1] = 1/k*P[:,:,1]/T[:,:,1]*(f[:,:,2]-f[:,:,1])/(Φ[:,:,2]-Φ[:,:,1])
        out[:,:,nz] = 1/k*P[:,:,nz]/T[:,:,nz]*(f[:,:,nz]-f[:,:,nz-1])/(Φ[:,:,nz]-Φ[:,:,nz-1])
        return out
    end

    println("Finished calculating divergence of Reynolds stress tensor")

    # Write to file
    ncwrite(UU,fout,"UU", start=[1,1,zm,1], count=[-1,-1,zb,1])
    ncwrite(UV,fout,"UV", start=[1,1,zm,1], count=[-1,-1,zb,1])
    ncwrite(UW,fout,"UW", start=[1,1,zm,1], count=[-1,-1,zb,1])
    ncwrite(VV,fout,"VV", start=[1,1,zm,1], count=[-1,-1,zb,1])
    ncwrite(VW,fout,"VW", start=[1,1,zm,1], count=[-1,-1,zb,1])
    ncwrite(WW,fout,"WW", start=[1,1,zm,1], count=[-1,-1,zb,1])
    ncwrite(uf,fout,"U", start=[1,1,zm,1], count=[-1,-1,zb,1])
    ncwrite(vf,fout,"V", start=[1,1,zm,1], count=[-1,-1,zb,1])
    ncwrite(vf,fout,"W", start=[1,1,zm,1], count=[-1,-1,zb,1])
    ncclose(fout)
end
