using NetCDF
using FFTW

filename = "foo.nc4"

# Read in data
u_raw = ncread(filename, "U")
v_raw = ncread(filename, "V")
w_raw = ncread(filename, "W")

# Linearly interpolate onto shared non-staggered grid
nx = shape(v_raw)[3]
ny = shape(w_raw)[2]
nz = shape(u_raw)[1]
xx = (0.5/nx):1/nx:(1-0.5/nx)
xs = 0:1/nx:1
yy = (0.5/ny):1/ny:(1-0.5/ny)
ys = 0:1/ny:1
ll = (0.5/nz):1/nz:(1-0.5/nz)
ls = 0:1/nz:1
U = 0.5*(u_raw[:,:,2:nx]+u_raw[:,:,1:nx-1])
V = 0.5*(v_raw[:,2:ny,:]+v_raw[:,1:ny-1,:])
W = 0.5*(w_raw[2:nz,:,:]+w_raw[1:nz-1,:,:])

# Windowing
Lw = 50/nx      # 50 grid cells ~ 150km
function fw(x, L)
    if x < 0
        return 0.0
    elseif x > L
        return 1.0
    else
        return exp(-(L/x)^2*exp(-L/(L-x)))
    end
end
wx = fw.(xx,Lw) .* fw.(1.-xx,Lw)
wy = fw.(yy,Lw) .* fw.(1.-yy,Lw)
U .*= wx.*wy
V .*= wx.*wy
W .*= wx.*wy

# High-pass filter
function high_pass!(f, Lx, Ly, Lhi)
    nz,ny,nx = size(f)
    for x = 1:nx
        for y = 1:ny/2  # Because rfft only stores half of the coefficients
            for l = 1:nz
                kx = mod(x+nx/2-1,nx)/Lx
                ky = mod(y+ny/2-1,ny)/Ly
                if sqrt(kx^2 + ky^2) < pi/Lhi
                    f[l,y,x] = 0.0
                end
            end
        end
    end
end
uf = rfft(U,dims=(2,3))
vf = rfft(V,dims=(2,3))
wf = rfft(W,dims=(2,3))
Lhi = 1000          # filter wavelength in km
Lx = nx*3           # width of domain in x
Ly = ny*3           # width of domain in y
high_pass!(uc,Lx,Ly,Lhi)
high_pass!(vc,Lx,Ly,Lhi)
high_pass!(wc,Lx,Ly,Lhi)
Uf = irfft(uf,dims=(2,3))
Vf = irfft(vf,dims=(2,3))
Wf = irfft(wf,dims=(2,3))

# Calculate Reynolds stresses
function low_pass!(f, Lx, Ly, Llo)
    nz,ny,nx = size(f)
    for x = 1:nx
        for y = 1:ny/2  # Because rfft only stores half of the coefficients
            for l = 1:nz
                kx = mod(x+nx/2-1,nx)/Lx
                ky = mod(y+ny/2-1,ny)/Ly
                if sqrt(kx^2 + ky^2) > pi/Llo
                    f[l,y,x] = 0.0
                end
            end
        end
    end
end
uuf = rfft(U.*U,dims=(2,3))
uvf = rfft(U.*V,dims=(2,3))
uwf = rfft(U.*W,dims=(2,3))
vvf = rfft(V.*V,dims=(2,3))
vwf = rfft(V.*W,dims=(2,3))
wwf = rfft(W.*W,dims=(2,3))
Llo = 100           # filter wavelength in km
low_pass!(uuf,Lx,Ly,Llo)
low_pass!(uvf,Lx,Ly,Llo)
low_pass!(uwf,Lx,Ly,Llo)
low_pass!(vvf,Lx,Ly,Llo)
low_pass!(vwf,Lx,Ly,Llo)
low_pass!(wwf,Lx,Ly,Llo)
UU = irfft(uuf,dims=(2,3))
UV = irfft(uvf,dims=(2,3))
UW = irfft(uwf,dims=(2,3))
VV = irfft(vvf,dims=(2,3))
VW = irfft(vwf,dims=(2,3))
WW = irfft(wwf,dims=(2,3))

# Extra! Extra! Calculating divergence of Reynolds stress tensor
T = ncread(filename, "T")
P = ncread(filename, "PB")
Φ_raw = ncread(filename, "PHB")
Φ = 0.5*(Φ_raw[2:nz,:,:].+Φ_raw[1:nz-1,:,:])
g = 9.81            # m/s^2
mm = 29970          # kg/mol
R = 8.315           # kg m^2/(s^2 mol K)
k = R/(g*mm)
Z = zeros((nz,ny,nx))
# dz = k T/P dΦ
for i = 2:nz
    # this is for debugging mostly
    @. Z[i,:,:] = Z[i-1,:,:] + k*(T[i,:,:]+T[i-1,:,:])/(P[i,:,:]+P[i-1,:,:])*(Φ[i,:,:]-Φ[i-1,:,:])
end
# d/dz = dl/dz d/dl = P/(k T) d/dΦ
function dZ(f, T, P, Φ, k)
    nz,ny,nx = shape(f)
    out = zeros((nz,ny,nx))
    for i = 2:nz-1
    @. out[i,:,:] = 1/k*P[i,:,:]/T[i,:,:]*(f[i+1,:,:]-f[i-1,:,:])/(Φ[i+1,:,:]-Φ[i-1,:,:])
    end
    out[1,:,:] = 1/k*P[1,:,:]/T[1,:,:]*(f[2,:,:]-f[1,:,:])/(Φ[2,:,:]-Φ[1,:,:])
    out[nz,:,:] = 1/k*P[nz,:,:]/T[nz,:,:]*(f[nz,:,:]-f[nz-1,:,:])/(Φ[nz,:,:]-Φ[nz-1,:,:])
    return out
end