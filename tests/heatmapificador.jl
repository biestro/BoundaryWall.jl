using Images
using CairoMakie
using FFTW

img = Images.load("./tests/img1.jpeg")

img_gray = Gray.(img)

img_arr = Float64.(img_gray)

let
  k = 3
  _c = :afmhot
  x = LinRange(0,1,10)
f=Figure(theme=theme_minimal())
ax=CairoMakie.Axis(f[1,1], yreversed=true)
hm=heatmap!(ax,img_arr[1:k:end, 1:k:end]', colormap=_c)
Colorbar(f[2,1], vertical=false,colormap=_c, colorrange=(0,1), flipaxis=false, ticks=0:0.2:1)
# , ticks=(x, ["min";fill(string(), length(x)-2);"max"]),
rowgap!(f.layout, -50)
# colsize!(f.layout, 1, Aspect(1,1))
ax.aspect=DataAspect()
hidedecorations!(ax)
hidespines!(ax)
# resize_to_layout!(f)
# save("./tests/hard_$_c.png", f, px_per_unit = 2)
f
end

img_fft = fft(img_arr)


filter(x::ComplexF64) = exp(-x^2)

sx,sy = size(img_fft)
w = 40


heatmap(log10.(abs2.(fftshift(img_fft)).+1))
filt = [filter(x/w)*filter(y/w) for x in -sx/2:sx/2-1, y in -sy/2:sy/2-1]
# img_new = img_fft .* filt
img_fft[125:end-125, 125:end-125] .= 0.0


let 
  fig = Figure()
  ax = Makie.Axis(fig[1,1],yreversed=true)
  heatmap!(ax,abs2.(ifft(img_fft)'))
  fig
end
