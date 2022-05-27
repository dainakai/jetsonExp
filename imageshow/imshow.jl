using GLMakie, Random
using Makie.Colors

cSize = 1024; rSize = 1024

function make_image(cSize, rSize)
    img = rand(cSize,rSize)
    return img
end

imgObservable = Observable(make_image(cSize, rSize))
imgObservable1 = Observable(make_image(cSize, rSize))

f = Figure(resolution = (2*cSize, rSize), figure_padding = 1, backgroundcolor = :gray80)
ax = Axis(f[1, 1], aspect = DataAspect(), yreversed = false, title = "Camera 1")
ax1 = Axis(f[1, 2], aspect = DataAspect(), yreversed = false, title = "Camera 2")

image!(ax, imgObservable)
image!(ax1, imgObservable1)

# show the plot here, just for this interactive example
display(f)
# f

for i in 1:100
    # update the image observable and therefore the plot
    imgObservable[] = make_image(cSize, rSize)
    imgObservable1[] = make_image(cSize, rSize)
    # wait to be able to see the changing plot
    sleep(0.01)
    # you could save the figure here if you wanted to
end
